//! Volume rendering via custom WGSL ray-march shader on a proxy cube.
//!
//! Uploads the 3D scalar field as a normalized R16Float 3D texture, then
//! renders a bounding-box cube with a custom material whose fragment shader
//! ray-marches through it in fractional-coordinate space (handles non-
//! orthogonal cells) with front-to-back compositing.

// The ShaderType derive emits per-field `fn check` helpers that are unused
// when the type is only written to uniform buffers via AsBindGroup.
#![allow(dead_code)]

use bevy::prelude::*;
use bevy::pbr::{Material, MaterialPlugin};
use bevy::render::render_asset::RenderAssetUsages;
use bevy::render::render_resource::{
    AsBindGroup, Extent3d, ShaderRef, ShaderType, TextureDimension, TextureFormat,
};
use crate::crystal::Lattice;
use crate::CubeResource;
use crate::Pes3dState;
use crate::VolumeProxy;

/// Uniform parameters shared between CPU state and the WGSL shader.
#[derive(ShaderType, Clone, Copy, PartialEq)]
pub struct VolumeParams {
    pub color_min: f32,
    pub color_max: f32,
    pub iso_value: f32,     // alpha falloff reference point (eV)
    pub alpha_scale: f32,   // opacity multiplier for low-energy regions
    pub alpha_falloff: f32, // gaussian decay rate above iso_value
    pub clip_x: Vec2,       // fractional clipping range [min, max]
    pub clip_y: Vec2,
    pub clip_z: Vec2,
    pub inv_lattice: Mat3,  // world → fractional coordinates
    pub steps: u32,         // ray-march step count
}

impl Default for VolumeParams {
    fn default() -> Self {
        Self {
            color_min: 0.0,
            color_max: 1.0,
            iso_value: 0.15,
            alpha_scale: 0.8,
            alpha_falloff: 0.25,
            clip_x: Vec2::new(0.0, 1.0),
            clip_y: Vec2::new(0.0, 1.0),
            clip_z: Vec2::new(0.0, 1.0),
            inv_lattice: Mat3::IDENTITY,
            steps: 128,
        }
    }
}

/// Custom material: 3D volume texture + ray-march fragment shader.
#[derive(Asset, AsBindGroup, TypePath, Clone)]
pub struct VolumeMaterial {
    #[texture(0, dimension = "3d")]
    #[sampler(1)]
    pub volume: Handle<Image>,
    #[uniform(2)]
    pub params: VolumeParams,
}

/// Weak handle for the inline ray-march shader (inserted into Assets<Shader>
/// at plugin build time — ShaderRef::from(&str) would treat it as an asset
/// path, so the source is embedded via include_str! + from_wgsl instead).
pub const VOLUME_SHADER_HANDLE: Handle<Shader> =
    Handle::weak_from_u128(0x56e8c11a9a4f4f01);

impl Material for VolumeMaterial {
    fn fragment_shader() -> ShaderRef {
        VOLUME_SHADER_HANDLE.into()
    }
    fn alpha_mode(&self) -> AlphaMode {
        AlphaMode::Blend
    }
}

pub struct VolumePlugin;

impl Plugin for VolumePlugin {
    fn build(&self, app: &mut App) {
        let shader = Shader::from_wgsl(
            include_str!("shaders/volume.wgsl"),
            "volume.wgsl",
        );
        app.world_mut()
            .resource_mut::<Assets<Shader>>()
            .insert(VOLUME_SHADER_HANDLE.id(), shader);
        app.add_plugins(MaterialPlugin::<VolumeMaterial> {
            prepass_enabled: false,
            shadows_enabled: false,
            ..default()
        });
    }
}

/// Upload the 3D energy field as a normalized R16Float 3D texture.
///
/// Values are normalized to [0, 1] using [color_min, color_max] on the CPU
/// (R16Float must store [0,1] data — R32Float linear filtering is not
/// guaranteed on all adapters). NaN / non-finite → 1.0 (high energy → faint).
pub fn build_volume_texture(
    field: &[f32], nx: usize, ny: usize, nz: usize,
    color_min: f32, color_max: f32,
) -> Image {
    let range = (color_max - color_min).max(1e-6);
    let mut halfs: Vec<u8> = Vec::with_capacity(field.len() * 2);
    for &v in field {
        let t = if v.is_finite() {
            ((v - color_min) / range).clamp(0.0, 1.0)
        } else {
            1.0  // NaN / inf → high energy → transparent
        };
        let h = half::f16::from_f32(t);
        halfs.extend_from_slice(&h.to_bits().to_le_bytes());
    }
    Image::new(
        Extent3d { width: nx as u32, height: ny as u32, depth_or_array_layers: nz as u32 },
        TextureDimension::D3,
        halfs,
        TextureFormat::R16Float,
        RenderAssetUsages::RENDER_WORLD,
    )
}

/// Build the shader uniform struct from Pes3dState + lattice.
pub fn volume_params_from_state(state: &Pes3dState, lattice: &Lattice) -> VolumeParams {
    let inv = lattice.inverse_vectors();
    // inverse_vectors returns column vectors; glam Mat3 is column-major, so
    // the columns map directly: Mat3::from_cols(col0, col1, col2).
    VolumeParams {
        color_min: state.color_min,
        color_max: state.color_max,
        iso_value: state.vol_iso_ref,
        alpha_scale: state.alpha_scale,
        alpha_falloff: state.alpha_falloff,
        clip_x: Vec2::new(state.clip_x[0], state.clip_x[1]),
        clip_y: Vec2::new(state.clip_y[0], state.clip_y[1]),
        clip_z: Vec2::new(state.clip_z[0], state.clip_z[1]),
        inv_lattice: Mat3::from_cols(inv[0], inv[1], inv[2]),
        steps: state.vol_steps,
    }
}

/// Keep the volume texture and shader params in sync with Pes3dState.
///
/// - color_min/max change → rebuild the normalized texture (new handle)
/// - any other param change → update the uniform struct only
/// - Uses Local caches (NOT is_changed(), which is true every frame because
///   egui constructs its widgets with &mut references → DerefMut) so that
///   texture/uniform updates happen only when parameters actually change.
pub fn update_volume(
    pes3d_state: Option<Res<Pes3dState>>,
    cube: Option<Res<CubeResource>>,
    mut images: ResMut<Assets<Image>>,
    mut materials: ResMut<Assets<VolumeMaterial>>,
    proxy_q: Query<&MeshMaterial3d<VolumeMaterial>, With<VolumeProxy>>,
    mut prev_color: Local<(f32, f32)>,
    mut prev_params: Local<VolumeParams>,
) {
    let Some(ps) = pes3d_state.as_ref() else { return };
    let Some(cube) = cube.as_ref() else { return };

    // Rebuild texture only when the color mapping range actually changes
    let color_key = (ps.color_min, ps.color_max);
    let color_changed = *prev_color != color_key;

    let new_params = volume_params_from_state(ps, &cube.0.to_lattice());
    if !color_changed && *prev_params == new_params {
        return;  // nothing changed — keep the material untouched
    }

    for handle in proxy_q.iter() {
        let Some(mat) = materials.get_mut(&handle.0) else { continue };
        if color_changed {
            let img = build_volume_texture(
                &cube.0.field, cube.0.nx, cube.0.ny, cube.0.nz,
                ps.color_min, ps.color_max,
            );
            mat.volume = images.add(img);
        }
        mat.params = new_params;
        // Mut<VolumeMaterial> automatically marks the asset changed on drop
    }
    *prev_params = new_params;
    if color_changed {
        *prev_color = color_key;
    }
}

/// Spawn a unit cube (proxy geometry) that fills the lattice volume.
pub fn volume_proxy_mesh(lattice: &Lattice) -> Mesh {
    let vecs = lattice.to_vectors();
    // 8 corners of the fractional cell mapped to Cartesian
    let corners: [Vec3; 8] = [
        Vec3::ZERO,                              // (0,0,0)
        vecs[0],                                 // (1,0,0)
        vecs[0] + vecs[1],                       // (1,1,0)
        vecs[1],                                 // (0,1,0)
        vecs[2],                                 // (0,0,1)
        vecs[0] + vecs[2],                       // (1,0,1)
        vecs[0] + vecs[1] + vecs[2],             // (1,1,1)
        vecs[1] + vecs[2],                       // (0,1,1)
    ];

    let positions: Vec<[f32; 3]> = corners.iter().map(|c| [c.x, c.y, c.z]).collect();
    let normals: Vec<[f32; 3]> = vec![[0.0, 1.0, 0.0]; 8];
    let uvs: Vec<[f32; 2]> = vec![[0.0, 0.0]; 8];

    // 6 faces, 2 triangles each, 12 triangles, both windings
    let faces: [(usize, usize, usize, usize); 6] = [
        (0, 1, 2, 3), // bottom  (-Z face)
        (4, 5, 6, 7), // top     (+Z face)
        (0, 3, 7, 4), // left    (-X face)
        (1, 2, 6, 5), // right   (+X face)
        (0, 1, 5, 4), // front   (-Y face)
        (2, 3, 7, 6), // back    (+Y face)
    ];
    let mut indices: Vec<u32> = Vec::new();
    for (a, b, c, d) in &faces {
        // Both windings for double-sided
        indices.extend_from_slice(&[*a as u32, *b as u32, *c as u32]);
        indices.extend_from_slice(&[*a as u32, *c as u32, *b as u32]);
        indices.extend_from_slice(&[*a as u32, *c as u32, *d as u32]);
        indices.extend_from_slice(&[*a as u32, *d as u32, *c as u32]);
    }

    let mut mesh = Mesh::new(
        bevy::render::mesh::PrimitiveTopology::TriangleList,
        RenderAssetUsages::RENDER_WORLD,
    );
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions);
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normals);
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, uvs);
    mesh.insert_indices(bevy::render::mesh::Indices::U32(indices));
    mesh
}
