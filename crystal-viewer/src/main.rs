//! Crystal Viewer — Interactive 3D crystal structure visualization

mod crystal; mod resources; mod picking; mod ui; mod pes;
mod cube_reader; mod marching_cubes; mod volume_render; mod slice_plane;
mod sphere_section;

use bevy::prelude::*;
use bevy::input::mouse::{MouseMotion, MouseScrollUnit, MouseWheel};
use bevy::render::mesh::{Mesh, Indices, PrimitiveTopology};
use bevy_egui::EguiContexts;
use crystal::{CrystalData, Lattice, PhononModesData};
use picking::{PickingState, click_pick, hover_pick, highlight_atoms};
use ui::{AtomInfo, CrystalMeta, ui_system};
use pes::PesData;
use cube_reader::{CubeData, is_cube, parse_cube};
use marching_cubes::marching_cubes_mesh;
use volume_render::{VolumeMaterial, VolumeParams, volume_proxy_mesh, build_volume_texture, update_volume};
use slice_plane::{generate_slice_texture, slice_plane_mesh};
use std::f32::consts::PI;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let json_path = args.get(1).map(|s| s.as_str()).unwrap_or("");

    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugins(bevy_egui::EguiPlugin)
        .add_plugins(volume_render::VolumePlugin)
        .insert_resource(CrystalPath(json_path.to_string()))
        .init_resource::<PanelRects>()
        .add_systems(Startup, setup)
        .insert_resource(RotateState { angle_deg: 45.0 })
        .add_systems(Update, (ui_system, orbit_camera).chain())
        .add_systems(Update, rotate_camera_keys.after(ui_system))
        .add_systems(Update, (click_pick, hover_pick).chain().after(ui_system))
        .add_systems(Update, highlight_atoms.after(ui_system))
        .add_systems(Update, move_selected_atom.after(highlight_atoms))
        .add_systems(Update, (add_atom_system, delete_atom_system))
        .add_systems(Update, (display_mode_system, sync_atom_radii).chain())
        .add_systems(Update, sync_axes_visibility.after(display_mode_system))
        .add_systems(Update, (toggle_projection, sync_arrow_visibility))
        .add_systems(Update, (toggle_surface_combined, update_color_clip, update_pes_surface).chain().after(ui_system))
        .add_systems(Update, toggle_pes3d_mode.after(ui_system))
        .add_systems(Update, (update_isosurface, update_isosurface_mesh).chain().after(ui_system))
        .add_systems(Update, update_surface_mesh.after(ui_system))
        .add_systems(Update, update_volume.after(ui_system))
        .add_systems(Update, update_slices_inner.after(ui_system))
        .run();
}

#[derive(Resource)] struct CrystalPath(String);
#[derive(Resource, Clone)] struct CrystalStore { data: CrystalData, json_path: String }
#[derive(Resource)] struct LatticeData { vecs: [Vec3; 3], inv: [Vec3; 3] }
#[derive(Resource)] struct ImageOffsets(pub Vec<Vec3>);  // fractional offset for each expanded atom

/// egui panel screen rects — used to block camera zoom/pick when mouse is over UI
#[derive(Resource, Default)]
pub struct PanelRects {
    pub left:   Option<bevy_egui::egui::Rect>,
    pub right:  Option<bevy_egui::egui::Rect>,
    pub bottom: Option<bevy_egui::egui::Rect>,
}

/// State machine for Add Atom UI flow
#[derive(Resource, Default)]
struct AddAtomState {
    show_table: bool,
    selected_element: Option<String>,
    coord_x: String,
    coord_y: String,
    coord_z: String,
    /// (element, frac_x, frac_y, frac_z) — populated by UI, consumed by add_atom_system
    pending: Option<(String, f32, f32, f32)>,
    next_id: u32,
}

#[derive(Resource)]
struct CachedSphere(Handle<Mesh>);
#[derive(Resource)] struct MoveState { step: f32 }
#[derive(Resource)] pub struct RotateState { pub angle_deg: f32 }
#[derive(Component)] pub struct MainCamera;
#[derive(Component)] struct FollowCamera;
#[derive(Component)] struct AtomMarker;
#[derive(Component)] struct BondMarker;
#[derive(Component)] struct CellMarker;
#[derive(Component)] struct DisplacementArrow;

#[derive(Resource, Clone)]
pub struct PhononState {
    pub frequency: f32,
    pub ir_intensity: f32,
    pub mode_charge_norm: f32,
    pub mode_index: u32,
    pub scale_factor: f32,
    pub show_arrows: bool,
}

#[derive(Resource)]
#[derive(Component)]
struct CellAxes;

#[derive(Resource, Clone)]
pub struct PesState {
    pub plane: String,
    pub nx: usize,
    pub ny: usize,
    pub e_min: f64,
    pub e_max: f64,
    pub scan_mode: String,
    pub has_energies: bool,
    pub show_surface: bool,
    pub color_step: i32,
}

/// Convert step counter to color_clip f32 (exponential, symmetric).
const COLOR_STEP_MIN: i32 = 0;
const COLOR_STEP_MAX: i32 = 13;
fn step_to_clip(step: i32) -> f32 {
    (0.8_f32.powi(step)).max(0.05)
}

#[derive(Component)]
struct PesSurface;

#[derive(Clone, PartialEq)]
pub enum VisMode { Isosurface, Volume, Slice, Sphere, Migration }

#[derive(Clone, Copy, PartialEq, Eq, Default)]
pub enum IsoMaterial {
    #[default]
    Opaque,           // alpha = 1.0, no transparency
    SemiTransparent,  // alpha = 0.7, light transparency
    Transparent,      // alpha = 0.3, heavy transparency
}

impl IsoMaterial {
    pub fn alpha(&self) -> f32 {
        match self {
            IsoMaterial::Opaque => 1.0,
            IsoMaterial::SemiTransparent => 0.7,
            IsoMaterial::Transparent => 0.3,
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            IsoMaterial::Opaque => "Opaque",
            IsoMaterial::SemiTransparent => "Semi-transparent",
            IsoMaterial::Transparent => "Transparent",
        }
    }
}

#[derive(Resource, Clone)]
#[allow(dead_code)]
pub struct Pes3dState {
    pub nx: usize, pub ny: usize, pub nz: usize,
    pub e_min: f32, pub e_max: f32,
    pub color_min: f32, pub color_max: f32,  // user-configurable color mapping range
    pub clip_x: [f32; 2], pub clip_y: [f32; 2], pub clip_z: [f32; 2],  // fractional clipping ranges [0,1]
    pub has_energies: bool, pub has_expanded: bool,
    pub show_surface: bool, pub vis_mode: VisMode,
    pub color_clip: f32, pub iso_value: f32,
    pub iso_step: f32,  // user-configurable isosurface step size (eV)
    pub iso_material: IsoMaterial,  // material preset
    pub sphere_center_idx: usize,     // mode 7: cube atom index (usize::MAX = custom frac)
    pub sphere_center_custom: [f32; 3],  // mode 7: custom center (fractional)
    pub sphere_radius: f32,           // mode 7: sphere radius (Å)
    pub mig_e_cap: f32,               // mode 8: max relative energy on the surface (eV)
    pub mig_show_shell: bool,         // mode 8: show cage shells
    pub mig_show_cap: bool,           // mode 8: show window caps
    pub alpha_scale: f32,   // volume render: opacity multiplier
    pub alpha_falloff: f32, // volume render: gaussian window width (energy layer selector)
    pub vol_iso_ref: f32,   // volume render: band-pass window center (eV)
    pub vol_steps: u32,     // volume render: ray-march step count
    pub slice_axis: u8, pub slice_pos: f32,
}
#[derive(Resource)] pub struct CubeResource(pub CubeData);
#[derive(Component)] struct IsoSurface;
#[derive(Component)] pub struct VolumeProxy;
#[derive(Component)] struct SlicePlaneMesh;
#[derive(Component)] struct SphereSurface;  // mode 7 sphere + mode 8 migration mesh

#[derive(Resource)]
struct DisplayMode {
    mode: u8,          // 1=ball-stick, 2=space-filling, 3=wireframe
    show_bonds: bool,
    show_cell: bool,
    show_axes: bool,
}

/// Saved initial camera state for R-key reset
#[derive(Resource)]
struct CameraInit(CameraState);

#[derive(Resource, Clone)]
pub struct CameraState { pub focus: Vec3, pub radius: f32, pub rot: Quat }

#[derive(Default)] struct InputState { rotating: bool }

#[derive(Resource, PartialEq)]
pub enum ProjMode { Perspective, Orthographic }

fn ortho_projection(scale: f32) -> Projection {
    Projection::Orthographic(OrthographicProjection {
        scaling_mode: bevy::render::camera::ScalingMode::AutoMin {
            min_width: scale, min_height: scale,
        },
        near: -1000.0, far: 1000.0,
        ..OrthographicProjection::default_3d()
    })
}

fn toggle_projection(
    keys: Res<ButtonInput<KeyCode>>,
    mut proj_mode: ResMut<ProjMode>,
    mut ortho_scale: ResMut<OrthoScale>,
    mut cam_state: ResMut<CameraState>,
    mut camera_q: Query<&mut Projection, With<MainCamera>>,
    mut contexts: EguiContexts,
) {
    if contexts.ctx_mut().wants_keyboard_input() { return; }
    if keys.just_pressed(KeyCode::KeyP) {
        *proj_mode = match *proj_mode {
            ProjMode::Perspective => {
                // Sync: derive ortho_scale from perspective radius
                ortho_scale.0 = cam_state.radius * 1.09;
                ProjMode::Orthographic
            }
            ProjMode::Orthographic => {
                // Sync: derive perspective radius from ortho_scale
                cam_state.radius = ortho_scale.0 / 1.09;
                ProjMode::Perspective
            }
        };
        let Ok(mut proj) = camera_q.get_single_mut() else { return };
        *proj = match *proj_mode {
            ProjMode::Perspective => Projection::Perspective(PerspectiveProjection {
                fov: 1.0, ..default()
            }),
            ProjMode::Orthographic => ortho_projection(ortho_scale.0),
        };
    }
}

#[derive(Resource)]
pub struct OrthoScale(pub f32);  // vertical world units visible in ortho mode

fn orbit_camera(
    mut camera_q: Query<&mut Transform, (With<MainCamera>, Without<FollowCamera>)>,
    mut proj_q: Query<&mut Projection, With<MainCamera>>,
    mut light_q: Query<&mut Transform, With<FollowCamera>>,
    mut cam_state: ResMut<CameraState>,
    cam_init: Res<CameraInit>,
    proj_mode: Res<ProjMode>,
    mut ortho_scale: ResMut<OrthoScale>,
    mut contexts: EguiContexts,
    mut input: Local<InputState>,
    mouse_btn: Res<ButtonInput<MouseButton>>,
    mut mouse_motion: EventReader<MouseMotion>,
    mut mouse_wheel: EventReader<MouseWheel>,
    keys: Res<ButtonInput<KeyCode>>,
    panel_rects: Res<PanelRects>,
) {
    let Ok(mut cam) = camera_q.get_single_mut() else { return };

    // Right mouse → rotate using local axes (no gimbal lock)
    if mouse_btn.just_pressed(MouseButton::Right) { input.rotating = true; }
    if mouse_btn.just_released(MouseButton::Right) { input.rotating = false; }

    for motion in mouse_motion.read() {
        let d = motion.delta;
        if input.rotating {
            let right = cam_state.rot * Vec3::X;
            let up    = cam_state.rot * Vec3::Y;
            let delta = Quat::from_axis_angle(right, -d.y * 0.005)
                      * Quat::from_axis_angle(up,    -d.x * 0.005);
            cam_state.rot = (delta * cam_state.rot).normalize();
        }
    }

    let ctx = contexts.ctx_mut();
    let over_panel = ctx.input(|i| i.pointer.interact_pos()).map_or(false, |pos| {
        panel_rects.left.map_or(false, |r| r.contains(pos))
            || panel_rects.right.map_or(false, |r| r.contains(pos))
            || panel_rects.bottom.map_or(false, |r| r.contains(pos))
    });
    for ev in mouse_wheel.read() {
        if over_panel { continue; }
        let dy = match ev.unit {
            MouseScrollUnit::Line => ev.y * 0.1,
            MouseScrollUnit::Pixel => ev.y * 0.001,
        };
        if *proj_mode == ProjMode::Orthographic {
            ortho_scale.0 = (ortho_scale.0 - dy * ortho_scale.0 * 0.1).clamp(1.0, 200.0);
            cam_state.radius = ortho_scale.0 / 1.09;  // keep in sync
            // Rebuild ortho projection with new scale
            if let Ok(mut proj) = proj_q.get_single_mut() {
                *proj = ortho_projection(ortho_scale.0);
            }
        } else {
            cam_state.radius = (cam_state.radius - dy * cam_state.radius * 0.1).clamp(1.0, 100.0);
            ortho_scale.0 = cam_state.radius * 1.09;  // keep in sync
        }
    }

    if keys.just_pressed(KeyCode::KeyR) {
        *cam_state = cam_init.0.clone();
        ortho_scale.0 = cam_init.0.radius * 1.09;
        if let Ok(mut proj) = proj_q.get_single_mut() {
            if *proj_mode == ProjMode::Orthographic {
                *proj = ortho_projection(ortho_scale.0);
            }
        }
    }

    // Camera position: from focus, step back radius units along forward (-Z) direction
    let forward = cam_state.rot * Vec3::NEG_Z;
    let pos = cam_state.focus - forward * cam_state.radius;
    cam.translation = pos;
    cam.rotation = cam_state.rot;

    for mut lt in light_q.iter_mut() {
        let cam_right = cam.rotation * Vec3::X;
        let cam_up = cam.rotation * Vec3::Y;
        let lp = pos + cam_right * cam_state.radius * 0.8 + cam_up * cam_state.radius * 0.6;
        lt.translation = lp;
        *lt = lt.looking_at(cam_state.focus, Vec3::Y);
    }
}

// ── Arrow key camera rotation ──

fn rotate_camera_keys(
    keys: Res<ButtonInput<KeyCode>>,
    mut cam_state: ResMut<CameraState>,
    rotate_state: Res<RotateState>,
    mut contexts: EguiContexts,
) {
    if contexts.ctx_mut().wants_keyboard_input() { return; }

    let angle_rad = rotate_state.angle_deg.clamp(1.0, 90.0).to_radians();
    let mut delta = Quat::IDENTITY;

    // ← → : yaw around camera's up axis (not world Y)
    let up = cam_state.rot * Vec3::Y;
    if keys.just_pressed(KeyCode::ArrowRight) {
        delta = Quat::from_axis_angle(up, -angle_rad);
    }
    if keys.just_pressed(KeyCode::ArrowLeft) {
        delta = Quat::from_axis_angle(up, angle_rad);
    }

    // ↑ ↓ : pitch around camera's right axis
    if keys.just_pressed(KeyCode::ArrowUp) {
        let right = cam_state.rot * Vec3::X;
        delta = Quat::from_axis_angle(right, angle_rad);
    }
    if keys.just_pressed(KeyCode::ArrowDown) {
        let right = cam_state.rot * Vec3::X;
        delta = Quat::from_axis_angle(right, -angle_rad);
    }

    if delta != Quat::IDENTITY {
        cam_state.rot = (delta * cam_state.rot).normalize();
    }
}

// ── Atom movement ──

fn move_selected_atom(
    keys: Res<ButtonInput<KeyCode>>,
    mut picking: ResMut<PickingState>,
    mut atoms: Query<&mut Transform, With<AtomMarker>>,
    mut contexts: EguiContexts,
    move_state: Res<MoveState>,
    mut crystal: ResMut<CrystalStore>,
    lattice: Res<LatticeData>,
    offsets: Res<ImageOffsets>,
) {
    if contexts.ctx_mut().wants_keyboard_input() { return; }
    if picking.selected < 0 { return; }
    let i = picking.selected as usize;
    if i >= picking.parent_indices.len() { return; }

    let step = move_state.step.clamp(0.01, 10.0);
    let mut dx = 0.0_f32;
    let mut dy = 0.0_f32;
    let mut dz = 0.0_f32;

    if keys.just_pressed(KeyCode::KeyH) { dx = step; }
    if keys.just_pressed(KeyCode::KeyK) { dx = -step; }
    if keys.just_pressed(KeyCode::KeyU) { dy = step; }
    if keys.just_pressed(KeyCode::KeyM) { dy = -step; }
    if keys.just_pressed(KeyCode::KeyI) { dz = step; }
    if keys.just_pressed(KeyCode::KeyN) { dz = -step; }

    if dx != 0.0 || dy != 0.0 || dz != 0.0 {
        let parent = picking.parent_indices[i];
        let cart_delta = Vec3::new(dx, dy, dz);
        let frac_delta = Lattice::apply_inverse(&lattice.inv, cart_delta);

        if parent < crystal.data.atoms.len() {
            // Get current parent coords in fractional space
            let (mut pf_x, mut pf_y, mut pf_z) = if crystal.data.positions_fractional {
                (crystal.data.atoms[parent].x, crystal.data.atoms[parent].y, crystal.data.atoms[parent].z)
            } else {
                let cart = Vec3::new(crystal.data.atoms[parent].x, crystal.data.atoms[parent].y, crystal.data.atoms[parent].z);
                let frac = Lattice::apply_inverse(&lattice.inv, cart);
                (frac.x, frac.y, frac.z)
            };
            // Apply fractional delta
            pf_x += frac_delta.x;
            pf_y += frac_delta.y;
            pf_z += frac_delta.z;
            // Write back in original storage format
            if crystal.data.positions_fractional {
                crystal.data.atoms[parent].x = pf_x;
                crystal.data.atoms[parent].y = pf_y;
                crystal.data.atoms[parent].z = pf_z;
            } else {
                crystal.data.atoms[parent].x = lattice.vecs[0].x * pf_x + lattice.vecs[1].x * pf_y + lattice.vecs[2].x * pf_z;
                crystal.data.atoms[parent].y = lattice.vecs[0].y * pf_x + lattice.vecs[1].y * pf_y + lattice.vecs[2].y * pf_z;
                crystal.data.atoms[parent].z = lattice.vecs[0].z * pf_x + lattice.vecs[1].z * pf_y + lattice.vecs[2].z * pf_z;
            }
            crystal.data.modified = true;

            // Regenerate all images in fractional space, convert to Cartesian
            let parent_frac = Vec3::new(pf_x, pf_y, pf_z);
            let vecs = &lattice.vecs;
            let siblings: Vec<usize> = picking.parent_indices.iter().enumerate()
                .filter(|(_, &p)| p == parent)
                .map(|(j, _)| j)
                .collect();
            for &j in &siblings {
                let img_frac = parent_frac + offsets.0[j];
                let cart = img_frac.x * vecs[0] + img_frac.y * vecs[1] + img_frac.z * vecs[2];
                picking.atom_positions[j] = cart;
                let entity = picking.atom_entities[j];
                if let Ok(mut transform) = atoms.get_mut(entity) {
                    transform.translation = cart;
                }
            }
        }

        picking.modified = true;

        // Auto-save modified positions to JSON
        if let Err(e) = crystal.data.write_to_file(&crystal.json_path) {
            eprintln!("  Failed to save modified positions: {}", e);
        }
    }
}

fn add_atom_system(
    mut add_state: ResMut<AddAtomState>,
    mut crystal: ResMut<CrystalStore>,
    mut picking: ResMut<PickingState>,
    mut offsets: ResMut<ImageOffsets>,
    mut atom_info: ResMut<ui::AtomInfo>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    sphere: Res<CachedSphere>,
    mut commands: Commands,
) {
    let Some((el, fx, fy, fz)) = add_state.pending.take() else { return };

    let label = format!("new_{}", add_state.next_id);
    add_state.next_id += 1;
    atom_info.labels.push(label.clone());

    // Append to CrystalData
    crystal.data.atoms.push(crystal::AtomData {
        element: el.clone(), x: fx, y: fy, z: fz, label,
    });
    crystal.data.modified = true;

    let parent = crystal.data.atoms.len() - 1;
    let frac = Vec3::new(fx, fy, fz);
    let (positions, image_offsets) = crystal.data.expand_single_atom(frac);

    // Spawn entities for each image
    let color = resources::element_color(&el);
    let mat_handle = materials.add(StandardMaterial {
        base_color: color, metallic: 0.2, perceptual_roughness: 0.5, ..default()
    });
    let mut handles = Vec::new();
    let mut entities = Vec::new();
    for pos in &positions {
        let entity = commands.spawn((
            Mesh3d(sphere.0.clone()),
            MeshMaterial3d(mat_handle.clone()),
            Transform::from_translation(*pos),
            AtomMarker,
        )).id();
        handles.push(mat_handle.clone());
        entities.push(entity);
    }

    // Update PickingState
    picking.add_images(positions, handles, entities, parent);

    // Update ImageOffsets
    for off in image_offsets {
        offsets.0.push(off);
    }

    // Update AtomInfo
    atom_info.elements.push(el.clone());
    atom_info.radii.push(resources::covalent_radius(&el));

    // Auto-save
    let _ = crystal.data.write_to_file(&crystal.json_path);
}

fn sync_atom_radii(
    atom_info: Res<ui::AtomInfo>,
    picking: Res<PickingState>,
    mut atoms: Query<&mut Transform, (With<AtomMarker>, Without<BondMarker>)>,
    display: Res<DisplayMode>,
    mut initialized: Local<bool>,
) {
    // Skip first-frame change detection (resource insertion triggers is_changed)
    if !*initialized { *initialized = true; return; }
    if display.mode != 1 || !atom_info.is_changed() { return; }
    for (i, mut t) in atoms.iter_mut().enumerate() {
        if i < picking.parent_indices.len() {
            let p = picking.parent_indices[i];
            if p < atom_info.elements.len() && p < atom_info.radii.len() {
                let default_r = resources::covalent_radius(&atom_info.elements[p]);
                t.scale = Vec3::splat(atom_info.radii[p] / default_r);
            }
        }
    }
}

fn delete_atom_system(
    keys: Res<ButtonInput<KeyCode>>,
    mut picking: ResMut<PickingState>,
    mut crystal: ResMut<CrystalStore>,
    mut offsets: ResMut<ImageOffsets>,
    mut atom_info: ResMut<ui::AtomInfo>,
    mut commands: Commands,
    mut contexts: EguiContexts,
) {
    if contexts.ctx_mut().wants_keyboard_input() { return; }
    if !keys.just_pressed(KeyCode::KeyD) || picking.selected < 0 { return; }
    let i = picking.selected as usize;
    if i >= picking.parent_indices.len() { return; }
    let parent = picking.parent_indices[i];
    if parent >= crystal.data.atoms.len() { return; }

    // Collect entities and image-offset indices before mutation
    let to_despawn: Vec<Entity> = (0..picking.parent_indices.len())
        .filter(|&j| picking.parent_indices[j] == parent)
        .map(|j| picking.atom_entities[j])
        .collect();
    let off_indices: Vec<usize> = (0..picking.parent_indices.len())
        .filter(|&j| picking.parent_indices[j] == parent)
        .collect();

    // Despawn entities
    for entity in &to_despawn {
        commands.entity(*entity).despawn();
    }

    // Remove from ImageOffsets (descending)
    let mut sorted = off_indices.clone();
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    for j in &sorted {
        offsets.0.remove(*j);
    }

    // Remove from PickingState (handles renumbering)
    picking.remove_images(parent);

    // Remove from canonical data
    crystal.data.atoms.remove(parent);
    if parent < atom_info.elements.len() {
        atom_info.elements.remove(parent);
        atom_info.labels.remove(parent);
        atom_info.radii.remove(parent);
    }

    picking.selected = -1;
    picking.modified = true;
    crystal.data.modified = true;
    let _ = crystal.data.write_to_file(&crystal.json_path);
}

// ── Geometry ──

fn uv_sphere(radius: f32, sec: u32, stk: u32) -> Mesh {
    let mut pos = Vec::new(); let mut nrm = Vec::new(); let mut idx = Vec::new();
    for i in 0..=stk {
        let phi = PI * i as f32 / stk as f32;
        let y = -radius * phi.cos(); let rr = radius * phi.sin();
        for j in 0..=sec {
            let t = 2.0 * PI * j as f32 / sec as f32;
            let x = rr * t.cos(); let z = rr * t.sin();
            pos.push([x, y, z]);
            let len = (x*x + y*y + z*z).sqrt().max(0.0001);
            nrm.push([x/len, y/len, z/len]);
        }
    }
    for i in 0..stk { for j in 0..sec {
        let f = i * (sec + 1) + j; let s = f + sec + 1;
        idx.extend_from_slice(&[f, s, f+1, s, s+1, f+1]);
    }}
    let mut m = Mesh::new(PrimitiveTopology::TriangleList, Default::default());
    m.insert_attribute(Mesh::ATTRIBUTE_POSITION, pos);
    m.insert_attribute(Mesh::ATTRIBUTE_NORMAL, nrm);
    m.insert_indices(Indices::U32(idx));
    m
}

// ── Scene ──

/// Percentile of the valid (finite) field values — used for the default
/// color-mapping range so that outlier points don't wash out the colormap.
fn energy_percentile(field: &[f32], p: f64) -> f32 {
    let mut vals: Vec<f32> = field.iter().filter(|v| v.is_finite()).copied().collect();
    if vals.is_empty() { return 0.0; }
    vals.sort_by(|a, b| a.total_cmp(b));
    let i = p * (vals.len() - 1) as f64;
    let lo = i.floor() as usize;
    let hi = i.ceil() as usize;
    let frac = (i - lo as f64) as f32;
    vals[lo] + (vals[hi] - vals[lo]) * frac
}

fn setup(
    mut commands: Commands, mut meshes: ResMut<Assets<Mesh>>,
    mut images: ResMut<Assets<Image>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut vol_materials: ResMut<Assets<VolumeMaterial>>, crystal_path: Res<CrystalPath>,
) {
    // ── Cube file detection ──
    let cube_opt: Option<CubeData> = if is_cube(&crystal_path.0) {
        match parse_cube(&crystal_path.0) {
            Ok(cd) => { println!("Loaded cube: {}x{}x{} grid", cd.nx, cd.ny, cd.nz); Some(cd) }
            Err(e) => { eprintln!("Cube parse error: {}", e); None }
        }
    } else { None };

    // ── 2D PES from cube (nz==1) → convert to PesData ──
    let pes_from_cube: Option<PesData> = cube_opt.as_ref()
        .and_then(|c| if c.is_pes_2d() { PesData::from_cube(c).ok() } else { None });

    let data = if let Some(ref cube) = cube_opt {
        if cube.is_pes_2d() {
            // 2D PES cube: exclude mobile atom from structure rendering
            let mobile_idx = cube.pes_meta.as_ref().and_then(|m| m.mobile_idx).unwrap_or(0);
            crystal_data_from_pes_cube(cube, mobile_idx)
        } else {
            crystal_data_from_cube(cube)
        }
    } else if crystal_path.0.is_empty() {
        default_cu_fcc()
    } else {
        CrystalData::from_json(&crystal_path.0).unwrap_or_else(|e| {
            eprintln!("Failed to load {}: {}. Using Cu FCC.", crystal_path.0, e);
            default_cu_fcc()
        })
    };

    let (positions, parent_indices, image_offsets) = data.expand_to_cell();
    let center = data.center();
    let corners = data.cell_corners();
    let edges = CrystalData::cell_edges();
    let n = positions.len();

    let lattice = LatticeData {
        vecs: data.lattice.to_vectors(),
        inv: data.lattice.inverse_vectors(),
    };
    spawn_cell_axes(&mut commands, &mut meshes, &mut materials, &lattice);
    commands.insert_resource(lattice);

    let sphere = meshes.add(uv_sphere(0.5, 32, 32));
    commands.insert_resource(CachedSphere(sphere.clone()));
    let mut handles = Vec::with_capacity(n);
    let mut entities = Vec::with_capacity(n);

    for (i, pos) in positions.iter().enumerate() {
        let parent = parent_indices[i];
        let el = &data.atoms[parent].element;
        let color = resources::element_color(el);
        let mat_handle = materials.add(StandardMaterial {
            base_color: color, metallic: 0.2, perceptual_roughness: 0.5, ..default()
        });
        let entity = commands.spawn((
            Mesh3d(sphere.clone()), MeshMaterial3d(mat_handle.clone()),
            Transform::from_translation(*pos), AtomMarker,
        )).id();
        handles.push(mat_handle);
        entities.push(entity);
    }

    // Bonds
    let bond_mat = materials.add(StandardMaterial {
        base_color: Color::srgba(0.65, 0.65, 0.65, 0.5),
        alpha_mode: AlphaMode::Blend, ..default()
    });
    for i in 0..n { for j in (i+1)..n {
        // Skip bonds between images of the same parent atom
        if parent_indices[i] == parent_indices[j] { continue; }
        let pi = parent_indices[i];
        let pj = parent_indices[j];
        if resources::has_bond(&data.atoms[pi].element, &data.atoms[pj].element,
                                positions[i].distance(positions[j]), 1.2) {
            spawn_bond(&mut commands, &mut meshes, &bond_mat, positions[i], positions[j], BondMarker);
        }
    }}

    // Cell edges
    let edge_mat = materials.add(StandardMaterial {
        base_color: Color::srgba(0.4, 0.4, 0.4, 0.35),
        alpha_mode: AlphaMode::Blend, ..default()
    });
    for (a, b) in edges {
        let s = corners[a]; let e = corners[b];
        let mid = (s + e) * 0.5; let dir = (e - s).normalize_or_zero(); let len = s.distance(e);
        commands.spawn((
            Mesh3d(meshes.add(Cuboid::new(0.04, 0.04, len))),
            MeshMaterial3d(edge_mat.clone()),
            Transform::from_translation(mid).with_rotation(Quat::from_rotation_arc(Vec3::Z, dir)),
            CellMarker,
        ));
    }

    commands.insert_resource(PickingState::new(positions.clone(), handles, entities, parent_indices.clone()));
    commands.insert_resource(ImageOffsets(image_offsets));

    // Store crystal data + path for auto-save on modification
    commands.insert_resource(CrystalStore {
        data: data.clone(),
        json_path: crystal_path.0.clone(),
    });
    commands.insert_resource(MoveState { step: 0.5 });
    commands.insert_resource(AddAtomState::default());

    // Atom metadata for UI
    let elements: Vec<String> = data.atoms.iter().map(|a| a.element.clone()).collect();
    let labels: Vec<String> = data.atoms.iter().enumerate()
        .map(|(_i, a)| if a.label.is_empty() { a.element.clone() } else { a.label.clone() })
        .collect();
    let radii: Vec<f32> = elements.iter().map(|el| resources::covalent_radius(el)).collect();
    commands.insert_resource(AtomInfo::new(elements, labels, radii));
    commands.insert_resource(DisplayMode { mode: 1, show_bonds: false, show_cell: true, show_axes: false });

    // Crystal metadata for UI (strip .json extension from filename)
    let fname = if crystal_path.0.is_empty() {
        "Cu FCC (demo)".to_string()
    } else {
        let raw = std::path::Path::new(&crystal_path.0)
            .file_name().map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| crystal_path.0.clone());
        raw.strip_suffix(".json").map(|s| s.to_string()).unwrap_or(raw)
    };
    commands.insert_resource(CrystalMeta {
        filename: fname,
        a: data.lattice.a, b: data.lattice.b, c: data.lattice.c,
        alpha: data.lattice.alpha, beta: data.lattice.beta, gamma: data.lattice.gamma,
    });

    // Compute default ortho scale from cell bounding box
    let cell_corners = data.cell_corners();
    let mut cmin = Vec3::splat(f32::MAX);
    let mut cmax = Vec3::splat(f32::MIN);
    for c in &cell_corners {
        cmin = cmin.min(*c);
        cmax = cmax.max(*c);
    }
    let cell_diag = (cmax - cmin).length();
    let ortho_scale = (cell_diag * 1.3).max(5.0);

    // Init camera (orthographic by default, looking down Z at XY plane)
    let radius = ortho_scale;
    let initial_rot = Quat::IDENTITY;  // camera looks along -Z (top-down XY view)
    commands.insert_resource(CameraState { focus: center, radius, rot: initial_rot });
    commands.insert_resource(CameraInit(CameraState { focus: center, radius, rot: initial_rot }));
    commands.insert_resource(ProjMode::Orthographic);
    commands.insert_resource(OrthoScale(ortho_scale));

    // ── Phonon mode data: store state and spawn displacement arrows ──
    let phonon_state = if let Some(ref pm) = data.phonon_modes {
        println!("Phonon mode {}: freq={:.2} cm⁻¹, IR={:.4}, |p_m|={:.4}",
            pm.mode_index, pm.frequency, pm.ir_intensity, pm.mode_charge_norm);
        let state = PhononState {
            frequency: pm.frequency,
            ir_intensity: pm.ir_intensity,
            mode_charge_norm: pm.mode_charge_norm,
            mode_index: pm.mode_index,
            scale_factor: 1.0,
            show_arrows: true,
        };
        // Spawn displacement arrows
        spawn_phonon_arrows(&mut commands, &mut meshes, &mut materials, &data, pm, &positions, &parent_indices);
        Some(state)
    } else {
        None
    };
    // Always insert the resource (None if no phonon data)
    if let Some(ps) = phonon_state {
        commands.insert_resource(ps);
    }

    commands.spawn((
        DirectionalLight { illuminance: 8000.0, shadows_enabled: false, ..default() },
        Transform::default(), FollowCamera,
    ));
    commands.spawn((
        DirectionalLight { illuminance: 2000.0, ..default() },
        Transform::from_xyz(0.0, -5.0, 0.0).looking_at(center, Vec3::Y),
    ));

    commands.spawn((
        Camera3d::default(),
        ortho_projection(ortho_scale),
        Transform::default(),
        MainCamera,
    ));

    // ── PES 2D surface mesh (from cube) ──
    if let Some(pes) = pes_from_cube.as_ref() {
        let e_range = pes.energy_range();
        if let Some((mesh, tex_image)) = pes.generate_surface(1.0) {
            let tex_handle = materials.add(StandardMaterial {
                base_color_texture: Some(images.add(tex_image)),
                alpha_mode: AlphaMode::Blend,
                unlit: true,
                ..default()
            });
            commands.spawn((
                Mesh3d(meshes.add(mesh)),
                MeshMaterial3d(tex_handle),
                PesSurface,
            ));
            println!("PES surface: {}×{} grid, {} energies, plane={}",
                pes.nx, pes.ny,
                if pes.has_energies { "with" } else { "no" },
                pes.plane);
        }

        commands.insert_resource(PesState {
            plane: pes.plane.clone(),
            nx: pes.nx,
            ny: pes.ny,
            e_min: e_range.map(|(min, _)| min).unwrap_or(0.0),
            e_max: e_range.map(|(_, max)| max).unwrap_or(0.0),
            scan_mode: pes.scan_mode.clone(),
            has_energies: pes.has_energies,
            show_surface: pes.has_energies,
            color_step: 0,
        });
        commands.insert_resource(PesDataResource((*pes).clone()));
    }

    // ── PES 3D surface mesh (isosurface + slice planes) — 3D cubes only ──
    if let Some(ref cube) = cube_opt {
        if cube.is_pes_2d() { /* handled by 2D surface above */ }
        else {
        let lattice = cube.to_lattice();
        let e_min = cube.field.iter().cloned().fold(f32::MAX, f32::min);
        let e_max = cube.field.iter().cloned().fold(f32::MIN, f32::max);
        // Default color mapping range: 5th–95th percentile of valid data.
        // The raw min/max range is dominated by a few outlier points (e.g. 0–3223 eV
        // with 91% of data in 400–500 eV), which makes the jet colormap render
        // almost everything blue. Percentiles give visible color variation.
        let c_min = energy_percentile(&cube.field, 0.05);
        let c_max = energy_percentile(&cube.field, 0.95);
        // Extract actual fractional range from cube metadata (fallback to [0,1])
        let frac_range = cube.pes_meta.as_ref()
            .and_then(|m| {
                let fx = m.fx_range.unwrap_or([0.0, 1.0]);
                let fy = m.fy_range.unwrap_or([0.0, 1.0]);
                let fz = m.fz_range.unwrap_or([0.0, 1.0]);
                Some([fx, fy, fz])
            })
            .unwrap_or([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]);
        println!("PES 3D: {}x{}x{} grid, E=[{:.4}, {:.4}], {} atoms",
            cube.nx, cube.ny, cube.nz, e_min, e_max, cube.atoms.len());

        // Spawn volume proxy (mode 5) — custom ray-march material
        let vol_mesh = meshes.add(volume_proxy_mesh(&lattice));
        let vol_img = images.add(build_volume_texture(
            &cube.field, cube.nx, cube.ny, cube.nz, c_min, c_max));
        let vol_mat = vol_materials.add(VolumeMaterial {
            volume: vol_img,
            params: VolumeParams::default(),
        });
        commands.spawn((Mesh3d(vol_mesh), MeshMaterial3d(vol_mat), VolumeProxy, Visibility::Hidden));

        // Volume-mode ISO ref: median energy (window center in the dense band).
        // Isosurface mode gets the classic 15% point of the full range.
        let iso_value_vol = energy_percentile(&cube.field, 0.5);
        // Spawn MC isosurface
        let iso_value = e_min + 0.15 * (e_max - e_min);
        if let Some(mc_mesh) = marching_cubes_mesh(
            &cube.field, cube.nx, cube.ny, cube.nz,
            &frac_range, &lattice,
            iso_value, e_min, e_max,
            [0.0, 1.0], [0.0, 1.0], [0.0, 1.0],  // no clipping at startup
        ) {
            let solid_mat = materials.add(StandardMaterial {
                base_color: Color::srgba(0.3, 0.6, 1.0, 0.5),
                alpha_mode: AlphaMode::Blend,
                unlit: true,
                cull_mode: None,
                ..default()
            });
            commands.spawn((
                Mesh3d(meshes.add(mc_mesh)),
                MeshMaterial3d(solid_mat),
                IsoSurface,
            ));
        }

        // Spawn slice planes (hidden until VisMode::Slice)
        for (axis, pos) in [(2u8, 0.5f32), (1, 0.5), (0, 0.5)] {
            let sm = meshes.add(slice_plane_mesh(&lattice, axis, pos));
            let tex = generate_slice_texture(&cube.field, cube.nx, cube.ny, cube.nz, axis, pos, e_min, e_max, 1.0);
            let mat = materials.add(StandardMaterial {
                base_color_texture: Some(images.add(tex)),
                alpha_mode: AlphaMode::Blend, unlit: true, ..default()
            });
            commands.spawn((Mesh3d(sm), MeshMaterial3d(mat), SlicePlaneMesh, Visibility::Hidden));
        }

        commands.insert_resource(CubeResource(cube.clone()));
        let default_iso_step = (e_max - e_min) * 0.02;  // 2% of range as default
        commands.insert_resource(Pes3dState {
            nx: cube.nx, ny: cube.ny, nz: cube.nz,
            e_min, e_max,
            color_min: c_min, color_max: c_max,  // percentile-based default range
            clip_x: [0.0, 1.0], clip_y: [0.0, 1.0], clip_z: [0.0, 1.0],  // default: no clipping
            has_energies: true, has_expanded: false,
            show_surface: true, vis_mode: VisMode::Isosurface,
            color_clip: 1.0, iso_value,
            iso_step: default_iso_step,
            iso_material: IsoMaterial::SemiTransparent,
            // Mode 7: default center = first S atom (cage center if present)
            sphere_center_idx: cube.atoms.iter().position(|a| a.z == 16).unwrap_or(0),
            sphere_center_custom: [0.25, 0.25, 0.25],
            sphere_radius: 1.4,
            // Mode 8: migration window (relative eV), shells + caps on
            mig_e_cap: 3.0, mig_show_shell: true, mig_show_cap: true,
            alpha_scale: 0.8, alpha_falloff: 0.25, vol_steps: 128,
            // Volume ISO ref starts at the data median so the band-pass window
            // sits in the dense energy band by default.
            vol_iso_ref: iso_value_vol,
            slice_axis: 2, slice_pos: 0.5,
        });
        } // end else (3D path)
    }

    if cube_opt.is_some() {
        println!("PES 3D mode: {} atoms. 4:iso 5:volume 6:slice 7:sphere 8:migration S:toggle -/+:iso [ ]:clip", n);
    } else if pes_from_cube.is_some() {
        println!("PES 2D mode: {} atoms. Right-drag: rotate | Scroll: zoom | S: surface toggle.", n);
    } else {
        println!("Loaded {} atoms. Right-drag to rotate, scroll to zoom, click to select.", n);
    }
}

/// Build a cone mesh pointing along +Y (apex at top, base at origin)
fn cone_mesh(base_radius: f32, height: f32, segments: usize) -> Mesh {
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut indices = Vec::new();
    let apex = [0.0, height, 0.0];
    let base_y = 0.0;

    // Base center + ring vertices
    let base_center_idx = positions.len();
    positions.push([0.0, base_y, 0.0]);
    normals.push([0.0, -1.0, 0.0]);

    for s in 0..=segments {
        let angle = (s as f32 / segments as f32) * 2.0 * std::f32::consts::PI;
        let x = base_radius * angle.cos();
        let z = base_radius * angle.sin();
        positions.push([x, base_y, z]);
        // Side normal: outward perpendicular to cone surface
        let slant = base_radius.atan2(height);
        normals.push([angle.cos() * slant.cos(), slant.sin(), angle.sin() * slant.cos()]);
    }

    let apex_idx = positions.len();
    positions.push(apex);
    normals.push([0.0, 1.0, 0.0]);

    // Base disc (fan)
    for s in 0..segments {
        indices.push((base_center_idx + 1 + s) as u32);
        indices.push((base_center_idx + 1 + (s + 1) % segments) as u32);
        indices.push(base_center_idx as u32);
    }

    // Side triangles
    for s in 0..segments {
        indices.push((base_center_idx + 1 + s) as u32);
        indices.push((base_center_idx + 1 + (s + 1) % segments) as u32);
        indices.push(apex_idx as u32);
    }

    let mut mesh = Mesh::new(PrimitiveTopology::TriangleList, bevy::asset::RenderAssetUsages::default());
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions);
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normals);
    mesh.insert_indices(Indices::U32(indices));
    mesh
}

fn spawn_phonon_arrows(
    commands: &mut Commands,
    meshes: &mut ResMut<Assets<Mesh>>,
    materials: &mut ResMut<Assets<StandardMaterial>>,
    data: &CrystalData,
    pm: &PhononModesData,
    positions: &[Vec3],
    parent_indices: &[usize],
) {
    let shaft_r = 0.06;
    let head_r = 0.14;       // cone base radius
    let head_h = 0.35;       // cone height

    // Map each asymmetric atom to its first cell-image position
    let n_asym = data.atoms.len();
    let first_pos: Vec<Vec3> = (0..n_asym)
        .map(|p| parent_indices.iter().position(|&x| x == p)
            .map(|idx| positions[idx]).unwrap_or(Vec3::ZERO))
        .collect();

    let cone = meshes.add(cone_mesh(head_r, head_h, 16));

    for (i, disp) in pm.atom_displacements.iter().enumerate() {
        if i >= n_asym { break; }
        let mag = (disp.dx * disp.dx + disp.dy * disp.dy + disp.dz * disp.dz).sqrt();
        if mag < 1e-8 { continue; }

        let base_pos = first_pos[i];
        let dir = Vec3::new(disp.dx, disp.dy, disp.dz) / mag;
        // Shaft extends from atom to (displacement - cone_height)
        let shaft_len = (mag - head_h).max(0.02);
        let rot = Quat::from_rotation_arc(Vec3::Y, dir);

        // Color: green (low contribution) → red (high contribution)
        let contrib = disp.contribution.clamp(0.0, 1.0);
        let color = if pm.mode_charge_norm > 1e-8 {
            Color::srgb(0.3 + 0.7 * contrib, 0.8 * (1.0 - contrib), 0.2)
        } else {
            Color::srgb(0.8, 0.8, 0.8)
        };

        let mat = materials.add(StandardMaterial {
            base_color: color,
            emissive: LinearRgba::rgb(
                color.to_linear().red * 1.5,
                color.to_linear().green * 1.5,
                color.to_linear().blue * 0.8,
            ),
            unlit: true, ..default()
        });

        if shaft_len > 0.02 {
            let shaft = meshes.add(Cylinder::new(shaft_r, shaft_len));
            let shaft_mid = base_pos + dir * shaft_len * 0.5;
            commands.spawn((
                Mesh3d(shaft),
                MeshMaterial3d(mat.clone()),
                Transform::from_translation(shaft_mid).with_rotation(rot),
                DisplacementArrow,
            ));
        }

        // Cone head: base at end of shaft, pointing outward
        let head_base = base_pos + dir * shaft_len;
        commands.spawn((
            Mesh3d(cone.clone()),
            MeshMaterial3d(mat.clone()),
            Transform::from_translation(head_base).with_rotation(rot),
            DisplacementArrow,
        ));
    }
}

fn spawn_cell_axes(
    commands: &mut Commands,
    meshes: &mut ResMut<Assets<Mesh>>,
    materials: &mut ResMut<Assets<StandardMaterial>>,
    lattice: &LatticeData,
) {
    let shaft_r = 0.04;
    let head_r = 0.12;
    let head_h = 0.3;
    let cone = meshes.add(cone_mesh(head_r, head_h, 16));

    let colors = [
        Color::srgb(0.9, 0.2, 0.2), // X red
        Color::srgb(0.2, 0.9, 0.2), // Y green
        Color::srgb(0.2, 0.4, 1.0), // Z blue
    ];

    for axis in 0..3 {
        let dir = lattice.vecs[axis].normalize();
        let total_len = lattice.vecs[axis].length() * 1.5;
        let shaft_len = (total_len - head_h).max(0.02);
        let rot = Quat::from_rotation_arc(Vec3::Y, dir);

        let mat = materials.add(StandardMaterial {
            base_color: colors[axis],
            unlit: true,
            depth_bias: -10.0,
            ..default()
        });

        let shaft_mid = dir * shaft_len * 0.5;
        if shaft_len > 0.02 {
            let shaft = meshes.add(Cylinder::new(shaft_r, shaft_len));
            commands.spawn((
                Mesh3d(shaft),
                MeshMaterial3d(mat.clone()),
                Transform::from_translation(shaft_mid).with_rotation(rot),
                CellAxes,
                Visibility::Hidden,
            ));
        }

        let head_pos = dir * shaft_len;
        commands.spawn((
            Mesh3d(cone.clone()),
            MeshMaterial3d(mat.clone()),
            Transform::from_translation(head_pos).with_rotation(rot),
            CellAxes,
            Visibility::Hidden,
        ));
    }
}

fn spawn_bond(
    commands: &mut Commands, meshes: &mut ResMut<Assets<Mesh>>,
    material: &Handle<StandardMaterial>, a: Vec3, b: Vec3, marker: BondMarker,
) {
    let mid = (a + b) * 0.5; let dir = (b - a).normalize_or_zero(); let len = a.distance(b);
    commands.spawn((
        Mesh3d(meshes.add(Cylinder::new(0.08, len))),
        MeshMaterial3d(material.clone()),
        Transform::from_translation(mid).with_rotation(Quat::from_rotation_arc(Vec3::Y, dir)),
        marker,
        Visibility::Hidden,  // hidden by default, press B to show
    ));
}

fn display_mode_system(
    keys: Res<ButtonInput<KeyCode>>,
    mut display: ResMut<DisplayMode>,
    mut atoms: Query<&mut Transform, (With<AtomMarker>, Without<BondMarker>)>,
    picking: Res<PickingState>,
    atom_info: Res<ui::AtomInfo>,
    bonds: Query<Entity, With<BondMarker>>,
    cells: Query<Entity, With<CellMarker>>,
    mut commands: Commands,
    mut contexts: EguiContexts,
) {
    if contexts.ctx_mut().wants_keyboard_input() { return; }
    let update_scales = |mode: &mut u8| -> bool {
        let changed = matches!(
            (keys.just_pressed(KeyCode::Digit1), keys.just_pressed(KeyCode::Digit2), keys.just_pressed(KeyCode::Digit3)),
            (true, _, _) | (_, true, _) | (_, _, true)
        );
        if !changed { return false; }
        if keys.just_pressed(KeyCode::Digit1) { *mode = 1; }
        if keys.just_pressed(KeyCode::Digit2) { *mode = 2; }
        if keys.just_pressed(KeyCode::Digit3) { *mode = 3; }
        true
    };

    if update_scales(&mut display.mode) {
        // Collect parent atom scales
        let n = atom_info.elements.len();
        let mut scales: Vec<f32> = vec![1.0; n];
        for i in 0..n {
            scales[i] = match display.mode {
                1 => 1.0,   // ball-stick: uniform viewing size
                2 => resources::ionic_radius(&atom_info.elements[i]) / 0.5,
                3 => resources::atomic_radius(&atom_info.elements[i]) / 0.5,
                _ => 1.0,
            };
        }
        for (i, mut t) in atoms.iter_mut().enumerate() {
            if i < picking.parent_indices.len() {
                let p = picking.parent_indices[i];
                if p < scales.len() {
                    t.scale = Vec3::splat(scales[p]);
                }
            }
        }
        // Show/hide bonds per mode
        display.show_bonds = display.mode != 2;
        let vis = if display.show_bonds { Visibility::Visible } else { Visibility::Hidden };
        for e in bonds.iter() { commands.entity(e).insert(vis); }
    }
    if keys.just_pressed(KeyCode::KeyB) {
        display.show_bonds = !display.show_bonds;
        let vis = if display.show_bonds { Visibility::Visible } else { Visibility::Hidden };
        for e in bonds.iter() { commands.entity(e).insert(vis); }
    }
    if keys.just_pressed(KeyCode::KeyC) {
        display.show_cell = !display.show_cell;
        let vis = if display.show_cell { Visibility::Visible } else { Visibility::Hidden };
        for e in cells.iter() { commands.entity(e).insert(vis); }
    }
    if keys.just_pressed(KeyCode::KeyA) {
        display.show_axes = !display.show_axes;
    }
}

fn default_cu_fcc() -> CrystalData {
    crystal::CrystalData {
        lattice: crystal::Lattice { a: 3.615, b: 3.615, c: 3.615, alpha: 90.0, beta: 90.0, gamma: 90.0 },
        atoms: vec![
            crystal::AtomData { element: "Cu".into(), x: 0.0, y: 0.0, z: 0.0, label: "1".into() },
            crystal::AtomData { element: "Cu".into(), x: 1.8075, y: 1.8075, z: 0.0, label: "2".into() },
            crystal::AtomData { element: "Cu".into(), x: 1.8075, y: 0.0, z: 1.8075, label: "3".into() },
            crystal::AtomData { element: "Cu".into(), x: 0.0, y: 1.8075, z: 1.8075, label: "4".into() },
        ],
        positions_fractional: false,
        modified: false,
        phonon_modes: None,
    }
}

/// Convert CubeData atoms (excluding mobile atom) to CrystalData for 2D PES rendering.
fn crystal_data_from_pes_cube(cube: &CubeData, mobile_idx: usize) -> CrystalData {
    let lattice = cube.to_lattice();
    let inv = lattice.inverse_vectors();
    crystal::CrystalData {
        lattice,
        atoms: cube.atoms.iter().enumerate()
            .filter(|(i, _)| *i != mobile_idx)
            .map(|(_, a)| {
                let cart_x = a.x as f32; let cart_y = a.y as f32; let cart_z = a.z_coord as f32;
                let fx = inv[0][0]*cart_x + inv[0][1]*cart_y + inv[0][2]*cart_z;
                let fy = inv[1][0]*cart_x + inv[1][1]*cart_y + inv[1][2]*cart_z;
                let fz = inv[2][0]*cart_x + inv[2][1]*cart_y + inv[2][2]*cart_z;
                crystal::AtomData {
                    element: atom_z_to_symbol(a.z),
                    x: fx, y: fy, z: fz,
                    label: String::new(),
                }
            }).collect(),
        positions_fractional: true,
        modified: false,
        phonon_modes: None,
    }
}

/// Convert CubeData atoms to CrystalData for shared rendering.
fn crystal_data_from_cube(cube: &CubeData) -> CrystalData {
    let lattice = cube.to_lattice();
    crystal::CrystalData {
        lattice,
        atoms: cube.atoms.iter().map(|a| crystal::AtomData {
            element: atom_z_to_symbol(a.z),
            x: a.x, y: a.y, z: a.z_coord,
            label: String::new(),
        }).collect(),
        positions_fractional: false,
        modified: false,
        phonon_modes: None,
    }
}

/// Atomic number → element symbol (supports up to Xe, 54).
pub(crate) fn atom_z_to_symbol(z: i32) -> String {
    const SYM: &[&str] = &[
        "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe",
    ];
    if z >= 0 && (z as usize) < SYM.len() { SYM[z as usize].to_string() } else { format!("Z{}", z) }
}

/// Toggle surface visibility with S key (handles both 2D and 3D PES).
fn toggle_surface_combined(
    keys: Res<ButtonInput<KeyCode>>,
    mut pes_state: Option<ResMut<PesState>>,
    mut pes3d_state: Option<ResMut<Pes3dState>>,
    mut vis_q: ParamSet<(
        Query<&mut Visibility, With<PesSurface>>,
        Query<&mut Visibility, With<IsoSurface>>,
        Query<&mut Visibility, With<SlicePlaneMesh>>,
        Query<&mut Visibility, With<VolumeProxy>>,
        Query<&mut Visibility, With<SphereSurface>>,
    )>,
    mut contexts: EguiContexts,
) {
    if contexts.ctx_mut().wants_keyboard_input() { return; }
    if !keys.just_pressed(KeyCode::KeyS) { return; }
    // 3D takes priority
    if let Some(ref mut ps) = pes3d_state {
        ps.show_surface = !ps.show_surface;
        let vis = if ps.show_surface { Visibility::Inherited } else { Visibility::Hidden };
        if ps.vis_mode == VisMode::Isosurface {
            for mut v in vis_q.p1().iter_mut() { *v = vis; }
        } else if ps.vis_mode == VisMode::Volume {
            for mut v in vis_q.p3().iter_mut() { *v = vis; }
        } else if ps.vis_mode == VisMode::Slice {
            for mut v in vis_q.p2().iter_mut() { *v = vis; }
        } else if matches!(ps.vis_mode, VisMode::Sphere | VisMode::Migration) {
            for mut v in vis_q.p4().iter_mut() { *v = vis; }
        }
        return;
    }
    if let Some(ref mut ps) = pes_state {
        ps.show_surface = !ps.show_surface;
        let vis = if ps.show_surface { Visibility::Inherited } else { Visibility::Hidden };
        for mut v in vis_q.p0().iter_mut() { *v = vis; }
    }
}

/// Adjust color clip range with -/+ keys for 2D PES.
/// Integer step counter (0..14) → symmetric exponential scale via 0.8^step.
/// just_pressed → instant response; held → throttle-repeat at ~6.7 Hz.
fn update_color_clip(
    keys: Res<ButtonInput<KeyCode>>,
    time: Res<Time>,
    mut pes_state: Option<ResMut<PesState>>,
    mut contexts: EguiContexts,
    mut repeat_timer: Local<f32>,
) {
    if contexts.ctx_mut().wants_keyboard_input() { return; }
    let ps = match pes_state.as_mut() { Some(ps) => ps, None => return };

    if keys.just_pressed(KeyCode::KeyR) {
        ps.color_step = 0;
        *repeat_timer = 0.0;
        return;
    }

    let just_plus  = keys.just_pressed(KeyCode::Equal) || keys.just_pressed(KeyCode::NumpadAdd);
    let just_minus = keys.just_pressed(KeyCode::Minus) || keys.just_pressed(KeyCode::NumpadSubtract);
    let held_plus  = keys.pressed(KeyCode::Equal) || keys.pressed(KeyCode::NumpadAdd);
    let held_minus = keys.pressed(KeyCode::Minus) || keys.pressed(KeyCode::NumpadSubtract);

    let active = if just_plus || just_minus {
        *repeat_timer = 0.0;
        true
    } else if held_plus || held_minus {
        *repeat_timer += time.delta_secs();
        if *repeat_timer >= 0.15 {
            *repeat_timer = 0.0;
            true
        } else {
            false
        }
    } else {
        *repeat_timer = 0.0;
        false
    };

    if !active { return; }

    let new_step = if just_plus || held_plus {
        (ps.color_step + 1).min(COLOR_STEP_MAX)
    } else {
        (ps.color_step - 1).max(COLOR_STEP_MIN)
    };
    if new_step != ps.color_step {
        ps.color_step = new_step;
    }
}

/// Stored 2D PES data for mesh rebuild on color_clip change.
#[derive(Resource)]
struct PesDataResource(pes::PesData);

/// Rebuild 2D PES surface mesh when color_clip changes.
fn update_pes_surface(
    pes_state: Option<Res<PesState>>,
    pes_data: Option<Res<PesDataResource>>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut images: ResMut<Assets<Image>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    surface_q: Query<Entity, With<PesSurface>>,
) {
    let ps = match pes_state.as_ref() { Some(ps) => ps, None => return };
    if !ps.is_changed() { return; }
    let pd = match pes_data.as_ref() { Some(pd) => pd, None => return };

    for e in surface_q.iter() { commands.entity(e).despawn(); }

    if let Some((mesh, tex)) = pd.0.generate_surface(step_to_clip(ps.color_step)) {
        let mat = materials.add(StandardMaterial {
            base_color_texture: Some(images.add(tex)),
            alpha_mode: AlphaMode::Blend, unlit: true, ..default()
        });
        let vis = if ps.show_surface { Visibility::Inherited } else { Visibility::Hidden };
        commands.spawn((Mesh3d(meshes.add(mesh)), MeshMaterial3d(mat), PesSurface, vis));
    }
}

/// Switch 3D PES visualization mode (4=isosurface, 5=volume, 6=slice,
/// 7=fixed-radius sphere, 8=radial-stationary migration surface).
fn toggle_pes3d_mode(
    keys: Res<ButtonInput<KeyCode>>,
    mut pes3d_state: Option<ResMut<Pes3dState>>,
    mut vis_q: ParamSet<(
        Query<&mut Visibility, With<IsoSurface>>,
        Query<&mut Visibility, With<VolumeProxy>>,
        Query<&mut Visibility, With<SlicePlaneMesh>>,
        Query<&mut Visibility, With<SphereSurface>>,
    )>,
    mut contexts: EguiContexts,
) {
    if contexts.ctx_mut().wants_keyboard_input() { return; }
    let ps = match pes3d_state.as_mut() { Some(ps) => ps, None => return };

    let new_mode = if keys.just_pressed(KeyCode::Digit4) {
        Some(VisMode::Isosurface)
    } else if keys.just_pressed(KeyCode::Digit5) {
        Some(VisMode::Volume)
    } else if keys.just_pressed(KeyCode::Digit6) {
        Some(VisMode::Slice)
    } else if keys.just_pressed(KeyCode::Digit7) {
        Some(VisMode::Sphere)
    } else if keys.just_pressed(KeyCode::Digit8) {
        Some(VisMode::Migration)
    } else { None };

    if let Some(mode) = new_mode {
        ps.vis_mode = mode;
        let show = ps.show_surface;
        let vis_on = if show { Visibility::Inherited } else { Visibility::Hidden };
        for mut v in vis_q.p0().iter_mut() { *v = if ps.vis_mode == VisMode::Isosurface { vis_on } else { Visibility::Hidden }; }
        for mut v in vis_q.p1().iter_mut() { *v = if ps.vis_mode == VisMode::Volume { vis_on } else { Visibility::Hidden }; }
        for mut v in vis_q.p2().iter_mut() { *v = if ps.vis_mode == VisMode::Slice { vis_on } else { Visibility::Hidden }; }
        for mut v in vis_q.p3().iter_mut() {
            *v = if matches!(ps.vis_mode, VisMode::Sphere | VisMode::Migration) { vis_on } else { Visibility::Hidden };
        }
    }
}

/// Adjust isosurface isovalue with -/+ keys.
/// Tap → single step; hold → continuous after 400ms initial delay, then 80ms repeat.
fn update_isosurface(
    keys: Res<ButtonInput<KeyCode>>,
    time: Res<Time>,
    mut pes3d_state: Option<ResMut<Pes3dState>>,
    mut contexts: EguiContexts,
    mut hold_time: Local<f32>,
) {
    if contexts.ctx_mut().wants_keyboard_input() { return; }
    let ps = match pes3d_state.as_mut() { Some(ps) => ps, None => return };
    // +/- adjusts iso_value (Isosurface mode) or vol_iso_ref (Volume mode —
    // the energy-layer window center). update_isosurface_mesh has its own
    // Isosurface guard so the MC mesh is not rebuilt in Volume mode.
    if ps.vis_mode == VisMode::Slice { return; }

    let just_plus  = keys.just_pressed(KeyCode::Equal) || keys.just_pressed(KeyCode::NumpadAdd);
    let just_minus = keys.just_pressed(KeyCode::Minus) || keys.just_pressed(KeyCode::NumpadSubtract);
    let held_plus  = keys.pressed(KeyCode::Equal) || keys.pressed(KeyCode::NumpadAdd);
    let held_minus = keys.pressed(KeyCode::Minus) || keys.pressed(KeyCode::NumpadSubtract);

    let delta = ps.iso_step;

    // Instant single-step on tap
    if just_plus {
        if ps.vis_mode == VisMode::Volume {
            ps.vol_iso_ref = (ps.vol_iso_ref + delta).min(ps.e_max);
        } else {
            ps.iso_value = (ps.iso_value + delta).min(ps.e_max);
        }
        *hold_time = 0.0;
        return;
    }
    if just_minus {
        if ps.vis_mode == VisMode::Volume {
            ps.vol_iso_ref = (ps.vol_iso_ref - delta).max(ps.e_min);
        } else {
            ps.iso_value = (ps.iso_value - delta).max(ps.e_min);
        }
        *hold_time = 0.0;
        return;
    }

    // Continuous repeat on hold (400ms initial delay, then every 80ms)
    if !held_plus && !held_minus { *hold_time = 0.0; return; }

    *hold_time += time.delta_secs();
    let initial_delay = 0.4;  // wait 400ms before first auto-repeat
    let repeat_interval = 0.08;  // then every 80ms

    if *hold_time < initial_delay { return; }
    let elapsed = *hold_time - initial_delay;
    let steps = (elapsed / repeat_interval).floor() as i32;
    if steps < 1 { return; }
    *hold_time = initial_delay + elapsed - steps as f32 * repeat_interval;

    let step_delta = delta * steps as f32;
    if held_plus {
        if ps.vis_mode == VisMode::Volume {
            ps.vol_iso_ref = (ps.vol_iso_ref + step_delta).min(ps.e_max);
        } else {
            ps.iso_value = (ps.iso_value + step_delta).min(ps.e_max);
        }
    } else {
        if ps.vis_mode == VisMode::Volume {
            ps.vol_iso_ref = (ps.vol_iso_ref - step_delta).max(ps.e_min);
        } else {
            ps.iso_value = (ps.iso_value - step_delta).max(ps.e_min);
        }
    }
}

/// Parameters of the last-built isosurface — the mesh is rebuilt ONLY when
/// these actually change. (is_changed() cannot be used: egui constructs its
/// DragValue/Slider widgets with &mut references every frame, which marks
/// Pes3dState as changed every frame via DerefMut → per-frame despawn/respawn
/// → visible flickering.)
#[derive(Default, Clone, Copy, PartialEq)]
struct IsoBuildCache {
    iso_value: f32,
    color_min: f32,
    color_max: f32,
    clip_x: [f32; 2],
    clip_y: [f32; 2],
    clip_z: [f32; 2],
    iso_material: IsoMaterial,  // material preset (alpha is baked into the mesh material)
}

/// Rebuild isosurface mesh when isovalue/color/clip actually change.
fn update_isosurface_mesh(
    pes3d_state: Option<Res<Pes3dState>>,
    cube: Option<Res<CubeResource>>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut images: ResMut<Assets<Image>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    iso_q: Query<Entity, With<IsoSurface>>,
    mut cache: Local<IsoBuildCache>,
) {
    let ps = match pes3d_state.as_ref() { Some(ps) => ps, None => return };
    if ps.vis_mode != VisMode::Isosurface { return; }

    let cur = IsoBuildCache {
        iso_value: ps.iso_value,
        color_min: ps.color_min,
        color_max: ps.color_max,
        clip_x: ps.clip_x,
        clip_y: ps.clip_y,
        clip_z: ps.clip_z,
        iso_material: ps.iso_material,
    };
    if *cache == cur { return; }  // nothing actually changed — keep the mesh
    *cache = cur;

    let cube = match cube.as_ref() { Some(c) => c, None => return };
    let lattice = cube.0.to_lattice();
    let frac_range = cube.0.pes_meta.as_ref()
        .and_then(|m| {
            let fx = m.fx_range.unwrap_or([0.0, 1.0]);
            let fy = m.fy_range.unwrap_or([0.0, 1.0]);
            let fz = m.fz_range.unwrap_or([0.0, 1.0]);
            Some([fx, fy, fz])
        })
        .unwrap_or([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]);

    // Despawn old isosurface
    for e in iso_q.iter() { commands.entity(e).despawn(); }

    // Build new mesh
    if let Some(mc_mesh) = marching_cubes_mesh(
        &cube.0.field, cube.0.nx, cube.0.ny, cube.0.nz,
        &frac_range, &lattice,
        ps.iso_value, ps.color_min, ps.color_max,
        ps.clip_x, ps.clip_y, ps.clip_z,
    ) {
        let colormap = PesData::generate_surface_static(ps.color_min, ps.color_max);
        let alpha = ps.iso_material.alpha();
        let tex = materials.add(StandardMaterial {
            base_color_texture: Some(images.add(colormap)),
            base_color: Color::srgba(1.0, 1.0, 1.0, alpha),  // apply alpha from material preset
            alpha_mode: if alpha < 1.0 { AlphaMode::Blend } else { AlphaMode::Opaque },
            unlit: true,
            cull_mode: None,
            ..default()
        });
        let vis = if ps.show_surface { Visibility::Inherited } else { Visibility::Hidden };
        commands.spawn((Mesh3d(meshes.add(mc_mesh)), MeshMaterial3d(tex), IsoSurface, vis));
    }
}

/// Parameters of the last-built shared 7/8 surface mesh (SphereSurface
/// entity) — rebuilt only when these actually change (same flicker rationale
/// as IsoBuildCache). `built_mode` records which mode the current mesh was
/// built for: without it, switching 8→7 sees identical sphere params and
/// skips the rebuild, leaving the migration mesh on screen (and vice versa).
#[derive(Default, Clone, PartialEq)]
struct SurfaceBuildCache {
    built_mode: Option<VisMode>,
    center_frac: [f32; 3],
    radius: f32,
    e_cap: f32,
    show_shell: bool,
    show_cap: bool,
    centers: Vec<[f32; 3]>,
    color_min: f32,
    color_max: f32,
    material: IsoMaterial,
}

/// Rebuild the shared mode-7 (sphere) / mode-8 (migration) surface mesh when
/// its parameters — or the active mode — actually change. Modes 4/5/6 use
/// their own entities and are toggled by visibility only.
fn update_surface_mesh(
    pes3d_state: Option<Res<Pes3dState>>,
    cube: Option<Res<CubeResource>>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    surf_q: Query<Entity, With<SphereSurface>>,
    mut cache: Local<SurfaceBuildCache>,
) {
    let ps = match pes3d_state.as_ref() { Some(ps) => ps, None => return };
    let cube = match cube.as_ref() { Some(c) => c, None => return };
    let lattice = cube.0.to_lattice();

    let centers: Vec<Vec3> = if ps.vis_mode == VisMode::Migration {
        sphere_section::detect_cage_centers(&cube.0)
    } else { Vec::new() };

    let cur = match ps.vis_mode {
        VisMode::Sphere => {
            // Resolve the sphere center to fractional coordinates
            let center_frac = if ps.sphere_center_idx == usize::MAX {
                Vec3::from_array(ps.sphere_center_custom)
            } else {
                let i = ps.sphere_center_idx.min(cube.0.atoms.len().saturating_sub(1));
                let a = &cube.0.atoms[i];
                let inv = lattice.inverse_vectors();
                Lattice::apply_inverse(&inv, Vec3::new(a.x, a.y, a.z_coord))
            };
            SurfaceBuildCache {
                built_mode: Some(VisMode::Sphere),
                center_frac: center_frac.to_array(),
                radius: ps.sphere_radius,
                e_cap: 0.0, show_shell: false, show_cap: false,
                centers: Vec::new(),
                color_min: ps.color_min, color_max: ps.color_max,
                material: ps.iso_material,
            }
        }
        VisMode::Migration => SurfaceBuildCache {
            built_mode: Some(VisMode::Migration),
            center_frac: [0.0; 3],
            radius: 0.0,
            e_cap: ps.mig_e_cap,
            show_shell: ps.mig_show_shell,
            show_cap: ps.mig_show_cap,
            centers: centers.iter().map(|c| c.to_array()).collect(),
            color_min: ps.color_min, color_max: ps.color_max,
            material: ps.iso_material,
        },
        _ => return,   // modes 4/5/6 have their own entities
    };

    if *cache == cur { return; }  // nothing actually changed — keep the mesh

    for e in surf_q.iter() { commands.entity(e).despawn(); }

    let mesh = match ps.vis_mode {
        VisMode::Sphere => sphere_section::sphere_section_mesh(
            &cube.0.field, cube.0.nx, &lattice,
            Vec3::from_array(cur.center_frac), ps.sphere_radius,
            ps.color_min, ps.color_max, ps.iso_material),
        VisMode::Migration => match sphere_section::migration_surface_mesh(
            &cube.0.field, cube.0.nx, &lattice,
            &centers, ps.mig_e_cap, ps.mig_show_shell, ps.mig_show_cap,
            ps.color_min, ps.color_max, ps.iso_material) {
            Some(m) => m, None => return,
        },
        _ => return,
    };
    let alpha = ps.iso_material.alpha();
    let mat = materials.add(StandardMaterial {
        base_color: Color::WHITE,   // per-vertex colors replace base_color
        alpha_mode: if alpha < 1.0 { AlphaMode::Blend } else { AlphaMode::Opaque },
        unlit: true,
        cull_mode: None,
        ..default()
    });
    let vis = if ps.show_surface { Visibility::Inherited } else { Visibility::Hidden };
    commands.spawn((Mesh3d(meshes.add(mesh)), MeshMaterial3d(mat), SphereSurface, vis));
    *cache = cur;   // only after a successful build — a failed build retries
}

/// Parameters of the last-built slice planes — rebuild only when these
/// actually change (same flicker rationale as IsoBuildCache).
#[derive(Default, Clone, Copy, PartialEq)]
struct SliceBuildCache {
    slice_axis: u8,
    slice_pos: f32,
    color_min: f32,
    color_max: f32,
    color_clip: f32,
}

/// Update slice plane textures when axis/position/color_clip changes.
fn update_slices_inner(
    pes3d_state: Option<Res<Pes3dState>>,
    cube: Option<Res<CubeResource>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut images: ResMut<Assets<Image>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    slice_q: Query<Entity, With<SlicePlaneMesh>>,
    mut commands: Commands,
    mut cache: Local<SliceBuildCache>,
) {
    let ps = match pes3d_state.as_ref() { Some(ps) => ps, None => return };
    let cur = SliceBuildCache {
        slice_axis: ps.slice_axis,
        slice_pos: ps.slice_pos,
        color_min: ps.color_min,
        color_max: ps.color_max,
        color_clip: ps.color_clip,
    };
    if *cache == cur { return; }  // nothing actually changed — keep planes
    *cache = cur;
    let cube = match cube.as_ref() { Some(c) => c, None => return };

    // Regenerate all slice textures (despawn old, spawn new)
    for entity in slice_q.iter() { commands.entity(entity).despawn(); }
    let lattice = cube.0.to_lattice();
    for axis in 0u8..3 {
        let pos = if axis == ps.slice_axis { ps.slice_pos } else { 0.5 };
        let sm = meshes.add(slice_plane_mesh(&lattice, axis, pos));
        let tex = generate_slice_texture(
            &cube.0.field, cube.0.nx, cube.0.ny, cube.0.nz,
            axis, pos, ps.color_min, ps.color_max, ps.color_clip,
        );
        let mat = materials.add(StandardMaterial {
            base_color_texture: Some(images.add(tex)),
            alpha_mode: AlphaMode::Blend, unlit: true, ..default()
        });
        let vis = if ps.show_surface && ps.vis_mode == VisMode::Slice { Visibility::Inherited } else { Visibility::Hidden };
        commands.spawn((Mesh3d(sm), MeshMaterial3d(mat), SlicePlaneMesh, vis));
    }
}

/// Toggle arrow visibility based on PhononState
fn sync_arrow_visibility(
    phonon_state: Option<Res<PhononState>>,
    mut arrow_q: Query<&mut Visibility, With<DisplacementArrow>>,
) {
    let Some(state) = phonon_state else { return };
    if !state.is_changed() { return; }
    let vis = if state.show_arrows { Visibility::Inherited } else { Visibility::Hidden };
    for mut v in arrow_q.iter_mut() {
        *v = vis;
    }
}

/// Toggle cell axes visibility based on DisplayMode.show_axes
fn sync_axes_visibility(
    display: Res<DisplayMode>,
    mut axes_q: Query<&mut Visibility, With<CellAxes>>,
) {
    if !display.is_changed() { return; }
    let vis = if display.show_axes { Visibility::Inherited } else { Visibility::Hidden };
    for mut v in axes_q.iter_mut() {
        *v = vis;
    }
}
