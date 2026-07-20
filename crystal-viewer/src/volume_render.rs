//! Volume rendering via custom WGSL ray-march shader on a proxy cube.
//!
//! Uploads the 3D scalar field as a 3D texture (R32Float), then renders
//! a bounding-box cube with a custom shader that ray-marches through it.

use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};
use bevy::render::render_asset::RenderAssetUsages;
use crate::crystal::Lattice;

/// Upload the 3D energy field as a GPU 3D texture.
pub fn upload_volume_texture(
    field: &[f32], nx: usize, ny: usize, nz: usize,
) -> Image {
    let size = Extent3d {
        width: nx as u32,
        height: ny as u32,
        depth_or_array_layers: nz as u32,
    };
    let bytes: Vec<u8> = field.iter()
        .flat_map(|v| v.to_le_bytes().to_vec())
        .collect();
    Image::new(
        size,
        TextureDimension::D3,
        bytes,
        TextureFormat::R32Float,
        RenderAssetUsages::RENDER_WORLD,
    )
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
