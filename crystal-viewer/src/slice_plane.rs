//! Interactive orthogonal slicing planes through a 3D scalar field.
//!
//! Three planes (XY, XZ, YZ) positioned at a user-controlled fractional coordinate,
//! each with dynamically-updated jet-colormap textures.

use bevy::prelude::*;
use bevy::render::mesh::{Indices, PrimitiveTopology};
use bevy::render::render_asset::RenderAssetUsages;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};
use crate::crystal::Lattice;
use super::pes::jet_rgb;

/// Generate a 2D RGBA slice texture from a 3D field at the given fractional position.
///
/// `axis`: 0=YZ (slice at fx), 1=XZ (slice at fy), 2=XY (slice at fz).
/// `position`: fractional coordinate along the slice axis (0.0..1.0).
pub fn generate_slice_texture(
    field: &[f32],
    nx: usize, ny: usize, nz: usize,
    axis: u8,
    position: f32,
    e_min: f32,
    e_max: f32,
    color_clip: f32,
) -> Image {
    let e_range = (e_max - e_min).max(1e-10);
    let (w, h, src_stride) = match axis {
        0 => (ny, nz, nx * ny),         // YZ slice: width=ny, height=nz
        1 => (nx, nz, nx * ny),         // XZ slice: width=nx, height=nz
        _ => (nx, ny, nx * ny),         // XY slice: width=nx, height=ny
    };

    let res = 128u32; // texture resolution (max dimension)
    let (tw, th) = (res.min(w as u32 * 4), res.min(h as u32 * 4));
    let mut pixels: Vec<u8> = Vec::with_capacity(tw as usize * th as usize * 4);

    for py in 0..th {
        let fy = py as f32 / (th - 1).max(1) as f32;
        for px in 0..tw {
            let fx = px as f32 / (tw - 1).max(1) as f32;

            // Map pixel (fx, fy) → field index
            let (ei, ej, ek) = match axis {
                0 => {
                    let idx = (position * (nx - 1) as f32).round() as usize;
                    (idx.min(nx - 1),
                     (fy * (ny - 1) as f32).round() as usize,
                     (fx * (nz - 1) as f32).round() as usize)
                }
                1 => {
                    ((fx * (nx - 1) as f32).round() as usize,
                     (position * (ny - 1) as f32).round() as usize,
                     (fy * (nz - 1) as f32).round() as usize)
                }
                _ => {
                    ((fx * (nx - 1) as f32).round() as usize,
                     (fy * (ny - 1) as f32).round() as usize,
                     (position * (nz - 1) as f32).round() as usize)
                }
            };

            let idx = ek * src_stride + ej.min(ny - 1) * nx + ei.min(nx - 1);
            let val = if idx < field.len() { field[idx] } else { e_max };
            // NaN (rejected atom-overlap holes) → transparent; jet_rgb(NaN)
            // would produce garbage colors at the hole edges.
            if !val.is_finite() {
                pixels.extend_from_slice(&[0, 0, 0, 0]);
                continue;
            }
            let t = ((val - e_min) / (e_range * color_clip)).clamp(0.0, 1.0);
            let (r, g, b) = jet_rgb(t);
            pixels.extend_from_slice(&[r, g, b, 220]);
        }
    }

    Image::new(
        Extent3d { width: tw, height: th, depth_or_array_layers: 1 },
        TextureDimension::D2,
        pixels,
        TextureFormat::Rgba8UnormSrgb,
        RenderAssetUsages::RENDER_WORLD,
    )
}

/// Build a quad mesh for a slice plane in the given axis-aligned orientation.
/// `axis`: 0=YZ, 1=XZ, 2=XY.
/// `position`: fractional coordinate along the axis.
pub fn slice_plane_mesh(
    lattice: &Lattice,
    axis: u8,
    position: f32,
) -> Mesh {
    let vecs = lattice.to_vectors();

    let (u_vec, v_vec, base) = match axis {
        0 => (vecs[1], vecs[2], vecs[0] * position),   // YZ: base on X
        1 => (vecs[0], vecs[2], vecs[1] * position),   // XZ: base on Y
        _ => (vecs[0], vecs[1], vecs[2] * position),   // XY: base on Z
    };

    let normal = u_vec.cross(v_vec).normalize();
    let corners = [
        base,
        base + u_vec,
        base + u_vec + v_vec,
        base + v_vec,
    ];
    let positions: Vec<[f32; 3]> = corners.iter().map(|c| [c.x, c.y, c.z]).collect();
    let normals = vec![[normal.x, normal.y, normal.z]; 4];
    let uvs: Vec<[f32; 2]> = vec![[0.0, 1.0], [1.0, 1.0], [1.0, 0.0], [0.0, 0.0]];
    let indices = vec![0u32, 1, 2, 0, 2, 3, 0, 2, 1, 0, 3, 2]; // double-sided

    let mut mesh = Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::RENDER_WORLD);
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions);
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normals);
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, uvs);
    mesh.insert_indices(Indices::U32(indices));
    mesh
}
