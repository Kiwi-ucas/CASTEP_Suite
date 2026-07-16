//! PES (Potential Energy Surface) data types and surface mesh generation.
//!
//! Parses `pes_metadata.json` produced by Fortran PosCASTEP PES scan module.
//! Generates a 3D surface mesh with vertex colors from a jet colormap.

use bevy::prelude::*;
use bevy::render::mesh::{Indices, PrimitiveTopology};
use bevy::render::render_asset::RenderAssetUsages;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};
use serde::{Deserialize, Serialize};
use crate::crystal::Lattice;

// ── PES data types (match Fortran write_pes_metadata_json output) ──

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PesData {
    #[serde(rename = "type", default)]
    pub data_type: Option<String>,
    pub plane: String,
    pub nx: usize,
    pub ny: usize,
    pub fx_range: [f64; 2],
    pub fy_range: [f64; 2],
    pub mobile_atom: MobileAtom,
    pub scan_mode: String,
    pub lattice: Lattice,
    pub structure_atoms: Vec<PesAtom>,
    pub energies: Vec<Option<f64>>,
    pub has_energies: bool,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct MobileAtom {
    pub index: i32,
    pub element: String,
}

/// PES atom with fractional coordinates (fx/fy/fz vs AtomData's x/y/z).
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PesAtom {
    pub element: String,
    pub fx: f64,
    pub fy: f64,
    pub fz: f64,
}

impl PesData {
    /// Load PES metadata from a JSON file.
    pub fn from_json(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let content = std::fs::read_to_string(path)?;
        Ok(serde_json::from_str(&content)?)
    }

    /// Quick check whether JSON content is a PES scan file.
    pub fn detect(json_str: &str) -> bool {
        json_str.contains("\"pes_scan\"")
    }

    /// Return (plane_axis_0, plane_axis_1) as indices into lattice vectors [a,b,c].
    pub fn plane_axes(&self) -> (usize, usize) {
        match self.plane.as_str() {
            "xz" => (0, 2),
            "yz" => (1, 2),
            _     => (0, 1), // "xy" or default
        }
    }

    /// Generate surface mesh + colormap texture, or None if no energies.
    ///
    /// Returns `(Mesh, Image)` where:
    /// - Mesh: subdivided triangle mesh (bilinear interpolation, passes through data)
    /// - Image: 256×1 colormap lookup texture (jet: blue→cyan→green→yellow→red)
    ///
    /// Null/missing grid points are treated as very-high energy, clipped at a fixed
    /// threshold above the valid range → flat red plateau ("volcano crater rim").
    pub fn generate_surface(&self) -> Option<(Mesh, Image)> {
        if !self.has_energies || self.nx < 2 || self.ny < 2 {
            return None;
        }

        // ── Subdivide grid for smoothness ──
        let (sub_energies, nx, ny) = self.subdivide_energies();
        let n = nx * ny;

        let lattice_vecs = self.lattice.to_vectors();
        let (pa0, pa1) = self.plane_axes();
        let normal = lattice_vecs[pa0].cross(lattice_vecs[pa1]).normalize();

        // ── Energy range from valid subdivided points ──
        let mut e_min = f64::MAX;
        let mut e_max_valid = f64::MIN;
        for e in &sub_energies {
            if let Some(v) = e {
                e_min = e_min.min(*v);
                e_max_valid = e_max_valid.max(*v);
            }
        }
        if e_min >= e_max_valid {
            return None;
        }
        let e_range = e_max_valid - e_min;

        // Clip threshold: 30% above valid max. Null points assigned 5× range above valid max.
        let e_clip = e_max_valid + e_range * 0.3;
        let e_fill = e_max_valid + e_range * 5.0;
        let clip_range = e_clip - e_min;

        // Height scale: ~30% of the third lattice vector length
        let third_idx = 3 - pa0 - pa1; // 0+1→2, 0+2→1, 1+2→0
        let max_height = lattice_vecs[third_idx].length() * 0.3;

        // Base plane = mobile atom's out-of-plane fractional coordinate
        let mobile = &self.structure_atoms[self.mobile_atom.index as usize];
        let frac_perp: f32 = match third_idx {
            1 => mobile.fy as f32,
            2 => mobile.fz as f32,
            _ => mobile.fx as f32,
        };

        // ── Build vertex arrays ──
        let mut positions: Vec<[f32; 3]> = Vec::with_capacity(n);
        let mut uvs: Vec<[f32; 2]> = Vec::with_capacity(n);
        // Track whether each vertex exceeds clip threshold
        let mut over_clip: Vec<bool> = Vec::with_capacity(n);

        for j in 0..ny {
            let fy = self.fy_range[0]
                + (self.fy_range[1] - self.fy_range[0]) * (j as f64) / ((ny - 1).max(1) as f64);
            for i in 0..nx {
                let fx = self.fx_range[0]
                    + (self.fx_range[1] - self.fx_range[0]) * (i as f64) / ((nx - 1).max(1) as f64);

                let mut frac = Vec3::ZERO;
                frac[pa0] = fx as f32;
                frac[pa1] = fy as f32;
                frac[third_idx] = frac_perp;
                let cart = frac.x * lattice_vecs[0]
                         + frac.y * lattice_vecs[1]
                         + frac.z * lattice_vecs[2];

                // Energy → capped height + colormap
                let raw_e = sub_energies[j * nx + i].unwrap_or(e_fill);
                let capped = raw_e.min(e_clip);
                let t = ((capped - e_min) / clip_range) as f32;
                let height = t.max(0.0) * max_height;

                let pos = cart + normal * height;
                positions.push([pos.x, pos.y, pos.z]);
                uvs.push([t, 0.0]);
                over_clip.push(raw_e > e_clip);
            }
        }

        // ── Build triangle indices: skip only where ALL 3 vertices exceed e_clip ──
        let mut indices: Vec<u32> = Vec::new();
        for j in 0..ny - 1 {
            for i in 0..nx - 1 {
                let a = (j * nx + i) as u32;
                let b = a + 1;
                let c = a + nx as u32;
                let d = c + 1;

                let skip_abc = over_clip[a as usize] && over_clip[b as usize] && over_clip[c as usize];
                let skip_bdc = over_clip[b as usize] && over_clip[d as usize] && over_clip[c as usize];

                if !skip_abc {
                    indices.extend_from_slice(&[a, b, c, a, c, b]);
                }
                if !skip_bdc {
                    indices.extend_from_slice(&[b, d, c, b, c, d]);
                }
            }
        }

        // ── Compute smooth vertex normals ──
        let mut normals: Vec<[f32; 3]> = vec![[0.0; 3]; n];
        for tri in indices.chunks(3) {
            let i0 = tri[0] as usize;
            let i1 = tri[1] as usize;
            let i2 = tri[2] as usize;
            let p0 = Vec3::from_array(positions[i0]);
            let p1 = Vec3::from_array(positions[i1]);
            let p2 = Vec3::from_array(positions[i2]);
            let face_n = (p1 - p0).cross(p2 - p0);
            let len = face_n.length();
            if len < 1e-10 { continue; }
            let face_n = face_n / len;
            normals[i0] = [normals[i0][0] + face_n.x, normals[i0][1] + face_n.y, normals[i0][2] + face_n.z];
            normals[i1] = [normals[i1][0] + face_n.x, normals[i1][1] + face_n.y, normals[i1][2] + face_n.z];
            normals[i2] = [normals[i2][0] + face_n.x, normals[i2][1] + face_n.y, normals[i2][2] + face_n.z];
        }
        for nrm in normals.iter_mut() {
            let len = (nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]).sqrt();
            if len > 1e-10 {
                nrm[0] /= len; nrm[1] /= len; nrm[2] /= len;
            } else {
                nrm[0] = normal.x; nrm[1] = normal.y; nrm[2] = normal.z;
            }
        }

        // ── Build mesh ──
        let mut mesh = Mesh::new(
            PrimitiveTopology::TriangleList,
            RenderAssetUsages::RENDER_WORLD,
        );
        mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions);
        mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normals);
        mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, uvs);
        mesh.insert_indices(Indices::U32(indices));

        // ── Build colormap texture (256×1, jet: blue→cyan→green→yellow→red) ──
        let tex_w = 256u32;
        let tex_size = Extent3d {
            width: tex_w,
            height: 1,
            depth_or_array_layers: 1,
        };
        let mut pixels: Vec<u8> = Vec::with_capacity(tex_w as usize * 4);
        for x in 0..tex_w {
            let t = x as f32 / (tex_w - 1).max(1) as f32;
            let (r, g, b) = jet_rgb(t);
            pixels.extend_from_slice(&[r, g, b, 200]);
        }
        let image = Image::new(
            tex_size,
            TextureDimension::D2,
            pixels,
            TextureFormat::Rgba8UnormSrgb,
            RenderAssetUsages::RENDER_WORLD,
        );

        Some((mesh, image))
    }

    /// Compute (min, max) energy across the grid.
    pub fn energy_range(&self) -> Option<(f64, f64)> {
        let mut e_min = f64::MAX;
        let mut e_max = f64::MIN;
        let mut any = false;
        for e in &self.energies {
            if let Some(val) = e {
                e_min = e_min.min(*val);
                e_max = e_max.max(*val);
                any = true;
            }
        }
        if any { Some((e_min, e_max)) } else { None }
    }

    /// Subdivide the energy grid using bilinear interpolation.
    ///
    /// Each original cell is split into `N×N` sub-cells (N = SUBDIVISION).
    /// Original data points at `(i*N, j*N)` are preserved exactly.
    /// Interior points use bilinear interpolation from the 4 surrounding corners.
    /// Null corners → sub-points near them also become null (conservative).
    const SUBDIVISION: usize = 4;

    fn subdivide_energies(&self) -> (Vec<Option<f64>>, usize, usize) {
        let n = Self::SUBDIVISION;
        let snx = (self.nx - 1) * n + 1;
        let sny = (self.ny - 1) * n + 1;
        let mut sub = vec![None; snx * sny];

        for sj in 0..sny {
            let j0 = sj / n;                    // original row index (floor)
            let j1 = ((sj / n) + 1).min(self.ny - 1); // next original row
            let fj = (sj % n) as f64 / n as f64; // fractional within cell (0..1)

            for si in 0..snx {
                let i0 = si / n;
                let i1 = ((si / n) + 1).min(self.nx - 1);
                let fi = (si % n) as f64 / n as f64;

                let v00 = self.energies[j0 * self.nx + i0];
                let v10 = self.energies[j0 * self.nx + i1];
                let v01 = self.energies[j1 * self.nx + i0];
                let v11 = self.energies[j1 * self.nx + i1];

                sub[sj * snx + si] = if si % n == 0 && sj % n == 0 {
                    // Original grid point — preserved exactly
                    v00
                } else if let (Some(e00), Some(e10), Some(e01), Some(e11)) =
                    (v00, v10, v01, v11)
                {
                    // All 4 corners valid → bilinear interpolation
                    Some(
                        e00 * (1.0 - fi) * (1.0 - fj)
                            + e10 * fi * (1.0 - fj)
                            + e01 * (1.0 - fi) * fj
                            + e11 * fi * fj,
                    )
                } else {
                    // Any corner is null → sub-point stays null
                    None
                };
            }
        }

        (sub, snx, sny)
    }
}

// ── Colormap helpers ──

/// Jet colormap (blue → cyan → green → yellow → red).
/// Returns (r, g, b) as u8 values 0–255.
fn jet_rgb(t: f32) -> (u8, u8, u8) {
    let t = t.clamp(0.0, 1.0);
    let r = if t < 0.375 {
        0.0
    } else if t < 0.625 {
        (t - 0.375) / 0.25
    } else if t < 0.875 {
        1.0
    } else {
        1.0 - (t - 0.875) / 0.125 * 0.5
    };
    let g = if t < 0.125 {
        0.0
    } else if t < 0.375 {
        (t - 0.125) / 0.25
    } else if t < 0.625 {
        1.0
    } else if t < 0.875 {
        1.0 - (t - 0.625) / 0.25
    } else {
        0.0
    };
    let b = if t < 0.125 {
        0.5 + t / 0.125 * 0.5
    } else if t < 0.375 {
        1.0
    } else if t < 0.625 {
        1.0 - (t - 0.375) / 0.25
    } else {
        0.0
    };
    ((r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8)
}

/// Map energy value to a jet colormap Color (for vertex/UI display).
#[allow(dead_code)]
pub fn energy_color(e: f64, e_min: f64, e_max: f64) -> Color {
    let t = if e_max > e_min {
        ((e - e_min) / (e_max - e_min)).clamp(0.0, 1.0) as f32
    } else {
        0.5
    };
    let (r, g, b) = jet_rgb(t);
    Color::srgb_u8(r, g, b)
}
