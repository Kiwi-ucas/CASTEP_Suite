//! PES-on-surface visualization: fixed-radius sphere (mode 7) and the
//! radial-stationary migration surface (mode 8).
//!
//! Mode 7 (sphere): user picks a center (cube atom or custom fractional
//! coords) and a radius in Å. Every mesh vertex samples the cube field
//! trilinearly (periodic wrap) and is colored by energy (jet, per-vertex
//! colors via `Mesh::ATTRIBUTE_COLOR` — Bevy's pbr shader replaces
//! base_color with vertex color, alpha included).
//!
//! Mode 8 (migration): for every cage-center atom (auto-detected: an atom
//! with ≥6 Li within 3.0 Å), walk radial energy profiles along each
//! direction. The FIRST radial MINIMUM forms the cage shell (in-cage
//! migration surface); the first radial MAXIMUM beyond it forms the window
//! cap (inter-cage crossing surface, where the Li's energy peaks while
//! hopping between cages — its color shows the barrier landscape, darkest
//! spot = true saddle). Energy-capped: only points with E ≤ e_cap survive
//! (directions blocked by framework atoms become holes).

use bevy::prelude::*;
use bevy::render::mesh::{Indices, PrimitiveTopology};
use bevy::render::render_asset::RenderAssetUsages;
use crate::cube_reader::CubeData;
use crate::crystal::Lattice;
use crate::IsoMaterial;
use super::pes::jet_rgb;

/// Trilinear sample of the cubic (n)³ field at fractional coords with
/// periodic wrap. Returns None when any corner is NaN (hole point).
pub fn sample_trilinear(field: &[f32], n: usize, frac: Vec3) -> Option<f32> {
    let f = frac - frac.floor();            // wrap to [0,1)
    let g = f * (n - 1) as f32;             // grid coords in [0, n-1]
    let i0f = g.floor();
    let w = g - i0f;                        // fractional weights
    let i0 = i0f.as_uvec3() % UVec3::splat(n as u32);
    let i1 = (i0f + Vec3::ONE).as_uvec3() % UVec3::splat(n as u32);
    let at = |ix: u32, iy: u32, iz: u32| {
        field[(iz as usize * n + iy as usize) * n + ix as usize] as f32
    };
    let c000 = at(i0.x, i0.y, i0.z);
    let c100 = at(i1.x, i0.y, i0.z);
    let c010 = at(i0.x, i1.y, i0.z);
    let c110 = at(i1.x, i1.y, i0.z);
    let c001 = at(i0.x, i0.y, i1.z);
    let c101 = at(i1.x, i0.y, i1.z);
    let c011 = at(i0.x, i1.y, i1.z);
    let c111 = at(i1.x, i1.y, i1.z);
    let corners = [c000, c100, c010, c110, c001, c101, c011, c111];
    if !corners.iter().all(|v| v.is_finite()) { return None; }
    let wx = w.x; let wy = w.y; let wz = w.z;
    let x0 = c000 * (1.0 - wx) + c100 * wx;
    let x1 = c010 * (1.0 - wx) + c110 * wx;
    let y0 = x0 * (1.0 - wy) + x1 * wy;
    let x2 = c001 * (1.0 - wx) + c101 * wx;
    let x3 = c011 * (1.0 - wx) + c111 * wx;
    let y1 = x2 * (1.0 - wy) + x3 * wy;
    Some(y0 * (1.0 - wz) + y1 * wz)
}

/// Cartesian → fractional via the inverse lattice matrix.
fn cart_to_frac(inv: &[Vec3; 3], pos: Vec3) -> Vec3 {
    Lattice::apply_inverse(inv, pos)
}

fn jet_f32(t: f32) -> [f32; 3] {
    let (r, g, b) = jet_rgb(t);
    [r as f32 / 255.0, g as f32 / 255.0, b as f32 / 255.0]
}

/// Direction vector of the UV-sphere grid (theta rows, phi cols).
fn sphere_dir(it: usize, ip: usize, nt: usize, np: usize) -> Vec3 {
    let theta = std::f32::consts::PI * it as f32 / nt as f32;
    let phi = std::f32::consts::TAU * ip as f32 / np as f32;
    Vec3::new(theta.sin() * phi.cos(), theta.sin() * phi.sin(), theta.cos())
}

// ─────────────────────────────────────────────────────────────────────────
// Mode 7: fixed-radius sphere section
// ─────────────────────────────────────────────────────────────────────────

/// Build a UV sphere of radius `radius` (Å) centered at `center_frac`
/// (fractional), colored by the trilinearly sampled energy (jet, vertex
/// colors). NaN samples (holes) get a neutral gray — visible "no data".
pub fn sphere_section_mesh(
    field: &[f32], n: usize, lattice: &Lattice,
    center_frac: Vec3, radius: f32,
    color_min: f32, color_max: f32, material: IsoMaterial,
) -> Mesh {
    const NT: usize = 48;
    const NP: usize = 96;
    let vecs = lattice.to_vectors();
    let inv = lattice.inverse_vectors();
    let center_cart = center_frac.x * vecs[0] + center_frac.y * vecs[1] + center_frac.z * vecs[2];
    let range = (color_max - color_min).max(1e-6);
    let alpha = material.alpha();

    let nv = (NT + 1) * NP;
    let mut positions: Vec<[f32; 3]> = Vec::with_capacity(nv);
    let mut normals: Vec<[f32; 3]> = Vec::with_capacity(nv);
    let mut colors: Vec<[f32; 4]> = Vec::with_capacity(nv);
    for it in 0..=NT {
        for ip in 0..NP {
            let dir = sphere_dir(it, ip, NT, NP);
            let pos = center_cart + dir * radius;
            let frac = cart_to_frac(&inv, pos);
            let e = sample_trilinear(field, n, frac);
            let rgb = match e {
                Some(e) => jet_f32(((e - color_min) / range).clamp(0.0, 1.0)),
                None => [0.30, 0.30, 0.30],   // NaN hole → neutral gray
            };
            positions.push([pos.x, pos.y, pos.z]);
            normals.push([dir.x, dir.y, dir.z]);
            colors.push([rgb[0], rgb[1], rgb[2], alpha]);
        }
    }
    let mut indices: Vec<u32> = Vec::with_capacity(NT * NP * 6);
    for it in 0..NT {
        for ip in 0..NP {
            let a = (it * NP + ip) as u32;
            let b = (it * NP + (ip + 1) % NP) as u32;
            let c = ((it + 1) * NP + (ip + 1) % NP) as u32;
            let d = ((it + 1) * NP + ip) as u32;
            indices.extend_from_slice(&[a, b, c, a, c, d]);
        }
    }

    let mut mesh = Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::RENDER_WORLD);
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions);
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normals);
    mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR, colors);
    mesh.insert_indices(Indices::U32(indices));
    mesh
}

// ─────────────────────────────────────────────────────────────────────────
// Mode 8: radial-stationary migration surface (cage shells + window caps)
// ─────────────────────────────────────────────────────────────────────────

/// Auto-detect cage centers: atoms with ≥6 Li neighbors within 3.0 Å
/// (periodic distance). Returns fractional coordinates.
pub fn detect_cage_centers(cube: &CubeData) -> Vec<Vec3> {
    let lattice = cube.to_lattice();
    let vecs = lattice.to_vectors();
    let inv = lattice.inverse_vectors();
    let mut centers = Vec::new();
    for (i, a) in cube.atoms.iter().enumerate() {
        if a.z != 16 { continue; }          // S candidates
        let frac_s = cart_to_frac(&inv, Vec3::new(a.x, a.y, a.z_coord));
        let mut n_li = 0usize;
        for b in cube.atoms.iter() {
            if b.z != 3 { continue; }       // Li
            let frac_l = cart_to_frac(&inv, Vec3::new(b.x, b.y, b.z_coord));
            let mut df = frac_l - frac_s;
            df -= df.round();               // periodic wrap
            let d = (df.x * vecs[0] + df.y * vecs[1] + df.z * vecs[2]).length();
            if d < 3.0 { n_li += 1; }
        }
        if n_li >= 6 {
            let mut f = frac_s;
            f -= f.floor();
            centers.push(f);
        }
        let _ = i;
    }
    centers
}

/// Build the merged migration surface: per cage, the first radial minimum
/// (cage shell) and the first radial maximum beyond it (window cap), both
/// energy-capped at `e_cap` (cube values are E−E_min, so this is a relative
/// energy). Cap points are deduplicated across cages (spatial hash, 0.25 Å)
/// to avoid z-fighting between coincident caps of neighboring cages.
///
/// Returns None when no cage center yields any valid surface point.
pub fn migration_surface_mesh(
    field: &[f32], n: usize, lattice: &Lattice,
    centers: &[Vec3],
    e_cap: f32, show_shell: bool, show_cap: bool,
    color_min: f32, color_max: f32, material: IsoMaterial,
) -> Option<Mesh> {
    const NT: usize = 48;
    const NP: usize = 96;
    const R_MIN: f32 = 0.60;
    const R_MAX: f32 = 4.00;
    const R_STEP: f32 = 0.05;
    const SHELL_R_MIN: f32 = 0.75;   // ignore minima inside the S-core repulsion

    let vecs = lattice.to_vectors();
    let inv = lattice.inverse_vectors();
    let range = (color_max - color_min).max(1e-6);
    let alpha = material.alpha();

    let mut positions: Vec<[f32; 3]> = Vec::new();
    let mut normals: Vec<[f32; 3]> = Vec::new();
    let mut colors: Vec<[f32; 4]> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();

    // Accepted cap positions (cartesian) in a spatial hash — cross-cage dedup.
    let cell = 0.25_f32;
    let mut cap_hash: std::collections::HashMap<(i32, i32, i32), Vec<Vec3>> =
        std::collections::HashMap::new();
    let hash_key = |p: Vec3| {
        (
            (p.x / cell).floor() as i32,
            (p.y / cell).floor() as i32,
            (p.z / cell).floor() as i32,
        )
    };
    let cap_taken = |p: Vec3, hash: &std::collections::HashMap<(i32, i32, i32), Vec<Vec3>>| -> bool {
        let k = hash_key(p);
        for dx in -1i32..=1 {
            for dy in -1i32..=1 {
                for dz in -1i32..=1 {
                    if let Some(bucket) = hash.get(&(k.0 + dx, k.1 + dy, k.2 + dz)) {
                        for q in bucket {
                            if (p - *q).length() < cell { return true; }
                        }
                    }
                }
            }
        }
        false
    };

    for center_frac in centers {
        let center_cart = center_frac.x * vecs[0] + center_frac.y * vecs[1] + center_frac.z * vecs[2];
        let mut shell_valid = vec![false; (NT + 1) * NP];
        let mut cap_valid = vec![false; (NT + 1) * NP];
        let mut shell_vert: Vec<Option<usize>> = vec![None; (NT + 1) * NP];
        let mut cap_vert: Vec<Option<usize>> = vec![None; (NT + 1) * NP];

        for it in 0..=NT {
            for ip in 0..NP {
                let dir = sphere_dir(it, ip, NT, NP);
                // Radial energy profile
                let mut prev: Option<(f32, f32)> = None;   // (r, E) previous sample
                let mut prev2: Option<(f32, f32)> = None;  // (r, E) two samples back
                let mut shell: Option<(f32, f32, Vec3)> = None;
                let mut cap: Option<(f32, f32, Vec3)> = None;
                let mut r = R_MIN;
                while r <= R_MAX {
                    let pos = center_cart + dir * r;
                    let frac = cart_to_frac(&inv, pos);
                    let e = sample_trilinear(field, n, frac);
                    if let (Some(pp), Some(p), Some(ev)) = (prev2, prev, e) {
                        let de_prev = p.1 - pp.1;
                        let de_cur = ev - p.1;
                        let s_prev = if de_prev > 1e-4 { 1 } else if de_prev < -1e-4 { -1 } else { 0 };
                        let s_cur = if de_cur > 1e-4 { 1 } else if de_cur < -1e-4 { -1 } else { 0 };
                        if s_prev != 0 && s_cur != 0 && s_prev != s_cur {
                            // Extremum between pp and r (parabolic refinement)
                            let denom = pp.1 - 2.0 * p.1 + ev;
                            let r_star = if denom.abs() > 1e-9 {
                                p.0 + 0.5 * R_STEP * (pp.1 - ev) / denom
                            } else { p.0 };
                            let e_star = p.1 + (r_star - p.0) / R_STEP * 0.5 * (ev - pp.1);
                            let is_min = s_prev < 0;   // − → + = minimum
                            if is_min {
                                if shell.is_none() && r_star >= SHELL_R_MIN && e_star <= e_cap {
                                    shell = Some((r_star, e_star, center_cart + dir * r_star));
                                }
                            } else if shell.is_some() && cap.is_none() && e_star <= e_cap {
                                cap = Some((r_star, e_star, center_cart + dir * r_star));
                            }
                        }
                    }
                    if e.is_none() {
                        prev = None; prev2 = None;   // NaN gap — restart the profile
                    } else {
                        prev2 = prev;
                        prev = Some((r, e.unwrap()));
                    }
                    r += R_STEP;
                }
                let idx = it * NP + ip;
                if let Some((_, e, pos)) = shell {
                    if show_shell {
                        shell_valid[idx] = true;
                        let rgb = jet_f32(((e - color_min) / range).clamp(0.0, 1.0));
                        shell_vert[idx] = Some(positions.len());
                        positions.push([pos.x, pos.y, pos.z]);
                        normals.push([dir.x, dir.y, dir.z]);
                        colors.push([rgb[0], rgb[1], rgb[2], alpha]);
                    }
                }
                if let Some((_, e, pos)) = cap {
                    if show_cap && !cap_taken(pos, &cap_hash) {
                        cap_valid[idx] = true;
                        cap_hash.entry(hash_key(pos)).or_default().push(pos);
                        let rgb = jet_f32(((e - color_min) / range).clamp(0.0, 1.0));
                        cap_vert[idx] = Some(positions.len());
                        positions.push([pos.x, pos.y, pos.z]);
                        normals.push([dir.x, dir.y, dir.z]);
                        colors.push([rgb[0], rgb[1], rgb[2], alpha]);
                    }
                }
            }
        }

        // Triangulate quads (skip any quad with missing corners → holes)
        let push_sheet = |valid: &[bool], vert: &[Option<usize>], indices: &mut Vec<u32>| {
            for it in 0..NT {
                for ip in 0..NP {
                    let a = it * NP + ip;
                    let b = it * NP + (ip + 1) % NP;
                    let c = (it + 1) * NP + (ip + 1) % NP;
                    let d = (it + 1) * NP + ip;
                    if valid[a] && valid[b] && valid[c] && valid[d] {
                        let (va, vb, vc, vd) = (vert[a].unwrap() as u32, vert[b].unwrap() as u32,
                                                vert[c].unwrap() as u32, vert[d].unwrap() as u32);
                        indices.extend_from_slice(&[va, vb, vc, va, vc, vd]);
                    }
                }
            }
        };
        push_sheet(&shell_valid, &shell_vert, &mut indices);
        push_sheet(&cap_valid, &cap_vert, &mut indices);
    }

    if positions.is_empty() { return None; }
    let mut mesh = Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::RENDER_WORLD);
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions);
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normals);
    mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR, colors);
    mesh.insert_indices(Indices::U32(indices));
    Some(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_cube() -> CubeData {
        crate::cube_reader::parse_cube("/tmp/sphere_test.cube")
            .expect("synthetic cube must parse")
    }

    fn colors_of(mesh: &Mesh) -> &Vec<[f32; 4]> {
        match mesh.attribute(Mesh::ATTRIBUTE_COLOR).unwrap() {
            bevy::render::mesh::VertexAttributeValues::Float32x4(v) => v,
            other => panic!("unexpected color format: {:?}", other),
        }
    }

    /// Auto-detection must find the four 4b cage-center S atoms.
    #[test]
    fn detects_cage_centers() {
        let cube = test_cube();
        let centers = detect_cage_centers(&cube);
        assert_eq!(centers.len(), 4, "expected 4 cage centers, got {}", centers.len());
        // All centers must be 4b-type: coordinates 0.25/0.75 (float tolerance)
        for c in &centers {
            let on_quarter = (c.x - 0.25).abs() < 1e-3 || (c.x - 0.75).abs() < 1e-3;
            assert!(on_quarter, "unexpected center x: {}", c.x);
        }
    }

    /// Mode 7: sphere at the C1–C2 midpoint (0.5,0.5,0.25) r=2.4 Å — passes
    /// through the C1 well (blue) and the barrier side (red); the NaN box on
    /// its surface renders neutral gray.
    #[test]
    fn sphere_section_colors() {
        let cube = test_cube();
        let lattice = cube.to_lattice();
        let mesh = sphere_section_mesh(
            &cube.field, cube.nx, &lattice,
            Vec3::new(0.5, 0.5, 0.25), 2.4,
            0.0, 0.5, IsoMaterial::Opaque,
        );
        assert_eq!(mesh.count_vertices(), (48 + 1) * 96, "UV sphere vertex count");
        let colors = colors_of(&mesh);
        let max_r = colors.iter().map(|c| c[0]).fold(0.0f32, f32::max);
        let max_b = colors.iter().map(|c| c[2]).fold(0.0f32, f32::max);
        assert!(max_r >= 0.45, "no high-energy (red) vertices, max_r={}", max_r);
        assert!(max_b >= 0.45, "no low-energy (blue) vertices, max_b={}", max_b);
        assert!(colors.iter().any(|c| c[0] > c[2] + 0.2), "no red-dominant vertices");
        assert!(colors.iter().any(|c| c[2] > c[0] + 0.2), "no blue-dominant vertices");
        // NaN box on the sphere → gray vertices (0.30, 0.30, 0.30)
        let gray = colors.iter().filter(|c| (c[0] - 0.30).abs() < 1e-4).count();
        assert!(gray >= 4, "expected gray NaN vertices, got {}", gray);
        // Opaque preset → alpha 1.0
        assert!(colors.iter().all(|c| c[3] == 1.0));
    }

    /// Mode 8: shell (first radial min ≈ well at 2.4 Å) + cap (barrier
    /// between the two wells) both present; cap dedup keeps < 4×4704×2 pts.
    #[test]
    fn migration_surface_shell_and_cap() {
        let cube = test_cube();
        let lattice = cube.to_lattice();
        let centers = detect_cage_centers(&cube);
        let mesh = migration_surface_mesh(
            &cube.field, cube.nx, &lattice,
            &centers, 3.0, true, true,
            0.0, 0.5, IsoMaterial::SemiTransparent,
        ).expect("migration surface must build");
        let nv = mesh.count_vertices();
        // 4 cages × 4704 grid points — shells + caps; dedup must cut the caps
        assert!(nv > 5000, "expected shells, got only {} vertices", nv);
        assert!(nv < 4 * 4704 * 2, "suspicious vertex count {}", nv);
        let colors = colors_of(&mesh);
        // Shell wells ≈ E_min (blue) and the inter-cage barrier (red) both present
        assert!(colors.iter().any(|c| c[2] > c[0] + 0.2), "no blue shell wells");
        assert!(colors.iter().any(|c| c[0] > c[2] + 0.2), "no barrier (red) cap points");
        assert!(colors.iter().all(|c| c[3] == 0.7), "SemiTransparent alpha");
    }
}
