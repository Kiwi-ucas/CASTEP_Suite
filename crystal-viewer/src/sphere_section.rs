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
//! direction. The FIRST radial MINIMUM forms the cage shell — the migration
//! surface. Shells are made hole-free (spike suppression without rejection,
//! gap interpolation, pole merging) and are WELDED across cages where the
//! window bulges near-touch (~0.3 Å), so adjacent cages' shells merge into
//! one connected surface — the inter-cage migration path.

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
// Mode 8: radial-stationary migration surface (cage shells, welded across
// the inter-cage window bulges)
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

// ─────────────────────────────────────────────────────────────────────────
// Global vertex welding — adjacent cages' shells near-touch at the shared
// window bulges (~0.3 Å apart). Welding MERGES the near-coincident vertices
// into shared ones, so the two bulges connect into one continuous migration
// surface between cages.
// ─────────────────────────────────────────────────────────────────────────

/// Merge distance: adjacent shells at a window bulge are within a few tenths
/// of an Ångstrom of each other. Below this, two vertices become one.
/// Real-position (same image) welding, so it can be raised generously without
/// the periodic-image stretch the user saw as lines through the cell.
const WELD_TOL: f32 = 0.6;
/// Spatial-hash cell size for the weld lookup (fine enough that only the 27
/// neighbor cells need checking).
const WELD_CELL: f32 = 0.15;

/// Push a vertex, or reuse an existing one from a DIFFERENT cage within
/// WELD_TOL. Returns the global vertex index. Same-cage vertices never weld
/// (the direction grid spacing ~0.16 Å is < WELD_TOL, so welding intra-cage
/// would collapse a shell to points). Cross-cage welding connects the
/// near-touching window bulges of adjacent shells.
///
/// The distance is measured in REAL cartesian space (no periodic folding):
/// welding must only merge vertices that genuinely coincide in the rendered
/// cell. Folding the coordinates first makes two vertices that are the same
/// PHYSICAL point but different periodic IMAGES (e.g. one at cart 0.05, the
/// other at cart 10.33) weld, and then the re-used vertex sits 10 Å from the
/// welder's neighbours — a triangle stretched across the whole cell (the
/// "lines through the cell" the user saw). Real-position welding connects
/// the in-cell windows (both cages reach the saddle in the same image) and
/// leaves the cell-face wraps separate, exactly like v0.3.3.
fn weld_or_push(
    pos: Vec3, normal: Vec3, color: [f32; 4], cage: usize,
    positions: &mut Vec<[f32; 3]>, normals: &mut Vec<[f32; 3]>, colors: &mut Vec<[f32; 4]>,
    vertex_cage: &mut Vec<usize>,
    hash: &mut std::collections::HashMap<(i32, i32, i32), Vec<(usize, Vec3, usize)>>,
    cross_pairs: &mut std::collections::HashSet<(usize, usize)>,
    stats: &mut (usize, usize),
) -> usize {
    let key = (
        (pos.x / WELD_CELL).floor() as i32,
        (pos.y / WELD_CELL).floor() as i32,
        (pos.z / WELD_CELL).floor() as i32,
    );
    for dx in -1i32..=1 {
        for dy in -1i32..=1 {
            for dz in -1i32..=1 {
                if let Some(bucket) = hash.get(&(key.0 + dx, key.1 + dy, key.2 + dz)) {
                    for (vi, q, c) in bucket {
                        // Cross-cage weld (real positions only).
                        if *c != cage && (*q - pos).length() < WELD_TOL {
                            cross_pairs.insert(((*c).min(cage), (*c).max(cage)));
                            stats.1 += 1;
                            return *vi;
                        }
                    }
                }
            }
        }
    }
    stats.0 += 1;
    let vi = positions.len();
    positions.push([pos.x, pos.y, pos.z]);
    normals.push([normal.x, normal.y, normal.z]);
    colors.push(color);
    vertex_cage.push(cage);
    hash.entry(key).or_default().push((vi, pos, cage));
    vi
}

/// Build the merged migration surface: per cage, the first radial minimum
/// (cage shell) forms the migration surface, energy-capped at `e_cap` (cube
/// values are E−E_min, so this is a relative energy).
///
/// The shell gets the full spike-suppression pipeline (opt 1 pre-smoothing,
/// opt 2 prominence, opt 3 Laplacian), but none of it is allowed to create a
/// hole — prominence only picks WHICH minimum (a rejected one falls back to
/// the first valid minimum), a smoothed radius that crosses above e_cap falls
/// back to the pre-smoothing radius, remaining gaps are interpolated from the
/// grid neighbours, and the coincident pole vertices are merged (with
/// degenerate triangles dropped). Adjacent cages' shells near-touch at the
/// window bulges (~0.3 Å apart) and are WELDED (real-position) into one
/// connected surface — the inter-cage migration path.
///
/// Returns None when no cage center yields any valid surface point.
pub fn migration_surface_mesh(
    field: &[f32], n: usize, lattice: &Lattice,
    centers: &[Vec3],
    e_cap: f32, show_shell: bool,
    color_min: f32, color_max: f32, material: IsoMaterial,
) -> Option<Mesh> {
    const NT: usize = 48;
    const NP: usize = 96;
    const R_MIN: f32 = 0.60;
    const R_MAX: f32 = 4.00;
    const R_STEP: f32 = 0.05;
    const SHELL_R_MIN: f32 = 0.75;   // ignore minima inside the S-core repulsion
    // Radial pre-smoothing (opt 1): moving average before extremum detection.
    const SMOOTH_WINDOW: usize = 5;   // window (odd); NaN breaks the run
    // Extremum validation (opt 2): a shell minimum must be prominent —
    // (window max − min) ≥ MIN_DEPTH — or the scan continues (spike gate).
    const MIN_DEPTH: f32 = 0.01;      // eV
    const DEPTH_STEPS: usize = 3;     // depth-check window (3 × 0.05 Å = 0.15 Å)
    // Shell-radius Laplacian smoothing (opt 3): removes directional jumps.
    const SMOOTH_ITERS: usize = 3;    // iterations (φ wraps, θ clamps at poles)

    let n_r = 1 + ((R_MAX - R_MIN) / R_STEP).round() as usize;
    let half_w = SMOOTH_WINDOW / 2;

    let vecs = lattice.to_vectors();
    let inv = lattice.inverse_vectors();
    let range = (color_max - color_min).max(1e-6);
    let alpha = material.alpha();

    let mut positions: Vec<[f32; 3]> = Vec::new();
    let mut normals: Vec<[f32; 3]> = Vec::new();
    let mut colors: Vec<[f32; 4]> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();
    // Creator cage of each vertex — used to suppress duplicate quads: a shell
    // quad whose four vertices were ALL created by another cage (welded, so
    // that cage already renders the merged window region) is not emitted.
    let mut vertex_cage: Vec<usize> = Vec::new();

    // Global weld hash (cross-cage shared vertices) + weld counters.
    let mut weld_hash: std::collections::HashMap<(i32, i32, i32), Vec<(usize, Vec3, usize)>> =
        std::collections::HashMap::new();
    let mut weld_stats: (usize, usize) = (0, 0);   // (new, welded)
    // Distinct cage pairs sharing at least one welded vertex — direct measure
    // of inter-cage connectivity via the shared window.
    let mut cross_pairs: std::collections::HashSet<(usize, usize)> =
        std::collections::HashSet::new();
    let mut n_shell_dirs = 0usize;
    // Quads skipped because all four vertices came from another cage (a
    // fully-welded duplicate region) — the window merge de-duplication.
    let mut n_suppressed = 0usize;
    // Duplicate-face removal (belt-and-suspenders on top of quad suppression).
    let mut seen_tris: std::collections::HashSet<(u32, u32, u32)> =
        std::collections::HashSet::new();

    for (cage_idx, &center_frac) in centers.iter().enumerate() {
        let center_cart = center_frac.x * vecs[0] + center_frac.y * vecs[1] + center_frac.z * vecs[2];
        let mut shell_valid = vec![false; (NT + 1) * NP];
        let mut shell_vert: Vec<Option<usize>> = vec![None; (NT + 1) * NP];
        // Prominence-validated shell radius (opt 2), Laplacian-smoothed (opt 3),
        // then emitted.
        let mut r_shell: Vec<Option<f32>> = vec![None; (NT + 1) * NP];
        // First valid minimum radius per direction — the fallback for opt-2
        // holes, so prominence picks WHICH minimum but never creates a hole.
        let mut first_min_r: Vec<Option<f32>> = vec![None; (NT + 1) * NP];

        // ── radial scan: shell = first prominent minimum (≤ e_cap), with opt 2
        // (prominence gate) + opt 3 (Laplacian) spike suppression. ──
        for it in 0..=NT {
            for ip in 0..NP {
                let dir = sphere_dir(it, ip, NT, NP);
                let idx = it * NP + ip;
                // 1) Sample the full radial profile, then 2) moving-average
                //    smooth it (NaN breaks the run) before extremum detection.
                let mut profile = Vec::with_capacity(n_r);
                for i in 0..n_r {
                    let r = R_MIN + i as f32 * R_STEP;
                    let pos = center_cart + dir * r;
                    let frac = cart_to_frac(&inv, pos);
                    profile.push(sample_trilinear(field, n, frac));
                }
                let mut sm = vec![None; n_r];
                for i in 0..n_r {
                    if profile[i].is_none() { continue; }
                    let lo = i.saturating_sub(half_w);
                    let hi = (i + half_w).min(n_r - 1);
                    let (mut sum, mut cnt) = (0.0f32, 0u32);
                    for j in lo..=hi {
                        if let Some(v) = profile[j] { sum += v; cnt += 1; }
                    }
                    if cnt > 0 { sm[i] = Some(sum / cnt as f32); }
                }
                // Extremum scan on the smoothed profile (slope-sign threshold
                // 1e-5 so flat extrema are still detected).
                let mut prev: Option<(f32, f32)> = None;   // (r, E) previous sample
                let mut prev2: Option<(f32, f32)> = None;  // (r, E) two samples back
                let mut shell: Option<(f32, f32, Vec3)> = None;
                // First valid minimum (r ≥ SHELL_R_MIN, e ≤ e_cap) — the shell
                // falls back to THIS when prominence (opt 2) rejects the min, so
                // opt 2 picks WHICH minimum but never creates a hole.
                let mut first_min_dir: Option<(f32, f32, Vec3)> = None;
                for i in 0..n_r {
                    let r = R_MIN + i as f32 * R_STEP;
                    let e = sm[i];
                    if let (Some(pp), Some(p), Some(ev)) = (prev2, prev, e) {
                        let de_prev = p.1 - pp.1;
                        let de_cur = ev - p.1;
                        let s_prev = if de_prev > 1e-5 { 1 } else if de_prev < -1e-5 { -1 } else { 0 };
                        let s_cur = if de_cur > 1e-5 { 1 } else if de_cur < -1e-5 { -1 } else { 0 };
                        if s_prev != 0 && s_cur != 0 && s_prev != s_cur {
                            // Extremum between pp and r (parabolic refinement)
                            let denom = pp.1 - 2.0 * p.1 + ev;
                            let r_star = if denom.abs() > 1e-9 {
                                p.0 + 0.5 * R_STEP * (pp.1 - ev) / denom
                            } else { p.0 };
                            let e_star = p.1 + (r_star - p.0) / R_STEP * 0.5 * (ev - pp.1);
                            let is_min = s_prev < 0;   // − → + = minimum
                            if is_min {
                                if first_min_dir.is_none() && r_star >= SHELL_R_MIN && e_star <= e_cap {
                                    first_min_dir = Some((r_star, e_star, center_cart + dir * r_star));
                                }
                                if shell.is_none() {
                                    // Prominence (opt 2): window max − min over
                                    // ±DEPTH_STEPS must be ≥ MIN_DEPTH — rejects
                                    // spurious shallow minima (shell spikes).
                                    // (No `continue` here: skipping the prev
                                    // update would leave stale profile state and
                                    // shift later extremum positions.)
                                    if r_star >= SHELL_R_MIN && e_star <= e_cap {
                                        let lo = i.saturating_sub(DEPTH_STEPS);
                                        let hi = (i + DEPTH_STEPS).min(n_r - 1);
                                        let mut w_max = f32::NEG_INFINITY;
                                        let mut w_ok = true;
                                        for j in lo..=hi {
                                            match sm[j] {
                                                Some(v) => w_max = w_max.max(v),
                                                None => { w_ok = false; break; }
                                            }
                                        }
                                        if w_ok && w_max - e_star >= MIN_DEPTH {
                                            shell = Some((r_star, e_star, center_cart + dir * r_star));
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if e.is_none() {
                        prev = None; prev2 = None;   // NaN gap — restart the profile
                    } else {
                        prev2 = prev;
                        prev = Some((r, e.unwrap()));
                    }
                }
                if let Some((r, _, _)) = shell {
                    if show_shell { r_shell[idx] = Some(r); }   // emitted after smoothing
                }
                first_min_r[idx] = first_min_dir.map(|(r, _, _)| r);
            }
        }

        // Opt 2 must never create a hole: where prominence rejected the shell,
        // fall back to the first valid minimum. Then interpolate any remaining
        // gaps (monotonic / NaN directions) from the 4 grid neighbours so the
        // shell closes into a sphere (χ=2) instead of a holey surface.
        for i in 0..(NT + 1) * NP {
            if r_shell[i].is_none() { r_shell[i] = first_min_r[i]; }
        }
        for _ in 0..4 {
            let mut changed = false;
            let mut next = r_shell.clone();
            for it in 0..=NT {
                for ip in 0..NP {
                    let i0 = it * NP + ip;
                    if next[i0].is_some() { continue; }
                    let ipm = (ip + NP - 1) % NP;
                    let ipp = (ip + 1) % NP;
                    let (mut sum, mut cnt) = (0.0f32, 0u32);
                    if it > 0 { if let Some(r) = r_shell[(it - 1) * NP + ip] { sum += r; cnt += 1; } }
                    if it < NT { if let Some(r) = r_shell[(it + 1) * NP + ip] { sum += r; cnt += 1; } }
                    if let Some(r) = r_shell[it * NP + ipm] { sum += r; cnt += 1; }
                    if let Some(r) = r_shell[it * NP + ipp] { sum += r; cnt += 1; }
                    if cnt >= 2 { next[i0] = Some(sum / cnt as f32); changed = true; }
                }
            }
            r_shell = next;
            if !changed { break; }
        }

        // Opt 3: Laplacian smoothing of the shell radius field (NaN-aware,
        // φ wraps, θ clamps at the poles) — removes directional jumps that
        // read as spikes on the mesh — then emit shell vertices with energy
        // re-sampled at the smoothed radius (keep the e_cap cutoff).
        if show_shell {
            let r_orig = r_shell.clone();   // pre-smoothing radii, for the clamp
            let mut r_field = r_shell;
            for _ in 0..SMOOTH_ITERS {
                let mut next = vec![None; (NT + 1) * NP];
                for it in 0..=NT {
                    for ip in 0..NP {
                        let i0 = it * NP + ip;
                        let Some(r) = r_field[i0] else { continue };
                        let (mut sum, mut cnt) = (r, 1);
                        let ipm = (ip + NP - 1) % NP;
                        let ipp = (ip + 1) % NP;
                        if let Some(v) = r_field[it * NP + ipm] { sum += v; cnt += 1; }
                        if let Some(v) = r_field[it * NP + ipp] { sum += v; cnt += 1; }
                        if it > 0 { if let Some(v) = r_field[(it - 1) * NP + ip] { sum += v; cnt += 1; } }
                        if it < NT { if let Some(v) = r_field[(it + 1) * NP + ip] { sum += v; cnt += 1; } }
                        next[i0] = Some(sum / cnt as f32);
                    }
                }
                r_field = next;
            }
            for it in 0..=NT {
                for ip in 0..NP {
                    let i0 = it * NP + ip;
                    let Some(r) = r_field[i0] else { continue };
                    let dir = sphere_dir(it, ip, NT, NP);
                    let frac = cart_to_frac(&inv, center_cart + dir * r);
                    let e = sample_trilinear(field, n, frac);
                    // Opt 3 must never create a hole: if the smoothed radius's
                    // energy crosses above e_cap, fall back to the pre-smoothing
                    // radius (which is ≤ e_cap) instead of dropping the point.
                    let (r_use, e_use) = match e {
                        Some(e) if e <= e_cap => (r, Some(e)),
                        _ => {
                            let ro = r_orig[i0].unwrap_or(r);
                            (ro, sample_trilinear(field, n, cart_to_frac(&inv, center_cart + dir * ro)))
                        }
                    };
                    if let Some(e) = e_use {
                        if e <= e_cap {
                            shell_valid[i0] = true;
                            let rgb = jet_f32(((e - color_min) / range).clamp(0.0, 1.0));
                            let pos = center_cart + dir * r_use;
                            shell_vert[i0] = Some(weld_or_push(
                                pos, dir, [rgb[0], rgb[1], rgb[2], alpha], cage_idx,
                                &mut positions, &mut normals, &mut colors,
                                &mut vertex_cage, &mut weld_hash, &mut cross_pairs, &mut weld_stats));
                        }
                    }
                }
            }
        }

        // Merge the coincident pole vertices (θ=0 and θ=NT: every φ maps to
        // the same direction → the same position) into a single index, so the
        // shell closes into a true sphere (χ=2) instead of a cylinder with a
        // degenerate fan of zero-area triangles at each pole.
        let verts = &mut shell_vert;
        for pole_it in [0usize, NT] {
            let pole_idx = (0..NP).find_map(|ip| verts[pole_it * NP + ip]);
            if let Some(vi) = pole_idx {
                for ip in 0..NP {
                    if verts[pole_it * NP + ip].is_some() {
                        verts[pole_it * NP + ip] = Some(vi);
                    }
                }
            }
        }

        // Triangulate quads. A quad whose four vertices were ALL created by
        // other cages is SKIPPED — it is a fully-welded duplicate region that
        // the owning cage already renders, so the merged window keeps one
        // clean surface instead of two cross-hatched sheets. Missing corners →
        // holes.
        let mut push_sheet = |valid: &[bool], vert: &[Option<usize>], indices: &mut Vec<u32>| {
            for it in 0..NT {
                for ip in 0..NP {
                    let a = it * NP + ip;
                    let b = it * NP + (ip + 1) % NP;
                    let c = (it + 1) * NP + (ip + 1) % NP;
                    let d = (it + 1) * NP + ip;
                    if valid[a] && valid[b] && valid[c] && valid[d] {
                        let (va, vb, vc, vd) = (vert[a].unwrap() as u32, vert[b].unwrap() as u32,
                                                vert[c].unwrap() as u32, vert[d].unwrap() as u32);
                        if ![va, vb, vc, vd].iter().any(|&v| vertex_cage[v as usize] == cage_idx) {
                            n_suppressed += 1;   // fully welded — the owning cage renders it
                            continue;
                        }
                        let sorted = |x: u32, y: u32, z: u32| {
                            let mut t = [x, y, z]; t.sort(); (t[0], t[1], t[2])
                        };
                        // Drop degenerate triangles (two coincident verts, e.g.
                        // the pole fans after merging) — they're zero-area.
                        if va != vb && vb != vc && va != vc && seen_tris.insert(sorted(va, vb, vc)) {
                            indices.extend_from_slice(&[va, vb, vc]);
                        }
                        if va != vc && vc != vd && va != vd && seen_tris.insert(sorted(va, vc, vd)) {
                            indices.extend_from_slice(&[va, vc, vd]);
                        }
                    }
                }
            }
        };
        n_shell_dirs += shell_valid.iter().filter(|&&v| v).count();
        push_sheet(&shell_valid, &shell_vert, &mut indices);
    }

    eprintln!(
        "[mode8] cages={} verts={} tris={} shell_dirs={} weld=new{} merged{} pairs={} suppressed={}",
        centers.len(), positions.len(), indices.len() / 3,
        n_shell_dirs, weld_stats.0, weld_stats.1, cross_pairs.len(), n_suppressed);

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

    /// Mode 8: the shell (first radial min ≈ well at 2.4 Å) is present and
    /// energy-capped; the shell wells are low-energy (blue).
    #[test]
    fn migration_surface_shell() {
        let cube = test_cube();
        let lattice = cube.to_lattice();
        let centers = detect_cage_centers(&cube);
        let mesh = migration_surface_mesh(
            &cube.field, cube.nx, &lattice,
            &centers, 3.0, true,
            0.0, 0.5, IsoMaterial::SemiTransparent,
        ).expect("migration surface must build");
        let nv = mesh.count_vertices();
        // 4 cages × 4704 grid points, minus pole merging and welded dupes
        assert!(nv > 5000, "expected shells, got only {} vertices", nv);
        assert!(nv < 4 * 4704 * 2, "suspicious vertex count {}", nv);
        let colors = colors_of(&mesh);
        // Shell wells ≈ E_min (blue) present
        assert!(colors.iter().any(|c| c[2] > c[0] + 0.2), "no blue shell wells");
        assert!(colors.iter().all(|c| c[3] == 0.7), "SemiTransparent alpha");
    }
}
