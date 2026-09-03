//! Slab cross-section (general (hkl) Miller indices) and vacuum layer builders.
//!
//! Both builders work on the viewer's canonical frame: input atoms are
//! Cartesian w.r.t. `Lattice::to_vectors()`, and the returned `CrystalData`
//! is likewise Cartesian w.r.t. the new lattice's standard-setting frame.
//!
//! The builders are independent: vacuum can be applied to any structure
//! (bulk or already-cut slab). After a slab cut, the new c axis is the
//! surface normal, so "vacuum along c" adds vacuum along the slab normal.
//!
//! In-plane cells:
//! - `InPlaneBasis::Primitive` — primitive cell of the *true* 2D Bravais
//!   lattice of the layer (translation stabilizer of the first layer's
//!   in-plane atom set; e.g. fcc (111) → hexagonal, 1 atom/layer;
//!   fcc (001) → square a/√2, 1 atom/layer).
//! - `InPlaneBasis::Orthogonal` — a right-angle (≈90°) in-plane cell,
//!   searched among small-integer combinations of the primitive vectors,
//!   yielding a 90°/90°/90° 3D cell when the layer lattice allows.
//! - `InPlaneBasis::Conventional` — in-plane cell of the conventional
//!   lattice only (legacy behaviour).
//! The in-plane cell is expanded by user factors U (horizontal) × V
//! (vertical) of the chosen basis pair.

use crate::crystal::{AtomData, CrystalData, InPlaneBasis, Lattice, SlabMeta, VacuumMeta};
use bevy::prelude::Vec3;

/// Slab cross-section parameters.
#[derive(Debug, Clone, Copy)]
pub struct SlabParams {
    pub h: i32,
    pub k: i32,
    pub l: i32,
    /// Slab start position along the surface normal, in Angstrom.
    pub start_ang: f32,
    /// Slab thickness along the normal, in Angstrom. 0.0 → 3·d_hkl.
    pub thickness_ang: f32,
    /// In-plane horizontal / vertical expansion: the in-plane cell is the
    /// U×V supercell of the chosen 2D basis pair (1 = no expansion).
    pub u: u8,
    pub v: u8,
    /// Which in-plane 2D basis to build from.
    pub basis: InPlaneBasis,
    /// Orthogonal candidate index (0 = smallest area) used when
    /// `basis == InPlaneBasis::Orthogonal`.
    pub orth_idx: u8,
    /// MS-style explicit in-plane vectors: the new in-plane cell vectors
    /// as integer combinations of the conventional lattice vectors
    /// (`ap = u_vec·(a,b,c)`, Materials Studio supercell-matrix rows,
    /// e.g. `(0,1,0)` = the b axis). When both are set they override
    /// `u`/`v` and `basis`; both must be strictly in-plane w.r.t. the
    /// (hkl) normal and non-collinear.
    pub u_vec: Option<[i32; 3]>,
    pub v_vec: Option<[i32; 3]>,
}

impl Default for SlabParams {
    fn default() -> Self {
        Self {
            h: 0,
            k: 0,
            l: 0,
            start_ang: 0.0,
            thickness_ang: 0.0,
            u: 1,
            v: 1,
            basis: InPlaneBasis::Primitive,
            orth_idx: 0,
            u_vec: None,
            v_vec: None,
        }
    }
}

/// Vacuum layer parameters.
#[derive(Debug, Clone, Copy)]
pub struct VacuumParams {
    /// 1, 2, 3 → extend the a / b / c axis.
    pub axis: u8,
    /// Vacuum thickness in Angstrom (0 = no vacuum; cell thickness =
    /// structure thickness).
    pub thickness_ang: f32,
    /// Legacy placement: false = whole vacuum on top; true = split V/2 on
    /// each side. Superseded by `position` whenever `position > 0`.
    pub both_sides: bool,
    /// Vertical position of the structure inside the vacuum-extended cell,
    /// 0..=1: 0 = structure at the bottom (vacuum on top), 1 = at the top
    /// (vacuum below), 0.5 = centred. When 0.0, `both_sides` decides.
    pub position: f32,
}

impl Default for VacuumParams {
    fn default() -> Self {
        Self {
            axis: 3,
            thickness_ang: 15.0,
            both_sides: false,
            position: 0.0,
        }
    }
}

/// Inverse of a 3×3 matrix given as column vectors. `None` when singular.
fn inv3(cols: [Vec3; 3]) -> Option<[Vec3; 3]> {
    let m = [
        [cols[0].x, cols[1].x, cols[2].x],
        [cols[0].y, cols[1].y, cols[2].y],
        [cols[0].z, cols[1].z, cols[2].z],
    ];
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() < 1e-12 {
        return None;
    }
    let i = 1.0 / det;
    Some([
        Vec3::new(
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * i,
            -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) * i,
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * i,
        ),
        Vec3::new(
            -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) * i,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * i,
            -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) * i,
        ),
        Vec3::new(
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * i,
            -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) * i,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * i,
        ),
    ])
}

/// Multiply inverse-matrix rows by vector: `M⁻¹ · v`.
fn apply_inv3(inv: &[Vec3; 3], v: Vec3) -> Vec3 {
    Vec3::new(
        inv[0].x * v.x + inv[1].x * v.y + inv[2].x * v.z,
        inv[0].y * v.x + inv[1].y * v.y + inv[2].y * v.z,
        inv[0].z * v.x + inv[1].z * v.y + inv[2].z * v.z,
    )
}

/// Interplanar spacing d_hkl of the (h,k,l) planes of this cell, in
/// Angstrom. `None` for the (0,0,0) index.
pub fn slab_d_hkl(data: &CrystalData, h: i32, k: i32, l: i32) -> Option<f32> {
    let [ia, ib, ic] = data.lattice.inverse_vectors();
    let g = (h as f32) * ia + (k as f32) * ib + (l as f32) * ic;
    let len = g.length();
    if len < 1e-9 {
        return None;
    }
    Some(1.0 / len)
}

fn gcd(a: i32, b: i32) -> i32 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// Extended GCD: returns `(g, x, y)` with `a·x + b·y = g` (|g| = gcd(a,b);
/// the sign of `g` may follow `a`, callers normalise).

/// (hkl) layer period `gidx·d_hkl` (gidx = gcd(h,k,l)): the normal span of
/// the repeating Bézout vector. The UI's fractional s / T inputs are
/// fractions of this period.
pub fn slab_period(data: &CrystalData, h: i32, k: i32, l: i32) -> Option<f32> {
    let d = slab_d_hkl(data, h, k, l)?;
    let gidx = gcd(gcd(h, k), l).max(1);
    Some(gidx as f32 * d)
}

/// Nearest atomic layer at or above (`up = true`) / at or below
/// (`up = false`) the cut position `s`, measured along the surface normal,
/// in Angstrom. Snaps the cut to an exact atomic plane so the slab bottom
/// face passes through an atom layer instead of landing between layers.
pub fn slab_layer_snap(data: &CrystalData, h: i32, k: i32, l: i32, s: f32, up: bool) -> Option<f32> {
    if data.atoms.is_empty() {
        return None;
    }
    let [ia, ib, ic] = data.lattice.inverse_vectors();
    let n = (h as f32) * ia + (k as f32) * ib + (l as f32) * ic;
    let gl = n.length();
    if gl < 1e-9 {
        return None;
    }
    let n_hat = n / gl;
    let [va, vb, vc] = data.lattice.to_vectors();
    let gidx = gcd(gcd(h, k), l).max(1);
    let span = gidx as f32 / gl;
    if span < 1e-9 {
        return None;
    }

    let inv = data.lattice.inverse_vectors();
    let search_w = 2.0 * span; // window around s (covers ≥ 4 layer periods)
    let tol = (span * 1e-6).max(1e-4);
    let mut best: Option<f32> = None;
    for at in &data.atoms {
        let f = if data.positions_fractional {
            Vec3::new(at.x, at.y, at.z)
        } else {
            apply_inv3(&inv, Vec3::new(at.x, at.y, at.z))
        };
        let p0 = (f.x * va + f.y * vb + f.z * vc).dot(n_hat);
        let q_lo = (((s - search_w - p0 - tol) / span).floor()) as i32;
        let q_hi = (((s + search_w - p0 + tol) / span).ceil()) as i32;
        for q in q_lo..=q_hi {
            let z = p0 + (q as f32) * span;
            if up {
                if z > s + tol && (best.is_none() || z < best.unwrap()) {
                    best = Some(z);
                }
            } else if z < s - tol && (best.is_none() || z > best.unwrap()) {
                best = Some(z);
            }
        }
    }
    best
}

/// Effective vertical position fraction (0..=1): the fraction of the
/// vacuum thickness that lies *below* the structure in the
/// vacuum-extended cell.
pub fn vacuum_position(p: &VacuumParams) -> f32 {
    if p.position > 0.0 {
        p.position.clamp(0.0, 1.0)
    } else if p.both_sides {
        0.5
    } else {
        0.0
    }
}

/// Search for the in-plane primitive basis of the (h,k,l) planes from
/// *conventional-cell translations only*. Returns `(ap, bp, found_area)`
/// with `ap×bp` right-handed w.r.t. `n_hat`.
///
/// Candidate translations are the *exactly in-plane* integer combinations
/// `(i·a + j·b + k·c)` (n̂·m ≈ 0), |i|,|j|,|k| ≤ 2 (expanding to ≤ 4 when
/// nothing is found), and the pair with the smallest area is kept. This
/// is the in-plane primitive cell of the *conventional* lattice:
/// - primitive (P) input cells → exact primitive in-plane cell;
/// - centered (F/I) conventional cells → possibly a supercell of the
///   true Bravais-primitive in-plane cell (e.g. fcc (111) → 4-atom layer).
fn find_inplane_basis(
    data: &CrystalData,
    h: i32,
    k: i32,
    l: i32,
    n_hat: Vec3,
    va: Vec3,
    vb: Vec3,
    vc: Vec3,
    g: f32,
) -> Option<(Vec3, Vec3, f32)> {
    let ref_len = va.length().max(vb.length()).max(vc.length()).max(1.0);
    let inplane_tol = 1e-5 * ref_len;

    // In-plane lattice vectors: n̂·(i·a + j·b + k·c) ≈ 0.
    let collect = |limit: i32| -> Vec<Vec3> {
        let mut out = Vec::new();
        for i in -limit..=limit {
            for j in -limit..=limit {
                for kk in -limit..=limit {
                    if i == 0 && j == 0 && kk == 0 {
                        continue;
                    }
                    let m = (i as f32) * va + (j as f32) * vb + (kk as f32) * vc;
                    if m.dot(n_hat).abs() < inplane_tol && m.length() > 0.02 * ref_len {
                        out.push(m);
                    }
                }
            }
        }
        out
    };
    let mut cand = collect(2);
    if cand.len() < 2 {
        cand = collect(4); // larger-index (hkl) on exotic cells
    }
    if cand.len() < 2 {
        return None;
    }

    let mut best: Option<(usize, usize, f32)> = None;
    let area_eps = 1e-6 * ref_len * ref_len;
    for i in 0..cand.len() {
        for j in (i + 1)..cand.len() {
            let area = cand[i].cross(cand[j]).dot(n_hat).abs();
            if area < area_eps {
                continue; // collinear
            }
            if best.is_none() || area < best.unwrap().2 {
                best = Some((i, j, area));
            }
        }
    }
    let (bi, bj, area_found) = best?;
    // Theoretical kernel area of this cell: V·|n|/q (q = gcd of h,k,l).
    // A smaller found area is a finer primitive division (good); a larger
    // one means a supercell was found — warn but keep it.
    let q = gcd(gcd(h, k), l).max(1);
    let a_expected = data.lattice.cell_volume() * g / q as f32;
    if area_found > a_expected * 1.05 {
        eprintln!(
            "  [slab] in-plane cell is a supercell (found {:.3} Å² vs primitive ≤ {:.3} Å²)",
            area_found, a_expected
        );
    }
    let mut ap = cand[bi];
    let mut bp = cand[bj];
    // 2D Gaussian reduction: keep the lattice (unimodular ops) but turn the
    // sheared search basis into a short, near-orthogonal one so the new
    // cell parameters stay intuitive.
    for _ in 0..8 {
        if ap.length() > bp.length() {
            std::mem::swap(&mut ap, &mut bp);
        }
        if bp.length() < 1e-9 {
            break;
        }
        let t = (bp.dot(ap) / ap.dot(ap)).round();
        bp -= t * ap;
        if bp.length() > ap.length() {
            std::mem::swap(&mut ap, &mut bp);
        }
        if ap.length() < 1e-9 {
            break;
        }
        let t = (ap.dot(bp) / bp.dot(bp)).round();
        ap -= t * bp;
    }
    // Right-handed with the normal: a'×b'·n̂ > 0
    if ap.cross(bp).dot(n_hat) < 0.0 {
        std::mem::swap(&mut ap, &mut bp);
    }
    Some((ap, bp, area_found))
}

// ────────────────────────── true 2D Bravais lattice ──────────────────────────

/// The *true* in-plane 2D Bravais lattice of the slab: the in-plane part
/// of the crystal's translation group — every structure-preserving
/// translation `t` with n̂·t = 0 (conventional-cell combinations and
/// centering vectors, half-integer window |q|,|r|,|s| ≤ 2), reduced to a
/// primitive right-handed pair `(u2, w2)` w.r.t. `n_hat`
/// (fcc (111) → hexagonal a/√2, 1 atom/layer; fcc (001) → square a/√2).
/// Falls back to the conventional in-plane basis `(ap0, bp0)`.
fn find_layer_bravais(
    data: &CrystalData,
    n_hat: Vec3,
    ap0: Vec3,
    bp0: Vec3,
    ref_len: f32,
) -> (Vec3, Vec3, Vec3, Vec3) {
    let [va, vb, vc] = data.lattice.to_vectors();
    let inv = data.lattice.inverse_vectors();
    let frac: Vec<Vec3> = data
        .atoms
        .iter()
        .map(|at| {
            if data.positions_fractional {
                Vec3::new(at.x, at.y, at.z)
            } else {
                apply_inv3(&inv, Vec3::new(at.x, at.y, at.z))
            }
        })
        .collect();
    let elems: Vec<&str> = data.atoms.iter().map(|at| at.element.as_str()).collect();

    // Candidate in-plane translations: fractional (q/2, r/2, s/2) with
    // q, r, s ∈ [−4, 4] (conventional + all Bravais centerings).
    let in_plane_tol = 1e-4 * ref_len;
    let frac_tol = 5e-3; // ≈ 0.018 Å at a = 3.6
    let mut cand: Vec<Vec3> = Vec::new();
    for q in -4..=4 {
        for r in -4..=4 {
            for s in -4..=4 {
                if q == 0 && r == 0 && s == 0 {
                    continue;
                }
                let tf = Vec3::new(q as f32 * 0.5, r as f32 * 0.5, s as f32 * 0.5);
                let t = tf.x * va + tf.y * vb + tf.z * vc;
                if n_hat.dot(t).abs() > in_plane_tol {
                    continue; // not in-plane
                }
                // Structure preservation: every atom maps onto an atom
                // of the same element (mod conventional translations).
                let mod1 = |x: f32| -> f32 {
                    let m = x.rem_euclid(1.0);
                    m.min(1.0 - m)
                };
                let mut ok = true;
                for i in 0..frac.len() {
                    let tx = frac[i].x + tf.x;
                    let ty = frac[i].y + tf.y;
                    let tz = frac[i].z + tf.z;
                    let mut found = false;
                    for j in 0..frac.len() {
                        if elems[i] != elems[j] {
                            continue;
                        }
                        let dx = mod1(frac[j].x - tx);
                        let dy = mod1(frac[j].y - ty);
                        let dz = mod1(frac[j].z - tz);
                        if dx < frac_tol && dy < frac_tol && dz < frac_tol {
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        ok = false;
                        break;
                    }
                }
                if ok {
                    cand.push(t);
                }
            }
        }
    }
    if cand.is_empty() {
        return (ap0, bp0, ap0, bp0);
    }

    // Primitive pair: shortest candidate, then shortest non-collinear.
    let mut sorted = cand.clone();
    sorted.sort_by(|a, b| a.length().partial_cmp(&b.length()).unwrap());
    let u2 = sorted[0];
    let mut w2: Option<Vec3> = None;
    for t in sorted.iter().skip(1) {
        let c = t.dot(u2) / u2.length().powi(2);
        let resid = (t - c * u2).length();
        if resid > 0.01 * t.length().max(1e-6) {
            w2 = Some(*t);
            break;
        }
    }
    let Some(w2) = w2 else {
        return (ap0, bp0, ap0, bp0); // rank-1 in-plane lattice (degenerate)
    };

    // Closure check: every candidate must lie in the Z-span of (u2, w2)
    // (2D Cramer in the orthonormal in-plane frame u = ap0̂, v = n̂×u).
    let u = ap0 / ap0.length();
    let v = n_hat.cross(u);
    let u2u = u2.dot(u);
    let u2v = u2.dot(v);
    let w2u = w2.dot(u);
    let w2v = w2.dot(v);
    let d = u2u * w2v - w2u * u2v;
    if d.abs() < 1e-12 {
        return (ap0, bp0, ap0, bp0);
    }
    let mut closed = true;
    for c in &cand {
        let cu = c.dot(u);
        let cv = c.dot(v);
        let alpha = (w2v * cu - w2u * cv) / d;
        let beta = (u2u * cv - u2v * cu) / d;
        if alpha.round() - alpha.abs() > 0.02 || beta.round() - beta.abs() > 0.02 {
            closed = false;
            break;
        }
    }
    if !closed {
        eprintln!(
            "  [slab] in-plane lattice closure check failed — using the conventional in-plane basis"
        );
        return (ap0, bp0, ap0, bp0);
    }
    if u2.cross(w2).dot(n_hat) < 0.0 {
        return (u2, -w2, ap0, bp0);
    }
    (u2, w2, ap0, bp0)
}

/// Orthogonal (≈90°) in-plane candidate pairs of the L2 primitive basis
/// `(u2, w2)`: small-integer combinations with |i|,|j| ≤ 4, p·q ≈ 0,
/// right-handed, sorted by area, multiples of an earlier candidate
/// removed. Empty when the layer lattice admits no right-angle pair
/// within the search bound.
pub fn ortho_candidates(u2: &Vec3, w2: &Vec3, n_hat: &Vec3) -> Vec<(Vec3, Vec3)> {
    let mut out: Vec<(Vec3, Vec3, f32)> = Vec::new();
    for i in -4..=4 {
        for j in -4..=4 {
            if i == 0 && j == 0 {
                continue;
            }
            let p = (i as f32) * u2 + (j as f32) * w2;
            let pl = p.length();
            if pl < 1e-6 {
                continue;
            }
            for i2 in -4..=4 {
                for j2 in -4..=4 {
                    if i2 == 0 && j2 == 0 {
                        continue;
                    }
                    let q = (i2 as f32) * u2 + (j2 as f32) * w2;
                    let ql = q.length();
                    if ql < 1e-6 {
                        continue;
                    }
                    // |cos θ| < 1e-3 → within 0.06° of 90°
                    if p.dot(q).abs() > 1e-3 * pl * ql {
                        continue;
                    }
                    if p.cross(q).dot(*n_hat) < 0.0 {
                        continue; // right-handed w.r.t. the normal
                    }
                    let area = p.cross(q).length();
                    let dup = out.iter().any(|(po, qo, _)| {
                        po.distance(p) < 1e-3 * pl && qo.distance(q) < 1e-3 * ql
                    });
                    if dup {
                        continue;
                    }
                    out.push((p, q, area));
                }
            }
        }
    }
    out.sort_by(|x, y| x.2.partial_cmp(&y.2).unwrap());
    // Drop multiples of an earlier (smaller) candidate in the same direction.
    let mut pruned: Vec<(Vec3, Vec3, f32)> = Vec::new();
    for c in out {
        let redundant = pruned.iter().any(|(po, qo, _)| {
            let pl = c.0.length();
            let ql = c.1.length();
            let dup_p = po.dot(c.0).abs() > 0.999 * pl * po.length()
                && pl >= 2.0 * po.length() - 1e-6;
            let dup_q = qo.dot(c.1).abs() > 0.999 * ql * qo.length()
                && ql >= 2.0 * qo.length() - 1e-6;
            dup_p || dup_q
        });
        if !redundant {
            pruned.push(c);
        }
    }
    pruned.into_iter().map(|(p, q, _)| (p, q)).collect()
}

/// Live in-plane cell info for the UI: primitive L2 dimensions, the layer
/// atom density and the orthogonal candidates.
pub struct SlabLayerInfo {
    /// L2 primitive basis dimensions (Å) and included angle (deg).
    pub a_prim: f32,
    pub b_prim: f32,
    pub gamma_prim: f32,
    /// Atoms per layer in the L2 primitive cell.
    pub atoms_per_layer: u32,
    /// Orthogonal candidates: (|p|, |q|, area) in Å, Å, Å².
    pub ortho: Vec<(f32, f32, f32)>,
}

/// Compute the in-plane basis actually used by `build_slab` for `p`
/// (honoring U/V, basis choice and orthogonal candidate) — for preview.
/// Returns (ap, bp, n_hat, d_hkl, basis_used).
pub fn inplane_basis_params(
    data: &CrystalData,
    p: &SlabParams,
) -> Option<(Vec3, Vec3, Vec3, f32, InPlaneBasis)> {
    if data.atoms.is_empty() {
        return None;
    }
    let [ia, ib, ic] = data.lattice.inverse_vectors();
    let n_recip = (p.h as f32) * ia + (p.k as f32) * ib + (p.l as f32) * ic;
    let g = n_recip.length();
    if g < 1e-9 {
        return None;
    }
    let n_hat = n_recip / g;
    let (ap, bp, basis_used) = layer_bravais_pipeline(data, p, n_hat, g)?;
    Some((ap, bp, n_hat, 1.0 / g, basis_used))
}

/// MS-style explicit in-plane vectors: `p.u_vec`/`p.v_vec` are integer
/// combinations of the conventional lattice vectors (the `(i j k)` rows
/// of a Materials Studio supercell matrix, e.g. `(0,1,0)` = the b
/// axis). When both are given they override the chosen 2D basis and the
/// U/V expansion. Both must be strictly in-plane w.r.t. the surface
/// normal and non-collinear; otherwise `None`.
fn explicit_inplane_vectors(
    data: &CrystalData,
    p: &SlabParams,
    n_hat: Vec3,
    ref_len: f32,
) -> Option<(Vec3, Vec3)> {
    let u = p.u_vec?;
    let v = p.v_vec?;
    let [va, vb, vc] = data.lattice.to_vectors();
    let make = |c: [i32; 3]| c[0] as f32 * va + c[1] as f32 * vb + c[2] as f32 * vc;
    let ap = make(u);
    let bp = make(v);
    let in_plane_tol = 1e-4 * ref_len;
    if n_hat.dot(ap).abs() > in_plane_tol || n_hat.dot(bp).abs() > in_plane_tol {
        return None; // not in-plane for this (hkl)
    }
    if ap.cross(bp).dot(n_hat).abs() < 1e-9 * ref_len * ref_len {
        return None; // collinear
    }
    Some((ap, bp))
}

/// Full layer-lattice pipeline for `p` (basis choice, orthogonal
/// candidate, U/V factors) — the in-plane pair `build_slab` uses.
/// Returns (ap, bp, basis actually used) — `None` when the (hkl) is
/// degenerate or no slab atoms exist.
fn layer_bravais_pipeline(
    data: &CrystalData,
    p: &SlabParams,
    n_hat: Vec3,
    g: f32,
) -> Option<(Vec3, Vec3, InPlaneBasis)> {
    let [va, vb, vc] = data.lattice.to_vectors();
    // Explicit (i j k) in-plane vectors override basis choice & U/V.
    if p.u_vec.is_some() || p.v_vec.is_some() {
        let ref_len = va.length().max(vb.length()).max(vc.length()).max(1.0);
        return explicit_inplane_vectors(data, p, n_hat, ref_len)
            .map(|(ap, bp)| (ap, bp, InPlaneBasis::Conventional));
    }
    let (ap0, bp0, _) =
        find_inplane_basis(data, p.h, p.k, p.l, n_hat, va, vb, vc, g)?;
    let ref_len = va.length().max(vb.length()).max(vc.length()).max(1.0);

    let (u2, w2, _, _) = find_layer_bravais(data, n_hat, ap0, bp0, ref_len);

    let (ap, bp, basis_used) = match p.basis {
        InPlaneBasis::Conventional => (ap0, bp0, InPlaneBasis::Conventional),
        InPlaneBasis::Orthogonal => {
            let cands = ortho_candidates(&u2, &w2, &n_hat);
            match cands.get(p.orth_idx as usize % cands.len().max(1)) {
                Some((pp, qq)) => (*pp, *qq, InPlaneBasis::Orthogonal),
                None => {
                    eprintln!(
                        "  [slab] no right-angle in-plane cell for this layer lattice — falling back to the primitive basis"
                    );
                    (u2, w2, InPlaneBasis::Primitive)
                }
            }
        }
        InPlaneBasis::Primitive => (u2, w2, InPlaneBasis::Primitive),
    };
    let u = p.u.max(1) as f32;
    let v = p.v.max(1) as f32;
    Some((u * ap, v * bp, basis_used))
}

/// UI-facing layer info (see `SlabLayerInfo`).
pub fn slab_layer_info(data: &CrystalData, p: &SlabParams) -> Option<SlabLayerInfo> {
    let [ia, ib, ic] = data.lattice.inverse_vectors();
    let n_recip = (p.h as f32) * ia + (p.k as f32) * ib + (p.l as f32) * ic;
    let g = n_recip.length();
    if g < 1e-9 {
        return None;
    }
    let n_hat = n_recip / g;
    let d_hkl = 1.0 / g;
    let [va, vb, vc] = data.lattice.to_vectors();
    let (ap0, bp0, _) = find_inplane_basis(data, p.h, p.k, p.l, n_hat, va, vb, vc, g)?;
    let ref_len = va.length().max(vb.length()).max(vc.length()).max(1.0);
    let gidx = gcd(gcd(p.h, p.k), p.l).max(1);
    let span = gidx as f32 / g;
    // One-layer sample: images in [0, 0.49·span)
    let kept = collect_slab_images(data, n_hat, 0.0, 0.49 * span, d_hkl);
    if kept.is_empty() {
        return None;
    }
    let (u2, w2, _, _) = find_layer_bravais(data, n_hat, ap0, bp0, ref_len);
    let a_prim = u2.length();
    let b_prim = w2.length();
    let prod = (a_prim * b_prim).max(1e-30);
    let sin_g = (u2.cross(w2).dot(n_hat).abs() / prod).min(1.0);
    let cos_g = (u2.dot(w2) / prod).clamp(-1.0, 1.0);
    let gamma = sin_g.atan2(cos_g).to_degrees();
    let a_l2 = u2.cross(w2).dot(n_hat).abs();
    let a_b0 = ap0.cross(bp0).dot(n_hat).abs();

    // Distinct in-plane points of the first layer, wrapped into B0.
    let dedup_tol = 1e-3 * ref_len;
    let z_min = kept.iter().map(|k| k.2).fold(f32::INFINITY, f32::min);
    let layer_tol = (span * 0.3).max(1e-3);
    let u = ap0 / ap0.length();
    let v = n_hat.cross(u);
    let ap_u = ap0.dot(u);
    let ap_v = ap0.dot(v);
    let bp_u = bp0.dot(u);
    let bp_v = bp0.dot(v);
    let det0 = ap_u * bp_v - ap_v * bp_u;
    let mut pts: Vec<Vec3> = Vec::new();
    for (_i, pos, zl) in &kept {
        if *zl > z_min + layer_tol {
            continue;
        }
        let r_in = pos - pos.dot(n_hat) * n_hat;
        let fa = (r_in.dot(u) * bp_v - r_in.dot(v) * bp_u) / det0;
        let fb = (ap_u * r_in.dot(v) - r_in.dot(u) * ap_v) / det0;
        let pp = fa.rem_euclid(1.0) * ap0 + fb.rem_euclid(1.0) * bp0;
        if !pts.iter().any(|q| q.distance(pp) < dedup_tol) {
            pts.push(pp);
        }
    }
    // Atoms per layer = distinct in-plane points of the first z-cluster,
    // wrapped into the L2 primitive cell.
    let u2_u = u2.dot(u);
    let u2_v = u2.dot(v);
    let w2_u = w2.dot(u);
    let w2_v = w2.dot(v);
    let d2 = u2_u * w2_v - w2_u * u2_v;
    let mut l2_pts: Vec<Vec3> = Vec::new();
    for p in &pts {
        let pp = if d2.abs() > 1e-12 {
            let x = p.dot(u);
            let y = p.dot(v);
            let alpha = (w2_v * x - w2_u * y) / d2;
            let beta = (u2_u * y - u2_v * x) / d2;
            let snap = |f: f32| {
                let w = f.rem_euclid(1.0);
                if w > 1.0 - 1e-4 {
                    w - 1.0
                } else {
                    w
                }
            };
            snap(alpha) * u2 + snap(beta) * w2
        } else {
            *p
        };
        if !l2_pts.iter().any(|q| q.distance(pp) < dedup_tol) {
            l2_pts.push(pp);
        }
    }
    let mut atoms_per_layer = l2_pts.len() as u32;
    let mut a_prim = a_prim;
    let mut b_prim = b_prim;
    let mut gamma_prim = gamma;
    let _ = a_l2;
    let _ = a_b0;
    // Explicit MS-style (i j k) vectors: report their cell dimensions;
    // atoms per layer scale with the area ratio vs the L2 primitive.
    if p.u_vec.is_some() || p.v_vec.is_some() {
        if let Some((ap, bp)) = explicit_inplane_vectors(data, p, n_hat, ref_len) {
            a_prim = ap.length();
            b_prim = bp.length();
            let prod = (a_prim * b_prim).max(1e-30);
            let sin_g = (ap.cross(bp).dot(n_hat).abs() / prod).min(1.0);
            let cos_g = (ap.dot(bp) / prod).clamp(-1.0, 1.0);
            gamma_prim = sin_g.atan2(cos_g).to_degrees();
            let ratio = ap.cross(bp).dot(n_hat).abs() / a_l2.max(1e-30);
            atoms_per_layer = ((atoms_per_layer as f32 * ratio).round()).max(1.0) as u32;
        }
    }
    let ortho = ortho_candidates(&u2, &w2, &n_hat)
        .into_iter()
        .map(|(pp, qq)| (pp.length(), qq.length(), pp.cross(qq).length()))
        .collect();
    Some(SlabLayerInfo {
        a_prim,
        b_prim,
        gamma_prim,
        atoms_per_layer,
        ortho,
    })
}

/// Every image (under the conventional translation group) of every atom
/// whose normal projection falls in `[s, s+T)`:
/// `(source index, Cartesian position, z_local = proj − s)`.
///
/// Images are enumerated over the full 3D translation group
/// `(m, n, p) ∈ [−M, M]³` — not just the 1-D Bézout w-family — so that
/// in-plane cells finer than any single w-strip (e.g. the FCC (111)
/// hexagonal L2 cell) still collect the complete slab content.
fn collect_slab_images(
    data: &CrystalData,
    n_hat: Vec3,
    s: f32,
    thickness: f32,
    d_hkl: f32,
) -> Vec<(usize, Vec3, f32)> {
    let [va, vb, vc] = data.lattice.to_vectors();
    let inv = data.lattice.inverse_vectors();
    let na = va.dot(n_hat);
    let nb = vb.dot(n_hat);
    let nc = vc.dot(n_hat);
    let step = (na.abs() + nb.abs() + nc.abs()) / 3.0;
    let m_min = 2 + (thickness / d_hkl.max(1e-9)).ceil() as i32;
    let m_max = if step > 1e-9 {
        ((thickness / step).ceil() as i32 + 1).max(0)
    } else {
        0
    };
    let m = m_min.max(m_max);
    let tol = (1e-4_f32 * d_hkl.max(1e-6)).max(1e-6_f32);
    let mut kept: Vec<(usize, Vec3, f32)> = Vec::new();
    for i in 0..data.atoms.len() {
        let f = if data.positions_fractional {
            Vec3::new(data.atoms[i].x, data.atoms[i].y, data.atoms[i].z)
        } else {
            apply_inv3(&inv, Vec3::new(data.atoms[i].x, data.atoms[i].y, data.atoms[i].z))
        };
        let p0 = (f.x * va + f.y * vb + f.z * vc).dot(n_hat);
        for mm in -m..=m {
            for nn in -m..=m {
                for pp in -m..=m {
                    let z = p0
                        + (mm as f32) * na
                        + (nn as f32) * nb
                        + (pp as f32) * nc;
                    if z > s - tol && z < s + thickness - tol {
                        let pos = f.x * va
                            + f.y * vb
                            + f.z * vc
                            + (mm as f32) * va
                            + (nn as f32) * vb
                            + (pp as f32) * vc;
                        kept.push((i, pos, z - s));
                    }
                }
            }
        }
    }
    kept
}

/// Build the slab cross-section of `data` along the (h,k,l) planes:
/// keeps atoms whose normal projection falls in `[s, s+T)`, with the slab
/// thickness T in Angstrom (0 → 3·d_hkl) and start position s in
/// Angstrom. The in-plane cell is the U×V supercell of the chosen 2D
/// basis (`p.basis`: primitive / orthogonal / conventional, see
/// `InPlaneBasis`), c' = T·n̂, atoms re-expressed in the new
/// standard-setting frame.
pub fn build_slab(data: &CrystalData, p: &SlabParams) -> Result<CrystalData, String> {
    if p.h == 0 && p.k == 0 && p.l == 0 {
        return Err("Miller indices are all zero".into());
    }
    let t = if p.thickness_ang > 0.0 { p.thickness_ang } else { 0.0 };
    if data.atoms.is_empty() {
        return Err("No atoms in structure".into());
    }

    let [va, vb, vc] = data.lattice.to_vectors();
    let [ia, ib, ic] = data.lattice.inverse_vectors();
    let n_recip = (p.h as f32) * ia + (p.k as f32) * ib + (p.l as f32) * ic;
    let g = n_recip.length();
    if g < 1e-9 {
        return Err("Invalid Miller indices for this lattice".into());
    }
    let n_hat = n_recip / g;
    let d_hkl = 1.0 / g;
    let thickness = if t > 0.0 { t } else { 3.0 * d_hkl };
    let s = p.start_ang;

    // ── Candidate images: full conventional translation group over
    // [s, s+T) — (m,n,p) cell translations, not just the 1-D Bézout
    // w-family, so in-plane cells finer than any single w-strip
    // (e.g. the FCC (111) L2 hexagonal cell) still collect complete
    // content. ──
    let kept = collect_slab_images(data, n_hat, s, thickness, d_hkl);
    if kept.is_empty() {
        return Err("No atoms fall inside the slab region [s, s+T) — check position/thickness".into());
    }

    // ── In-plane cell: chosen 2D basis × U/V ──
    let ref_len = va.length().max(vb.length()).max(vc.length()).max(1.0);
    let (ap0, bp0, _) = find_inplane_basis(data, p.h, p.k, p.l, n_hat, va, vb, vc, g)
        .ok_or("Degenerate in-plane lattice for this (hkl) — try a different Miller index or cell setting")?;

    let (ap, bp, basis_used) = if p.u_vec.is_some() || p.v_vec.is_some() {
        // Explicit MS-style (i j k) vectors override basis choice & U/V.
        let (ap, bp) = explicit_inplane_vectors(data, p, n_hat, ref_len).ok_or_else(|| {
            format!(
                "U/V in-plane vectors ({:?} / {:?}) must be in-plane for ({},{},{}) and non-collinear — e.g. for (100): U=(0 1 0), V=(0 0 1)",
                p.u_vec.unwrap_or([0; 3]),
                p.v_vec.unwrap_or([0; 3]),
                p.h, p.k, p.l
            )
        })?;
        (ap, bp, InPlaneBasis::Conventional)
    } else {
        match p.basis {
            InPlaneBasis::Conventional => (ap0, bp0, InPlaneBasis::Conventional),
            InPlaneBasis::Orthogonal | InPlaneBasis::Primitive => {
                let (u2, w2, _, _) = find_layer_bravais(data, n_hat, ap0, bp0, ref_len);
                match p.basis {
                    InPlaneBasis::Orthogonal => {
                        let cands = ortho_candidates(&u2, &w2, &n_hat);
                        match cands.get(p.orth_idx as usize % cands.len().max(1)) {
                            Some((pp, qq)) => (*pp, *qq, InPlaneBasis::Orthogonal),
                            None => {
                                eprintln!(
                                    "  [slab] no right-angle in-plane cell for this layer lattice — falling back to the primitive basis"
                                );
                                (u2, w2, InPlaneBasis::Primitive)
                            }
                        }
                    }
                    _ => (u2, w2, InPlaneBasis::Primitive),
                }
            }
        }
    };
    let uf = if p.u_vec.is_some() || p.v_vec.is_some() {
        1.0 // explicit vectors already carry the full in-plane size
    } else {
        p.u.max(1) as f32
    };
    let vf = if p.u_vec.is_some() || p.v_vec.is_some() {
        1.0
    } else {
        p.v.max(1) as f32
    };
    let ap = uf * ap;
    let bp = vf * bp;

    // ── Place atoms in the slab cell: wrap in-plane, keep z_local ──
    // Orthonormal in-plane basis: u = ap̂, v = n̂×u (then u×v = n̂)
    let u = ap / ap.length();
    let v = n_hat.cross(u);
    let ap_u = ap.dot(u);
    let ap_v = ap.dot(v);
    let bp_u = bp.dot(u);
    let bp_v = bp.dot(v);
    let det2 = ap_u * bp_v - ap_v * bp_u; // > 0 by construction

    let cp = thickness * n_hat;
    let l_new = Lattice::from_cartesian_vectors(ap, bp, cp);
    let std_basis = l_new.to_vectors();

    let mut atoms: Vec<AtomData> = Vec::with_capacity(kept.len());
    for (src_idx, pos, z_local) in &kept {
        let r_in = pos - z_local * n_hat;
        let x = r_in.dot(u);
        let y = r_in.dot(v);
        // Cramer: [ap_u bp_u; ap_v bp_v]·[fa;fb] = [x;y], D = det2
        let fa = (x * bp_v - y * bp_u) / det2;
        let fb = (ap_u * y - x * ap_v) / det2;
        // Wrap to [0, 1), snapping f32 noise near 1.0 back to 0 so
        // boundary images of the same point do not survive as duplicates.
        let wrap = |f: f32| {
            let w = f.rem_euclid(1.0);
            if w > 1.0 - 1e-4 {
                w - 1.0
            } else {
                w
            }
        };
        let fa = wrap(fa);
        let fb = wrap(fb);
        let out = fa * std_basis[0] + fb * std_basis[1] + (z_local / thickness) * std_basis[2];
        let at = &data.atoms[*src_idx];
        let new_at = AtomData {
            element: at.element.clone(),
            x: out.x,
            y: out.y,
            z: out.z,
            label: at.label.clone(),
        };
        // Bézout q-images can map several slab copies onto the same point
        // of the in-plane cell — keep only one.
        let dup = atoms.iter().any(|o| {
            o.element == new_at.element
                && (o.x - new_at.x).abs() < 5e-3
                && (o.y - new_at.y).abs() < 5e-3
                && (o.z - new_at.z).abs() < 5e-3
        });
        if !dup {
            atoms.push(new_at);
        }
    }

    Ok(CrystalData {
        lattice: l_new,
        atoms,
        positions_fractional: false,
        modified: true,
        phonon_modes: None,
        slab: Some(SlabMeta {
            h: p.h,
            k: p.k,
            l: p.l,
            start_ang: s,
            thickness_ang: thickness,
            u: p.u,
            v: p.v,
            basis: basis_used,
            u_vec: p.u_vec,
            v_vec: p.v_vec,
        }),
        vacuum: None,
        supercell: None,
    })
}

/// Build a vacuum layer of `thickness_ang` Angstrom by extending the chosen
/// cell axis (1=a, 2=b, 3=c) along its own direction: the new cell
/// thickness is `c_old + V`. Atom positions are placed according to
/// `vacuum_position(p)`: 0 = whole vacuum on top, 0.5 = split V/2 on each
/// side, 1 = whole vacuum below. `V = 0` is a no-op thickness-wise.
pub fn build_vacuum(data: &CrystalData, p: &VacuumParams) -> Result<CrystalData, String> {
    if p.axis < 1 || p.axis > 3 {
        return Err("Vacuum axis must be 1 (a), 2 (b) or 3 (c)".into());
    }
    if !(p.thickness_ang >= 0.0) || p.thickness_ang.is_nan() {
        return Err("Vacuum thickness must be non-negative".into());
    }
    if data.atoms.is_empty() {
        return Err("No atoms in structure".into());
    }

    let idx = (p.axis - 1) as usize;
    let [ca, cb, cc] = data.lattice.to_vectors();
    let mut axes = [ca, cb, cc];
    let axis_len = axes[idx].length();
    if axis_len < 1e-9 {
        return Err("Degenerate cell axis".into());
    }
    let a_hat = axes[idx] / axis_len;
    axes[idx] = axes[idx] + p.thickness_ang * a_hat;

    let inv_new = inv3(axes).ok_or("Degenerate new lattice")?;
    let l_new = Lattice::from_cartesian_vectors(axes[0], axes[1], axes[2]);
    let std_basis = l_new.to_vectors(); // standard-setting frame for output

    let shift = vacuum_position(p) * p.thickness_ang * a_hat;

    let mut atoms = Vec::with_capacity(data.atoms.len());
    for at in &data.atoms {
        // Atom position in the current standard-setting Cartesian frame
        let r = if data.positions_fractional {
            at.x * ca + at.y * cb + at.z * cc
        } else {
            Vec3::new(at.x, at.y, at.z)
        };
        let f_new = apply_inv3(&inv_new, r + shift);
        let out = f_new.x * std_basis[0] + f_new.y * std_basis[1] + f_new.z * std_basis[2];
        atoms.push(AtomData {
            element: at.element.clone(),
            x: out.x,
            y: out.y,
            z: out.z,
            label: at.label.clone(),
        });
    }

    Ok(CrystalData {
        lattice: l_new,
        atoms,
        positions_fractional: false,
        modified: true,
        phonon_modes: None,
        slab: data.slab,
        vacuum: Some(VacuumMeta {
            axis: p.axis,
            thickness_ang: p.thickness_ang,
            both_sides: p.both_sides,
            position: p.position,
        }),
        supercell: data.supercell,
    })
}

/// Supercell expansion parameters: integer multipliers along a, b, c
/// (1×1×1 = the current cell itself; 2×1×1 = one extra cell along a).
#[derive(Debug, Clone, Copy, Default)]
pub struct SupercellParams {
    pub x: i32,
    pub y: i32,
    pub z: i32,
}

/// Expand the structure into a [x y z] supercell: the lattice vectors are
/// multiplied by (x, y, z) and every atom is replicated into the
/// x·y·z copies at the integer translations (i, j, k) with i < x,
/// j < y, k < z, re-expressed in the enlarged cell. The supercell
/// provenance is recorded on the result (slab/vacuum provenance is kept
/// — a c-superlaced slab+vacuum stack is the intended catalysis use).
/// Multipliers < 1 are clamped to 1 (a no-op axis).
pub fn build_supercell(data: &CrystalData, p: &SupercellParams) -> Result<CrystalData, String> {
    let x = p.x.max(1);
    let y = p.y.max(1);
    let z = p.z.max(1);
    if data.atoms.is_empty() {
        return Err("No atoms in structure".into());
    }

    let [ca, cb, cc] = data.lattice.to_vectors();
    let a_new = ca * x as f32;
    let b_new = cb * y as f32;
    let c_new = cc * z as f32;
    let inv_new = inv3([a_new, b_new, c_new]).ok_or("Degenerate new lattice")?;
    let l_new = Lattice::from_cartesian_vectors(a_new, b_new, c_new);
    let basis = l_new.to_vectors(); // standard-setting frame for output

    let mut atoms = Vec::with_capacity(data.atoms.len() * x as usize * y as usize * z as usize);
    for at in &data.atoms {
        // Atom position in the current standard-setting Cartesian frame
        let r = if data.positions_fractional {
            at.x * ca + at.y * cb + at.z * cc
        } else {
            Vec3::new(at.x, at.y, at.z)
        };
        for i in 0..x {
            for j in 0..y {
                for k in 0..z {
                    let pos = r + ca * i as f32 + cb * j as f32 + cc * k as f32;
                    let f_new = apply_inv3(&inv_new, pos);
                    let out = f_new.x * basis[0] + f_new.y * basis[1] + f_new.z * basis[2];
                    atoms.push(AtomData {
                        element: at.element.clone(),
                        x: out.x,
                        y: out.y,
                        z: out.z,
                        label: at.label.clone(),
                    });
                }
            }
        }
    }

    Ok(CrystalData {
        lattice: l_new,
        atoms,
        positions_fractional: false,
        modified: true,
        phonon_modes: None,
        slab: data.slab,
        vacuum: data.vacuum,
        supercell: Some([x, y, z]),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cubic(a: f32) -> CrystalData {
        CrystalData {
            lattice: Lattice { a, b: a, c: a, alpha: 90.0, beta: 90.0, gamma: 90.0 },
            atoms: vec![AtomData {
                element: "Cu".into(),
                x: 0.0,
                y: 0.0,
                z: 0.0,
                label: String::new(),
            }],
            positions_fractional: false,
            modified: false,
            phonon_modes: None,
            slab: None,
            vacuum: None,
            supercell: None,
        }
    }

    fn fcc(a: f32) -> CrystalData {
        CrystalData {
            lattice: Lattice { a, b: a, c: a, alpha: 90.0, beta: 90.0, gamma: 90.0 },
            atoms: vec![
                AtomData { element: "Cu".into(), x: 0.0, y: 0.0, z: 0.0, label: String::new() },
                AtomData { element: "Cu".into(), x: 0.0, y: a / 2.0, z: a / 2.0, label: String::new() },
                AtomData { element: "Cu".into(), x: a / 2.0, y: 0.0, z: a / 2.0, label: String::new() },
                AtomData { element: "Cu".into(), x: a / 2.0, y: a / 2.0, z: 0.0, label: String::new() },
            ],
            positions_fractional: false,
            modified: false,
            phonon_modes: None,
            slab: None,
            vacuum: None,
            supercell: None,
        }
    }

    #[test]
    fn d_hkl_cubic() {
        let d = cubic(4.0);
        assert!((slab_d_hkl(&d, 1, 0, 0).unwrap() - 4.0).abs() < 1e-4);
        assert!(slab_d_hkl(&d, 0, 0, 0).is_none());
    }

    #[test]
    fn slab_period_cubic() {
        let d = cubic(3.6);
        assert!((slab_period(&d, 0, 0, 1).unwrap() - 3.6).abs() < 1e-4);
        let f = fcc(3.6);
        // (210): gidx = 1 → period = d_210
        let p = slab_period(&f, 2, 1, 0).unwrap();
        assert!((p - slab_d_hkl(&f, 2, 1, 0).unwrap()).abs() < 1e-4);
    }

    #[test]
    fn vacuum_along_c_keeps_atoms_and_extends_cell() {
        let mut d = cubic(3.6);
        d.atoms.push(AtomData {
            element: "Cu".into(),
            x: 1.8,
            y: 1.8,
            z: 3.0,
            label: String::new(),
        });
        let out = build_vacuum(&d, &VacuumParams { axis: 3, thickness_ang: 15.0, both_sides: false, position: 0.0 }).unwrap();
        assert!((out.lattice.c - 18.6).abs() < 1e-3, "c: {}", out.lattice.c);
        // Top-only: atom Cartesian z is unchanged
        let z0 = out.atoms[1].z;
        assert!((z0 - 3.0).abs() < 1e-3, "z: {}", z0);
        assert!(out.vacuum.is_some());
        assert!(out.slab.is_none());
    }

    #[test]
    fn vacuum_both_sides_shifts_atoms_by_half() {
        let mut d = cubic(3.6);
        d.atoms.push(AtomData {
            element: "Cu".into(),
            x: 1.8,
            y: 1.8,
            z: 3.0,
            label: String::new(),
        });
        let out = build_vacuum(&d, &VacuumParams { axis: 3, thickness_ang: 15.0, both_sides: true, position: 0.0 }).unwrap();
        let z0 = out.atoms[1].z;
        // Atom moves up by V/2 = 7.5 in Cartesian (top-only z was 3.0 → 10.5)
        assert!((z0 - 10.5).abs() < 1e-3, "z: {}", z0);
    }

    #[test]
    fn vacuum_position_top_and_bottom() {
        let mut d = cubic(3.6);
        d.atoms.push(AtomData {
            element: "Cu".into(),
            x: 0.0,
            y: 0.0,
            z: 3.0,
            label: String::new(),
        });
        // position = 1.0: whole vacuum below the structure
        let out = build_vacuum(&d, &VacuumParams { axis: 3, thickness_ang: 15.0, both_sides: false, position: 1.0 }).unwrap();
        assert!((out.lattice.c - 18.6).abs() < 1e-3);
        assert!((out.atoms[1].z - 18.0).abs() < 1e-3, "z: {}", out.atoms[1].z);
        // position = 0.5 (explicit) == both_sides
        let out = build_vacuum(&d, &VacuumParams { axis: 3, thickness_ang: 15.0, both_sides: false, position: 0.5 }).unwrap();
        assert!((out.atoms[1].z - 10.5).abs() < 1e-3, "z: {}", out.atoms[1].z);
    }

    #[test]
    fn vacuum_zero_thickness_is_noop_thickness() {
        let mut d = cubic(3.6);
        d.atoms.push(AtomData {
            element: "Cu".into(),
            x: 0.0,
            y: 0.0,
            z: 1.8,
            label: String::new(),
        });
        let out = build_vacuum(&d, &VacuumParams { axis: 3, thickness_ang: 0.0, both_sides: false, position: 0.0 }).unwrap();
        assert!((out.lattice.c - 3.6).abs() < 1e-6, "c: {}", out.lattice.c);
        assert!((out.atoms[1].z - 1.8).abs() < 1e-3, "z: {}", out.atoms[1].z);
    }

    #[test]
    fn vacuum_rejects_bad_input() {
        let d = cubic(3.6);
        assert!(build_vacuum(&d, &VacuumParams { axis: 0, thickness_ang: 5.0, both_sides: false, position: 0.0 }).is_err());
        assert!(build_vacuum(&d, &VacuumParams { axis: 3, thickness_ang: -1.0, both_sides: false, position: 0.0 }).is_err());
    }

    #[test]
    fn slab_001_cubic_three_layers() {
        // One atom at the origin; (001) cut, 3 layers of c=3.6
        let d = cubic(3.6);
        let out = build_slab(&d, &SlabParams {
            h: 0, k: 0, l: 1,
            start_ang: 0.0,
            thickness_ang: 3.0 * 3.6, // explicit thickness
            ..Default::default()
        })
        .unwrap();
        assert!((out.lattice.c - 10.8).abs() < 1e-3, "c: {}", out.lattice.c);
        assert!((out.lattice.a - 3.6).abs() < 1e-3);
        assert!((out.lattice.b - 3.6).abs() < 1e-3);
        // Images at z = 0, 3.6, 7.2 → 3 atoms
        assert_eq!(out.atoms.len(), 3, "atom count");
        let zs: Vec<f32> = out.atoms.iter().map(|a| a.z).collect();
        for z in &zs {
            assert!(*z < 10.8 + 1e-3, "atom z {} in vacuum region", z);
        }
        assert!(out.slab.is_some());
        assert!((out.slab.unwrap().thickness_ang - 10.8).abs() < 1e-4);
    }

    #[test]
    fn slab_001_offset() {
        // s = 1.8: slab [1.8, 12.6) → images at 3.6, 7.2, 10.8
        let d = cubic(3.6);
        let out = build_slab(&d, &SlabParams {
            h: 0, k: 0, l: 1,
            start_ang: 1.8,
            thickness_ang: 10.8,
            ..Default::default()
        })
        .unwrap();
        assert_eq!(out.atoms.len(), 3);
        let zs: Vec<f32> = out.atoms.iter().map(|a| a.z).collect();
        for z in &zs {
            assert!((*z - 1.8) > -1e-3 && (*z - 1.8) < 10.8 + 1e-3, "z {} outside slab", z);
        }
    }

    #[test]
    fn slab_210_mixed_indices_no_missing_layers() {
        // Regression: stepping images along a single cell vector (the one
        // with the largest normal span, span = 2·d_210 here) skips every
        // other (210) layer. The Bézout stepping must yield all 3 layers.
        let a = 3.6;
        let d = cubic(a);
        let d210 = slab_d_hkl(&d, 2, 1, 0).unwrap();
        assert!((d210 - a / 5.0_f32.sqrt()).abs() < 1e-3, "d210: {}", d210);
        let out = build_slab(&d, &SlabParams {
            h: 2, k: 1, l: 0,
            start_ang: 0.0,
            thickness_ang: 3.0 * d210,
            ..Default::default()
        })
        .unwrap();
        assert_eq!(out.atoms.len(), 3, "all 3 layers of the (210) slab");
        assert!((out.lattice.c - 3.0 * d210).abs() < 1e-3, "c: {}", out.lattice.c);
        // Layer z-local positions: 0, d, 2d (coordinate along the new c axis)
        let c_std = out.lattice.to_vectors()[2];
        let c_hat = c_std / c_std.length();
        let mut zl: Vec<f32> = out
            .atoms
            .iter()
            .map(|at| Vec3::new(at.x, at.y, at.z).dot(c_hat))
            .collect();
        zl.sort_by(|x, y| x.partial_cmp(y).unwrap());
        assert!((zl[0]).abs() < 1e-3 && (zl[1] - d210).abs() < 1e-3 && (zl[2] - 2.0 * d210).abs() < 1e-3,
            "layer z-locals: {:?}", zl);
    }

    #[test]
    fn slab_211_mixed_indices_no_missing_layers() {
        // (211) on cubic: g = 1, Bézout vector w with G·w = 1.
        let a = 3.6;
        let d = cubic(a);
        let d211 = slab_d_hkl(&d, 2, 1, 1).unwrap();
        let out = build_slab(&d, &SlabParams {
            h: 2, k: 1, l: 1,
            start_ang: 0.0,
            thickness_ang: 3.0 * d211,
            ..Default::default()
        })
        .unwrap();
        assert_eq!(out.atoms.len(), 3, "all 3 layers of the (211) slab");
        assert!((out.lattice.c - 3.0 * d211).abs() < 1e-3, "c: {}", out.lattice.c);
    }

    #[test]
    fn slab_001_fcc_primitive_square() {
        // fcc conventional cell, (001): true 2D Bravais lattice of the
        // layer is the square a/√2 cell (1 atom/layer), finer than the
        // conventional a×b in-plane cell (2 atoms/layer).
        let a = 3.6;
        let d = fcc(a);
        let out = build_slab(&d, &SlabParams {
            h: 0, k: 0, l: 1,
            start_ang: 0.0,
            thickness_ang: 0.0, // → 3·d_001 = 3a
            ..Default::default()
        })
        .unwrap();
        // d_001 = a (conventional reciprocal), layers at 0, a/2, a, ...
        // 6 physical layers × 1 atom = 6 atoms
        let c_exp = 3.0 * a;
        assert!((out.lattice.c - c_exp).abs() < 1e-3, "c: {}", out.lattice.c);
        assert_eq!(out.atoms.len(), 6, "fcc (001) primitive: 6 layers × 1 atom");
        // In-plane cell: a/√2 square, right angles
        let s2 = a / 2.0_f32.sqrt();
        assert!((out.lattice.a - s2).abs() < 1e-2, "a': {}", out.lattice.a);
        assert!((out.lattice.b - s2).abs() < 1e-2, "b': {}", out.lattice.b);
        assert!((out.lattice.alpha - 90.0).abs() < 1e-2);
        assert!((out.lattice.gamma - 90.0).abs() < 1e-2, "gamma: {}", out.lattice.gamma);
    }

    #[test]
    fn slab_111_fcc_primitive_hex() {
        // fcc (111): true 2D Bravais lattice is hexagonal with
        // a_2d = a/√2 (1 atom per L2 cell; the B point of the FCC
        // layer is itself a Bravais point, ≡ 0 mod L2), versus the
        // conventional √3a² in-plane cell (4 atoms per plane).
        let a = 3.6;
        let d = fcc(a);
        let d_hkl = slab_d_hkl(&d, 1, 1, 1).unwrap();
        let out = build_slab(&d, &SlabParams {
            h: 1, k: 1, l: 1,
            start_ang: 0.0,
            thickness_ang: 3.0 * d_hkl,
            ..Default::default()
        })
        .unwrap();
        assert!((out.lattice.c - 3.0 * d_hkl).abs() < 1e-3, "c: {}", out.lattice.c);
        assert_eq!(out.atoms.len(), 3, "fcc (111) primitive: 3 planes × 1 atom (B point ≡ 0 mod L2)");
        let s2 = a / 2.0_f32.sqrt();
        assert!((out.lattice.a - s2).abs() < 1e-2, "a': {}", out.lattice.a);
        assert!((out.lattice.b - s2).abs() < 1e-2, "b': {}", out.lattice.b);
        // Hexagonal: γ = 60° (or 120°, Gauss-reduced angle ∈ (0°,90°])
        let g = out.lattice.gamma;
        assert!((g - 60.0).abs() < 1.0 || (g - 120.0).abs() < 1.0, "gamma: {}", g);
        let m = out.slab.unwrap();
        assert_eq!(m.basis, InPlaneBasis::Primitive);
        assert_eq!(m.u, 1);
        assert_eq!(m.v, 1);
    }

    #[test]
    fn slab_111_fcc_orthogonal_90deg() {
        // fcc (111) with the right-angle in-plane cell: a'×b' =
        // (a/√2) × (√3·a/√2), 90°/90°/90° — twice the L2 primitive
        // area: 2 L2 cells per z-plane.
        let a = 3.6;
        let d = fcc(a);
        let d_hkl = slab_d_hkl(&d, 1, 1, 1).unwrap();
        let out = build_slab(&d, &SlabParams {
            h: 1, k: 1, l: 1,
            start_ang: 0.0,
            thickness_ang: 3.0 * d_hkl,
            basis: InPlaneBasis::Orthogonal,
            ..Default::default()
        })
        .unwrap();
        assert!((out.lattice.c - 3.0 * d_hkl).abs() < 1e-3, "c: {}", out.lattice.c);
        assert_eq!(out.atoms.len(), 6, "fcc (111) orthogonal: 3 planes × 2 L2 cells");
        let s2 = a / 2.0_f32.sqrt();
        let (a1, a2) = (out.lattice.a, out.lattice.b);
        let lo = a1.min(a2);
        let hi = a1.max(a2);
        assert!((lo - s2).abs() < 1e-2, "short side: {}", lo);
        assert!((hi - (3.0_f32.sqrt() * s2)).abs() < 1e-2, "long side: {}", hi);
        assert!((out.lattice.alpha - 90.0).abs() < 1e-2, "alpha: {}", out.lattice.alpha);
        assert!((out.lattice.beta - 90.0).abs() < 1e-2, "beta: {}", out.lattice.beta);
        assert!((out.lattice.gamma - 90.0).abs() < 1e-2, "gamma: {}", out.lattice.gamma);
        assert_eq!(out.slab.unwrap().basis, InPlaneBasis::Orthogonal);
    }

    #[test]
    fn slab_111_fcc_conventional_supercell() {
        // Legacy conventional in-plane basis: √3a² cell, 4 atoms/layer.
        let a = 3.6;
        let d = fcc(a);
        let d_hkl = slab_d_hkl(&d, 1, 1, 1).unwrap();
        let out = build_slab(&d, &SlabParams {
            h: 1, k: 1, l: 1,
            start_ang: 0.0,
            thickness_ang: 3.0 * d_hkl,
            basis: InPlaneBasis::Conventional,
            ..Default::default()
        })
        .unwrap();
        assert_eq!(out.atoms.len(), 12, "fcc (111) conventional: 3 layers × 4 atoms");
        assert!((out.lattice.c - 3.0 * d_hkl).abs() < 1e-3);
        let a_new = out.lattice.a;
        assert!((a_new - 2.0_f32.sqrt() * a).abs() < 0.05, "a': {}", a_new);

        // Regression: the in-plane 2x2 Cramer solve must place the FCC
        // layer's 4 atoms at the conventional-cell in-plane positions
        // (fa,fb) = (0,0), (1/2,0), (0,1/2), (1/2,1/2) — wrong Cramer
        // terms scramble the in-plane positions for non-orthogonal (ap,bp).
        let b_new = out.lattice.b;
        let g_new = out.lattice.gamma.to_radians();
        let bottom: Vec<(f32, f32)> = out
            .atoms
            .iter()
            .filter(|at| at.z.abs() < 1e-3)
            .map(|at| (at.x, at.y))
            .collect();
        assert_eq!(bottom.len(), 4, "bottom layer atom count");
        for (fx, fy) in [(0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5)] {
            let ex = fx * a_new + fy * b_new * g_new.cos();
            let ey = fy * b_new * g_new.sin();
            assert!(
                bottom.iter().any(|(x, y)| (x - ex).abs() < 1e-2 && (y - ey).abs() < 1e-2),
                "missing in-plane position (fa,fb)=({}, {}) -> cartesian ({}, {}); got {:?}",
                fx, fy, ex, ey, bottom
            );
        }
    }

    #[test]
    fn slab_supercell_uv_expansion() {
        // U×V expansion of the (001) primitive in-plane cell.
        let a = 3.6;
        let d = cubic(a);
        let out = build_slab(&d, &SlabParams {
            h: 0, k: 0, l: 1,
            start_ang: 0.0,
            thickness_ang: 3.0 * a,
            u: 2,
            v: 3,
            ..Default::default()
        })
        .unwrap();
        assert!((out.lattice.a - 2.0 * a).abs() < 1e-3, "a: {}", out.lattice.a);
        assert!((out.lattice.b - 3.0 * a).abs() < 1e-3, "b: {}", out.lattice.b);
        assert_eq!(out.atoms.len(), 18, "2×3 in-plane images × 3 layers");
        assert_eq!(out.slab.unwrap().u, 2);
        assert_eq!(out.slab.unwrap().v, 3);
    }

    #[test]
    fn slab_layer_snap_cubic() {
        let d = cubic(3.6); // (001) layers at z = 0, 3.6, 7.2, ...
        assert!((slab_layer_snap(&d, 0, 0, 1, 1.0, true).unwrap() - 3.6).abs() < 1e-4);
        assert!((slab_layer_snap(&d, 0, 0, 1, 1.0, false).unwrap() - 0.0).abs() < 1e-4);
        assert!((slab_layer_snap(&d, 0, 0, 1, 5.0, true).unwrap() - 7.2).abs() < 1e-4);
        assert!((slab_layer_snap(&d, 0, 0, 1, 5.0, false).unwrap() - 3.6).abs() < 1e-4);
    }

    #[test]
    fn slab_layer_snap_fcc_half_layers() {
        // fcc (001): layers at z = 0, a/2, a, ... (span = a, atoms of the
        // two species interleave at a/2)
        let a = 3.6;
        let d = fcc(a);
        assert!((slab_layer_snap(&d, 0, 0, 1, 1.0, true).unwrap() - a / 2.0).abs() < 1e-4,
            "snap up: {}", slab_layer_snap(&d, 0, 0, 1, 1.0, true).unwrap());
        assert!((slab_layer_snap(&d, 0, 0, 1, 2.5, true).unwrap() - a).abs() < 1e-4);
        assert!((slab_layer_snap(&d, 0, 0, 1, 2.5, false).unwrap() - a / 2.0).abs() < 1e-4);
    }

    #[test]
    fn ortho_candidates_fcc111() {
        // fcc (111) L2 = hexagonal s × s (60°): the right-angle pair is
        // u × (u − 2w) = s × √3·s.
        let a = 3.6;
        let d = fcc(a);
        let info = slab_layer_info(&d, &SlabParams {
            h: 1, k: 1, l: 1,
            ..Default::default()
        })
        .unwrap();
        let s2 = a / 2.0_f32.sqrt();
        assert!((info.a_prim - s2).abs() < 1e-2, "a_prim: {}", info.a_prim);
        assert_eq!(info.atoms_per_layer, 1, "fcc (111): 1 atom per hexagonal cell");
        assert!(!info.ortho.is_empty(), "expected an orthogonal candidate");
        let (pl, ql, area) = info.ortho[0];
        let (lo, hi) = (pl.min(ql), pl.max(ql));
        assert!((lo - s2).abs() < 1e-2 && (hi - 3.0_f32.sqrt() * s2).abs() < 1e-2,
            "orthogonal pair: {} x {}", lo, hi);
        assert!((area - 3.0_f32.sqrt() * s2 * s2).abs() < 1e-2, "area: {}", area);
    }

    #[test]
    fn layer_info_cubic_001() {
        let a = 3.6;
        let d = cubic(a);
        let info = slab_layer_info(&d, &SlabParams { h: 0, k: 0, l: 1, ..Default::default() }).unwrap();
        assert!((info.a_prim - a).abs() < 1e-3, "a_prim: {}", info.a_prim);
        assert_eq!(info.atoms_per_layer, 1);
        // Square lattice: u × w already orthogonal → candidate available
        assert!(!info.ortho.is_empty());
    }

    #[test]
    fn slab_0001_hexagonal() {
        // Hexagonal: a=b=3.2, c=5.2, gamma=120°; 2-atom basis, fractional coords
        let d = CrystalData {
            lattice: Lattice { a: 3.2, b: 3.2, c: 5.2, alpha: 90.0, beta: 90.0, gamma: 120.0 },
            atoms: vec![
                AtomData { element: "Mg".into(), x: 0.0, y: 0.0, z: 0.0, label: String::new() },
                AtomData { element: "Mg".into(), x: 2.0 / 3.0, y: 1.0 / 3.0, z: 2.0 / 3.0, label: String::new() },
            ],
            positions_fractional: true,
            modified: false,
            phonon_modes: None,
            slab: None,
            vacuum: None,
            supercell: None,
        };
        let out = build_slab(&d, &SlabParams {
            h: 0, k: 0, l: 1,
            start_ang: 0.0,
            thickness_ang: 3.0 * 5.2,
            ..Default::default()
        })
        .unwrap();
        assert!((out.lattice.c - 15.6).abs() < 1e-2, "c: {}", out.lattice.c);
        // 6 layers (c/3 spacing × 3c) × 1 atom = 6
        assert_eq!(out.atoms.len(), 6, "hex (0001) slab atom count");
    }

    #[test]
    fn slab_zero_indices_rejected() {
        let d = cubic(3.6);
        let e = build_slab(&d, &SlabParams { h: 0, k: 0, l: 0, start_ang: 0.0, thickness_ang: 5.0, ..Default::default() }).unwrap_err();
        assert!(e.contains("all zero"), "{}", e);
    }

    #[test]
    fn slab_100_explicit_uv_vectors() {
        let a = 3.615;
        let d = fcc(a);
        // MS-style rows: U = 2b, V = c (in-plane for (100)).
        let p = SlabParams {
            h: 1, k: 0, l: 0,
            u_vec: Some([0, 2, 0]),
            v_vec: Some([0, 0, 1]),
            ..Default::default()
        };
        let out = build_slab(&d, &p).unwrap();
        let l = out.lattice;
        let mut ab = [l.a, l.b];
        ab.sort_by(|x, y| x.partial_cmp(y).unwrap());
        assert!((ab[0] - a).abs() < 0.01, "short side {} != a", ab[0]);
        assert!((ab[1] - 2.0 * a).abs() < 0.01, "long side {} != 2a", ab[1]);
        assert!((l.c - 3.0 * a).abs() < 0.02, "c' {} != 3a", l.c);
        // 6 (100) layers × 4 atoms per 2b×c cell
        assert_eq!(out.atoms.len(), 24, "{}", out.atoms.len());
        let meta = out.slab.as_ref().unwrap();
        assert_eq!(meta.u_vec, Some([0, 2, 0]));
        assert_eq!(meta.v_vec, Some([0, 0, 1]));
    }

    #[test]
    fn slab_explicit_uv_rejects_out_of_plane() {
        let d = fcc(3.615);
        let p = SlabParams {
            h: 1, k: 0, l: 0,
            u_vec: Some([1, 0, 0]), // a axis is the normal, not in-plane
            v_vec: Some([0, 0, 1]),
            ..Default::default()
        };
        assert!(build_slab(&d, &p).is_err());
    }

    #[test]
    fn slab_explicit_uv_rejects_collinear() {
        let d = fcc(3.615);
        let p = SlabParams {
            h: 1, k: 0, l: 0,
            u_vec: Some([0, 1, 0]),
            v_vec: Some([0, 2, 0]),
            ..Default::default()
        };
        assert!(build_slab(&d, &p).is_err());
    }

    #[test]
    fn inplane_explicit_uv_override() {
        let a = 3.615;
        let d = fcc(a);
        // (010): normal = b, so in-plane rows use a and c.
        let p = SlabParams {
            h: 0, k: 1, l: 0,
            u_vec: Some([2, 0, 0]),
            v_vec: Some([0, 0, 1]),
            ..Default::default()
        };
        let (ap, bp, _n, _d, _used) = inplane_basis_params(&d, &p).unwrap();
        assert!((ap.length() - 2.0 * a).abs() < 1e-3, "{}", ap.length());
        assert!((bp.length() - a).abs() < 1e-3, "{}", bp.length());
        // U along the normal → not in-plane → no in-plane basis at all.
        let bad = SlabParams {
            h: 0, k: 1, l: 0,
            u_vec: Some([0, 1, 0]),
            v_vec: Some([0, 0, 1]),
            ..Default::default()
        };
        assert!(inplane_basis_params(&d, &bad).is_none());
    }

    #[test]
    fn inplane_basis_params_honors_choice() {
        let a = 3.6;
        let d = fcc(a);
        let (ap, _bp, _, _, used) = inplane_basis_params(&d, &SlabParams {
            h: 1, k: 1, l: 1,
            u: 2,
            v: 1,
            ..Default::default()
        })
        .unwrap();
        let s2 = a / 2.0_f32.sqrt();
        assert!((ap.length() - 2.0 * s2).abs() < 1e-2, "U-expanded a': {}", ap.length());
        assert_eq!(used, InPlaneBasis::Primitive);
    }

    #[test]
    fn supercell_identity_is_noop() {
        let a = 4.0;
        let d = fcc(a);
        let out = build_supercell(&d, &SupercellParams { x: 1, y: 1, z: 1 }).unwrap();
        assert_eq!(out.atoms.len(), 4);
        assert!((out.lattice.a - a).abs() < 1e-5);
        assert_eq!(out.supercell, Some([1, 1, 1]));
        // positions unchanged (same frame, same atoms)
        for (i, o) in out.atoms.iter().enumerate() {
            let orig = &d.atoms[i];
            assert!((o.x - orig.x).abs() < 1e-4 && (o.y - orig.y).abs() < 1e-4
                && (o.z - orig.z).abs() < 1e-4, "atom {} moved", i);
        }
    }

    #[test]
    fn supercell_222_doubles_cell_and_replicates_atoms() {
        let a = 4.0;
        let d = fcc(a);
        let out = build_supercell(&d, &SupercellParams { x: 2, y: 2, z: 2 }).unwrap();
        assert_eq!(out.atoms.len(), 4 * 8);
        assert!((out.lattice.a - 2.0 * a).abs() < 1e-4, "a: {}", out.lattice.a);
        assert!((out.lattice.c - 2.0 * a).abs() < 1e-4, "c: {}", out.lattice.c);
        // corner atom (0,0,0) copies: x ∈ {0, a} etc. — check two of them
        let found = out.atoms.iter().any(|at| {
            (at.x - a).abs() < 1e-4 && at.y.abs() < 1e-4 && at.z.abs() < 1e-4
        });
        assert!(found, "missing corner-atom copy at (a, 0, 0)");
        // face atom (0, a/2, a/2) copy with i=1: (a, a/2, a/2)
        let found2 = out.atoms.iter().any(|at| {
            (at.x - a).abs() < 1e-4 && (at.y - a / 2.0).abs() < 1e-4
                && (at.z - a / 2.0).abs() < 1e-4
        });
        assert!(found2, "missing face-atom copy at (a, a/2, a/2)");
        assert_eq!(out.supercell, Some([2, 2, 2]));
    }

    #[test]
    fn supercell_211_doubles_only_a() {
        let a = 4.0;
        let d = cubic(a);
        let out = build_supercell(&d, &SupercellParams { x: 2, y: 1, z: 1 }).unwrap();
        assert_eq!(out.atoms.len(), 2);
        assert!((out.lattice.a - 2.0 * a).abs() < 1e-4);
        assert!((out.lattice.b - a).abs() < 1e-5);
        assert!((out.lattice.c - a).abs() < 1e-5);
        // the two copies sit at x = 0 and x = a
        let xs: Vec<f32> = out.atoms.iter().map(|at| at.x).collect();
        assert!(xs.iter().any(|x| x.abs() < 1e-4));
        assert!(xs.iter().any(|x| (x - a).abs() < 1e-4));
        assert_eq!(out.supercell, Some([2, 1, 1]));
    }

    #[test]
    fn supercell_keeps_slab_vacuum_provenance() {
        let a = 4.0;
        let d = fcc(a);
        let slab = build_slab(&d, &SlabParams { h: 0, k: 0, l: 1, thickness_ang: 2.0, ..Default::default() }).unwrap();
        assert!(slab.slab.is_some());
        let out = build_supercell(&slab, &SupercellParams { x: 1, y: 1, z: 2 }).unwrap();
        assert!(out.slab.is_some(), "slab provenance must survive a c-supercell");
        assert_eq!(out.supercell, Some([1, 1, 2]));
        // c doubled: the slab stack repeats with the vacuum gap intact
        let slab_c = slab.lattice.c;
        assert!((out.lattice.c - 2.0 * slab_c).abs() < 1e-4, "c: {}", out.lattice.c);
        assert_eq!(out.atoms.len(), slab.atoms.len() * 2);
    }

    #[test]
    fn supercell_clamps_below_one() {
        let d = cubic(4.0);
        let out = build_supercell(&d, &SupercellParams { x: 0, y: -3, z: 3 }).unwrap();
        assert_eq!(out.supercell, Some([1, 1, 3]));
        assert_eq!(out.atoms.len(), 3);
    }
}
