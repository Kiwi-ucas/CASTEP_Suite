//! PES 3D volume rendering — GPU ray marching in fractional-coordinate space.
//!
//! The proxy cube spans the unit cell (a parallelpiped in world space).
//! Ray origins/directions are transformed to fractional space via inv_lattice,
//! where the volume is the unit cube [0,1]^3 and sampling is trivially correct
//! even for non-orthogonal cells. Front-to-back compositing with early
//! termination; straight-alpha output paired with AlphaMode::Blend.

#define_import_path volume::volume_material

#import bevy_pbr::mesh_view_bindings::view
#import bevy_pbr::forward_io::VertexOutput

struct VolumeParams {
    color_min: f32,
    color_max: f32,
    iso_value: f32,
    alpha_scale: f32,
    alpha_falloff: f32,
    clip_x: vec2<f32>,
    clip_y: vec2<f32>,
    clip_z: vec2<f32>,
    inv_lattice: mat3x3<f32>,
    steps: u32,
}

@group(2) @binding(0) var volume_texture: texture_3d<f32>;
@group(2) @binding(1) var volume_sampler: sampler;
@group(2) @binding(2) var<uniform> p: VolumeParams;

// ── Jet colormap (mirrors pes.rs jet_rgb: blue → cyan → green → yellow → red) ──

fn jet(x_in: f32) -> vec3<f32> {
    let t = clamp(x_in, 0.0, 1.0);
    var r = 0.0;
    var g = 0.0;
    var b = 0.0;
    if (t < 0.375) {
        r = 0.0;
    } else if (t < 0.625) {
        r = (t - 0.375) / 0.25;
    } else if (t < 0.875) {
        r = 1.0;
    } else {
        r = 1.0 - (t - 0.875) / 0.125 * 0.5;
    }
    if (t < 0.125) {
        g = 0.0;
    } else if (t < 0.375) {
        g = (t - 0.125) / 0.25;
    } else if (t < 0.625) {
        g = 1.0;
    } else if (t < 0.875) {
        g = 1.0 - (t - 0.625) / 0.25;
    } else {
        g = 0.0;
    }
    if (t < 0.125) {
        b = 0.5 + t / 0.125 * 0.5;
    } else if (t < 0.375) {
        b = 1.0;
    } else if (t < 0.625) {
        b = 1.0 - (t - 0.375) / 0.25;
    } else {
        b = 0.0;
    }
    return vec3<f32>(r, g, b);
}

// ── Transfer function: energy → (RGB, alpha) ──
// Band-pass window centered on iso_value: alpha peaks at the ISO ref energy
// and falls off as a gaussian on BOTH sides. Dragging the ISO ref slides the
// energy window through the volume, revealing the spatial structure of each
// energy layer. (A one-sided "solid below / fade above" transfer function was
// unusable here: the 91% low-energy bulk saturates front-to-back compositing
// before the ISO-controlled region is ever reached.)
//
// Falloff controls the window width: smaller = wider window, larger = sharper
// layer selection.

fn transfer(v: f32, dt: f32) -> vec4<f32> {
    // v is already normalized to [0, 1] on the CPU by build_volume_texture
    // using [color_min, color_max] — do NOT re-normalize here (double
    // normalization would clamp every sample to 0).
    let t = clamp(v, 0.0, 1.0);
    let range = max(p.color_max - p.color_min, 1e-6);
    let t_iso = clamp((p.iso_value - p.color_min) / range, 0.0, 1.0);
    let x = t - t_iso;
    let alpha_raw = p.alpha_scale * exp(-p.alpha_falloff * x * x);
    // Beer-Lambert step scaling: per-sample alpha ∝ step size, so total
    // absorption along a ray is independent of the step count and depends
    // only on path length. Without this, dense sampling (128/192 steps)
    // saturates at the surface and hides the interior.
    //
    // Normalization: ×64 keeps the 64-step visual identical to the old
    // per-sample-alpha behavior (the user's preferred quality setting);
    // 128/192 steps then accumulate the SAME total absorption as 64 steps.
    let alpha = clamp(alpha_raw * dt * 64.0, 0.0, 1.0);
    return vec4<f32>(jet(t), alpha);
}

// ── Ray / axis-aligned box intersection (slab method), fractional space ──
// Returns (t_near, t_far); t_near > t_far means no intersection.

fn ray_box(ro: vec3<f32>, rd: vec3<f32>, lo: vec3<f32>, hi: vec3<f32>) -> vec2<f32> {
    var tn = -1e6;
    var tf = 1e6;
    for (var a = 0u; a < 3u; a++) {
        if (abs(rd[a]) < 1e-6) {
            if (ro[a] < lo[a] || ro[a] > hi[a]) {
                return vec2<f32>(1.0, -1.0);
            }
        } else {
            let t1 = (lo[a] - ro[a]) / rd[a];
            let t2 = (hi[a] - ro[a]) / rd[a];
            tn = max(tn, min(t1, t2));
            tf = min(tf, max(t1, t2));
        }
    }
    return vec2<f32>(tn, tf);
}

@fragment
fn fragment(mesh: VertexOutput) -> @location(0) vec4<f32> {
    // Orthographic vs perspective: Bevy's clip_from_view w-component is 1.0
    // for orthographic projections.
    let ortho = view.clip_from_view[3].w == 1.0;
    var rd_w = normalize(mesh.world_position.xyz - view.world_position);
    if (ortho) {
        rd_w = -normalize(view.world_from_view[2].xyz);
    }
    var ro_w = view.world_position;
    if (ortho) {
        ro_w = mesh.world_position.xyz - rd_w * 1e4;
    }

    // World → fractional space (parallelpiped → unit cube).
    // rd is NORMALIZED so that the ray parameter t is the actual fractional
    // distance. Without this, |rd| = 1/a (e.g. 1/10.28 for a 10.28 Å cubic
    // cell) and the step size dt is inflated by a factor of a — which makes
    // the Beer-Lambert per-sample alpha ~10× too large → surface saturation.
    let ro = p.inv_lattice * ro_w;
    let rd = normalize(p.inv_lattice * rd_w);

    // Intersect the unit cell and the user clipping box
    let tt = ray_box(ro, rd, vec3<f32>(0.0), vec3<f32>(1.0));
    let ct = ray_box(ro, rd,
        vec3<f32>(p.clip_x.x, p.clip_y.x, p.clip_z.x),
        vec3<f32>(p.clip_x.y, p.clip_y.y, p.clip_z.y));
    var tn = max(tt.x, ct.x);
    var tf = min(tt.y, ct.y);
    tn = max(tn, 0.0);  // camera inside the box: start from the camera
    if (tn > tf) {
        discard;
    }

    let dt = (tf - tn) / f32(p.steps);
    var t = tn + dt * 0.5;
    var rgb = vec3<f32>(0.0);
    var a = 0.0;
    for (var i = 0u; i < p.steps; i++) {
        if (a >= 0.98) {
            break;  // early termination — fully opaque
        }
        let s = textureSample(volume_texture, volume_sampler, ro + rd * t).r;
        let src = transfer(s, dt);
        if (src.a > 0.001) {
            rgb += (1.0 - a) * src.a * src.rgb;
            a += (1.0 - a) * src.a;
        }
        t += dt;
    }
    return vec4<f32>(rgb, a);  // straight alpha; AlphaMode::Blend handles the rest
}
