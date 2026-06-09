//! Crystal Viewer — Interactive 3D crystal structure visualization

mod crystal; mod resources; mod picking; mod ui;

use bevy::prelude::*;
use bevy::input::mouse::{MouseMotion, MouseScrollUnit, MouseWheel};
use bevy::render::mesh::{Mesh, Indices, PrimitiveTopology};
use crystal::CrystalData;
use picking::{PickingState, click_pick, hover_pick, highlight_atoms};
use ui::{AtomInfo, CrystalMeta, ui_system};
use std::f32::consts::PI;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let json_path = args.get(1).map(|s| s.as_str()).unwrap_or("");

    App::new()
        .add_plugins(
            DefaultPlugins
                .build()
                .disable::<bevy::audio::AudioPlugin>()
                .disable::<bevy::gilrs::GilrsPlugin>(),
        )
        .add_plugins(bevy_egui::EguiPlugin)
        .insert_resource(CrystalPath(json_path.to_string()))
        .add_systems(Startup, setup)
        .add_systems(Update, orbit_camera)
        .add_systems(Update, (click_pick, hover_pick).chain())
        .add_systems(Update, highlight_atoms)
        .add_systems(Update, ui_system.after(highlight_atoms))
        .add_systems(Update, display_mode_system)
        .add_systems(Update, toggle_projection)
        .insert_resource(ProjMode::Perspective)
        .insert_resource(OrthoScale(15.0))
        .run();
}

#[derive(Resource)] struct CrystalPath(String);
#[derive(Resource)] struct CrystalScene { center: Vec3 }
#[derive(Component)] struct MainCamera;
#[derive(Component)] struct FollowCamera;
#[derive(Component)] struct AtomMarker;
#[derive(Component)] struct BondMarker;
#[derive(Component)] struct CellMarker;

#[derive(Resource)]
struct DisplayMode {
    mode: u8,          // 1=ball-stick, 2=space-filling, 3=wireframe
    show_bonds: bool,
    show_cell: bool,
}

/// Saved initial camera state for R-key reset
#[derive(Resource)]
struct CameraInit { focus: Vec3, radius: f32, yaw: f32, pitch: f32 }

#[derive(Resource)]
pub struct CameraState { pub focus: Vec3, radius: f32, yaw: f32, pitch: f32 }

#[derive(Default)] struct InputState { rotating: bool }

#[derive(Resource, PartialEq)]
enum ProjMode { Perspective, Orthographic }

fn toggle_projection(
    keys: Res<ButtonInput<KeyCode>>,
    windows: Query<&Window, With<bevy::window::PrimaryWindow>>,
    mut proj_mode: ResMut<ProjMode>,
    ortho_scale: ResMut<OrthoScale>,
    mut camera_q: Query<&mut Projection, With<MainCamera>>,
) {
    if keys.just_pressed(KeyCode::KeyP) {
        *proj_mode = match *proj_mode {
            ProjMode::Perspective => ProjMode::Orthographic,
            ProjMode::Orthographic => ProjMode::Perspective,
        };
        let Ok(mut proj) = camera_q.get_single_mut() else { return };
        *proj = match *proj_mode {
            ProjMode::Perspective => Projection::Perspective(PerspectiveProjection {
                fov: 1.0, ..default()
            }),
            ProjMode::Orthographic => {
                let Ok(window) = windows.get_single() else {
                    return;
                };
                let w = window.physical_width() as f32;
                let h = window.physical_height() as f32;
                let aspect = w / h;
                let vh = ortho_scale.0;
                let vw = vh * aspect;
                Projection::Orthographic(OrthographicProjection {
                    scaling_mode: bevy::render::camera::ScalingMode::Fixed {
                        width: vw, height: vh,
                    },
                    near: 0.0,
                    far: 1000.0,
                    ..OrthographicProjection::default_3d()
                })
            }
        };
    }
}

#[derive(Resource)]
struct OrthoScale(f32);  // vertical world units visible in ortho mode

fn orbit_camera(
    mut camera_q: Query<&mut Transform, (With<MainCamera>, Without<FollowCamera>)>,
    mut proj_q: Query<&mut Projection, With<MainCamera>>,
    mut light_q: Query<&mut Transform, With<FollowCamera>>,
    mut cam_state: ResMut<CameraState>,
    cam_init: Res<CameraInit>,
    proj_mode: Res<ProjMode>,
    mut ortho_scale: ResMut<OrthoScale>,
    scene: Option<Res<CrystalScene>>,
    windows: Query<&Window, With<bevy::window::PrimaryWindow>>,
    mut input: Local<InputState>,
    mouse_btn: Res<ButtonInput<MouseButton>>,
    mut mouse_motion: EventReader<MouseMotion>,
    mut mouse_wheel: EventReader<MouseWheel>,
    keys: Res<ButtonInput<KeyCode>>,
) {
    let Ok(mut cam) = camera_q.get_single_mut() else { return };

    // Right mouse → rotate (no limit on pitch for infinite rotation)
    if mouse_btn.just_pressed(MouseButton::Right) { input.rotating = true; }
    if mouse_btn.just_released(MouseButton::Right) { input.rotating = false; }

    for motion in mouse_motion.read() {
        let d = motion.delta;
        if input.rotating {
            cam_state.yaw -= d.x * 0.005;
            cam_state.pitch += d.y * 0.005;   // no clamp — infinite rotation
        }
    }

    for ev in mouse_wheel.read() {
        let dy = match ev.unit {
            MouseScrollUnit::Line => ev.y * 0.1,
            MouseScrollUnit::Pixel => ev.y * 0.001,
        };
        if *proj_mode == ProjMode::Orthographic {
            ortho_scale.0 = (ortho_scale.0 - dy * ortho_scale.0 * 0.1).clamp(1.0, 200.0);
            // Rebuild ortho projection with new scale
            if let Ok(mut proj) = proj_q.get_single_mut() {
                if let Ok(window) = windows.get_single() {
                    let w = window.physical_width() as f32;
                    let h = window.physical_height() as f32;
                    let aspect = w / h;
                    let vh = ortho_scale.0;
                    let vw = vh * aspect;
                    *proj = Projection::Orthographic(OrthographicProjection {
                        scaling_mode: bevy::render::camera::ScalingMode::Fixed { width: vw, height: vh },
                        near: 0.0, far: 1000.0,
                        ..OrthographicProjection::default_3d()
                    });
                }
            }
        } else {
            cam_state.radius = (cam_state.radius - dy * cam_state.radius * 0.1).clamp(1.0, 100.0);
        }
    }

    if keys.just_pressed(KeyCode::KeyF) {
        if let Some(s) = scene.as_ref() { cam_state.focus = s.center; cam_state.radius = 10.0; }
    }
    if keys.just_pressed(KeyCode::KeyR) {
        *cam_state = CameraState {
            focus: cam_init.focus, radius: cam_init.radius,
            yaw: cam_init.yaw, pitch: cam_init.pitch,
        };
    }

    let pos = cam_state.focus + spherical(cam_state.radius, cam_state.yaw, cam_state.pitch);
    cam.translation = pos;
    cam.look_at(cam_state.focus, Vec3::Y);

    for mut lt in light_q.iter_mut() {
        let lp = pos + *cam.right() * cam_state.radius * 0.8 + *cam.up() * cam_state.radius * 0.6;
        lt.translation = lp;
        lt.look_at(cam_state.focus, Vec3::Y);
    }
}

fn spherical(r: f32, yaw: f32, pitch: f32) -> Vec3 {
    Vec3::new(r * pitch.cos() * yaw.sin(), r * pitch.sin(), r * pitch.cos() * yaw.cos())
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

fn setup(
    mut commands: Commands, mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>, crystal_path: Res<CrystalPath>,
) {
    let data = if crystal_path.0.is_empty() { default_cu_fcc() } else {
        CrystalData::from_json(&crystal_path.0).unwrap_or_else(|e| {
            eprintln!("Failed to load {}: {}. Using Cu FCC.", crystal_path.0, e);
            default_cu_fcc()
        })
    };

    let positions = data.cartesian_positions();
    let center = data.center();
    let corners = data.cell_corners();
    let edges = CrystalData::cell_edges();
    let n = positions.len();

    let sphere = meshes.add(uv_sphere(0.5, 32, 32));
    let mut handles = Vec::with_capacity(n);

    for (i, pos) in positions.iter().enumerate() {
        let el = &data.atoms[i].element;
        let color = resources::element_color(el);
        let mat_handle = materials.add(StandardMaterial {
            base_color: color, metallic: 0.2, perceptual_roughness: 0.5, ..default()
        });
        commands.spawn((
            Mesh3d(sphere.clone()), MeshMaterial3d(mat_handle.clone()),
            Transform::from_translation(*pos), AtomMarker,
        ));
        handles.push(mat_handle);
    }

    // Bonds
    let bond_mat = materials.add(StandardMaterial {
        base_color: Color::srgba(0.65, 0.65, 0.65, 0.5),
        alpha_mode: AlphaMode::Blend, ..default()
    });
    for i in 0..n { for j in (i+1)..n {
        if resources::has_bond(&data.atoms[i].element, &data.atoms[j].element,
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
        let mid = (s + e) * 0.5; let dir = (e - s).normalize(); let len = s.distance(e);
        commands.spawn((
            Mesh3d(meshes.add(Cuboid::new(0.04, 0.04, len))),
            MeshMaterial3d(edge_mat.clone()),
            Transform::from_translation(mid).with_rotation(Quat::from_rotation_arc(Vec3::Z, dir)),
            CellMarker,
        ));
    }

    commands.insert_resource(PickingState::new(positions, handles));

    // Atom metadata for UI
    let elements: Vec<String> = data.atoms.iter().map(|a| a.element.clone()).collect();
    let labels: Vec<String> = data.atoms.iter().enumerate()
        .map(|(_i, a)| if a.label.is_empty() { a.element.clone() } else { a.label.clone() })
        .collect();
    commands.insert_resource(AtomInfo::new(elements, labels));
    commands.insert_resource(DisplayMode { mode: 1, show_bonds: true, show_cell: true });

    // Crystal metadata for UI
    let fname = if crystal_path.0.is_empty() {
        "Cu FCC (demo)".to_string()
    } else {
        std::path::Path::new(&crystal_path.0)
            .file_name().map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| crystal_path.0.clone())
    };
    commands.insert_resource(CrystalMeta {
        filename: fname,
        a: data.lattice.a, b: data.lattice.b, c: data.lattice.c,
        alpha: data.lattice.alpha, beta: data.lattice.beta, gamma: data.lattice.gamma,
    });

    // Init camera
    let radius = 10.0; let yaw = -0.75; let pitch = 0.5;
    commands.insert_resource(CameraState { focus: center, radius, yaw, pitch });
    commands.insert_resource(CameraInit { focus: center, radius, yaw, pitch });
    commands.insert_resource(CrystalScene { center });

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
        Transform::from_xyz(10.0, 8.0, 10.0).looking_at(center, Vec3::Y),
        MainCamera,
    ));

    println!("Loaded {} atoms. Right-drag to rotate, scroll to zoom, click to select.", n);
}

fn spawn_bond(
    commands: &mut Commands, meshes: &mut ResMut<Assets<Mesh>>,
    material: &Handle<StandardMaterial>, a: Vec3, b: Vec3, marker: BondMarker,
) {
    let mid = (a + b) * 0.5; let dir = (b - a).normalize(); let len = a.distance(b);
    commands.spawn((
        Mesh3d(meshes.add(Cylinder::new(0.08, len))),
        MeshMaterial3d(material.clone()),
        Transform::from_translation(mid).with_rotation(Quat::from_rotation_arc(Vec3::Y, dir)),
        marker,
    ));
}

fn display_mode_system(
    keys: Res<ButtonInput<KeyCode>>,
    mut display: ResMut<DisplayMode>,
    mut atoms: Query<&mut Transform, (With<AtomMarker>, Without<BondMarker>)>,
    bonds: Query<Entity, With<BondMarker>>,
    cells: Query<Entity, With<CellMarker>>,
    mut commands: Commands,
) {
    if keys.just_pressed(KeyCode::Digit1) {
        display.mode = 1; display.show_bonds = true;
        // Scale atoms back to normal, show bonds & cell
        for mut t in atoms.iter_mut() { t.scale = Vec3::ONE; }
        for e in bonds.iter() { commands.entity(e).insert(Visibility::Visible); }
        for e in cells.iter() { commands.entity(e).insert(Visibility::Visible); }
    }
    if keys.just_pressed(KeyCode::Digit2) {
        display.mode = 2; display.show_bonds = false;
        // Space-filling: larger atoms, hide bonds
        for mut t in atoms.iter_mut() { t.scale = Vec3::splat(2.0); }
        for e in bonds.iter() { commands.entity(e).insert(Visibility::Hidden); }
        for e in cells.iter() { commands.entity(e).insert(Visibility::Visible); }
    }
    if keys.just_pressed(KeyCode::Digit3) {
        display.mode = 3;
        // Wireframe: tiny atoms, show bonds & cell
        for mut t in atoms.iter_mut() { t.scale = Vec3::splat(0.3); }
        for e in bonds.iter() { commands.entity(e).insert(Visibility::Visible); }
        for e in cells.iter() { commands.entity(e).insert(Visibility::Visible); }
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
    }
}
