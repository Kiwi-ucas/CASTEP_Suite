//! egui UI panels: atom list, info panel, toolbar

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};
use crate::picking::PickingState;
use crate::MoveState;
use crate::AddAtomState;
use crate::CrystalStore;
use crate::Lattice;
use crate::PhononState;
use crate::RotateState;
use crate::PanelRects;
use crate::PesState;
use crate::Pes3dState;
use crate::VisMode;
use crate::IsoMaterial;
use crate::resources;

#[derive(Resource)]
pub struct AtomInfo {
    pub elements: Vec<String>,
    pub labels: Vec<String>,
    pub radii: Vec<f32>,
}

#[derive(Resource)]
pub struct CrystalMeta {
    pub filename: String,
    pub a: f32, pub b: f32, pub c: f32,
    pub alpha: f32, pub beta: f32, pub gamma: f32,
}

impl AtomInfo {
    pub fn new(elements: Vec<String>, labels: Vec<String>, radii: Vec<f32>) -> Self {
        Self { elements, labels, radii }
    }
}

pub fn ui_system(
    mut contexts: EguiContexts,
    mut picking: ResMut<PickingState>,
    mut atom_info: ResMut<AtomInfo>,
    meta: Option<Res<CrystalMeta>>,
    mut move_state: ResMut<MoveState>,
    mut add_state: ResMut<AddAtomState>,
    crystal: Option<Res<CrystalStore>>,
    mut phonon_state: Option<ResMut<PhononState>>,
    mut rotate_state: ResMut<RotateState>,
    mut panel_rects: ResMut<PanelRects>,
    pes_state: Option<Res<PesState>>,
    mut pes3d_state: Option<ResMut<Pes3dState>>,
) {
    let ctx = contexts.ctx_mut();

    // Pre-compute: parent indices for selected/hovered images
    let sel_parent = if picking.selected >= 0 && (picking.selected as usize) < picking.parent_indices.len() {
        Some(picking.parent_indices[picking.selected as usize])
    } else { None };
    let hov_parent = if picking.hovered >= 0 && (picking.hovered as usize) < picking.parent_indices.len() {
        Some(picking.parent_indices[picking.hovered as usize])
    } else { None };

    // Pre-compute first expanded-atom index for each asymmetric-unit atom
    let count = atom_info.elements.len();
    let first_images: Vec<Option<usize>> = (0..count)
        .map(|i| picking.parent_indices.iter().position(|&p| p == i))
        .collect();

    // ── Bottom toolbar ──
    let bottom_resp = egui::TopBottomPanel::bottom("toolbar").show(ctx, |ui| {
        ui.horizontal(|ui| {
            if let Some(p) = sel_parent {
                let el = if p < atom_info.elements.len() {
                    atom_info.elements[p].as_str()
                } else { "?" };
                let img_count = picking.parent_indices.iter().filter(|&&x| x == p).count();
                ui.label(egui::RichText::new(
                    format!("Atom {} ({}) selected, {} images   H/K:\u{b1}X  U/M:\u{b1}Y  I/N:\u{b1}Z   D: delete",
                        p, el, img_count)
                ).strong());
            } else {
                let mut base = "\u{1f5b0} Right-drag: rotate | Scroll: zoom | Click: select | 1/2/3: mode | P: proj | A: axes | B: bonds | C: cell | R: reset".to_string();
                if let Some(ref pes) = pes_state {
                    if pes.has_energies {
                        base.push_str(" | +/-: color range");
                    }
                }
                ui.label(base);
            }
        });
    });
    panel_rects.bottom = Some(bottom_resp.response.rect);

    // ── Left panel: atom list (asymmetric unit only) ──
    let left_resp = egui::SidePanel::left("atom_list")
        .resizable(false).default_width(180.0)
        .show(ctx, |ui| {
            ui.heading("Atoms");
            ui.separator();
            if ui.button("\u{2795} Add Atom").clicked() {
                add_state.show_table = true;
            }
            ui.separator();
            egui::ScrollArea::vertical().auto_shrink([false; 2]).show(ui, |ui| {
                for i in 0..count {
                    let label = format!("{:2}. {:3}", i, atom_info.elements[i]);
                    let selected = sel_parent == Some(i);
                    let hovered = hov_parent == Some(i);

                    let mut text = egui::RichText::new(&label);
                    if selected {
                        text = text.color(egui::Color32::from_rgb(255, 170, 30)).strong();
                    } else if hovered {
                        text = text.color(egui::Color32::from_rgb(180, 180, 180));
                    }
                    ui.horizontal(|ui| {
                        let response = ui.selectable_label(selected, text);
                        if response.clicked() {
                            if let Some(img) = first_images[i] {
                                picking.selected = img as i32;
                            }
                        }
                        if response.hovered() {
                            if let Some(img) = first_images[i] {
                                picking.hovered = img as i32;
                            }
                        }
                        // Radius edit
                        if i < atom_info.radii.len() {
                            let mut r = atom_info.radii[i];
                            let resp = ui.add(
                                egui::DragValue::new(&mut r)
                                    .speed(0.01)
                                    .range(0.1..=5.0)
                                    .suffix(" \u{c5}")
                            );
                            if resp.changed() && i < atom_info.radii.len() {
                                atom_info.radii[i] = r;
                                picking.modified = true;
                            }
                        }
                    });
                }
            });
        });
    panel_rects.left = Some(left_resp.response.rect);

    // ── Right panel ──
    let right_resp = egui::SidePanel::right("info_panel")
        .resizable(false).default_width(220.0)
        .show(ctx, |ui| {
            // Crystal info header
            if let Some(m) = meta.as_ref() {
                ui.heading(&m.filename);
                ui.separator();
                ui.label("Cell Parameters:");
                ui.label(format!("a = {:.3} \u{c5}", m.a));
                ui.label(format!("b = {:.3} \u{c5}", m.b));
                ui.label(format!("c = {:.3} \u{c5}", m.c));
                ui.label(format!("\u{3b1} = {:.2}\u{b0}", m.alpha));
                ui.label(format!("\u{3b2} = {:.2}\u{b0}", m.beta));
                ui.label(format!("\u{3b3} = {:.2}\u{b0}", m.gamma));
                // Cell volume from lattice
                if let Some(c) = crystal.as_ref() {
                    let vol = c.data.lattice.cell_volume();
                    ui.label(format!("Volume: {:.2} \u{c5}\u{b3}", vol));
                }
                ui.label(format!("Asym. atoms: {}", atom_info.elements.len()));
                ui.label(format!("Displayed: {}", picking.atom_positions.len()));
                ui.separator();
            } else {
                ui.heading("Crystal Viewer");
            }

            // Selected atom details
            if let Some(p) = sel_parent {
                if p < atom_info.elements.len() {
                    ui.label(egui::RichText::new("Selected Atom").strong());
                    ui.label(format!("Element: {}", atom_info.elements[p]));
                    ui.label(format!("Asym. index: {}", p));
                    if picking.selected >= 0 && (picking.selected as usize) < picking.atom_positions.len() {
                        let cart = picking.atom_positions[picking.selected as usize];
                        ui.label("Cartesian (\u{c5}):");
                        ui.label(format!("  X: {:.4}", cart.x));
                        ui.label(format!("  Y: {:.4}", cart.y));
                        ui.label(format!("  Z: {:.4}", cart.z));
                        // Show fractional coords if lattice data available
                        if let Some(c) = crystal.as_ref() {
                            let inv = c.data.lattice.inverse_vectors();
                            let frac = Lattice::apply_inverse(&inv, cart);
                            ui.label("Fractional:");
                            ui.label(format!("  x: {:.4}", frac.x));
                            ui.label(format!("  y: {:.4}", frac.y));
                            ui.label(format!("  z: {:.4}", frac.z));
                        }
                    }
                }
            } else {
                ui.label("Click an atom\nto see details");
            }

            // ── PES scan info (2D) ──
            if let Some(ref pes) = pes_state {
                ui.separator();
                ui.label(egui::RichText::new("PES Scan").strong());
                ui.label(format!("Plane: {}", pes.plane.to_uppercase()));
                ui.label(format!("Grid: {}×{}", pes.nx, pes.ny));
                ui.label(format!("Mode: {}", pes.scan_mode));
                if pes.has_energies {
                    ui.label(format!("E range: {:.4} – {:.4} eV", pes.e_min, pes.e_max));
                    let pct = (crate::step_to_clip(pes.color_step) * 100.0 + 0.5) as u32;
                    ui.label(format!("Color range: {}% (+/-)", pct));
                    ui.label(format!("Surface: {}", if pes.show_surface { "ON (S)" } else { "OFF (S)" }));
                } else {
                    ui.label("No energies (run CASTEP first)");
                }
            }

            // ── PES 3D scan info ──
            if let Some(ref ps) = pes3d_state {
                ui.separator();
                ui.label(egui::RichText::new("PES 3D Scan").strong());
                ui.label(format!("Grid: {}×{}×{}", ps.nx, ps.ny, ps.nz));
                if ps.has_energies {
                    ui.label(format!("E range: {:.2} – {:.2} eV", ps.e_min, ps.e_max));
                    let mode_str = match ps.vis_mode {
                        VisMode::Isosurface => format!("Isosurface (4) iso={:.2}", ps.iso_value),
                        VisMode::Volume => "Volume (5)".to_string(),
                        VisMode::Slice => format!("Slice (6) axis={} pos={:.2}", ps.slice_axis, ps.slice_pos),
                    };
                    ui.label(format!("Mode: {}", mode_str));
                } else {
                    ui.label("No energies (run CASTEP first)");
                }
            }

            // ── Isosurface controls (only when in isosurface mode) ──
            if let Some(ref mut ps) = pes3d_state {
                if ps.has_energies && ps.vis_mode == VisMode::Isosurface {
                    ui.horizontal(|ui| {
                        ui.label("Isosurface step");
                        ui.add(egui::DragValue::new(&mut ps.iso_step)
                            .speed(0.1)
                            .range(0.01..=1000.0)
                            .suffix(" eV"));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Material");
                        egui::ComboBox::from_id_salt("iso_material")
                            .selected_text(ps.iso_material.name())
                            .show_ui(ui, |ui| {
                                ui.selectable_value(&mut ps.iso_material, IsoMaterial::Opaque, "Opaque");
                                ui.selectable_value(&mut ps.iso_material, IsoMaterial::SemiTransparent, "Semi-transparent");
                                ui.selectable_value(&mut ps.iso_material, IsoMaterial::Transparent, "Transparent");
                            });
                    });
                }
            }

            // ── Color mapping range (for all 3D PES modes) ──
            if let Some(ref mut ps) = pes3d_state {
                if ps.has_energies {
                    ui.separator();
                    ui.label(egui::RichText::new("Color Mapping").strong());
                    let e_min = ps.e_min;
                    let e_max = ps.e_max;
                    let color_min = ps.color_min;
                    let color_max = ps.color_max;
                    ui.horizontal(|ui| {
                        ui.label("Min");
                        ui.add(egui::DragValue::new(&mut ps.color_min)
                            .speed(1.0)
                            .range(e_min..=color_max)
                            .suffix(" eV"));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Max");
                        ui.add(egui::DragValue::new(&mut ps.color_max)
                            .speed(1.0)
                            .range(color_min..=e_max)
                            .suffix(" eV"));
                    });
                    if ui.button("Reset to full range").clicked() {
                        ps.color_min = ps.e_min;
                        ps.color_max = ps.e_max;
                    }
                }
            }

            // ── Spatial clipping (XYZ ranges) ──
            if let Some(ref mut ps) = pes3d_state {
                if ps.has_energies && ps.vis_mode == VisMode::Isosurface {
                    ui.separator();
                    ui.label(egui::RichText::new("Spatial Clipping").strong());

                    let clip_x = ps.clip_x;
                    let clip_y = ps.clip_y;
                    let clip_z = ps.clip_z;

                    // X axis
                    ui.horizontal(|ui| {
                        ui.label("X min");
                        ui.add(egui::DragValue::new(&mut ps.clip_x[0])
                            .speed(0.01)
                            .range(0.0..=clip_x[1]));
                    });
                    ui.horizontal(|ui| {
                        ui.label("X max");
                        ui.add(egui::DragValue::new(&mut ps.clip_x[1])
                            .speed(0.01)
                            .range(clip_x[0]..=1.0));
                    });

                    // Y axis
                    ui.horizontal(|ui| {
                        ui.label("Y min");
                        ui.add(egui::DragValue::new(&mut ps.clip_y[0])
                            .speed(0.01)
                            .range(0.0..=clip_y[1]));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Y max");
                        ui.add(egui::DragValue::new(&mut ps.clip_y[1])
                            .speed(0.01)
                            .range(clip_y[0]..=1.0));
                    });

                    // Z axis
                    ui.horizontal(|ui| {
                        ui.label("Z min");
                        ui.add(egui::DragValue::new(&mut ps.clip_z[0])
                            .speed(0.01)
                            .range(0.0..=clip_z[1]));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Z max");
                        ui.add(egui::DragValue::new(&mut ps.clip_z[1])
                            .speed(0.01)
                            .range(clip_z[0]..=1.0));
                    });

                    if ui.button("Reset clipping").clicked() {
                        ps.clip_x = [0.0, 1.0];
                        ps.clip_y = [0.0, 1.0];
                        ps.clip_z = [0.0, 1.0];
                    }
                }
            }

            ui.with_layout(egui::Layout::bottom_up(egui::Align::Min), |ui| {
                ui.horizontal(|ui| {
                    ui.label("Move Step");
                    ui.add(egui::DragValue::new(&mut move_state.step)
                        .speed(0.05)
                        .suffix(" \u{c5}"));
                });
                ui.horizontal(|ui| {
                    ui.label("Rotation Angle");
                    ui.add(egui::DragValue::new(&mut rotate_state.angle_deg)
                        .suffix("\u{b0}"));
                });
                ui.separator();
            });
        });
    panel_rects.right = Some(right_resp.response.rect);

    // ── Add Atom popup: periodic table ──
    if add_state.show_table {
        egui::Window::new("Periodic Table — Select Element")
            .collapsible(false).resizable(false)
            .default_width(800.0).default_height(400.0)
            .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
            .show(ctx, |ui| {
                let mut pt_map = std::collections::HashMap::new();
                for &(sym, r, c) in resources::PERIODIC_TABLE {
                    pt_map.insert((r, c), sym);
                }
                let rows = 9;
                let cols = 18;
                let cell_size = 36.0;
                ui.vertical_centered(|ui| {
                    ui.heading("Periodic Table of Elements");
                });
                ui.add_space(6.0);
                for r in 0..rows {
                    ui.horizontal(|ui| {
                        for c in 0..cols {
                            if let Some(sym) = pt_map.get(&(r, c)) {
                                let bg = resources::category_color(sym);
                                let rich = egui::RichText::new(*sym)
                                    .size(11.0)
                                    .color(egui::Color32::WHITE);
                                let btn = egui::Button::new(rich)
                                    .fill(bg)
                                    .min_size(egui::vec2(cell_size, cell_size));
                                if ui.add(btn).clicked() {
                                    add_state.selected_element = Some(sym.to_string());
                                    add_state.show_table = false;
                                    add_state.coord_x = "0.0".into();
                                    add_state.coord_y = "0.0".into();
                                    add_state.coord_z = "0.0".into();
                                }
                            } else {
                                ui.add_sized([cell_size, cell_size], egui::Label::new(""));
                            }
                        }
                    });
                    if r == 6 || r == 7 {
                        ui.add_space(20.0);
                    }
                }
                ui.add_space(8.0);
                // Legend
                ui.horizontal(|ui| {
                    for (label, color) in [
                        ("Alkali", egui::Color32::from_rgb(0x8B, 0x3A, 0x3A)),
                        ("Alk.earth", egui::Color32::from_rgb(0x8B, 0x69, 0x3A)),
                        ("Transition", egui::Color32::from_rgb(0x3A, 0x5A, 0x7A)),
                        ("Halogen", egui::Color32::from_rgb(0x3A, 0x6B, 0x3A)),
                        ("Noble gas", egui::Color32::from_rgb(0x5A, 0x3A, 0x7A)),
                        ("Non-metal", egui::Color32::from_rgb(0x3A, 0x6B, 0x6B)),
                        ("Lanthanide", egui::Color32::from_rgb(0x7A, 0x4A, 0x5A)),
                        ("Actinide", egui::Color32::from_rgb(0x6A, 0x3A, 0x4A)),
                    ] {
                        ui.add(egui::Button::new(
                            egui::RichText::new(label).size(9.0).color(egui::Color32::WHITE)
                        ).fill(color).min_size(egui::vec2(60.0, 18.0)));
                    }
                });
                ui.add_space(6.0);
                if ui.button("Cancel").clicked() {
                    add_state.show_table = false;
                    add_state.selected_element = None;
                }
            });
    }

    // ── Coordinate input popup ──
    if !add_state.show_table && add_state.selected_element.is_some() {
        let el = add_state.selected_element.as_ref().unwrap().clone();
        egui::Window::new(format!("Add {} atom — fractional coordinates", el))
            .collapsible(false).resizable(false)
            .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
            .show(ctx, |ui| {
                ui.label("Fractional coordinates:");
                ui.horizontal(|ui| {
                    ui.label("x:"); ui.text_edit_singleline(&mut add_state.coord_x);
                });
                ui.horizontal(|ui| {
                    ui.label("y:"); ui.text_edit_singleline(&mut add_state.coord_y);
                });
                ui.horizontal(|ui| {
                    ui.label("z:"); ui.text_edit_singleline(&mut add_state.coord_z);
                });
                ui.separator();
                ui.horizontal(|ui| {
                    if ui.button("Confirm").clicked() {
                        if let (Ok(x), Ok(y), Ok(z)) = (
                            add_state.coord_x.parse::<f32>(),
                            add_state.coord_y.parse::<f32>(),
                            add_state.coord_z.parse::<f32>(),
                        ) {
                            add_state.pending = Some((el, x, y, z));
                            add_state.selected_element = None;
                        }
                    }
                    if ui.button("Cancel").clicked() {
                        add_state.selected_element = None;
                    }
                });
            });

        // ── Phonon mode info panel ──
        if let Some(ref mut ps) = phonon_state {
            egui::SidePanel::right("phonon_panel")
                .resizable(false).default_width(200.0)
                .show(ctx, |ui| {
                    ui.heading("Phonon Mode");
                    ui.separator();
                    ui.label(format!("Mode: {}", ps.mode_index));
                    ui.label(format!("Freq: {:.2} cm⁻¹", ps.frequency));
                    ui.label(format!("IR: {:.4} (D/A)²/amu", ps.ir_intensity));
                    ui.label(format!("|p_m|: {:.4}", ps.mode_charge_norm));
                    ui.separator();
                    ui.label("Arrow scale:");
                    let mut sf = ps.scale_factor;
                    if ui.add(egui::Slider::new(&mut sf, 0.1..=10.0)).changed() {
                        ps.scale_factor = sf;
                        // Note: scale change requires re-spawning arrows;
                        // For now, the scale applies on next viewer relaunch
                    }
                    let mut show = ps.show_arrows;
                    if ui.checkbox(&mut show, "Show Arrows").changed() {
                        ps.show_arrows = show;
                    }
                });
        }
    }
}
