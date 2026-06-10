//! egui UI panels: atom list, info panel, toolbar

use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};
use crate::picking::PickingState;
use crate::MoveState;

#[derive(Resource)]
pub struct AtomInfo {
    pub elements: Vec<String>,
    pub labels: Vec<String>,
}

#[derive(Resource)]
pub struct CrystalMeta {
    pub filename: String,
    pub a: f32, pub b: f32, pub c: f32,
    pub alpha: f32, pub beta: f32, pub gamma: f32,
}

impl AtomInfo {
    pub fn new(elements: Vec<String>, labels: Vec<String>) -> Self {
        Self { elements, labels }
    }
}

pub fn ui_system(
    mut contexts: EguiContexts,
    mut picking: ResMut<PickingState>,
    atom_info: Res<AtomInfo>,
    meta: Option<Res<CrystalMeta>>,
    move_state: Option<Res<MoveState>>,
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
    egui::TopBottomPanel::bottom("toolbar").show(ctx, |ui| {
        ui.horizontal(|ui| {
            if let Some(p) = sel_parent {
                let step = move_state.as_ref().map(|m| m.step).unwrap_or(0.1);
                let el = if p < atom_info.elements.len() {
                    atom_info.elements[p].as_str()
                } else { "?" };
                // Count images for this parent
                let img_count = picking.parent_indices.iter().filter(|&&x| x == p).count();
                ui.label(egui::RichText::new(
                    format!("[Atom {} ({}) selected, {} images]  I/K:\u{b1}Z  J/L:\u{b1}X  U/O:\u{b1}Y  Step:{:.2}\u{c5}  [/]:adjust",
                        p, el, img_count, step)
                ).strong());
            } else {
                ui.label("\u{1f5b0} Right-drag: rotate | Scroll: zoom | Click: select | 1/2/3: mode | P: proj | B: bonds | C: cell | R: reset");
            }
        });
    });

    // ── Left panel: atom list (asymmetric unit only) ──
    egui::SidePanel::left("atom_list")
        .resizable(true).default_width(180.0)
        .show(ctx, |ui| {
            ui.heading("Atoms (asym.)");
            ui.separator();
            egui::ScrollArea::vertical().show(ui, |ui| {
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
                }
            });
        });

    // ── Right panel ──
    egui::SidePanel::right("info_panel")
        .resizable(true).default_width(220.0)
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
                        let pos = picking.atom_positions[picking.selected as usize];
                        ui.label(format!("Cart X: {:.4}", pos.x));
                        ui.label(format!("Cart Y: {:.4}", pos.y));
                        ui.label(format!("Cart Z: {:.4}", pos.z));
                    }
                }
                if ui.button("Deselect").clicked() { picking.selected = -1; }
            } else {
                ui.label("Click an atom\nto see details");
            }
        });
}
