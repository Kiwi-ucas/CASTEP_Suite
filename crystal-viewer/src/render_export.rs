//! Bevy-native offscreen render export (WYSIWYG).
//!
//! The viewer renders with the same PBR pipeline it displays with, so the
//! exported image is pixel-identical to what the user sees (minus egui
//! panels, which only exist on the window pass). Flow:
//!
//! 1. `RenderSettings.request` is set by the Render dialog.
//! 2. `start_offscreen_render` creates a large `Image` render target, a
//!    temporary camera that copies the main camera's transform/projection/
//!    tonemapping/exposure/MSAA, and a `Screenshot::image(handle)` component.
//! 3. After the next frame renders the target, Bevy triggers
//!    `ScreenshotCaptured`; `save_to_disk` writes the PNG/TIFF, and
//!    `finish_offscreen` cleans up the temporary camera + texture and
//!    reports the result in `RenderSettings.last_status`.

use bevy::prelude::*;
use bevy::core_pipeline::tonemapping::Tonemapping;
use bevy::render::camera::{Exposure, RenderTarget};
use bevy::render::render_asset::RenderAssetUsages;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat, TextureUsages};
use bevy::render::view::screenshot::{Screenshot, ScreenshotCaptured, save_to_disk};

use crate::RenderSettings;

/// Tonemapping operator applied by the render (default: TonyMcMapface,
/// same as Bevy's default camera).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum TonemapChoice {
    TonyMcMapface,
    Aces,
    AgX,
    Reinhard,
}

impl TonemapChoice {
    pub fn label(&self) -> &'static str {
        match self {
            Self::TonyMcMapface => "TonyMcMapface",
            Self::Aces => "ACES",
            Self::AgX => "AgX",
            Self::Reinhard => "Reinhard",
        }
    }
    pub fn to_bevy(&self) -> Tonemapping {
        match self {
            Self::TonyMcMapface => Tonemapping::TonyMcMapface,
            Self::Aces => Tonemapping::AcesFitted,
            Self::AgX => Tonemapping::AgX,
            Self::Reinhard => Tonemapping::Reinhard,
        }
    }
}

/// Output image formats (encoded by Bevy's built-in `image` support; the
/// `png` feature is always on, `tiff` is enabled in Cargo.toml).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ImgFormat {
    Png,
    Tiff,
}

impl ImgFormat {
    pub fn label(&self) -> &'static str {
        match self { Self::Png => "PNG", Self::Tiff => "TIFF" }
    }
    pub fn ext(&self) -> &'static str {
        match self { Self::Png => "png", Self::Tiff => "tiff" }
    }
}

/// Marker + bookkeeping on the temporary offscreen camera.
#[derive(Component)]
pub struct OffscreenRender {
    /// Absolute path the image was saved to (filled by `start_offscreen_render`).
    pub out_path: String,
}

/// Spawn the offscreen render camera when a request is pending.
pub fn start_offscreen_render(
    mut commands: Commands,
    mut settings: ResMut<RenderSettings>,
    cam_q: Query<(&Camera, &Transform, &Projection, &Tonemapping, &Exposure, &bevy::render::view::Msaa), With<crate::MainCamera>>,
    mut images: ResMut<Assets<Image>>,
) {
    if !settings.request {
        return;
    }
    settings.request = false;
    if settings.rendering.load(std::sync::atomic::Ordering::Relaxed) {
        eprintln!("[render] busy — ignoring request");
        return;
    }
    let Ok((cam, tf, proj, tonemapping, exposure, msaa)) = cam_q.get_single() else {
        *settings.last_status.lock().unwrap() = "[render] no main camera".to_string();
        return;
    };
    let w = settings.width.clamp(320, 8192);
    let h = settings.height.clamp(240, 8192);

    let mut image = Image::new_fill(
        Extent3d { width: w, height: h, depth_or_array_layers: 1 },
        TextureDimension::D2,
        &[0, 0, 0, 0],
        TextureFormat::Rgba8UnormSrgb,
        RenderAssetUsages::RENDER_WORLD,
    );
    image.texture_descriptor.usage =
        TextureUsages::TEXTURE_BINDING | TextureUsages::COPY_DST | TextureUsages::RENDER_ATTACHMENT;
    let handle = images.add(image);

    let out_path = std::env::current_dir()
        .map(|d| d.join(format!("render.{}", settings.format.ext())))
        .unwrap_or_else(|_| std::path::PathBuf::from(format!("render.{}", settings.format.ext())));
    let out_path_str = out_path.display().to_string();

    settings.rendering.store(true, std::sync::atomic::Ordering::Relaxed);
    *settings.last_status.lock().unwrap() = format!("rendering {w}×{h}…");

    let mut cmd = commands.spawn((
        Camera3d::default(),
        Camera {
            target: RenderTarget::Image(handle.clone()),
            clear_color: cam.clear_color,
            ..default()
        },
        tonemapping.clone(),
        exposure.clone(),
        msaa.clone(),
        proj.clone(),
        tf.clone(),
        Screenshot::image(handle),
        OffscreenRender { out_path: out_path_str.clone() },
    ));
    cmd.observe(save_to_disk(out_path_str))
        .observe(finish_offscreen);
}

/// Clean up after `ScreenshotCaptured`: despawn the temporary camera, free
/// the render-target texture, and publish the result.
fn finish_offscreen(
    trigger: Trigger<ScreenshotCaptured>,
    mut commands: Commands,
    q: Query<&OffscreenRender>,
    settings: Res<RenderSettings>,
) {
    let entity = trigger.entity();
    let Ok(or) = q.get(entity) else { return };
    let path = or.out_path.clone();
    let _ = trigger.event(); // image already consumed by save_to_disk
    commands.entity(entity).despawn();
    settings.rendering.store(false, std::sync::atomic::Ordering::Relaxed);
    *settings.last_status.lock().unwrap() = format!("saved {}", path);
    eprintln!("[render] {}", settings.last_status.lock().unwrap());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn img_format_labels() {
        assert_eq!(ImgFormat::Png.label(), "PNG");
        assert_eq!(ImgFormat::Png.ext(), "png");
        assert_eq!(ImgFormat::Tiff.label(), "TIFF");
        assert_eq!(ImgFormat::Tiff.ext(), "tiff");
    }

    #[test]
    fn render_settings_defaults() {
        let rs = RenderSettings::default();
        assert!(!rs.request && !rs.show_dialog);
        assert_eq!((rs.width, rs.height), (1920, 1080));
        assert_eq!(rs.format, ImgFormat::Png);
    }
}
