#include "raylib.h"
#include "../unified_sanctum_grammar.h"
#include "../unified_sanctum_rasterizer.h"
#include "../unified_sanctum_grammar.c"
#include <stdio.h>
#include <math.h>

#define MAX_PREVIEW_VOXELS 1048576

typedef struct {
    Vector3 pos;
    Color color;
} SimpleVoxel;

static SimpleVoxel preview_voxels[MAX_PREVIEW_VOXELS];
static int preview_voxel_count = 0;

static void plot_sanctum_preview(int gx, int gy, int gz, Color c) {
    if (preview_voxel_count >= MAX_PREVIEW_VOXELS) return;
    float voxel_size = 0.5f;
    float floor_size = 40.0f; // 80m arena
    float px = (gx + 0.5f) * voxel_size - floor_size;
    float py = (gy + 0.5f) * voxel_size;
    float pz = (gz + 0.5f) * voxel_size - floor_size;
    preview_voxels[preview_voxel_count++] = (SimpleVoxel){ (Vector3){ px, py, pz }, c };
}

static void render_sanctum_shot(uint32_t seed, int growth_steps,
                                Vector3 cam_pos, Vector3 cam_target, float fovy,
                                const char *out_path, const char *title, const char *subtitle) {
    preview_voxel_count = 0;
    SanctumCitadelPlan plan = generate_unified_sanctum(seed, growth_steps);
    int center_g = 80; // 40m / 0.5m = 80 voxels
    rasterize_sanctum_plan(&plan, center_g, center_g, 10, plot_sanctum_preview, NULL);

    Camera3D camera = { 0 };
    camera.position = cam_pos;
    camera.target = cam_target;
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = fovy;
    camera.projection = CAMERA_PERSPECTIVE;

    RenderTexture2D target = LoadRenderTexture(1280, 720);
    BeginTextureMode(target);
        // Atmospheric Twilight Horizon: Deep indigo twilight with subtle cosmic aurora
        ClearBackground((Color){ 28, 36, 48, 255 });

        BeginMode3D(camera);
            // Sacred terrain ground plane
            DrawCube((Vector3){ 0.0f, 4.5f, 0.0f }, 90.0f, 1.0f, 90.0f, (Color){ 44, 48, 54, 255 });
            DrawCubeWires((Vector3){ 0.0f, 4.5f, 0.0f }, 90.0f, 1.0f, 90.0f, (Color){ 32, 36, 42, 255 });

            float vs = 0.5f;
            for (int i = 0; i < preview_voxel_count; ++i) {
                DrawCube(preview_voxels[i].pos, vs, vs, vs, preview_voxels[i].color);
                DrawCubeWires(preview_voxels[i].pos, vs, vs, vs, ColorBrightness(preview_voxels[i].color, -0.18f));
            }
        EndMode3D();

        // Syncretic UI telemetry banner
        DrawRectangle(20, 20, 780, 85, Fade((Color){ 16, 22, 30, 255 }, 0.88f));
        DrawRectangleLines(20, 20, 780, 85, (Color){ 45, 225, 255, 255 });
        DrawText(title, 35, 30, 22, (Color){ 240, 248, 255, 255 });
        DrawText(subtitle, 35, 56, 16, (Color){ 255, 205, 90, 255 });
        DrawText(TextFormat("Nodes: %d | Voxels: %d | Sarsen: %d | Marble: %d | Titanium: %d | HardLight: %d | Steps: %d",
                            plan.node_count, preview_voxel_count, plan.count_stone, plan.count_marble,
                            plan.count_titanium, plan.count_hardlight, growth_steps),
                 35, 78, 14, (Color){ 120, 220, 255, 255 });

    EndTextureMode();

    Image img = LoadImageFromTexture(target.texture);
    ImageFlipVertical(&img);
    ExportImage(img, out_path);
    printf("Rendered [%s] -> %s (%d voxels, %d nodes)\n", title, out_path, preview_voxel_count, plan.node_count);

    UnloadImage(img);
    UnloadRenderTexture(target);
}

int main(void) {
    SetConfigFlags(FLAG_WINDOW_HIDDEN);
    InitWindow(1280, 720, "Unified Syncretic Citadel Gallery Renderer");
    if (!IsWindowReady()) return 1;

    const char *out_dir = "/Users/migr4479/.gemini/antigravity/brain/86453e72-1dcc-4b9c-ba5a-2222f2dbe7ff";

    // 1. Sprawling Aerial Panorama of the Monumental Megaron Acropolis (Step 20)
    render_sanctum_shot(2552, 20,
                        (Vector3){ 45.0f, 42.0f, 45.0f },
                        (Vector3){ 0.0f, 16.0f, 0.0f },
                        54.0f,
                        TextFormat("%s/unified_sanctum_overview.png", out_dir),
                        "The Monumental Precursor Citadel: Aerial Panorama",
                        "Colossal Megaron Sanctuary (Y=38), Perimeter Stoas, Obelisk Plazas & Stelae Avenues");

    // 2. Colossal Propylaea Portal & Great Hypostyle Hall
    render_sanctum_shot(777, 16,
                        (Vector3){ 0.0f, 18.0f, -36.0f },
                        (Vector3){ 0.0f, 16.0f, 0.0f },
                        52.0f,
                        TextFormat("%s/unified_sanctum_hypostyle_hall.png", out_dir),
                        "The Colossal Propylaea & Great Megaron Hall",
                        "Twin Canted Pylons, Inscribed Lintel, Enclosed Cella Walls, Coffered Ceiling & Braziers");

    // 3. West Obelisk Plaza & Processional Stelae Avenue
    render_sanctum_shot(117, 14,
                        (Vector3){ -28.0f, 18.0f, -14.0f },
                        (Vector3){ -12.0f, 10.0f, -2.0f },
                        52.0f,
                        TextFormat("%s/unified_sanctum_courtyard_plaza.png", out_dir),
                        "The Monolithic Obelisk Plaza & Stelae Avenue",
                        "Tapered Sarsen Needle, Sacrificial Hearth, PBF Water Rills, Megalith Rows & Braziers");

    // 4. Heroic Low-Angle Vista of the Towering Megastructure (Y=38)
    render_sanctum_shot(343, 18,
                        (Vector3){ 0.0f, 12.0f, -48.0f },
                        (Vector3){ 0.0f, 22.0f, 0.0f },
                        52.0f,
                        TextFormat("%s/unified_sanctum_acropolis_terrace.png", out_dir),
                        "Heroic Skyline: Apex Transmission Needle & Bastions",
                        "Multi-Tiered Ziggurat Podium, Canted Buttresses, Corner Bastions & Y=38 Oracle Core");

    CloseWindow();
    return 0;
}
