#include "raylib.h"
#include "../megalith_grammar.h"
#include "../megalith_rasterizer.h"
#include "../megalith_grammar.c"
#include <stdio.h>
#include <math.h>

#define MAX_PREVIEW_VOXELS 524288

typedef struct {
    Vector3 pos;
    Color color;
} SimpleVoxel;

static SimpleVoxel preview_voxels[MAX_PREVIEW_VOXELS];
static int preview_voxel_count = 0;

static void plot_megalith_preview(int gx, int gy, int gz, Color c) {
    if (preview_voxel_count >= MAX_PREVIEW_VOXELS) return;
    float voxel_size = 0.5f;
    float floor_size = 25.0f; // 50m arena
    float px = (gx + 0.5f) * voxel_size - floor_size;
    float py = (gy + 0.5f) * voxel_size;
    float pz = (gz + 0.5f) * voxel_size - floor_size;
    preview_voxels[preview_voxel_count++] = (SimpleVoxel){ (Vector3){ px, py, pz }, c };
}

static void render_megalith_shot(MegalithArchetype archetype, MegalithStage stage, uint32_t seed,
                                 Vector3 cam_pos, Vector3 cam_target, float fovy,
                                 const char *out_path, const char *title, const char *subtitle) {
    preview_voxel_count = 0;
    MegalithPlan plan = generate_megalith_structure(seed, archetype, stage);
    int center_g = 50; // Floor size 25m / 0.5m = 50 voxels
    rasterize_megalith_plan(&plan, center_g, center_g, 2, plot_megalith_preview);

    Camera3D camera = { 0 };
    camera.position = cam_pos;
    camera.target = cam_target;
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = fovy;
    camera.projection = CAMERA_PERSPECTIVE;

    RenderTexture2D target = LoadRenderTexture(1280, 720);
    BeginTextureMode(target);
        // Misty prehistoric moorland / neolithic dawn atmosphere
        ClearBackground((Color){ 165, 180, 185, 255 });

        BeginMode3D(camera);
            // Grassy windswept moorland floor
            DrawCube((Vector3){ 0.0f, 0.5f, 0.0f }, 60.0f, 1.0f, 60.0f, (Color){ 70, 95, 60, 255 });
            DrawCubeWires((Vector3){ 0.0f, 0.5f, 0.0f }, 60.0f, 1.0f, 60.0f, (Color){ 50, 75, 45, 255 });

            float vs = 0.5f;
            for (int i = 0; i < preview_voxel_count; ++i) {
                DrawCube(preview_voxels[i].pos, vs, vs, vs, preview_voxels[i].color);
                DrawCubeWires(preview_voxels[i].pos, vs, vs, vs, ColorBrightness(preview_voxels[i].color, -0.25f));
            }
        EndMode3D();

        // Atmospheric vignette banner
        DrawRectangle(20, 20, 620, 85, Fade(BLACK, 0.70f));
        DrawRectangleLines(20, 20, 620, 85, (Color){ 200, 195, 180, 255 });
        DrawText(title, 35, 30, 22, RAYWHITE);
        DrawText(subtitle, 35, 56, 16, (Color){ 180, 210, 175, 255 });
        DrawText(TextFormat("Nodes: %d | Voxels: %d | Seed: %u", plan.node_count, preview_voxel_count, seed), 35, 78, 15, (Color){ 220, 180, 80, 255 });

    EndTextureMode();

    Image img = LoadImageFromTexture(target.texture);
    ImageFlipVertical(&img);
    ExportImage(img, out_path);
    printf("Rendered [%s] -> %s (%d voxels)\n", title, out_path, preview_voxel_count);

    UnloadImage(img);
    UnloadRenderTexture(target);
}

int main(void) {
    SetConfigFlags(FLAG_WINDOW_HIDDEN);
    InitWindow(1280, 720, "Megalith Gallery Renderer");
    if (!IsWindowReady()) return 1;

    const char *out_dir = "/Users/migr4479/.gemini/antigravity/brain/86453e72-1dcc-4b9c-ba5a-2222f2dbe7ff";

    // 1. Dolmen Portal Tomb (Stage 1: Passage Grave Archetype)
    render_megalith_shot(MEGALITH_ARCHETYPE_PASSAGE_GRAVE, MEGALITH_STAGE_DOLMEN, 42,
                         (Vector3){ 15.0f, 10.5f, 15.0f },
                         (Vector3){ 0.0f, 2.5f, 0.0f },
                         35.0f,
                         TextFormat("%s/dolmen_portal_tomb.png", out_dir),
                         "Prehistoric Dolmen: Portal Tomb",
                         "2 Supporting Orthostats + Gravitationally Balanced Capstone Slab");

    // 2. Corbelled Burial Chamber Interior with Ritual Hearth
    render_megalith_shot(MEGALITH_ARCHETYPE_PASSAGE_GRAVE, MEGALITH_STAGE_CHAMBER, 1337,
                         (Vector3){ 0.0f, 2.2f, 5.2f },
                         (Vector3){ 0.0f, 2.0f, -0.5f },
                         52.0f,
                         TextFormat("%s/chamber_interior_hearth.png", out_dir),
                         "Tomb Sanctum: Corbelled Chamber & Ritual Hearth",
                         "Overhanging Drystone Corbels, Basal Orthostats & Glowing Embers");

    // 3. Elongated Solstice Dromos Corridor (Pre-Mound)
    render_megalith_shot(MEGALITH_ARCHETYPE_PASSAGE_GRAVE, MEGALITH_STAGE_PASSAGE_GRAVE, 1337,
                         (Vector3){ 12.0f, 8.5f, 16.0f },
                         (Vector3){ 0.0f, 2.5f, 4.0f },
                         45.0f,
                         TextFormat("%s/solstice_passage_corridor.png", out_dir),
                         "Neolithic Dromos: Solstice Passage Corridor",
                         "Paired Megalithic Wall Slabs & Transverse Roof Lintels");

    // 4. Earthen Tumulus Barrow with Kerb Ring & Solstice Passage
    render_megalith_shot(MEGALITH_ARCHETYPE_PASSAGE_GRAVE, MEGALITH_STAGE_TUMULUS, 1337,
                         (Vector3){ 18.0f, 15.0f, 24.0f },
                         (Vector3){ 0.0f, 4.0f, 6.0f },
                         50.0f,
                         TextFormat("%s/passage_grave_tumulus.png", out_dir),
                         "Neolithic Tumulus Barrow & Passage Grave",
                         "Earthen Mantle, Solstice Dromos Corridor & Peristalith Kerb Ring");

    // 4. Cromlech Stone Circle (Radial expansion)
    render_megalith_shot(MEGALITH_ARCHETYPE_STONE_CIRCLE, MEGALITH_STAGE_DOLMEN, 777,
                         (Vector3){ 17.0f, 12.0f, 17.0f },
                         (Vector3){ 0.0f, 3.0f, 0.0f },
                         48.0f,
                         TextFormat("%s/stone_circle_cromlech.png", out_dir),
                         "Prehistoric Cromlech: Sacred Stone Circle",
                         "Radial Affordance Boundary Expansion & Central Altar Orthostat");

    // 5. Monumental Henge with Central Trilithons and Processional Avenue
    render_megalith_shot(MEGALITH_ARCHETYPE_STONE_CIRCLE, MEGALITH_STAGE_HENGE, 999,
                         (Vector3){ 0.0f, 8.2f, 24.0f },
                         (Vector3){ 0.0f, 3.2f, 0.0f },
                         44.0f,
                         TextFormat("%s/monumental_henge_trilithons.png", out_dir),
                         "Monumental Henge: Trilithons & Solstice Avenue",
                         "Paired Orthostats with Horizontal Lintels & Processional Way");

    CloseWindow();
    return 0;
}
