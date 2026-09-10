#include "raylib.h"
#include "../greek_grammar.h"
#include "../greek_rasterizer.h"
#include "../greek_grammar.c"
#include <stdio.h>

#define MAX_PREVIEW_VOXELS 32768
typedef struct {
    Vector3 pos;
    Color color;
} SimpleVoxel;

static SimpleVoxel preview_voxels[MAX_PREVIEW_VOXELS];
static int preview_voxel_count = 0;

static void plot_preview_voxel(int gx, int gy, int gz, Color c) {
    if (preview_voxel_count >= MAX_PREVIEW_VOXELS) return;
    float voxel_size = 0.5f;
    float floor_size = 20.0f;
    float px = (gx + 0.5f) * voxel_size - floor_size;
    float py = (gy + 0.5f) * voxel_size;
    float pz = (gz + 0.5f) * voxel_size - floor_size;
    preview_voxels[preview_voxel_count++] = (SimpleVoxel){ (Vector3){ px, py, pz }, c };
}

int main(int argc, char **argv) {
    SetConfigFlags(FLAG_WINDOW_HIDDEN);
    InitWindow(1280, 720, "Temple Entrance Renderer");
    if (!IsWindowReady()) return 1;

    TemplePlan plan = generate_greek_temple(1337, TEMPLE_STAGE_PERIPTERAL, 2);
    int center_g = 40;
    rasterize_temple_plan(&plan, center_g, center_g, 3, plot_preview_voxel);

    Camera3D camera = { 0 };
    // Eye-level view approaching the front entrance steps
    camera.position = (Vector3){ 0.0f, 4.0f, 22.0f };
    camera.target = (Vector3){ 0.0f, 4.0f, 0.0f };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = 55.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    RenderTexture2D target = LoadRenderTexture(1280, 720);
    BeginTextureMode(target);
        ClearBackground((Color){ 140, 190, 230, 255 });
        BeginMode3D(camera);
            DrawCube((Vector3){ 0.0f, 0.5f, 0.0f }, 40.0f, 1.0f, 40.0f, (Color){ 75, 105, 70, 255 });
            DrawCubeWires((Vector3){ 0.0f, 0.5f, 0.0f }, 40.0f, 1.0f, 40.0f, (Color){ 55, 80, 50, 255 });

            float vs = 0.5f;
            for (int i = 0; i < preview_voxel_count; ++i) {
                DrawCube(preview_voxels[i].pos, vs, vs, vs, preview_voxels[i].color);
                DrawCubeWires(preview_voxels[i].pos, vs, vs, vs, ColorBrightness(preview_voxels[i].color, -0.2f));
            }
        EndMode3D();
    EndTextureMode();

    Image img = LoadImageFromTexture(target.texture);
    ImageFlipVertical(&img);
    const char *out_path = (argc > 1) ? argv[1] : "temple_entrance.png";
    ExportImage(img, out_path);

    UnloadImage(img);
    UnloadRenderTexture(target);
    CloseWindow();
    return 0;
}
