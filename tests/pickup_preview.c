#include "../weapon_renderer.h"
#include "rlgl.h"
#include <assert.h>

int main(void) {
    SetConfigFlags(FLAG_WINDOW_HIDDEN);
    InitWindow(960,720,"Pickup previews");
    assert(weapon_renderer_init());
    Camera3D camera={{2.6f,2.2f,3.6f},{0,.9f,0},{0,1,0},40,CAMERA_PERSPECTIVE};
    for(int small=0;small<2;small++)for(int glow=0;glow<2;glow++) {
        int width=small?480:960,height=small?360:720;
        RenderTexture2D scene=LoadRenderTexture(width,height);
        WeaponBloom bloom={0};
        Vector3 positions[]={{-.65f,1,0},{.65f,1,0}};
        BeginDrawing();BeginTextureMode(scene);ClearBackground((Color){18,25,35,255});
        BeginMode3D(camera);DrawPlane((Vector3){0},(Vector2){10,10},DARKGRAY);
        for(int n=0;n<2;n++)weapon_draw_pickup(n==0,positions[n],camera.position,1,WEAPON_COLOR);
        for(int n=0;n<2;n++)weapon_draw_pickup(n==0,positions[n],camera.position,1,WEAPON_GLASS);
        EndMode3D();
        if(glow) {
            EndTextureMode();assert(weapon_bloom_resize(&bloom,width,height));
            weapon_bloom_begin(&bloom,scene);BeginMode3D(camera);
            for(int n=0;n<2;n++)weapon_draw_pickup(n==0,positions[n],camera.position,1,WEAPON_EMISSION);
            EndMode3D();weapon_bloom_composite(&bloom,scene);
        }
        EndTextureMode();EndDrawing();
        Image capture=LoadImageFromTexture(scene.texture);ImageFlipVertical(&capture);
        assert(ExportImage(capture,TextFormat("artifacts/pickups-%s-%s.png",small?"split":"full",glow?"bloom":"plain")));
        UnloadImage(capture);weapon_bloom_unload(&bloom);UnloadRenderTexture(scene);
    }
    /* Opaque cover must suppress pickup emission, not just their color pass. */
    RenderTexture2D covered=LoadRenderTexture(480,360);WeaponBloom mask={0};
    camera.position=(Vector3){0,1,4};camera.target=(Vector3){0,1,0};
    BeginDrawing();BeginTextureMode(covered);ClearBackground(BLACK);BeginMode3D(camera);
    DrawCube((Vector3){0,1,1},4,4,.2f,GRAY);EndMode3D();EndTextureMode();
    assert(weapon_bloom_resize(&mask,480,360));weapon_bloom_begin(&mask,covered);
    BeginMode3D(camera);
    weapon_draw_pickup(true,(Vector3){-.6f,1,0},camera.position,1,WEAPON_EMISSION);
    weapon_draw_pickup(false,(Vector3){.6f,1,0},camera.position,1,WEAPON_EMISSION);
    EndMode3D();weapon_bloom_composite(&mask,covered);EndTextureMode();EndDrawing();
    Image hidden=LoadImageFromTexture(mask.emission.texture);Color *pixels=LoadImageColors(hidden);
    for(int n=0;n<hidden.width*hidden.height;n++)assert(!pixels[n].r && !pixels[n].g && !pixels[n].b);
    UnloadImageColors(pixels);UnloadImage(hidden);weapon_bloom_unload(&mask);UnloadRenderTexture(covered);
    weapon_renderer_shutdown();CloseWindow();return 0;
}
