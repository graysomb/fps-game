/* Uses the production mesh/shaders. Run from the repository root, with output
 * directory as argv[1]; optional argv[2] is --no-bloom. */
#include "../weapon_renderer.h"
#include "raymath.h"
#include "rlgl.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static WeaponBloom bloom;
static bool use_bloom = true;

static void save(RenderTexture2D scene, const char *directory, const char *name) {
    char path[1024]; snprintf(path,sizeof(path),"%s/%s.png",directory,name);
    Image image=LoadImageFromTexture(scene.texture); ImageFlipVertical(&image);
    assert(ExportImage(image,path)); UnloadImage(image);
}

static void render(RenderTexture2D scene, Camera3D camera, WeaponPose pose,
                   WeaponVisual *visual, bool gold, bool wall, bool first_person) {
    bool glow=use_bloom && weapon_bloom_resize(&bloom,scene.texture.width,scene.texture.height);
    BeginTextureMode(scene); ClearBackground((Color){29,55,94,255});
    BeginMode3D(camera);
    if (wall) DrawCube((Vector3){0,0,1},5,5,.3f,(Color){72,75,83,255});
    if (!first_person) weapon_draw(pose,visual,camera.position,1,gold,false,false);
    EndMode3D(); EndTextureMode();
    if (glow) {
        weapon_bloom_begin(&bloom,scene);
        if (first_person) weapon_clear_depth();
        BeginMode3D(camera); weapon_draw(pose,visual,camera.position,1,gold,false,true); EndMode3D();
        weapon_bloom_composite(&bloom,scene);
    } else BeginTextureMode(scene);
    if (first_person) {
        weapon_clear_depth(); BeginMode3D(camera);
        weapon_draw(pose,visual,camera.position,1,gold,false,false); EndMode3D();
    }
    EndTextureMode();
}

static void test_visual_state(void) {
    WeaponVisual v; weapon_visual_reset(&v);
    assert(weapon_shot_amount(&v,1)==0);
    weapon_visual_shot(&v,1); assert(v.shot_sequence==1);
    assert(weapon_shot_amount(&v,1)==1);
    assert(weapon_shot_amount(&v,1.2f)==0);
    weapon_visual_receive(&v,2,30,2); assert(fabsf(v.shot_time-1.97f)<.0001f);
    weapon_visual_receive(&v,2,0,3); assert(fabsf(v.shot_time-1.97f)<.0001f);
    weapon_visual_receive(&v,1,0,3); assert(v.shot_sequence==2);
    weapon_visual_receive(&v,3,255,3); assert(weapon_shot_amount(&v,3)==0);
    v.shot_sequence=UINT32_MAX; weapon_visual_receive(&v,0,0,4);
    assert(v.shot_sequence==0 && v.shot_time==4);
    weapon_visual_update(&v,true,.075f); assert(fabsf(v.claw_open-.5f)<.0001f);
    weapon_visual_update(&v,true,.075f); assert(v.claw_open==1);
    weapon_visual_update(&v,false,.15f); assert(v.claw_open==0);
}

static void test_melee_pose(void) {
    Vector3 forward={0,0,-1},right={1,0,0},up={0,1,0};
    for (int view=0;view<2;++view) {
        bool fp=view==0;
        WeaponPose idle=weapon_pose((Vector3){0},forward,right,up,fp,16.0f/9,NULL,0);
        for (int phase=0;phase<=100;++phase) {
            WeaponMeleePose melee={phase/100.0f,{-.29f,-.26f,-2.90f}};
            WeaponPose pose=weapon_pose((Vector3){0},forward,right,up,fp,16.0f/9,&melee,0);
            Matrix m=pose.transform;
            Vector3 x={m.m0,m.m1,m.m2},y={m.m4,m.m5,m.m6},z={m.m8,m.m9,m.m10};
            float scale=fp ? .36f : .56f;
            assert(fabsf(Vector3Length(x)-scale)<.00001f);
            assert(fabsf(Vector3Length(y)-scale)<.00001f);
            assert(fabsf(Vector3Length(z)-scale)<.00001f);
            assert(fabsf(Vector3DotProduct(x,y))<.00001f);
            assert(fabsf(Vector3DotProduct(y,z))<.00001f);
            if (melee.progress>=WEAPON_MELEE_ACTIVE_START && melee.progress<=WEAPON_MELEE_ACTIVE_END) {
                assert(Vector3Distance(pose.grip_heel,melee.strike_endpoint)<.00001f);
                WeaponPose firing=weapon_pose((Vector3){0},forward,right,up,fp,16.0f/9,&melee,1);
                assert(Vector3Distance(pose.muzzle,firing.muzzle)<.00001f);
            }
            if (phase==0 || phase==100) {
                assert(Vector3Distance(idle.muzzle,pose.muzzle)<.00001f);
                assert(Vector3Distance(idle.grip_heel,pose.grip_heel)<.00001f);
            }
        }
        const float knots[]={WEAPON_MELEE_ACTIVE_START,WEAPON_MELEE_PEAK,WEAPON_MELEE_ACTIVE_END};
        for (int i=0;i<3;++i) {
            WeaponMeleePose before={knots[i]-.00001f,{-.29f,-.26f,-2.90f}},after=before;
            after.progress+=.00002f;
            WeaponPose a=weapon_pose((Vector3){0},forward,right,up,fp,16.0f/9,&before,0);
            WeaponPose b=weapon_pose((Vector3){0},forward,right,up,fp,16.0f/9,&after,0);
            assert(Vector3Distance(a.muzzle,b.muzzle)<.0001f);
        }
    }
}

int main(int argc,char **argv) {
    const char *out=argc>1 ? argv[1] : "artifacts/weapon";
    use_bloom=argc<3 || strcmp(argv[2],"--no-bloom")!=0;
    test_visual_state(); test_melee_pose();
    SetConfigFlags(FLAG_WINDOW_HIDDEN); InitWindow(1280,720,"Weapon preview");
    assert(IsWindowReady()); assert(weapon_renderer_init());
    RenderTexture2D scene=LoadRenderTexture(1000,1000);
    WeaponVisual visual; weapon_visual_reset(&visual);
    Camera3D camera={.position={-2.0f,1.2f,-2.6f},.target={0,-.22f,-.10f},
                     .up={0,1,0},.fovy=42,.projection=CAMERA_PERSPECTIVE};
    WeaponPose pose={.transform=MatrixIdentity(),.muzzle={0,0,-.84f}};
    render(scene,camera,pose,&visual,false,false,false); save(scene,out,"reference-idle");
    Image visible=LoadImageFromTexture(scene.texture); Color *lit=LoadImageColors(visible);
    int cyan=0;
    for (int i=0;i<visible.width*visible.height;++i)
        if (lit[i].g>150 && lit[i].b>150 && lit[i].r<80) ++cyan;
    assert(cyan>100); UnloadImageColors(lit); UnloadImage(visible);
    visual.claw_open=1;
    render(scene,camera,pose,&visual,false,false,false); save(scene,out,"reference-tether");
    render(scene,camera,pose,&visual,true,false,false); save(scene,out,"reference-gold");
    visual.claw_open=0; weapon_visual_shot(&visual,1);
    render(scene,camera,pose,&visual,false,false,false); save(scene,out,"reference-firing");
    weapon_visual_reset(&visual);
    // Hidden gun must contribute no cyan pixels or bloom through an opaque wall.
    camera.position=(Vector3){0,0,3}; camera.target=(Vector3){0,0,0};
    render(scene,camera,pose,&visual,false,true,false); save(scene,out,"wall-occlusion");
    Image hidden=LoadImageFromTexture(scene.texture); Color *pixels=LoadImageColors(hidden);
    for (int i=0;i<hidden.width*hidden.height;++i) assert(!(pixels[i].g>150 && pixels[i].b>150 && pixels[i].r<80));
    UnloadImageColors(pixels); UnloadImage(hidden); UnloadRenderTexture(scene);
    // Aspect ratios used by 1, 2, and 3/4 local viewports, plus extreme pitch.
    for (int i=0;i<5;++i) {
        int width=i==1 ? 640 : 1280, height=720;
        if (i==2) width=640,height=360;
        scene=LoadRenderTexture(width,height);
        float pitch=i==3 ? 89*DEG2RAD : i==4 ? -89*DEG2RAD : 0;
        Vector3 forward={0,sinf(pitch),-cosf(pitch)},right={1,0,0};
        Vector3 up=Vector3CrossProduct(right,forward);
        camera=(Camera3D){.position={0},.target=forward,.up={0,1,0},.fovy=60,.projection=CAMERA_PERSPECTIVE};
        pose=weapon_pose((Vector3){0},forward,right,up,true,(float)width/height,NULL,0);
        render(scene,camera,pose,&visual,false,false,true);
        BeginTextureMode(scene);
        DrawLine(width/2-6,height/2,width/2+6,height/2,WHITE);
        DrawLine(width/2,height/2-6,width/2,height/2+6,WHITE);
        EndTextureMode();
        const char *names[]={"first-person","first-person-narrow","first-person-quarter","pitch-up","pitch-down"};
        save(scene,out,names[i]); UnloadRenderTexture(scene);
    }
    // Timings for four view-sized draws, independent of the physics backend.
    scene=LoadRenderTexture(640,360);
    double start=GetTime();
    for (int frame=0;frame<40;++frame) for (int view=0;view<4;++view)
        render(scene,camera,pose,&visual,false,false,true);
    rlDrawRenderBatchActive();
    Image sync=LoadImageFromTexture(scene.texture); UnloadImage(sync);
    printf("Four-view weapon rendering, bloom=%d: %.3f ms/frame\n",use_bloom,(GetTime()-start)*1000/40);
    UnloadRenderTexture(scene); weapon_bloom_unload(&bloom);
    weapon_renderer_shutdown(); CloseWindow(); return 0;
}
