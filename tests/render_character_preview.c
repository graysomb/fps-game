#include "../character_renderer.h"
#include "raymath.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

static void save(RenderTexture2D rt,const char *directory,const char *name) {
    char path[1024];snprintf(path,sizeof(path),"%s/%s.png",directory,name);
    Image im=LoadImageFromTexture(rt.texture);ImageFlipVertical(&im);
    assert(ExportImage(im,path));UnloadImage(im);
}
int main(int argc,char **argv) {
    const char *out=argc>1?argv[1]:"artifacts/rat";
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,720,"Rat preview");SetTargetFPS(0);
    assert(IsWindowReady());assert(character_renderer_init());assert(weapon_renderer_init());
    Color colors[]={SKYBLUE,RED,GREEN,GOLD};
    RenderTexture2D rt=LoadRenderTexture(960,720);
    Camera3D camera={.position={1,.65f,-1.4f},.target={0,.55f,0},.up={0,1,0},.fovy=38};
    for(int i=0;i<16;i++) {
        float angle=(i%4)*PI/2;
        if(i%4==3)angle=-PI/4;
        camera.position=(Vector3){sinf(angle)*3.0f,1.4f,-cosf(angle)*3.0f};
        CharacterPose pose=character_pose((Vector3){0,1,0},0,0,i/4==1?3:0,i/4!=2,i/4==3,colors[i/4],.72f);
        BeginTextureMode(rt);ClearBackground((Color){214,216,218,255});BeginMode3D(camera);
        DrawPlane((Vector3){0,0,0},(Vector2){8,8},(Color){190,193,196,255});
        character_draw(pose,false);
        if(!pose.downed) for(int s=-1;s<=1;s+=2)
            character_draw_arm(Vector3Transform((Vector3){s*.17f,.11f,0},character_body_transform(pose)),Vector3Transform((Vector3){s*.18f,.03f,-.16f},character_body_transform(pose)),pose.gear);
        character_draw_visor(pose);
        EndMode3D();DrawText("RAT / flat-shaded procedural mesh",20,20,22,DARKGRAY);EndTextureMode();
        char name[64];snprintf(name,sizeof(name),"pose-%d-view-%d",i/4,i%4);save(rt,out,name);
    }
    camera=(Camera3D){.target={0,0,-1},.up={0,1,0},.fovy=60};
    WeaponVisual visual={0};weapon_visual_reset(&visual);
    for(int quarter=0;quarter<2;quarter++) for(int kind=0;kind<2;kind++) for(int action=0;action<3;action++) {
        RenderTexture2D fp=LoadRenderTexture(quarter?480:960,quarter?360:720);
        WeaponMeleePose melee={.progress=action==1?.40f:-1,.strike_endpoint={.15f,-.2f,-2.7f}};
        WeaponPose pose=weapon_pose(kind,(Vector3){0},(Vector3){0,0,-1},(Vector3){1,0,0},(Vector3){0,1,0},true,4.0f/3,&melee,0);
        visual.claw_open=action==2?1:0;
        BeginTextureMode(fp);ClearBackground((Color){55,72,86,255});BeginMode3D(camera);
        weapon_draw(pose,&visual,camera.position,1,action==2,WEAPON_COLOR);
        character_draw_paws(pose,SKYBLUE,false);
        weapon_draw(pose,&visual,camera.position,1,action==2,WEAPON_GLASS);
        EndMode3D();DrawCircle(fp.texture.width/2,fp.texture.height/2,2,WHITE);EndTextureMode();
        char name[64];snprintf(name,sizeof(name),"paws-%d-weapon-%d-action-%d",quarter,kind,action);save(fp,out,name);
        UnloadRenderTexture(fp);
    }
    // A fully occluding wall must also hide the scouter rim in the emission mask.
    WeaponBloom bloom={0};assert(weapon_bloom_resize(&bloom,960,720));
    camera=(Camera3D){.position={0,.1f,-2},.target={0,.1f,0},.up={0,1,0},.fovy=45};
    CharacterPose rat=character_pose((Vector3){0},0,0,0,true,false,SKYBLUE,1);
    BeginTextureMode(rt);ClearBackground(BLACK);BeginMode3D(camera);
    character_draw(rat,false);DrawCube((Vector3){0,0,-.8f},2,2,.1f,GRAY);
    EndMode3D();EndTextureMode();weapon_bloom_begin(&bloom,rt);
    BeginMode3D(camera);character_draw(rat,true);EndMode3D();EndTextureMode();
    Image mask=LoadImageFromTexture(bloom.emission.texture);Color *pixels=LoadImageColors(mask);
    for(int i=0;i<mask.width*mask.height;i++)assert(pixels[i].r==0 && pixels[i].g==0 && pixels[i].b==0);
    UnloadImageColors(pixels);UnloadImage(mask);weapon_bloom_unload(&bloom);
    UnloadRenderTexture(rt);character_renderer_shutdown();character_renderer_shutdown();weapon_renderer_shutdown();CloseWindow();
    puts("Character previews passed");return 0;
}
