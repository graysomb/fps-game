#ifndef GOLEM_BODY_H
#define GOLEM_BODY_H
#include "golem_renderer.h"
#include "raymath.h"
/* One authoritative cell definition for rendering and physical construction. */
enum { GM_CHEST=GOLEM_JOINT_COUNT, GM_CORE, GM_EYE, GM_TURRET_L, GM_TURRET_R, GM_COUNT };
typedef void (*GolemCellCallback)(Vector3 center,Vector3 size,Color color);
typedef void (*GolemPartCallback)(int part);
static inline void go_body_grid(GolemCellCallback cube,int nx,int ny,int nz,Vector3 bottom,Color color){
 for(int y=ny-1;y>=0;y--)for(int z=0;z<nz;z++)for(int x=0;x<nx;x++)
  cube(Vector3Add(bottom,(Vector3){(x-(nx-1)*.5f)*.26f,(y+.5f)*.26f,(z-(nz-1)*.5f)*.26f}),(Vector3){.254f,.254f,.254f},color);
}
static inline void golem_body_generate(GolemCellCallback cube,GolemPartCallback finish){
 const Color grSteel={119,137,150,255},grDark={47,57,66,255},grOrange={204,121,44,255};
    go_body_grid(cube,4,2,2,(Vector3){0,-.25f,0},grDark);finish(GOLEM_PELVIS);
    for(int y=3;y>=0;y--)for(int z=0;z<3;z++)for(int x=0;x<6;x++){
        if(x>=1&&x<=4&&y>=1&&y<=2&&z<2)continue;
        cube((Vector3){(x-2.5f)*.26f,(y+.5f)*.26f,(z-1)*.26f},(Vector3){.254f,.254f,.254f},grSteel);
    }finish(GOLEM_TORSO);
    go_body_grid(cube,3,2,2,(Vector3){0,0,0},grSteel);cube((Vector3){0,-.13f,0},(Vector3){.35f,.25f,.35f},grDark);finish(GOLEM_HEAD);
    for(int s=0;s<2;s++){
        int a=s?GOLEM_ARM_R:GOLEM_ARM_L;
        go_body_grid(cube,2,2,2,(Vector3){0,-.72f,0},grSteel);go_body_grid(cube,3,2,3,(Vector3){0,-.19f,0},grOrange);finish(a);
        go_body_grid(cube,3,2,3,(Vector3){0,-.64f,0},grOrange);cube((Vector3){0,-.05f,0},(Vector3){.25f,.22f,.25f},grDark);finish(a+1);
        for(int y=0;y<2;y++)for(int z=0;z<2;z++)for(int x=0;x<3;x++)
            cube((Vector3){(x-1)*.26f,-(y+.5f)*.195f,(z-.5f)*.26f-.04f},(Vector3){.254f,.189f,.254f},grDark);
        finish(a+2);
        int l=s?GOLEM_THIGH_R:GOLEM_THIGH_L;
        go_body_grid(cube,2,3,2,(Vector3){0,-.78f,0},grSteel);finish(l);
        go_body_grid(cube,2,3,2,(Vector3){0,-.81f,0},grSteel);cube((Vector3){0,0,0},(Vector3){.54f,.18f,.54f},grOrange);finish(l+1);
        go_body_grid(cube,3,1,4,(Vector3){0,-.217f,-.15f},grDark);finish(l+2);
    }
    go_body_grid(cube,4,2,1,(Vector3){0,.26f,-.52f},grSteel);finish(GM_CHEST);
    cube((Vector3){0,.52f,-.06f},(Vector3){.40f,.40f,.40f},(Color){255,191,52,255});finish(GM_CORE);
    cube((Vector3){0,.30f,-.268f},(Vector3){.48f,.09f,.035f},(Color){255,170,31,255});finish(GM_EYE);
}
#endif
