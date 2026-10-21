/* Production renderer preview/regression, including existing weapon checks. */
#define main existing_weapon_preview_main
#include "render_weapon_preview.c"
#undef main
int main(int argc,char **argv){
    const char *out=argc>1?argv[1]:"artifacts/weapon-concepts/rendered";
    const char *names[]={"shotgun","cutter","adhesive","railgun","rocket","sword","builder"};
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,720,"Weapon concepts");
    assert(weapon_renderer_init());
    RenderTexture2D gallery=LoadRenderTexture(1440,780);
    BeginTextureMode(gallery);ClearBackground((Color){29,55,94,255});EndTextureMode();
    Camera3D camera={{-2.8f,1.5f,-3.5f},{0,-.1f,-.4f},{0,1,0},42,CAMERA_PERSPECTIVE};
    for(int k=WEAPON_IMPULSE_SHOTGUN;k<=WEAPON_REINFORCED_BUILDER;k++){
        WeaponVisual visual;weapon_visual_reset(&visual);visual.charge=.8f;visual.activity=1;visual.deployment=1;
        WeaponPose pose={.kind=k,.transform=MatrixIdentity()};
        RenderTexture2D scene=LoadRenderTexture(900,900);
        char name[80];
        for(int glow=0;glow<2;glow++){
            use_bloom=glow;render_time=.2f;render(scene,camera,pose,&visual,false,false,false);
            snprintf(name,sizeof(name),"%s-%s",names[k-2],glow?"bloom":"plain");save(scene,out,name);
            if(glow){
                int col=(k-2)%4,row=(k-2)/4;
                BeginTextureMode(gallery);
                DrawTexturePro(scene.texture,(Rectangle){0,0,900,-900},(Rectangle){col*360,row*390,360,360},(Vector2){0},0,WHITE);
                DrawText(names[k-2],col*360+20,row*390+350,22,WHITE);EndTextureMode();
            }
        }
        use_bloom=false;render_time=.2f;visual.charge=0;visual.activity=0;visual.deployment=0;
        render(scene,camera,pose,&visual,false,false,false);Image a=LoadImageFromTexture(scene.texture);Color *pa=LoadImageColors(a);
        render_time=.73f;visual.charge=1;visual.activity=1;visual.deployment=1;visual.attached=1;visual.cooldown=1;visual.shot_time=.72f;render(scene,camera,pose,&visual,false,false,false);
        Image b=LoadImageFromTexture(scene.texture);Color *pb=LoadImageColors(b);int changed=0;
        for(int p=0;p<a.width*a.height;p++)if(abs(pa[p].r-pb[p].r)+abs(pa[p].g-pb[p].g)+abs(pa[p].b-pb[p].b)>20)changed++;
        assert(changed>30);UnloadImageColors(pa);UnloadImageColors(pb);UnloadImage(a);UnloadImage(b);
        UnloadRenderTexture(scene);visual.cooldown=0;visual.shot_time=-1000;
        for(int small=0;small<2;small++){
            int w=small?640:1280,h=small?360:720;scene=LoadRenderTexture(w,h);
            Camera3D fp={{0},{0,0,-1},{0,1,0},60,CAMERA_PERSPECTIVE};
            pose=weapon_pose(k,(Vector3){0},(Vector3){0,0,-1},(Vector3){1,0,0},(Vector3){0,1,0},true,(float)w/h,NULL,0);
            assert(isfinite(pose.muzzle.x)&&isfinite(pose.muzzle.y)&&isfinite(pose.muzzle.z));
            render(scene,fp,pose,&visual,false,false,true);
            snprintf(name,sizeof(name),"%s-first-person%s",names[k-2],small?"-quarter":"");save(scene,out,name);
            UnloadRenderTexture(scene);
        }
    }
    save(gallery,out,"all-weapons");UnloadRenderTexture(gallery);
    weapon_bloom_unload(&bloom);weapon_renderer_shutdown();assert(weapon_renderer_init());weapon_renderer_shutdown();CloseWindow();
    puts("Seven weapon visuals: animation, bloom, first-person and cleanup passed");return 0;
}
