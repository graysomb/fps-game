static RenderTexture2D clarityShowcase;
static Model clarityShowcaseCube;
static Shader clarityShowcaseShader;
static bool clarityShowcaseReady;
static void clarity_menu_cube(Vector3 pos,float x,float y,float z,Color color) {
    DrawModelEx(clarityShowcaseCube,pos,(Vector3){0,1,0},0,(Vector3){x,y,z},color);
    DrawCubeWires(pos,x,y,z,(Color){(unsigned char)(color.r*.45f),(unsigned char)(color.g*.45f),(unsigned char)(color.b*.45f),color.a});
}
static int clarityShowcaseEntry;
static float clarityShowcaseTimer;
static void clarity_menu_shutdown(void) {
    feedback_shutdown();
    if(clarityShowcaseReady) {
        UnloadModel(clarityShowcaseCube);
        if(clarityShowcaseShader.id && clarityShowcaseShader.id!=rlGetShaderIdDefault())UnloadShader(clarityShowcaseShader);
        clarityShowcaseReady=false;
    }
    if(clarityShowcase.id){UnloadRenderTexture(clarityShowcase);clarityShowcase=(RenderTexture2D){0};}
}
static void clarity_menu_render(int type,float yaw,float now) {
    if(!clarityShowcase.id)clarityShowcase=LoadRenderTexture(720,720);
    if(!clarityShowcaseReady){
        clarityShowcaseCube=LoadModelFromMesh(GenMeshCube(1,1,1));
        clarityShowcaseShader=LoadShader("shaders/character.vert","shaders/menu_preview.frag");
        if(clarityShowcaseShader.id && clarityShowcaseShader.id!=rlGetShaderIdDefault())
            clarityShowcaseCube.materials[0].shader=clarityShowcaseShader;
        clarityShowcaseReady=true;
    }
    BeginTextureMode(clarityShowcase);ClearBackground((Color){12,23,32,255});
    Camera3D cam={{4,2.8f,4},{0,1,0},{0,1,0},38,CAMERA_PERSPECTIVE};if(type<0){cam.position=(Vector3){3,2,3};cam.target.y=.7f;}
    if(type==ENEMY_TYPE_LOBBER){cam.position.y=3.4f;cam.target.y=1.35f;cam.fovy=42;}
    BeginMode3D(cam);
    DrawCylinder((Vector3){0,-.1f,0},1.5f,1.5f,.15f,24,DARKGRAY);
    DrawCircle3D((Vector3){0,.02f,0},1.4f,(Vector3){1,0,0},90,SKYBLUE);
    if(type<0){CharacterPose pose=character_pose((Vector3){0,BASE_EYE_HEIGHT,0},yaw,0,0,true,false,SKYBLUE,now);character_draw(pose,false);character_draw_idle_arms(pose);character_draw_visor(pose);}
    else {
        Color color=type>=ENEMY_TYPE_SCRAPPER?enemy_color(type):type==ENEMY_TYPE_SWARMER?RED:SKYBLUE;
        float size=type==ENEMY_TYPE_SWARMER?.5f:type==ENEMY_TYPE_FOUNDRY?1.5f:1;
        float height=type==ENEMY_TYPE_FOUNDRY?1:type==ENEMY_TYPE_GOLIATH?.75f:size*.5f;
        robot_draw(type,(Vector3){0,height,0},yaw,size,color,now,1,0,true);
        rlPushMatrix();rlRotatef(yaw,0,1,0);
        if(type==ENEMY_TYPE_GOLIATH)for(int n=0;n<GOLIATH_MAX_ATTACHED;n++)clarity_menu_cube(v_add((Vector3){0,.75f,0},goliathArmorOffsets[n]),.48f,.48f,.48f,DARKGRAY);
        if(type==ENEMY_TYPE_SHIELDBEARER){rlDrawRenderBatchActive();rlDisableDepthMask();clarity_menu_cube((Vector3){0,1,-.8f},2,2,.1f,Fade(SKYBLUE,.35f));rlDrawRenderBatchActive();rlEnableDepthMask();}
        if(type==ENEMY_TYPE_FOUNDRY)for(int x=-1;x<=1;x+=2)for(int z=-1;z<=1;z+=2)clarity_menu_cube((Vector3){x*1.2f,.7f,z*1.2f},.4f,1.4f,.4f,ORANGE);
        rlPopMatrix();
    }
    DrawCircle3D((Vector3){0,.1f+(type<0?1.4f:2.6f)*(.5f+.5f*sinf(now)),0},1.25f,(Vector3){1,0,0},90,Fade(SKYBLUE,.45f));EndMode3D();EndTextureMode();
}
static void clarity_menu(void) {
    static const int roster[]={-1,ENEMY_TYPE_SWARMER,ENEMY_TYPE_STANDARD,ENEMY_TYPE_GOLIATH,ENEMY_TYPE_SCRAPPER,ENEMY_TYPE_SHIELDBEARER,ENEMY_TYPE_SAPPER,ENEMY_TYPE_LOBBER,ENEMY_TYPE_BUILDER,ENEMY_TYPE_ANCHOR,ENEMY_TYPE_FOUNDRY};
    static const char *names[]={"RAT OPERATOR","SWARMER","RIFLEMAN","GOLIATH","SCRAPPER","SHIELDBEARER","SAPPER","LOBBER","CONSTRUCTOR","ANCHOR","WALKING FOUNDRY"};
    static const char *ability[]={"Build. Break. Survive.","Fast contact attacker.","Keeps pressure with gunfire.","Harvests material to rebuild armor.","Collects and weaponizes scrap.","Turns its shield toward threats.","Plants timed demolition charges.","Throws large tethered structures.","Builds cover and access ramps.","Heals and braces three nearby allies.","Heavy weapons and support legs."};
    static const char *counter[]={"Use cover and the environment.","Listen for skittering; keep moving.","Break sightlines; close the gap.","Strip its armor; interrupt harvesting.","Strip its scavenged armor.","Melee shield sections or flank.","Shoot the charge or leave its radius.","Watch the marker and dodge.","Interrupt construction; harvest cover.","Break the external anchor.","Break two supports, then attack core."};
    if(!clarityShowcase.id)clarityShowcase=LoadRenderTexture(720,720);
    clarityShowcaseTimer+=GetFrameTime();int direction=couch_axis(true);
    if(direction || clarityShowcaseTimer>=10){clarityShowcaseEntry=(clarityShowcaseEntry+(direction?direction:1)+11)%11;clarityShowcaseTimer=0;}
    float yaw=fmodf((float)GetTime()*15,360);int type=roster[clarityShowcaseEntry];
    clarity_menu_render(type,yaw,(float)GetTime());
    int x=GetScreenWidth()/2+12,w=GetScreenWidth()-x-28,y=85,size=(int)fminf(w,GetScreenHeight()-230);
    DrawTexturePro(clarityShowcase.texture,(Rectangle){0,0,720,-720},(Rectangle){x+(w-size)/2,y,size,size},(Vector2){0},0,WHITE);
    DrawText(TextFormat("%s  |  %s",names[clarityShowcaseEntry],type<0?"NEW INHABITANT":robot_roles[type]),x,y+size+12,14,SKYBLUE);DrawText(ability[clarityShowcaseEntry],x,y+size+43,16,RAYWHITE);DrawText(counter[clarityShowcaseEntry],x,y+size+67,16,LIGHTGRAY);DrawText("Left / Right: browse scans",x,y+size+94,14,GRAY);
}
