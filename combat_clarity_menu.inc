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
        weapon_renderer_init();
        clarityShowcaseShader=LoadShader("shaders/menu_preview.vert","shaders/menu_preview.frag");
        if(clarityShowcaseShader.id && clarityShowcaseShader.id!=rlGetShaderIdDefault())
            clarityShowcaseCube.materials[0].shader=clarityShowcaseShader;
        clarityShowcaseReady=true;
    }
    BeginTextureMode(clarityShowcase);ClearBackground((Color){12,23,32,255});
    Camera3D cam={{4,2.8f,4},{0,1,0},{0,1,0},38,CAMERA_PERSPECTIVE};if(type<0){cam.position=(Vector3){3,2,3};cam.target.y=.7f;}
    if(type==ENEMY_TYPE_FOUNDRY){cam.position=(Vector3){5,3.6f,5};cam.target.y=1.2f;cam.fovy=44;}
    if(type==ENEMY_TYPE_LOBBER){cam.position.y=3.4f;cam.target.y=1.35f;cam.fovy=42;}
    BeginMode3D(cam);
    DrawCylinder((Vector3){0,-.1f,0},1.5f,1.5f,.15f,24,DARKGRAY);
    DrawCircle3D((Vector3){0,.02f,0},1.4f,(Vector3){1,0,0},90,SKYBLUE);
    if(type<0){CharacterPose pose=character_pose((Vector3){0,BASE_EYE_HEIGHT,0},yaw,0,0,true,false,SKYBLUE,now);character_draw(pose,false);character_draw_idle_arms(pose);character_draw_visor(pose);}
    else {rlPushMatrix();rlRotatef(yaw,0,1,0);Color color=type>=ENEMY_TYPE_SCRAPPER?enemy_color(type):type==ENEMY_TYPE_SWARMER?RED:GRAY;float size=type==ENEMY_TYPE_SWARMER?.5f:type==ENEMY_TYPE_FOUNDRY?1.5f:1;clarity_menu_cube((Vector3){0,1,0},size,size,size,color);
        Player preview={0};preview.enemyType=type;
        draw_player_face(&preview,(Vector3){0,1,0},size,color);
        enemy_draw_details(type,(Vector3){0,1,0},color);
        if(type==ENEMY_TYPE_STANDARD){
            WeaponMeleePose melee={.progress=-1};WeaponVisual visual={0};
            WeaponPose gun=weapon_pose(WEAPON_STANDARD,(Vector3){0,1,0},(Vector3){0,0,-1},(Vector3){1,0,0},(Vector3){0,1,0},false,1,&melee,0);
            weapon_draw(gun,&visual,cam.position,now,false,WEAPON_COLOR);
        }
        if(type==ENEMY_TYPE_SWARMER)clarity_menu_cube((Vector3){0,1.25f,0},.55f,.08f,.55f,YELLOW);
        if(type==ENEMY_TYPE_GOLIATH)for(int n=0;n<GOLIATH_MAX_ATTACHED;n++)clarity_menu_cube(v_add((Vector3){0,1,0},goliathArmorOffsets[n]),.48f,.48f,.48f,DARKGRAY);
        EnemyEquipment layout[ENEMY_EQUIPMENT_MAX];int equipmentCount=enemy_equipment_layout(type,layout);
        if(type==ENEMY_TYPE_SHIELDBEARER){rlDrawRenderBatchActive();rlDisableDepthMask();}
        for(int n=0;n<equipmentCount;n++)clarity_menu_cube(v_add((Vector3){0,1,0},layout[n].offset),.48f,.48f,.48f,
            type==ENEMY_TYPE_SHIELDBEARER?Fade(SKYBLUE,.35f):type==ENEMY_TYPE_FOUNDRY && layout[n].group<4?ORANGE:color);
        if(type==ENEMY_TYPE_SHIELDBEARER){rlDrawRenderBatchActive();rlEnableDepthMask();}
        if(type==ENEMY_TYPE_LOBBER)for(int x=0;x<2;x++)for(int z=0;z<2;z++)clarity_menu_cube((Vector3){x*.5f-.25f,2.5f,z*.5f-.25f},.49f,.49f,.49f,GOLD);
        if(type==ENEMY_TYPE_ANCHOR)DrawCylinder((Vector3){0,1.5f,0},.1f,.1f,1,8,SKYBLUE);
        rlPopMatrix();}
    DrawCircle3D((Vector3){0,.1f+(type<0?1.4f:2.6f)*(.5f+.5f*sinf(now)),0},1.25f,(Vector3){1,0,0},90,Fade(SKYBLUE,.45f));EndMode3D();EndTextureMode();
}
static void clarity_menu(void) {
    static const int roster[]={-1,ENEMY_TYPE_SWARMER,ENEMY_TYPE_STANDARD,ENEMY_TYPE_GOLIATH,ENEMY_TYPE_SHIELDBEARER,ENEMY_TYPE_SAPPER,ENEMY_TYPE_LOBBER,ENEMY_TYPE_BUILDER,ENEMY_TYPE_ANCHOR,ENEMY_TYPE_FOUNDRY};
    static const char *names[]={"RAT OPERATOR","SWARMER","RIFLEMAN","GOLIATH","SHIELDBEARER","SAPPER","LOBBER","CONSTRUCTOR","ANCHOR","WALKING FOUNDRY"};
    static const char *ability[]={"Build. Break. Survive.","Fast contact attacker.","Keeps pressure with gunfire.","Harvests material to rebuild armor.","Turns its shield toward threats.","Plants timed demolition charges.","Throws large tethered structures.","Builds cover and access ramps.","Heals and braces three nearby allies.","Heavy weapons and support legs."};
    static const char *counter[]={"Use cover and the environment.","Listen for skittering; keep moving.","Break sightlines; close the gap.","Melee or throw debris to strip armor; bullets bounce.","Melee shield sections or flank.","Shoot the charge or leave its radius.","Watch the marker and dodge.","Interrupt construction; harvest cover.","Break the external anchor.","Break two supports, then attack core."};
    if(!clarityShowcase.id)clarityShowcase=LoadRenderTexture(720,720);
    clarityShowcaseTimer+=GetFrameTime();int direction=couch_axis(true);
    if(direction || clarityShowcaseTimer>=10){clarityShowcaseEntry=(clarityShowcaseEntry+(direction?direction:1)+10)%10;clarityShowcaseTimer=0;}
    float yaw=fmodf((float)GetTime()*15,360);int type=roster[clarityShowcaseEntry];
    clarity_menu_render(type,yaw,(float)GetTime());
    int x=GetScreenWidth()/2+12,w=GetScreenWidth()-x-28,y=85,size=(int)fminf(w,GetScreenHeight()-230);
    DrawTexturePro(clarityShowcase.texture,(Rectangle){0,0,720,-720},(Rectangle){x+(w-size)/2,y,size,size},(Vector2){0},0,WHITE);
    DrawText(names[clarityShowcaseEntry],x,y+size+12,22,SKYBLUE);DrawText(ability[clarityShowcaseEntry],x,y+size+43,16,RAYWHITE);DrawText(counter[clarityShowcaseEntry],x,y+size+67,16,LIGHTGRAY);DrawText("Left / Right: browse scans",x,y+size+94,14,GRAY);
}
