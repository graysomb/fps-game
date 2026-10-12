static RenderTexture2D clarityShowcase;
static int clarityShowcaseEntry;
static float clarityShowcaseTimer;
static void clarity_menu_shutdown(void){if(clarityShowcase.id){UnloadRenderTexture(clarityShowcase);clarityShowcase=(RenderTexture2D){0};}}
static void clarity_menu(void) {
    static const int roster[]={-1,ENEMY_TYPE_SWARMER,ENEMY_TYPE_STANDARD,ENEMY_TYPE_GOLIATH,ENEMY_TYPE_SHIELDBEARER,ENEMY_TYPE_SAPPER,ENEMY_TYPE_LOBBER,ENEMY_TYPE_BUILDER,ENEMY_TYPE_ANCHOR,ENEMY_TYPE_FOUNDRY};
    static const char *names[]={"RAT OPERATOR","SWARMER","RIFLEMAN","GOLIATH","SHIELDBEARER","SAPPER","LOBBER","CONSTRUCTOR","ANCHOR","WALKING FOUNDRY"};
    static const char *ability[]={"Build. Break. Survive.","Fast contact attacker.","Keeps pressure with gunfire.","Harvests material to rebuild armor.","Turns its shield toward threats.","Plants timed demolition charges.","Throws large tethered structures.","Builds cover and access ramps.","Braces nearby allies.","Heavy weapons and support legs."};
    static const char *counter[]={"Use cover and the environment.","Listen for skittering; keep moving.","Break sightlines; close the gap.","Strip its armor; interrupt harvesting.","Melee shield sections or flank.","Shoot the charge or leave its radius.","Watch the marker and dodge.","Interrupt construction; harvest cover.","Break the external anchor.","Break two supports, then attack core."};
    if(!clarityShowcase.id)clarityShowcase=LoadRenderTexture(720,720);
    clarityShowcaseTimer+=GetFrameTime();int direction=couch_axis(true);
    if(direction || clarityShowcaseTimer>=10){clarityShowcaseEntry=(clarityShowcaseEntry+(direction?direction:1)+10)%10;clarityShowcaseTimer=0;}
    float yaw=fmodf((float)GetTime()*15,360);int type=roster[clarityShowcaseEntry];
    BeginTextureMode(clarityShowcase);ClearBackground((Color){12,23,32,255});
    Camera3D cam={{4,2.8f,4},{0,1,0},{0,1,0},38,CAMERA_PERSPECTIVE};if(type<0){cam.position=(Vector3){3,2,3};cam.target.y=.7f;}BeginMode3D(cam);
    DrawCylinder((Vector3){0,-.1f,0},1.5f,1.5f,.15f,24,DARKGRAY);
    DrawCircle3D((Vector3){0,.02f,0},1.4f,(Vector3){1,0,0},90,SKYBLUE);
    if(type<0){CharacterPose pose=character_pose((Vector3){0,BASE_EYE_HEIGHT,0},yaw,0,0,true,false,SKYBLUE,(float)GetTime());character_draw(pose,false);character_draw_visor(pose);}
    else {rlPushMatrix();rlRotatef(yaw,0,1,0);Color color=type>=ENEMY_TYPE_SCRAPPER?enemy_color(type):type==ENEMY_TYPE_SWARMER?RED:GRAY;float size=type==ENEMY_TYPE_SWARMER?.5f:type==ENEMY_TYPE_FOUNDRY?1.5f:1;DrawCube((Vector3){0,1,0},size,size,size,color);
        if(type==ENEMY_TYPE_GOLIATH)for(int n=0;n<GOLIATH_MAX_ATTACHED;n++)DrawCube(v_add((Vector3){0,1,0},goliathArmorOffsets[n]),.48f,.48f,.48f,DARKGRAY);
        if(type==ENEMY_TYPE_SHIELDBEARER){DrawCube((Vector3){0,1,-.8f},2,2,.1f,Fade(SKYBLUE,.35f));DrawCubeWires((Vector3){0,1,-.8f},2,2,.12f,BLUE);}
        if(type==ENEMY_TYPE_FOUNDRY)for(int x=-1;x<=1;x+=2)for(int z=-1;z<=1;z+=2)DrawCube((Vector3){x*1.2f,.7f,z*1.2f},.4f,1.4f,.4f,ORANGE);
        if(type==ENEMY_TYPE_BUILDER)DrawCube((Vector3){0,1.6f,0},1.2f,.2f,1.2f,YELLOW);
        if(type==ENEMY_TYPE_SAPPER)DrawCube((Vector3){0,1.5f,.5f},.5f,.5f,.5f,ORANGE);
        if(type==ENEMY_TYPE_LOBBER)DrawCubeWires((Vector3){0,2.5f,0},1,1,1,GOLD);
        if(type==ENEMY_TYPE_ANCHOR)DrawCylinder((Vector3){0,1.5f,0},.1f,.1f,1,8,SKYBLUE);
        rlPopMatrix();}
    DrawCircle3D((Vector3){0,.1f+(type<0?1.4f:2.6f)*(.5f+.5f*sinf((float)GetTime())),0},1.25f,(Vector3){1,0,0},90,Fade(SKYBLUE,.45f));EndMode3D();EndTextureMode();
    int x=GetScreenWidth()/2+12,w=GetScreenWidth()-x-28,y=85,size=(int)fminf(w,GetScreenHeight()-230);
    DrawTexturePro(clarityShowcase.texture,(Rectangle){0,0,720,-720},(Rectangle){x,y,w,size},(Vector2){0},0,WHITE);
    DrawText(names[clarityShowcaseEntry],x,y+size+12,22,SKYBLUE);DrawText(ability[clarityShowcaseEntry],x,y+size+43,16,RAYWHITE);DrawText(counter[clarityShowcaseEntry],x,y+size+67,16,LIGHTGRAY);DrawText("Left / Right: browse scans",x,y+size+94,14,GRAY);
}
