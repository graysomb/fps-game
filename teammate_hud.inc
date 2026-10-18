/* Screen-space markers deliberately remain visible through world geometry. */
static bool teammate_marker_visible(int viewer,int teammate) {
    return gameMode==GAME_MODE_FIREFIGHT && teammate!=viewer &&
        teammate<activePlayers && player_uses_rat(teammate) &&
        combatantTeam[viewer]==combatantTeam[teammate] &&
        (combatantLife[teammate]==LIFE_DOWNED ||
         (combatantLife[teammate]==LIFE_ALIVE && players[teammate].respawn_timer<=0));
}

static Vector2 teammate_marker_position(Camera3D camera,Vector3 target,int w,int h,bool *edge) {
    Vector3 forward=v_norm(v_sub(camera.target,camera.position));
    Vector3 right=v_norm(Vector3CrossProduct(forward,camera.up));
    Vector3 up=Vector3CrossProduct(right,forward);
    Vector3 delta=v_sub(target,camera.position);
    float depth=v_dot(delta,forward);
    float focal=h*.5f/tanf(camera.fovy*DEG2RAD*.5f);
    Vector2 direction={v_dot(delta,right)*focal/fmaxf(fabsf(depth),.01f),
                       -v_dot(delta,up)*focal/fmaxf(fabsf(depth),.01f)};
    float cx=w*.5f,cy=h*.5f;
    float rx=fmaxf(1,cx-65),ry=fmaxf(1,cy-85);
    *edge=depth<=0 || fabsf(direction.x)>rx || fabsf(direction.y)>ry;
    if(*edge) {
        // A target exactly behind the camera gets a stable downward arrow.
        if(fabsf(direction.x)+fabsf(direction.y)<.01f)direction.y=1;
        float scale=fmaxf(fabsf(direction.x)/rx,fabsf(direction.y)/ry);
        direction=Vector2Scale(direction,1/scale);
    }
    return (Vector2){cx+direction.x,cy+direction.y};
}

static void draw_teammate_markers(int viewer,Camera3D camera,int w,int h) {
    int font=h<500?12:16;
    for(int teammate=0;teammate<MAX_PLAYERS;teammate++) {
        if(!teammate_marker_visible(viewer,teammate))continue;
        bool downed=combatantLife[teammate]==LIFE_DOWNED,edge;
        Vector3 target=world_image(players[teammate].pos,camera.position);
        target.y+=downed?0:.65f;
        Vector2 point=teammate_marker_position(camera,target,w,h,&edge);
        Color color=player_palette_color(teammate);
        DrawCircleV(point,9,Fade(BLACK,.8f));
        DrawCircleLines((int)point.x,(int)point.y,8,color);
        if(edge) {
            Vector2 dir=Vector2Normalize(Vector2Subtract(point,(Vector2){w*.5f,h*.5f}));
            Vector2 tip=Vector2Add(point,Vector2Scale(dir,6));
            Vector2 base=Vector2Subtract(point,Vector2Scale(dir,3));
            Vector2 side={-dir.y*4,dir.x*4};
            DrawLineEx(tip,Vector2Add(base,side),2,color);
            DrawLineEx(tip,Vector2Subtract(base,side),2,color);
        } else DrawCircleV(point,3,color);
        float distance=v_length(world_delta(players[teammate].pos,players[viewer].pos));
        const char *label=TextFormat(downed?"P%d DOWN | %.0fm":"P%d | %.0fm",teammate+1,distance);
        int width=MeasureText(label,font),x=(int)point.x-width/2,y=(int)point.y+12;
        DrawRectangle(x-4,y-2,width+8,font+4,Fade(BLACK,.75f));
        DrawText(label,x,y,font,downed?ORANGE:color);
    }
}
