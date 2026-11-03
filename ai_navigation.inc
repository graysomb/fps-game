/* Bounded, cached ground routes. Only the host searches; no world navmesh rebuild.
 * Coordinates are local to the search origin so periodic seams are ordinary edges.
 * The collision sweep uses the actor's actual bounds, including large enemies. */
#define AI_NAV_SIDE 49
#define AI_NAV_CELLS (AI_NAV_SIDE*AI_NAV_SIDE)
#define AI_NAV_PATH 128
typedef struct {Vector3 points[AI_NAV_PATH],goal;int count,next;float retry;} AiNavigation;
static AiNavigation aiNavigation[MAX_COMBATANTS];
static int aiNavParent[AI_NAV_CELLS],aiNavQueue[AI_NAV_CELLS];
static Vector3 ai_nav_point(Vector3 origin,int cell) {
    origin.x+=cell%AI_NAV_SIDE-AI_NAV_SIDE/2;
    origin.z+=cell/AI_NAV_SIDE-AI_NAV_SIDE/2;
    return world_wrap(origin);
}
static bool ai_nav_clear(int slot,Vector3 from,Vector3 to) {
    Vector3 delta=world_delta(to,from);delta.y=0;
    int steps=(int)ceilf(v_length(delta)/.3f);if(steps<1)steps=1;
    for(int n=1;n<=steps;n++) {
        Vector3 p=world_wrap(v_add(from,v_mul(delta,(float)n/steps)));
        if(!worldPeriodic && (fabsf(p.x)>FLOOR_SIZE-.5f || fabsf(p.z)>FLOOR_SIZE-.5f))return false;
        if(player_aabb_hits_world_ex(p,player_collision_half(&players[slot]),slot))return false;
    }
    return true;
}
static void ai_nav_search(int slot,Vector3 goal) {
    AiNavigation *route=&aiNavigation[slot];Vector3 origin=players[slot].pos;
    route->count=route->next=0;route->goal=goal;route->retry=1.f+(slot%7)*.07f;
    for(int n=0;n<AI_NAV_CELLS;n++)aiNavParent[n]=-1;
    int center=AI_NAV_CELLS/2,head=0,tail=1,best=center;
    aiNavQueue[0]=center;aiNavParent[center]=center;
    float bestDistance=v_length(world_delta(goal,origin));
    const int dx[]={1,-1,0,0},dz[]={0,0,1,-1};
    while(head<tail && head<AI_NAV_CELLS) {
        int cell=aiNavQueue[head++],x=cell%AI_NAV_SIDE,z=cell/AI_NAV_SIDE;
        Vector3 pos=ai_nav_point(origin,cell);float d=v_length(world_delta(goal,pos));
        if(d<bestDistance){best=cell;bestDistance=d;}
        if(d<.8f)break;
        for(int k=0;k<4;k++) {
            int nx=x+dx[k],nz=z+dz[k];if(nx<0 || nx>=AI_NAV_SIDE || nz<0 || nz>=AI_NAV_SIDE)continue;
            int next=nz*AI_NAV_SIDE+nx;if(aiNavParent[next]!=-1)continue;
            if(!ai_nav_clear(slot,pos,ai_nav_point(origin,next)))continue;
            aiNavParent[next]=cell;aiNavQueue[tail++]=next;
        }
    }
    int reverse[AI_NAV_PATH],count=0;
    for(int c=best;c!=center && count<AI_NAV_PATH;c=aiNavParent[c])reverse[count++]=c;
    // Reject truncated routes rather than jumping to an unreachable waypoint.
    if(count==AI_NAV_PATH)return;
    while(count)route->points[route->count++]=ai_nav_point(origin,reverse[--count]);
}
static Vector3 ai_navigation_direction(Player *bot,Vector3 desired,float dt) {
    int slot=(int)(bot-players);if(slot<0 || slot>=MAX_COMBATANTS)return desired;
    AiNavigation *route=&aiNavigation[slot];route->retry=fmaxf(0,route->retry-dt);
    desired.y=0;float distance=v_length(desired);
    if(distance<.01f){route->count=route->next=0;return desired;}
    Vector3 goal=world_wrap(v_add(bot->pos,desired));
    // Some callers pass a unit steering vector. Extend it only for routing.
    if(distance<1.5f)goal=world_wrap(v_add(bot->pos,v_mul(v_norm(desired),4)));
    if(route->count && (fabsf(route->points[route->next].y-bot->pos.y)>1 ||
        (route->retry<=0 && v_length(world_delta(goal,route->goal))>3)))route->count=route->next=0;
    while(route->next<route->count && v_length(world_delta(route->points[route->next],bot->pos))<.25f)route->next++;
    if(route->next>=route->count)route->count=route->next=0;
    if(route->count && !ai_nav_clear(slot,bot->pos,route->points[route->next]))route->count=route->next=0;
    if(!route->count && route->retry<=0) {
        Vector3 probe=world_wrap(v_add(bot->pos,v_mul(v_norm(desired),fminf(distance,4))));
        if(!ai_nav_clear(slot,bot->pos,probe))ai_nav_search(slot,goal);
    }
    if(route->count)return world_delta(route->points[route->next],bot->pos);
    return desired;
}
