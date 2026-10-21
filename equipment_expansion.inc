/* Gameplay for the visual expansion. No simulation is performed by drawing. */
static uint32_t exPartTick[16];
static bool ex_suppressed(int i){return i>=0 && i<MAX_COMBATANTS && exSuppression[i]>0;}
static void ex_clear_status(int i){if(i<0||i>=MAX_COMBATANTS)return;exSuppression[i]=0;for(int n=0;n<EX_SPLATS;n++)if(exSplats[n].actor==i)exSplats[n].id=0;}
static void ex_reset_player(int i){if(i<0||i>=MAX_COMBATANTS)return;ex_clear_status(i);exBotPull[i]=exBotApplied[i]=(Vector3){0};++exIncarnation[i];if(i<MAX_COMBATANTS){memset(exGrenades[i],0,sizeof(exGrenades[i]));exSelectedGrenade[i]=0;exFailedBuild[i]=0;}}
static void ex_reset(void){memset(exGlobs,0,sizeof(exGlobs));memset(exFields,0,sizeof(exFields));memset(exSplats,0,sizeof(exSplats));memset(exPartTick,0,sizeof(exPartTick));memset(exLastCycleInput,0,sizeof(exLastCycleInput));for(int i=0;i<MAX_COMBATANTS;i++)ex_reset_player(i);exClock=0;ex_construction_reset();}
static int ex_grenade_count(int i,int kind){if(i<0||i>=MAX_COMBATANTS||kind<0||kind>=GRENADE_KIND_COUNT)return 0;return kind?exGrenades[i][kind-1]:grenadeInventory[i];}
static void ex_cycle(int i){if(i<0||i>=MAX_COMBATANTS)return;for(int n=1;n<=GRENADE_KIND_COUNT;n++){int k=(exSelectedGrenade[i]+n)%GRENADE_KIND_COUNT;if(ex_grenade_count(i,k)>0){exSelectedGrenade[i]=k;return;}}}
static bool ex_grenade_pickup(int i,int kind){if(i<0||i>=MAX_COMBATANTS||kind<1||kind>2||exGrenades[i][kind-1]>=2)return false;bool empty=ex_grenade_count(i,exSelectedGrenade[i])==0;exGrenades[i][kind-1]++;if(empty)exSelectedGrenade[i]=kind;return true;}
static bool ex_throw_available(int i,int kind){
    if(!ex_grenade_count(i,kind))return false;if(kind!=GRENADE_GRAVITY)return true;
    int count=0;for(int n=0;n<EX_FIELDS;n++)if(exFields[n].id){if(exFields[n].owner==i)return false;count++;}
    for(int n=0;n<CLARITY_GRENADES;n++)if(grenades[n].id && grenades[n].kind==GRENADE_GRAVITY){if(grenades[n].owner==i)return false;count++;}
    return count<EX_FIELDS;
}
static void ex_input(int seat,NetInputCommand *c){
    if(seat<0||seat>=MAX_PLAYERS||creativeModeActive||(enemyLab&&!equipmentLab))return;
    bool pressed=playerInput[seat]==INPUT_TYPE_GAMEPAD?IsGamepadButtonPressed(seat_device(seat),GAMEPAD_BUTTON_LEFT_FACE_DOWN):IsKeyPressed(seat==0?KEY_TAB:KEY_SEMICOLON);
    if(pressed)c->pressed|=NET_INPUT_GRENADE_CYCLE;
    if(fp_active() && firefightWaveState==FIREFIGHT_STATE_INTERMISSION && netTransport.role==NET_ROLE_CLIENT){
        bool pad=playerInput[seat]==INPUT_TYPE_GAMEPAD;int d=seat_device(seat);
        if(pad?IsGamepadButtonPressed(d,GAMEPAD_BUTTON_LEFT_FACE_LEFT):IsKeyPressed(seat==0?KEY_F:KEY_LEFT))c->pressed|=NET_INPUT_PREP_LEFT;
        if(pad?IsGamepadButtonPressed(d,GAMEPAD_BUTTON_LEFT_FACE_RIGHT):IsKeyPressed(seat==0?KEY_H:KEY_RIGHT))c->pressed|=NET_INPUT_PREP_RIGHT;
        if(pad?IsGamepadButtonPressed(d,GAMEPAD_BUTTON_RIGHT_FACE_DOWN):IsKeyPressed(seat==0?KEY_SPACE:KEY_RIGHT_SHIFT))c->pressed|=NET_INPUT_PREP_READY;
        if(pad?IsGamepadButtonPressed(d,GAMEPAD_BUTTON_RIGHT_FACE_UP):IsKeyPressed(seat==0?KEY_V:KEY_N))c->pressed|=NET_INPUT_PREP_MENU;
    }
}
static void ex_splat(Vector3 at,int actor,int voxel,bool impact){
    int free=-1,oldest=0,count=0,actorOld=-1;
    for(int n=0;n<EX_SPLATS;n++){ExpansionSplat *s=&exSplats[n];if(!s->id)free=n;if(s->age>exSplats[oldest].age)oldest=n;if(s->id && actor>=0 && s->actor==actor){count++;if(actorOld<0||s->age>exSplats[actorOld].age)actorOld=n;}}
    if(count>=3)free=actorOld;if(free<0)free=oldest;
    ExpansionSplat *s=&exSplats[free];*s=(ExpansionSplat){.id=exNextId++,.actor=actor,.pos=at,.impact=impact};
    if(actor>=0)s->local=Vector3RotateByAxisAngle(world_delta(at,players[actor].pos),(Vector3){0,1,0},-players[actor].yaw*DEG2RAD);
    else if(voxel>=0 && voxel<voxel_count){s->support=voxels[voxel].identity;s->local=Vector3Transform(world_image(at,voxels[voxel].pos),MatrixInvert(grenade_support_basis(voxel)));}
}
static void ex_dmr_cut(int owner,Ray ray,VoxelHit hit){
    int v=hit.id;if(v<0||v>=voxel_count)return;Vector3 at=world_wrap(v_add(ray.position,v_mul(ray.direction,hit.t)));
    if(!voxels[v].simulate){
        remove_static_voxels_in_region_recycle(hit.x,hit.x,hit.y,hit.y,hit.z,hit.z);
        activate_static_voxels_near_region(hit.x,hit.x,hit.y,hit.y,hit.z,hit.z,owner);
    }else {enemy_detach(v);detach_goliath_armor_voxel(v);remove_voxel_index(v);}
    mark_surface_neighbors(at);ex_splat(at,-1,-1,true);
    // Resolve each loose island once; never transmit impulse through cover.
    static bool applied[MAX_VOXELS];memset(applied,0,(size_t)voxel_count);int targets[32],count=0;float weights[32],total=0;
    for(int j=0;j<voxel_count && count<32;j++)if(voxels[j].simulate&&!voxels[j].fixed&&!voxels[j].isBullet && enemy_equipment_owner(j)<0 && get_goliath_owner_of_voxel(j)<0){
        Vector3 d=world_delta(voxels[j].pos,at);float dist=v_length(d);if(dist>1.5f)continue;
        if(applied[j])continue;
        VoxelHit block;if(dist>.3f&&first_voxel_hit_detailed((Ray){v_add(at,v_mul(ray.direction,.02f)),v_norm(d)},dist-.3f,-1,false,&block))continue;
        int nc=build_glue_cluster_indices(j,glueClusterIndices);for(int k=0;k<nc;k++)applied[glueClusterIndices[k]]=true;targets[count]=j;weights[count]=1-dist/1.5f;total+=weights[count++];
    }
    for(int n=0;n<count && total>0;n++)apply_tether_impact_to_dynamic_island(targets[n],0,v_mul(ray.direction,eqShotgunImpulse*.25f*weights[n]/total));
}
static bool ex_fire(int i,int gun,Ray ray){
    if(gun==EQ_GLOB){int free=-1;for(int n=0;n<EX_GLOBS;n++)if(!exGlobs[n].id){free=n;break;}if(free<0)return false;exGlobs[free]=(GlobProjectile){.id=exNextId++,.owner=i,.pos=ray.position,.velocity=v_mul(ray.direction,60)};return true;}
    float length=eqGuns[gun].range;VoxelHit hit;bool wall=first_voxel_hit_detailed(ray,length,-1,false,&hit);if(wall)length=hit.t;
    int target=eq_actor_ray(i,ray,&length);
    if(target>=0)eq_damage(i,target,eqGuns[gun].damage,gun==EQ_DMR,gun==EQ_DMR?DEATH_DMR:DEATH_GLOCK,ray.position);
    else if(wall && gun==EQ_DMR)ex_dmr_cut(i,ray,hit);
    eq_trace(player_weapon_pose(i,false).muzzle,v_add(ray.position,v_mul(ray.direction,length)),gun,i);return true;
}
#include "igloo_gameplay.inc"
#include "gravity_gameplay.inc"
static void ex_grenade_landed(GrenadeState *g){g->landed=true;if(g->kind==GRENADE_IGLOO)g->fuse=0;}
static void ex_grenade_deploy(GrenadeState g){
    if(g.kind==GRENADE_REGULAR){clarity_explode(g.pos,g.owner);return;}
    if(g.kind==GRENADE_IGLOO){
        if(!ex_igloo(g)){int i=g.owner;if(i>=0&&i<MAX_COMBATANTS&&exIncarnation[i]==g.incarnation&&players[i].respawn_timer<=0){exGrenades[i][1]=(int)fminf(2,exGrenades[i][1]+1);exSelectedGrenade[i]=GRENADE_IGLOO;exFailedBuild[i]=2;}}return;
    }
    for(int n=0;n<EX_FIELDS;n++)if(!exFields[n].id){exFields[n]=(GravityField){.id=exNextId++,.owner=g.owner,.pos=g.pos};play_sfx_at(SFX_SHIELD,g.pos,g.owner,2,.6f);return;}
}
static void ex_tick(float dt){
    exClock+=dt;bool client=netTransport.role==NET_ROLE_CLIENT;
    for(int i=0;i<MAX_COMBATANTS;i++)exSuppression[i]=fmaxf(0,exSuppression[i]-dt);
    for(int i=0;i<MAX_PLAYERS;i++)exFailedBuild[i]=fmaxf(0,exFailedBuild[i]-dt);
    for(int n=0;n<EX_SPLATS;n++){ExpansionSplat *s=&exSplats[n];if(!s->id)continue;s->age+=dt;if(s->age>=(s->impact?.65f:5)){s->id=0;continue;}
        if(s->actor>=0){if(players[s->actor].respawn_timer>0){s->id=0;continue;}s->pos=world_wrap(v_add(players[s->actor].pos,Vector3RotateByAxisAngle(s->local,(Vector3){0,1,0},players[s->actor].yaw*DEG2RAD)));}
        else if(s->support){int v=enemy_voxel(s->support);if(v>=0)s->pos=world_wrap(Vector3Transform(s->local,grenade_support_basis(v)));else if(!client)s->support=0;}
    }
    for(int n=0;n<EX_GLOBS;n++){GlobProjectile *g=&exGlobs[n];if(!g->id)continue;float step=fminf(dt,1-g->age);g->age+=dt;
        if(client){g->pos=world_wrap(v_add(g->pos,v_mul(g->velocity,step)));if(g->age>=1)g->id=0;continue;}
        float length=60*step;Ray ray={g->pos,v_norm(g->velocity)};VoxelHit hit;bool wall=first_voxel_hit_detailed(ray,length,-1,false,&hit);if(wall)length=hit.t;int actor=eq_actor_ray(g->owner,ray,&length);
        bool floor=ray.direction.y<0 && ray.position.y+ray.direction.y*length<0;if(floor){length=-ray.position.y/ray.direction.y;actor=-1;}
        g->pos=world_wrap(v_add(ray.position,v_mul(ray.direction,length)));
        if(actor>=0){if(eq_hostile(g->owner,actor)&&players[actor].invuln_timer<=0&&enemy_core_open(actor)){
                eq_damage(g->owner,actor,eqGuns[EQ_GLOB].damage,false,DEATH_GLOB,ray.position);
                if(players[actor].respawn_timer<=0){exSuppression[actor]=5;ex_splat(g->pos,actor,-1,false);}
            }g->id=0;
        }else if(wall||floor){ex_splat(g->pos,-1,wall&&!floor?hit.id:-1,false);g->id=0;}else if(g->age>=1)g->id=0;
    }
    int budget=32;for(int n=0;n<EX_FIELDS;n++)if(exFields[n].id){GravityField *f=&exFields[n];f->age+=dt;if(f->age>=5){f->id=0;continue;}if(!client){ex_gravity_activate(f,&budget);ex_gravity_pull(f,dt);}}
    ex_construction_tick(dt);
}
static void ex_draw(WeaponRenderPass pass){
    Vector3 eye=drawViewPlayerIndex>=0?players[drawViewPlayerIndex].pos:(Vector3){0};
    for(int n=0;n<EX_GLOBS;n++){GlobProjectile *g=&exGlobs[n];if(!g->id||!periodic_visible(g->pos,1))continue;WeaponEffect e={.transform=MatrixTranslate(g->pos.x,g->pos.y,g->pos.z),.age=g->age,.lifetime=1,.seed=g->id};weapon_draw_glob(&e,eye,exClock,pass);}
    for(int n=0;n<EX_SPLATS;n++){ExpansionSplat *s=&exSplats[n];if(!s->id||!periodic_visible(s->pos,1))continue;WeaponEffect e={.transform=MatrixTranslate(s->pos.x,s->pos.y,s->pos.z),.age=s->age,.lifetime=s->impact?.65f:5,.seed=s->id,.attached=true};if(s->impact)weapon_draw_dmr_impact(&e,eye,exClock,pass);else weapon_draw_glob(&e,eye,exClock,pass);}
    for(int n=0;n<EX_FIELDS;n++){GravityField *f=&exFields[n];if(!f->id||!periodic_visible(f->pos,6))continue;WeaponEffect e={.transform=MatrixMultiply(MatrixScale(2,2,2),MatrixTranslate(f->pos.x,f->pos.y,f->pos.z)),.age=f->age,.lifetime=5,.seed=f->id};weapon_draw_gravity(&e,eye,exClock,pass);WeaponVisual v={.shot_time=-1000,.deployment=clampf(f->age/.35f,0,1),.activity=1};weapon_draw((WeaponPose){.kind=WEAPON_BLACK_HOLE_BOMB,.transform=MatrixMultiply(MatrixScale(.6f,.6f,.6f),MatrixTranslate(f->pos.x,f->pos.y,f->pos.z))},&v,eye,exClock,false,pass);}
    float now=weaponFrameSampled?weaponFrameTime:(float)GetTime();
    for(int i=0;i<MAX_COMBATANTS;i++){int gun=eq_gun(i);if(gun<EQ_GLOB||players[i].respawn_timer>0)continue;float age=now-weaponVisuals[i].shot_time;if(age<0||age>.065f)continue;
        WeaponPose pose=player_weapon_pose(i,i==drawViewPlayerIndex);Vector3 end=v_add(pose.muzzle,v_mul(player_forward(&players[i]),.18f*(1-age/.065f)));
        Color tint=gun==EQ_GLOB?LIME:gun==EQ_DMR?SKYBLUE:(Color){255,220,150,255};DrawCylinderEx(pose.muzzle,end,.055f*(1-age/.065f),.001f,5,tint);
    }
    ex_construction_draw(pass);
}
static void ex_hud(int i,int w,int h){
    if(exSuppression[i]>0)DrawText(TextFormat("MATTER RESTORE BLOCKED %.1fs",exSuppression[i]),w/2-110,h/2+85,w<700?12:16,LIME);
    if(exFailedBuild[i]>0)DrawText("Shelter blocked - grenade returned",w/2-130,h/2+105,w<700?12:16,SKYBLUE);
    const char *name[]={"Frag","Gravity","Igloo"};int size=12;
    for(int k=0;k<3;k++){int x=w-252+k*84,y=h-60;DrawRectangle(x-3,y-2,81,16,Fade(BLACK,.7f));Color c=k==exSelectedGrenade[i]?GOLD:GRAY;DrawText(TextFormat("%s %d",name[k],ex_grenade_count(i,k)),x,y,size,c);}
}
static void ex_layout(void){
    const int types[]={PICKUP_GLOB,PICKUP_GLOCK,PICKUP_DMR,PICKUP_BLACK_HOLE,PICKUP_IGLOO};
    const Vector3 at[]={{26,0,20},{34,0,16},{9,0,5},{5,0,9},{32,0,29}};
    bool ff=gameMode==GAME_MODE_FIREFIGHT;
    for(int t=0;t<5;t++)for(int q=0;q<(ff?2:4);q++){Vector3 pos;if(!clarity_supply_position(crosscut_rotate(at[t],ff?q*2:q),0,&pos))continue;for(int n=0;n<MAX_PICKUPS;n++)if(!pickups[n].active&&pickups[n].respawnTimer<=0){pickups[n]=(Pickup){.pos=pos,.active=true,.type=types[t]};break;}}
}
