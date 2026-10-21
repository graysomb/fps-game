typedef struct {uint32_t id;Vector3 origin;float age;int rotation;} IglooReveal;
static IglooReveal exReveals[8];
static void ex_construction_reset(void){memset(exReveals,0,sizeof(exReveals));}
static Vector3 ex_igloo_cell(Vector3 origin,int x,int y,int z,int rotation){Vector3 p={(x-5)*VOXEL_SIZE,y*VOXEL_SIZE,(z-5)*VOXEL_SIZE};return world_wrap(v_add(origin,crosscut_rotate(p,rotation)));}
static bool ex_cell_free(Vector3 at,bool actors){
    int x=(int)floorf(at.x/VOXEL_SIZE),y=(int)floorf(at.y/VOXEL_SIZE),z=(int)floorf(at.z/VOXEL_SIZE);if(occupied(x,y,z))return false;
    Vector3 half={.249f,.249f,.249f};BoundingBox box={v_sub(at,half),v_add(at,half)};
    for(int v=0;v<voxel_count;v++)if(voxels[v].simulate){VoxelWorldBounds b;voxel_particle_world_bounds(&voxels[v],&b);Vector3 center={(b.minx+b.maxx)*.5f,(b.miny+b.maxy)*.5f,(b.minz+b.maxz)*.5f},shift=v_sub(world_image(center,at),center);if(CheckCollisionBoxes(box,(BoundingBox){v_add((Vector3){b.minx,b.miny,b.minz},shift),v_add((Vector3){b.maxx,b.maxy,b.maxz},shift)}))return false;}
    if(actors)for(int i=0;i<activePlayers;i++)if((players[i].respawn_timer<=0||combatantLife[i]==LIFE_DOWNED)&&CheckCollisionBoxes(box,combatant_bounds_at(i,world_image(players[i].pos,at),PLAYER_SIZE*.5f)))return false;
    return true;
}
static bool ex_igloo(GrenadeState g){
    if(g.owner<0||g.owner>=MAX_PLAYERS)return false;
    Vector3 ground=g.pos;VoxelHit hit;Ray down={v_add(g.pos,(Vector3){0,.2f,0}),{0,-1,0}};
    if(first_voxel_hit_detailed(down,3.2f,-1,false,&hit)){if(voxels[hit.id].simulate)return false;ground.y=down.position.y-hit.t;}
    else {if(g.pos.y>3)return false;ground.y=0;}
    if(fabsf(ground.y/VOXEL_SIZE-roundf(ground.y/VOXEL_SIZE))>.02f)return false;
    Vector3 origin={(floorf(ground.x/VOXEL_SIZE)+.5f)*VOXEL_SIZE,ground.y-.25f,(floorf(ground.z/VOXEL_SIZE)+.5f)*VOXEL_SIZE};
    Vector3 toward=world_delta(g.launchOrigin,origin);int rotation=0;float best=-1e9f;
    for(int r=0;r<4;r++){float score=v_dot(toward,crosscut_rotate((Vector3){0,0,-1},r));if(score>best){best=score;rotation=r;}}
    Vector3 cells[IGLOO_MAX_CELLS];int count=0;
    // Validate the entire foundation, entrance and shell before allocating anything.
    for(int z=0;z<IGLOO_SIDE;z++)for(int x=0;x<IGLOO_SIDE;x++)if(igloo_cell(x,0,z,false)){
        Vector3 at=ex_igloo_cell(origin,x,0,z,rotation);int gx=floorf(at.x/VOXEL_SIZE),gy=floorf(at.y/VOXEL_SIZE),gz=floorf(at.z/VOXEL_SIZE);
        if(!worldPeriodic && (fabsf(at.x)>worldFloorHalfSize-.25f||fabsf(at.z)>worldFloorHalfSize-.25f))return false;
        int support=table_get_static_only(gx,gy,gz);
        if(ground.y>0 && support<0)return false;
        if(support<0){if(!ex_cell_free(at,true))return false;cells[count++]=at;}
    }
    for(int y=1;y<IGLOO_HEIGHT;y++)for(int z=0;z<IGLOO_SIDE;z++)for(int x=0;x<IGLOO_SIDE;x++)if(igloo_cell(x,y,z,false)){
        Vector3 at=ex_igloo_cell(origin,x,y,z,rotation);if(!ex_cell_free(at,true)||count>=IGLOO_MAX_CELLS)return false;cells[count++]=at;
    }
    for(int y=1;y<=3;y++)for(int z=0;z<=3;z++)for(int x=4;x<=6;x++)if(!ex_cell_free(ex_igloo_cell(origin,x,y,z,rotation),false))return false;
    if(voxel_count+count>=MAX_VOXELS || MAX_PARTICLES-particle_pool_count+free_particle_count<count*8)return false;
    uint64_t made[IGLOO_MAX_CELLS];int used=0;
    for(int n=0;n<count;n++){Vector3 at=cells[n];int v=add_static_voxel_at_grid(floorf(at.x/VOXEL_SIZE),floorf(at.y/VOXEL_SIZE),floorf(at.z/VOXEL_SIZE),(Color){190,205,215,255},0);if(v<0)break;made[used++]=voxels[v].identity;voxels[v].reinforced=true;voxels[v].persistentConstruction=true;tag_owned_static_voxel(v,g.owner);mark_surface(v);mark_surface_neighbors(at);}
    if(used!=count){for(int n=0;n<used;n++){int v=enemy_voxel(made[n]);if(v>=0)remove_voxel_index(v);}return false;}
    int slot=0;for(int n=0;n<8;n++)if(!exReveals[n].id||exReveals[n].age>exReveals[slot].age)slot=n;
    exReveals[slot]=(IglooReveal){.id=exNextId++,.origin=origin,.rotation=rotation};play_sfx_at(SFX_SHIELD,origin,g.owner,2,.8f);return true;
}
static void ex_construction_tick(float dt){for(int n=0;n<8;n++)if(exReveals[n].id){exReveals[n].age+=dt;if(exReveals[n].age>=1.2f)exReveals[n].id=0;}}
static void ex_construction_draw(WeaponRenderPass pass){
    // Highlight real, surviving cells. The collision and visible shelter appear together.
    for(int n=0;n<8;n++){IglooReveal *r=&exReveals[n];if(!r->id||!periodic_visible(r->origin,4))continue;int course=(int)(r->age/1.2f*IGLOO_HEIGHT);
        Vector3 eye=drawViewPlayerIndex>=0?players[drawViewPlayerIndex].pos:(Vector3){0};WeaponVisual visual={.shot_time=-1000,.deployment=clampf(r->age/.3f,0,1),.build_pulse=1-r->age/1.2f};
        weapon_draw((WeaponPose){.kind=WEAPON_IGLOO_GRENADE,.transform=MatrixMultiply(MatrixScale(.45f,.45f,.45f),MatrixTranslate(r->origin.x,r->origin.y+.4f,r->origin.z))},&visual,eye,exClock,false,pass);
        Color c=Fade(SKYBLUE,1-r->age/1.2f);for(int z=0;z<IGLOO_SIDE;z++)for(int x=0;x<IGLOO_SIDE;x++)if(igloo_cell(x,course,z,false)){Vector3 at=ex_igloo_cell(r->origin,x,course,z,r->rotation);if(occupied(floorf(at.x/VOXEL_SIZE),floorf(at.y/VOXEL_SIZE),floorf(at.z/VOXEL_SIZE)))DrawCubeWires(at,.505f,.505f,.505f,c);}
        if(pass==WEAPON_COLOR)DrawCircle3D(v_add(r->origin,(Vector3){0,.26f,0}),3*r->age/1.2f,(Vector3){1,0,0},90,c);
    }
}
