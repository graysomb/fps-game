static void ex_gravity_activate(GravityField *f,int *budget){
    if(*budget<=0||f->activated>=256||f->age<.15f)return;
    static UnitVoxelBuffer units;units.count=0;
    int cx=floorf(f->pos.x/VOXEL_SIZE),cy=floorf(f->pos.y/VOXEL_SIZE),cz=floorf(f->pos.z/VOXEL_SIZE);
    // Batch removal prevents each cell from rebuilding support bookkeeping.
    // Coarse physical voxels are admitted whole, and counted by their unit volume.
    for(int scans=0;scans<512 && *budget>0 && f->activated<256 && f->cursor<25*25*25;scans++){
        int n=f->cursor++,x=world_grid(cx+n%25-12),y=cy+(n/25)%25-12,z=world_grid(cz+n/625-12);
        if(y<=0)continue;Vector3 at={(x+.5f)*VOXEL_SIZE,(y+.5f)*VOXEL_SIZE,(z+.5f)*VOXEL_SIZE};if(v_length(world_delta(at,f->pos))>6)continue;
        int v=table_get_static_only(x,y,z);if(v<0||voxel_is_fluid(&voxels[v])||voxels[v].pendingActivation||voxels[v].activationCooldownFrames>0)continue;
        Voxel *cell=&voxels[v];int x0,x1,y0,y1,z0,z1;voxel_grid_bounds(cell,&x0,&x1,&y0,&y1,&z0,&z1);
        int count=(x1-x0+1)*(y1-y0+1)*(z1-z0+1);if(y0<=0||count>*budget||count>256-f->activated)continue;
        if(voxel_count+units.count+count>=MAX_VOXELS || particle_pool_count-free_particle_count+8*(units.count+count)>=MAX_PARTICLES || recycleQueueCount+32>=MAX_VOXELS)break;
        bool inside=true;for(int z=z0;z<=z1&&inside;z++)for(int yy=y0;yy<=y1&&inside;yy++)for(int xx=x0;xx<=x1;xx++)if(v_length(world_delta((Vector3){(xx+.5f)*VOXEL_SIZE,(yy+.5f)*VOXEL_SIZE,(z+.5f)*VOXEL_SIZE},f->pos))>6)inside=false;
        if(!inside)continue;
        for(int z=z0;z<=z1;z++)for(int yy=y0;yy<=y1;yy++)for(int xx=x0;xx<=x1;xx++)unit_voxel_buffer_push(&units,xx,yy,z,cell->color,cell->type,false,v,cell->debugClusterTag,f->owner);
        cell->pendingActivation=true;f->activated+=count;*budget-=count;
    }
    if(units.count){remove_buffered_static_voxels(&units);int start=voxel_count;emit_unit_voxels_from_units(&units,false,true,-1);for(int v=start;v<voxel_count;v++){voxels[v].owner=f->owner;voxels[v].explosionUntil=clarityClock+5;voxels[v].activationCooldownFrames=120;}rebuild_voxel_hash();refresh_activation_surfaces(&units);glue_dynamic_voxel_to_static_neighbors();}
    units.count=0;
}
static Vector3 ex_pull_delta(Vector3 velocity,Vector3 delta,float acceleration,float dt,float ramp){
    float distance=v_length(delta);if(distance<.05f||distance>=6)return (Vector3){0};Vector3 direction=v_mul(delta,1/distance);
    float speed=fmaxf(0,12-v_dot(velocity,direction));return v_mul(direction,fminf(speed,acceleration*(1-distance/6)*dt*ramp));
}
static void ex_gravity_pull(GravityField *f,float dt){
    float ramp=clampf(f->age/.5f,0,1)*clampf((5-f->age)/.5f,0,1);
    // Wake first: waking assigns particle stamps and restores inverse masses.
    // Keep that separate from force deduplication across shared corners.
    for(int v=0;v<voxel_count;v++)if(voxels[v].simulate&&voxels[v].sleeping&&!voxels[v].fixed&&v_length(world_delta(f->pos,voxels[v].pos))<6&&enemy_equipment_owner(v)<0&&get_goliath_owner_of_voxel(v)<0)wake_sleeping_cluster_in_place(v);
    int stamp=++particle_sync_stamp;
    for(int v=0;v<voxel_count;v++){Voxel *x=&voxels[v];if(!x->simulate||x->fixed||x->isBullet||enemy_equipment_owner(v)>=0||get_goliath_owner_of_voxel(v)>=0)continue;
        Vector3 d=world_delta(f->pos,x->pos),dv=ex_pull_delta(x->vel,d,30,dt,ramp);if(v_length(dv)<1e-6f)continue;
        x->vel=v_add(x->vel,dv);x->sleeping=false;x->sleepFrames=0;x->wake_source=true;x->wake_timer=COARSENING_WAKE_FRAMES;x->owner=f->owner;x->activator=f->owner;x->explosionUntil=clarityClock+.5f;
        for(int c=0;c<8;c++){Particle *p=x->particles[c];if(p&&p->inv_mass>0&&p->sync_stamp!=stamp){p->sync_stamp=stamp;p->vel=v_add(p->vel,ex_pull_delta(p->vel,world_delta(f->pos,p->pos),30,dt,ramp));}}
    }
}

static bool ex_in_field(int i){
 if(i<0||i>=activePlayers||players[i].enemyType==ENEMY_TYPE_FOUNDRY)return false;
 for(int n=0;n<EX_FIELDS;n++)if(exFields[n].id){Vector3 d=world_delta(exFields[n].pos,players[i].pos);float length=v_length(d);VoxelHit hit;
  if(length<6 && (length<.3f||!first_voxel_hit_detailed((Ray){players[i].pos,v_norm(d)},length-.25f,-1,false,&hit)))return true;
 }return false;
}
static void ex_actor_pull(int i,float dt){
 if(i<0||i>=activePlayers||players[i].respawn_timer>0||players[i].enemyType==ENEMY_TYPE_FOUNDRY)return;
 bool bot=is_player_bot(i);if(bot)exBotPull[i]=v_mul(exBotPull[i],expf(-dt*1.5f));
 Vector3 velocity=v_add(players[i].vel,bot?exBotPull[i]:(Vector3){0}),delta={0};
 for(int n=0;n<EX_FIELDS;n++){GravityField *f=&exFields[n];if(!f->id)continue;
  float ramp=clampf(f->age/.5f,0,1)*clampf((5-f->age)/.5f,0,1);
  Vector3 d=world_delta(f->pos,players[i].pos);float distance=v_length(d);if(distance>=6)continue;VoxelHit block;
  if(distance>.3f&&first_voxel_hit_detailed((Ray){players[i].pos,v_norm(d)},distance-.25f,-1,false,&block))continue;
  Vector3 dv=ex_pull_delta(v_add(velocity,delta),d,players[i].enemyType==ENEMY_TYPE_GOLIATH?6:18,dt,ramp);delta=v_add(delta,dv);
 }
 if(bot){
  // AI chooses a desired walking velocity each frame. Keep pull momentum
  // separately so that decision cannot erase the field's acceleration.
  exBotPull[i].x+=delta.x;exBotPull[i].z+=delta.z;float speed=v_length(exBotPull[i]);if(speed>12)exBotPull[i]=v_mul(exBotPull[i],12/speed);
  players[i].vel=v_add(players[i].vel,v_add(exBotPull[i],(Vector3){0,delta.y,0}));exBotApplied[i]=exBotPull[i];
 }else players[i].vel=v_add(players[i].vel,delta);
 if(delta.y>.01f)players[i].onGround=false;
}
