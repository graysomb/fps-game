/* Diagnostic only: keep real attached armor still and exercise the production
 * projectile paths, with rendered cube bounds as an independent coverage check. */
static bool debug_setup_goliath_armor_shots(DebugRunState *state) {
    debug_reset_world();ex_reset();eq_reset_all();
    gameMode=GAME_MODE_FIREFIGHT;gameState=GAME_STATE_PLAYING;activePlayers=2;
    enemyLab=false;creativeModeActive=false;worldPeriodic=false;
    debug_init_enemy_player(0,(Vector3){0,2,8},ENEMY_TYPE_STANDARD,INPUT_TYPE_KEYBOARD);
    debug_init_enemy_player(1,(Vector3){0,2,0},ENEMY_TYPE_GOLIATH,INPUT_TYPE_BOT_MEDIUM);
    combatantTeam[0]=0;combatantTeam[1]=1;
    spawn_goliath_armor(1);update_goliath_armor_positions(1);rebuild_voxel_hash();
    state->camera=(Camera3D){{5,5,8},{0,2,0},{0,1,0},50,CAMERA_PERSPECTIVE};
    char path[DEBUG_PATH_CAPACITY];debug_join_path(path,sizeof(path),state->output_dir,"armor-shots.csv");
    FILE *out=fopen(path,"w");if(!out)return false;
    fprintf(out,"gun,position,hz,angle,height,lateral,covered,armor_distance,body_distance,damage\n");
    int covered[2]={0},leaks[2]={0},openHits[2]={0};
    const int rates[]={30,60,144};
    for(int gun=0;gun<2;gun++)for(int position=0;position<2;position++)for(int rate=0;rate<3;rate++)for(int angle=0;angle<8;angle++)for(int y=-1;y<=1;y++)for(int x=-1;x<=1;x++){
        players[1].pos=position?(Vector3){.13f,2.17f,.23f}:(Vector3){0,2,0};
        update_goliath_armor_positions(1);rebuild_voxel_hash();
        float a=angle*PI/4;Vector3 outward={sinf(a),0,cosf(a)},side={cosf(a),0,-sinf(a)};
        Vector3 at=v_add(players[1].pos,v_add(v_mul(side,x*.25f),(Vector3){0,y*.35f,0}));
        Ray ray={v_add(at,v_mul(outward,8)),v_mul(outward,-1)};
        players[0].pos=ray.position;players[1].matter=100;players[1].isExposed=false;players[1].invuln_timer=0;
        players[1].respawn_timer=0;exSuppression[1]=0;
        RayCollision body=GetRayCollisionBox(ray,combatant_bounds_at(1,players[1].pos,PLAYER_SIZE*.5f));
        float armor=1000;
        for(int n=0;n<goliathStates[1].count;n++){
            int v=enemy_voxel(goliathStates[1].voxels[n].voxelIdentity);if(v<0)continue;
            float h=voxels[v].rest_edge*.5f;Vector3 half={h,h,h};
            RayCollision hit=GetRayCollisionBox(ray,(BoundingBox){v_sub(voxels[v].pos,half),v_add(voxels[v].pos,half)});
            if(hit.hit)armor=fminf(armor,hit.distance);
        }
        bool blocked=body.hit && armor<body.distance;
        if(gun){ex_fire(0,EQ_GLOB,ray);}else{
            int v=addVoxel(ray.position.x,ray.position.y,ray.position.z,false,true,WHITE,0);
            voxels[v].isBullet=true;voxels[v].owner=0;voxels[v].activator=0;voxels[v].vel=v_mul(ray.direction,60);voxels[v].shotOrigin=ray.position;
        }
        rebuild_voxel_hash();
        for(int f=0;f<rates[rate]/2;f++){if(gun)ex_tick(1.f/rates[rate]);else update_projectiles(1.f/rates[rate]);}
        float damage=100-players[1].matter;
        if(blocked){covered[gun]++;if(damage>0)leaks[gun]++;}else if(damage>0)openHits[gun]++;
        fprintf(out,"%s,%d,%d,%d,%.2f,%.2f,%d,%.4f,%.4f,%.1f\n",gun?"glob":"default",position,rates[rate],angle*45,y*.35f,x*.25f,blocked,armor,body.hit?body.distance:-1,damage);
        for(int v=voxel_count-1;v>=0;v--)if(voxels[v].isBullet)remove_voxel_index(v);
        ex_reset();rebuild_voxel_hash();
    }
    fclose(out);
    for(int gun=0;gun<2;gun++){
        fprintf(stderr,"ARMOR %s: covered=%d damaged-through-coverage=%d open-body-hits=%d\n",gun?"glob":"default",covered[gun],leaks[gun],openHits[gun]);
        if(!covered[gun])debug_add_failure(state,"No covered %s shots exercised",gun?"glob":"default");
        if(leaks[gun])debug_add_failure(state,"%s: %d/%d shots damaged through visible intact armor",gun?"glob":"default",leaks[gun],covered[gun]);
    }
    // Precise armor queries must preserve periodic wrapping and filter semantics.
    worldPeriodic=true;worldFloorHalfSize=40;players[1].pos=(Vector3){39.8f,2,0};
    update_goliath_armor_positions(1);rebuild_voxel_hash();
    Ray seam={{-37,2,0},{-1,0,0}};VoxelHit hit;
    if(!first_voxel_hit_detailed(seam,8,-1,false,&hit) || get_goliath_owner_of_voxel(hit.id)!=1)
        debug_add_failure(state,"Armor did not block across periodic boundary");
    if(first_voxel_hit_detailed_ex(seam,8,-1,false,1,&hit))debug_add_failure(state,"Ignored armor owner still blocked ray");
    if(first_voxel_hit_detailed(seam,8,-1,true,&hit))debug_add_failure(state,"Static-only ray hit attached armor");
    worldPeriodic=false;players[1].pos=(Vector3){0,2,0};update_goliath_armor_positions(1);rebuild_voxel_hash();
    // Positive controls: an entirely unarmored Goliath must take both shots.
    for(int n=goliathStates[1].count-1;n>=0;n--){int v=enemy_voxel(goliathStates[1].voxels[n].voxelIdentity);if(v>=0)remove_voxel_index(v);}
    goliathStates[1].count=0;rebuild_voxel_hash();
    for(int gun=0;gun<2;gun++){
        players[1].matter=100;players[1].isExposed=false;
        Ray ray={v_add(players[1].pos,(Vector3){0,0,8}),(Vector3){0,0,-1}};players[0].pos=ray.position;
        if(gun)ex_fire(0,EQ_GLOB,ray);else{int v=addVoxel(ray.position.x,ray.position.y,ray.position.z,false,true,WHITE,0);voxels[v].isBullet=true;voxels[v].owner=0;voxels[v].activator=0;voxels[v].vel=v_mul(ray.direction,60);}
        rebuild_voxel_hash();for(int f=0;f<30;f++){if(gun)ex_tick(1.f/60);else update_projectiles(1.f/60);}
        fprintf(stderr,"UNARMORED %s damage=%.1f\n",gun?"glob":"default",100-players[1].matter);
        if(players[1].matter>=100)debug_add_failure(state,"Unarmored positive control did not take damage");
    }
    spawn_goliath_armor(1);update_goliath_armor_positions(1);rebuild_voxel_hash();
    return true;
}
