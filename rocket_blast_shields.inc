// The limb shield is the box drawn around all of its attached cells, including
// the gaps between them. Capture it before the explosion can detach any cells.
typedef struct { BoundingBox box[GM_COUNT]; bool present[GM_COUNT]; } RocketLimbBounds;
static void rocket_limb_bounds(RocketLimbBounds *bounds){
    if(netTransport.role==NET_ROLE_CLIENT){
        if(!goRemoteActive || goRemoteDead)return;
        for(int p=1;p<GM_COUNT;p++)if(p!=GM_CORE && goRemoteParts[p].alive && goRemoteParts[p].shield>0){
            bounds->present[p]=true;bounds->box[p]=goRemoteParts[p].box;
        }
    }else if(!goCombatEnabled || goCombatDead)return;
    else for(int c=0;c<goCount;c++){
        int p=goCells[c].part;
        if(p<=0 || p==GM_CORE || goShields[p].value<=0 || !go_part_alive(p))continue;
        int v=go_find(goCells[c].id,goCells[c].index);
        if(v<0)continue;
        BoundingBox *b=&bounds->box[p];
        if(!bounds->present[p]){
            bounds->present[p]=true;b->min=b->max=voxels[v].particles[0]->pos;
        }
        for(int k=0;k<8;k++){
            Vector3 corner=voxels[v].particles[k]->pos;
            b->min=Vector3Min(b->min,corner);b->max=Vector3Max(b->max,corner);
        }
    }
    for(int p=1;p<GM_COUNT;p++)if(bounds->present[p]){
        Vector3 pad={.025f,.025f,.025f};
        bounds->box[p].min=v_sub(bounds->box[p].min,pad);
        bounds->box[p].max=v_add(bounds->box[p].max,pad);
    }
}

// Return the entry fraction of a segment. Touching a face while travelling
// outwards leaves the front side exposed; travelling inwards is blocked.
static bool rocket_crosses_box(Vector3 from,Vector3 segment,BoundingBox box,float *entry){
    Vector3 center=v_mul(v_add(box.min,box.max),.5f),o=world_delta(from,center);
    Vector3 h=v_mul(v_sub(box.max,box.min),.5f);float near=0,far=1;
    for(int k=0;k<3;k++){
        float d=(&segment.x)[k],pos=(&o.x)[k],half=(&h.x)[k];
        if(fabsf(d)<1e-7f){if(fabsf(pos)>=half)return false;continue;}
        float a=(-half-pos)/d,b=(half-pos)/d;
        near=fmaxf(near,fminf(a,b));far=fminf(far,fmaxf(a,b));
        if(near>=far)return false;
    }
    if(entry)*entry=near;
    return true;
}
typedef struct {BoundingBox box;uint64_t id;} RocketPanel;
typedef struct RocketBlastCover {
    Vector3 center;
    RocketLimbBounds limbs;
    RocketPanel *panels;
    int count,capacity;
    bool failed;
} RocketBlastCover;

static int rocket_first_limb(const RocketLimbBounds *bounds,Vector3 from,Vector3 segment,float *fraction){
    int closest=-1;
    for(int p=1;p<GM_COUNT;p++)if(bounds->present[p]){
        float t;
        if(rocket_crosses_box(from,segment,bounds->box[p],&t) && t<*fraction){*fraction=t;closest=p;}
    }
    return closest;
}
static void rocket_snapshot_cover(RocketBlastCover *cover,float reach){
    rocket_limb_bounds(&cover->limbs);
    for(int i=0;i<voxel_count;i++){
        Voxel *v=&voxels[i];
        if(!spawn_glass(i) && !clarity_shield_voxel(i))continue;
        Vector3 d=world_delta(v->pos,cover->center);
        if(v_dot(d,d)>(reach+VOXEL_SIZE)*(reach+VOXEL_SIZE))continue;
        if(cover->count==cover->capacity){
            int capacity=cover->capacity?cover->capacity*2:16;
            RocketPanel *more=realloc(cover->panels,capacity*sizeof(*more));
            if(!more){cover->failed=true;return;}
            cover->panels=more;cover->capacity=capacity;
        }
        Vector3 h={VOXEL_SIZE*.5f,VOXEL_SIZE*.5f,VOXEL_SIZE*.5f};
        cover->panels[cover->count++]=(RocketPanel){{v_sub(v->pos,h),v_add(v->pos,h)},v->identity};
    }
}
static bool rocket_cover_blocks(const RocketBlastCover *cover,Vector3 pos,uint64_t target,int targetLimb){
    if(!cover)return false;
    if(cover->failed)return true;
    Vector3 delta=world_delta(pos,cover->center);
    for(int p=1;p<GM_COUNT;p++)if(cover->limbs.present[p] && p!=targetLimb &&
        rocket_crosses_box(cover->center,delta,cover->limbs.box[p],NULL))return true;
    for(int n=0;n<cover->count;n++)if(cover->panels[n].id!=target &&
        rocket_crosses_box(cover->center,delta,cover->panels[n].box,NULL))return true;
    return false;
}
