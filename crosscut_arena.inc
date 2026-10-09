/* Crosscut uses signed world coordinates; box bounds are exclusive at the top. */
static Vector3 crosscut_rotate(Vector3 p,int quarter) {
    for(int i=0;i<quarter;i++) {float x=p.x;p.x=-p.z;p.z=x;}
    return p;
}
static void crosscut_box(int q,float x0,float x1,float y0,float y1,float z0,float z1,Color c) {
    Vector3 a=crosscut_rotate((Vector3){x0,y0,z0},q),b=crosscut_rotate((Vector3){x1,y1,z1},q);
    add_static_box_at_grid((int)roundf(fminf(a.x,b.x)/VOXEL_SIZE),(int)roundf(fmaxf(a.x,b.x)/VOXEL_SIZE)-1,
        (int)roundf(y0/VOXEL_SIZE),(int)roundf(y1/VOXEL_SIZE)-1,
        (int)roundf(fminf(a.z,b.z)/VOXEL_SIZE),(int)roundf(fmaxf(a.z,b.z)/VOXEL_SIZE)-1,c);
}
static Vector3 crosscut_spawn_point(int index) {
    return crosscut_rotate((Vector3){index%2?-4.5f:4.5f,BASE_EYE_HEIGHT,15.5f},index/2);
}
static Vector3 crosscut_supply_point(int index) {
    return crosscut_rotate((Vector3){9,1,11},index);
}
static void buildCrosscut(void) {
    Color stone={190,184,159,255},support={236,132,43,255};
    const Color routes[4]={{61,156,182,255},{187,88,107,255},{137,158,73,255},{146,113,185,255}};
    for(int q=0;q<4;q++) {
        /* Terraces and two three-metre-wide stepped approaches. */
        crosscut_box(q,-4,4,0,2.5f,9,12,stone);
        for(int s=0;s<5;s++) {
            crosscut_box(q,4+s,5+s,0,2.5f-.5f*s,9,12,routes[q]);
            crosscut_box(q,-5-s,-4-s,0,2.5f-.5f*s,9,12,routes[q]);
        }
        crosscut_box(q,-1,1,2.5f,3.5f,11,12,routes[q]);
        /* Short independent slabs, separated by half-metre seams. */
        for(int s=0;s<3;s++) {
            float z=2+s*2.5f;
            crosscut_box(q,-1.5f,1.5f,2,2.5f,z,z+2,stone);
            crosscut_box(q,-1.5f,-1,0,2,z+.5f,z+1,support);
            crosscut_box(q,1,1.5f,0,2,z+.5f,z+1,support);
        }
        /* Spawn screens leave both ends open to the perimeter circuit. */
        crosscut_box(q,-6.5f,6.5f,0,3,13,13.5f,stone);
        crosscut_box(q,-.5f,.5f,0,3,13.5f,17,stone);
        crosscut_box(q,10,11,0,1,6,7,routes[q]);
        crosscut_box(q,11.5f,12.5f,0,1,6,7,stone);
    }
    crosscut_box(0,-1.5f,1.5f,2,2.5f,-1.5f,1.5f,stone);
    for(int q=0;q<4;q++) crosscut_box(q,1,1.5f,0,2,1,1.5f,support);
}
static int crosscutSpawnCursor=0;
static bool crosscutInitialSpawns=false;
