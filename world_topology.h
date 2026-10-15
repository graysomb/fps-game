#ifndef WORLD_TOPOLOGY_H
#define WORLD_TOPOLOGY_H
/* Horizontal torus. Geometry/constraints use ordinary subtraction; gameplay
 * distances and broadphase coordinates explicitly use these helpers. */
static bool worldPeriodic=false;
static bool worldPeriodicDisabled=false;
static float world_wrap_axis(float x) {
    return worldPeriodic?x-floorf((x+FLOOR_SIZE)/(2*FLOOR_SIZE))*(2*FLOOR_SIZE):x;
}
static Vector3 world_wrap(Vector3 p) {
    p.x=world_wrap_axis(p.x);p.z=world_wrap_axis(p.z);return p;
}
static Vector3 world_delta(Vector3 a,Vector3 b) {
    return world_wrap((Vector3){a.x-b.x,a.y-b.y,a.z-b.z});
}
static Vector3 world_image(Vector3 p,Vector3 reference) {
    Vector3 d=world_delta(p,reference);return (Vector3){reference.x+d.x,reference.y+d.y,reference.z+d.z};
}
static int world_grid(int x) {
    if(!worldPeriodic)return x;
    int half=(int)lroundf(FLOOR_SIZE/VOXEL_SIZE),period=2*half;
    int q=(x+half)%period;return (q<0?q+period:q)-half;
}
#endif
