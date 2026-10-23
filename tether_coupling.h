/* Weak tethers exchange momentum; a standard voxel has eight mass-one corners. */
#define TETHER_GRAPPLE_PLAYER_MASS (16.f*TETHER_REFERENCE_MASS)
#define TETHER_GRAPPLE_ACCELERATION 30.f
#define TETHER_GRAPPLE_VOXEL_ACCELERATION 60.f
#define TETHER_GRAPPLE_MAX_APPROACH 12.f
#define TETHER_GRAPPLE_STOP_DISTANCE 1.5f
typedef struct {
    bool active,anchored;
    uint64_t identity;
    int generation,voxelCount,face,count;
    float mass;
    Vector3 target,velocity,force,impulse;
    Particle *corner[4];
} WeakTetherCoupling;
static WeakTetherCoupling weakTethers[MAX_PLAYERS];
static bool resolve_tether_voxel(Player *player);
static Vector3 player_tether_hold_position(const Player *player);
static bool player_tether_visual_target(int player_index,Vector3 *target);
static bool weak_tether_enabled(int i);
static bool weak_tether_any(void);
static bool weak_tether_sample(int i,WeakTetherCoupling *coupling);
static Vector3 weak_tether_impulse(Vector3 error,Vector3 relative,float inverseMass,float dt);
static void weak_tether_prepare(float dt);
static void weak_tether_integrate(float dt);
