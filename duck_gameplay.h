#include "duck_renderer.h"
#define MAP_DUCKS 4
#define DUCK_RESPAWN_SECONDS 90.0f
#define DUCK_ROUTE_RADIUS 3.0f
#define DUCK_SWIM_SPEED .36f
typedef struct MapDuck {
    Vector3 home,pos,puff;
    float yaw,phase,respawn;
    uint32_t bursts;
    bool alive;
} MapDuck;
static MapDuck mapDucks[MAP_DUCKS];
static bool ducksInitialized;
static float duckClock;
static uint32_t duckReceiveTick;
static void ducks_reset(void);
static void ducks_tick(float dt);
static void ducks_draw(void);
static bool ducks_ray(Ray ray,float distance,float radius,int ignore,bool penetrate);
static void ducks_segment(Vector3 start,Vector3 end,float radius,int ignore);
static void ducks_blast(Vector3 pos,float radius);
static void ducks_snapshot(int peer);
static bool ducks_receive(NetMessageType type,NetReader *r,uint32_t tick);
