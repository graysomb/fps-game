/* Authored arena data: deliberately independent of physics particles and live actors. */
#define CM_SPAWNS 64
#define CM_NAME 64
#define CM_PATH 512
#define CM_MIN_WIDTH 40
#define CM_MAX_WIDTH 256
enum {CM_PLAYER,CM_ENTRANCE,CM_CAMP,CM_BOSS};
enum {CM_BODY=1,CM_REMOVE=2,CM_TRIGGER=4};
typedef struct {int x,y,z;Color color;unsigned char material,reinforced,trap,flags;} CmCell;
typedef struct {Vector3 pos;float yaw;int kind,group;} CmSpawn;
typedef struct {char name[CM_NAME];float width;bool periodic;int count,capacity;CmCell *cells;
 int pickupCount;Pickup supplies[MAX_PICKUPS];int spawnCount;CmSpawn spawns[CM_SPAWNS];
 int duckCount;Vector3 ducks[MAP_DUCKS];int trapCount;unsigned seeds[FOUNDRY_TRAPS];} CmMap;
static CmMap *cmDocument,*cmRuntime;
static char cmSelectedPath[CM_PATH],cmError[256];
static bool cmTesting;
static bool cm_apply(const CmMap *m,bool editing);
static bool cm_load_selected(void);
static bool cm_spawn(int kind,int ordinal,Vector3 *pos,float *yaw);
static void ce_begin(void);
static void ce_update(float dt);
static void ce_hud(int i,int w,int h);
static void ce_world(int i);
static void ce_overlay(void);
static bool ce_frame(void);
static void ce_library(bool lobby);
static bool ce_validate_selected(bool firefight);
static void ce_shutdown(void);
