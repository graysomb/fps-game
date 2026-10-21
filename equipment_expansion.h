enum { GRENADE_REGULAR, GRENADE_GRAVITY, GRENADE_IGLOO, GRENADE_KIND_COUNT };
enum { EX_GLOBS=128, EX_FIELDS=4, EX_SPLATS=64 };
typedef struct {uint32_t id;int owner;Vector3 pos,velocity;float age;} GlobProjectile;
typedef struct {uint32_t id;int owner,activated,cursor;Vector3 pos;float age;} GravityField;
typedef struct {uint32_t id;int actor;uint64_t support;Vector3 pos,local;float age;bool impact;} ExpansionSplat;
static GlobProjectile exGlobs[EX_GLOBS];
static GravityField exFields[EX_FIELDS];
static ExpansionSplat exSplats[EX_SPLATS];
static float exSuppression[MAX_COMBATANTS],exFailedBuild[MAX_COMBATANTS];
static uint32_t exIncarnation[MAX_COMBATANTS],exNextId=1;
static int exGrenades[MAX_COMBATANTS][2],exSelectedGrenade[MAX_COMBATANTS];
static float exClock;
static bool ex_suppressed(int i);
static void ex_reset(void);
static void ex_reset_player(int i);
static void ex_clear_status(int i);
static int ex_grenade_count(int i,int kind);
static void ex_cycle(int i);
static bool ex_grenade_pickup(int i,int kind);
static bool ex_throw_available(int i,int kind);
static void ex_grenade_landed(GrenadeState *g);
static void ex_grenade_deploy(GrenadeState g);
static void ex_tick(float dt);
static void ex_draw(WeaponRenderPass pass);
static void ex_hud(int i,int w,int h);
static void ex_snapshot(int peer);
static bool ex_receive(NetMessageType type,NetReader *r,uint32_t tick);
static void ex_input(int seat,NetInputCommand *c);
static void ex_layout(void);

static void ex_construction_reset(void);

static void ex_actor_pull(int i,float dt);
static bool ex_in_field(int i);

static uint32_t exLastCycleInput[MAX_PLAYERS];

static Vector3 exBotPull[MAX_COMBATANTS],exBotApplied[MAX_COMBATANTS];
