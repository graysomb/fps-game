/* Included after combat globals. Equipment references survive voxel compaction. */
#define ENEMY_EQUIPMENT_MAX 48
#define ENEMY_STRUCTURE_MAX 24
#define ENEMY_LOB_MAX_VOXELS 64
#define ENEMY_CHARGE_MAX 64
#define ENEMY_CONSTRUCTION_OWNER (-2) /* Temporary world material: never enqueue map restoration. */
typedef struct { uint64_t id; Vector3 offset; int group; } EnemyEquipment;
typedef struct {
    float cooldown, windup, melee, exposure, repair, spawnTimer, disabled;
    float searchTimer, stompCooldown, coreHealth, coreMax, repositionTimer, payloadFallback;
    int shieldHits[4];
    float shieldFlash[4], exposedAge, healingPulse, coreCue, blockedCue;
    bool coreOpen, supplied;
    int actions, work;
    const char *reason;
    Vector3 aim, buildOrigin;
    uint64_t selected;
    EnemyEquipment equipment[ENEMY_EQUIPMENT_MAX];
    int equipmentCount, target, state, repairs, spawned, lostMask, buildKind, buildQuarter;
    unsigned generation;
} EnemySpecialState;
typedef struct { bool active; int owner; uint64_t id,support; Vector3 pos; float fuse,beep; } EnemyCharge;
typedef struct { uint64_t id, ids[ENEMY_LOB_MAX_VOXELS]; int count, owner; float age, hitCooldown[MAX_COMBATANTS]; Vector3 target; } EnemyLob;
#define ENEMY_STRUCTURE_CELLS 512
typedef struct {
    int owner; unsigned generation; uint64_t ids[ENEMY_STRUCTURE_CELLS]; int count;
    bool tower,built; uint32_t tag; Vector3 origin; int quarter,occupant;
    float construction,repairClock,pulse,blocked,approach,bestDistance;
} EnemyStructure;
typedef struct {uint32_t tag;int waypoint,state,target,post;float stuck,blocked,lost,jump,retry;Vector3 last;} TowerRifle;
static TowerRifle towerRifles[MAX_COMBATANTS];
static uint32_t towerNextTag=1,towerNetTick;
static void tower_reset(void);
static void tower_release(int i);
static void tower_builder(int i,float dt);
static bool tower_rifle(int i,float dt);
static bool tower_stationed(int i);
static void tower_tick(float dt);
static void tower_draw(void);
static EnemySpecialState enemySpecial[MAX_COMBATANTS];
static EnemyCharge enemyCharges[ENEMY_CHARGE_MAX];
static EnemyLob enemyLobs[MAX_COMBATANTS*2];
static EnemyStructure enemyStructures[ENEMY_STRUCTURE_MAX];
static unsigned enemyGeneration;
static bool enemyCleanupPending;
static bool enemyLab=false,enemyLabPaused=false,enemyLabGod=true,enemyLabStep=false;
static int enemyLabCase=ENEMY_TYPE_SCRAPPER;
static bool enemyLabAll=false;
static float enemyLabClock;
static const char *enemy_names[ENEMY_TYPE_COUNT]={"normal","swarmer","goliath","scrapper","shieldbearer","sapper","lobber","builder","anchor","foundry"};
static const char *enemy_hints[7]={"Strip its scavenged armor before it grows","Bait its turning shield, then flank","Shoot the blinking charge before detonation","Dodge the marker or tether the incoming scrap","Builds a rifleman watchtower; kill the nearby builder to stop rapid repairs","Break the external anchor to release its allies","Break 2 orange supports; shoot the gold core until destroyed"};
static bool enemy_is_special(int i) {return i>=0 && i<MAX_COMBATANTS && players[i].enemyType>=ENEMY_TYPE_SCRAPPER;}
static float enemy_matter(EnemyType type);
static void enemy_init(int i);
static void enemy_release(int i);
static void enemy_on_death(int i);
static Color enemy_color(EnemyType t);
static void enemy_cleanup(bool removeEquipment);
static int enemy_equipment_owner(int v);
static void enemy_detach(int v);
static bool enemy_core_open(int i);
static float enemy_knockback_scale(int i);
static void enemy_update(int i,float dt);
static void enemy_tick(float dt);
static void enemy_update_equipment(int i);
static void enemy_draw(void);
static void enemy_lab_start(int which);
static void enemy_lab_exit(void);
static bool enemy_lab_controls(void);
static void enemy_lab_hud(int w,int h);
static bool enemy_run_smoke(const char *which);

static void enemy_stock_scrap(bool initial);
static void enemy_scrap_bundle(Vector3 origin);
static bool enemy_lob_hit(int voxel,int victim);
static void enemy_drop_payload(int i);

static const char *enemy_description(int type){return type==ENEMY_TYPE_GOLIATH?"Strips terrain to rebuild armor; interrupt its harvesting":type>=ENEMY_TYPE_SCRAPPER && type<ENEMY_TYPE_COUNT?enemy_hints[type-ENEMY_TYPE_SCRAPPER]:"Use cover and listen for attacks";}
