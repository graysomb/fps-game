/* Included after combat globals. Equipment references survive voxel compaction. */
#define ENEMY_EQUIPMENT_MAX 48
#define ENEMY_STRUCTURE_MAX 24
#define ENEMY_CONSTRUCTION_OWNER (-2) /* Temporary world material: never enqueue map restoration. */
typedef struct { uint64_t id; Vector3 offset; int group; } EnemyEquipment;
typedef struct {
    float cooldown, windup, melee, exposure, repair, spawnTimer, disabled;
    Vector3 aim, buildOrigin;
    uint64_t selected;
    EnemyEquipment equipment[ENEMY_EQUIPMENT_MAX];
    int equipmentCount, target, state, repairs, spawned, lostMask, buildKind, buildQuarter;
    unsigned generation;
} EnemySpecialState;
typedef struct { bool active; int owner; uint64_t id,support; Vector3 pos; float fuse,beep; } EnemyCharge;
typedef struct { uint64_t id; int owner; float age; Vector3 target; } EnemyLob;
typedef struct { int owner; unsigned generation; uint64_t ids[96]; int count; } EnemyStructure;
static EnemySpecialState enemySpecial[MAX_COMBATANTS];
static EnemyCharge enemyCharges[MAX_COMBATANTS];
static EnemyLob enemyLobs[MAX_COMBATANTS*2];
static EnemyStructure enemyStructures[ENEMY_STRUCTURE_MAX];
static unsigned enemyGeneration;
static bool enemyCleanupPending;
static bool enemyLab=false,enemyLabPaused=false,enemyLabGod=true,enemyLabStep=false;
static int enemyLabCase=ENEMY_TYPE_SCRAPPER;
static bool enemyLabAll=false;
static float enemyLabClock;
static const char *enemy_names[ENEMY_TYPE_COUNT]={"normal","swarmer","goliath","scrapper","shieldbearer","sapper","lobber","builder","anchor","foundry"};
static const char *enemy_hints[7]={"Strip its scavenged armor before it grows","Flank the slow-turning shield","Shoot the blinking charge before detonation","Dodge the marker or tether the incoming scrap","Interrupt construction; harvest its cover","Break the external anchor to release its allies","Break two supports, then expose and finish the core"};
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
