/* Run-local earned equipment. owner -1 is on the ground; -2 is consumed. */
typedef struct {uint32_t id;int type,owner,slot;Vector3 pos;float groundAge;} FfItem;
typedef struct {int type,wave,chassis;char tier;} FfEquipmentRule;
#define FP_DROP_SECONDS 60.0f
#define FP_DUPLICATE_MATTER 25.0f
static bool fpCacheReady;
static Vector3 fpCacheOrigin,fpCachePositions[16];
static void fp_cache_place(void);
static FfItem *fpIncoming;
static FfItem *fpItems;static size_t fpCount,fpCapacity;
static uint32_t fpNextId=1,fpRevision=1,fpGeneration=1,fpUnlocked;
static int fpCarrierType[MAX_COMBATANTS];
static float fpCarrierClock[MAX_COMBATANTS];
static uint64_t fpCarrierPanel[MAX_COMBATANTS];
static int fpPending[8],fpPendingCount,fpRepeatCursor;
typedef struct {int mate,target;Vector3 known,patrol;float memory,patrolTime;bool regroup;} FfRiflePair;
/* mate is actor index + 1, so zero-initialized state is unpaired. */
static FfRiflePair fpPairs[MAX_COMBATANTS];
static void fp_pair_clear(int i);
static void fp_pairs_tick(float dt);
static bool fp_pair_valid(int i);
static void fp_pair_link(int a,int b);
static bool fp_spawn_pair(int living,int cap);
static bool fp_pair_pending(void);
static bool firefightCheckpoints=false,fpRestoring=false;
static float fpRestartTimer,fpBannerTimer,fpSupplyTimer;
static int fpBannerQueue[32],fpBannerRead,fpBannerWrite;
typedef struct {
    bool valid;int wave,lives,score,kills,repeatCursor;uint32_t nextId,unlocked;
    int remaining[FF_TYPE_COUNT],pending[8],pendingCount,total;
    unsigned upgrades[MAX_PLAYERS][3];Player humans[MAX_PLAYERS];Equipment equipment[MAX_PLAYERS];
    int grenades[MAX_PLAYERS],special[MAX_PLAYERS][2],selected[MAX_PLAYERS],revives[MAX_PLAYERS];
    FfItem *items;size_t count;
} FfCheckpoint;
static FfCheckpoint fpCheckpoint;
/* Disk saves are armed only by an explicit offline New Game / Continue. */
static bool fpSaveEnabled;
static const char *fpSaveError;
static void fp_save_write(void);
static void fp_save_clear(void);
static void fp_start_new_game(void);
static bool fp_active(void);
static int fp_rule(int type);
static int fp_team_size(void);
static bool fp_carrier(int i);
static void fp_new_run(void);
static void fp_shutdown(void);
static void fp_prepare_wave(int wave);
static void fp_prepare_round(int wave,int round);
static void fp_spawn(int i,int chassis);
static void fp_death(int i,Vector3 position);
static int fp_near(int i,bool automatic);
static int fp_cache_near(int i);
static bool fp_interact(int i);
static void fp_tick(float dt);
static bool fp_retry(void);
static void fp_capture(void);
static void fp_restore(void);
static void fp_draw(WeaponRenderPass pass);
static void fp_hud(int i,int w,int h);
static void fp_radar(int i,int cx,int cy,int radius,float range);
static void fp_update_carrier(int i,float dt);
static int fp_carrier_visual(int i);
static void fp_snapshot(int peer);
static bool fp_receive(NetMessageType type,NetReader *reader,uint32_t tick);
static void fp_network_reset(void);
static void net_send_world_to(int slot);
