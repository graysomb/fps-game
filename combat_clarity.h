#ifndef COMBAT_CLARITY_H
#define COMBAT_CLARITY_H
#define CLARITY_GRENADES 32
#define CLARITY_EFFECTS 64
typedef struct {uint32_t id; int owner; Vector3 pos,vel; float fuse;} GrenadeState;
typedef struct {uint32_t id; int kind,type,attacker; Vector3 pos; float age,yaw;} ClarityEffect;
static GrenadeState grenades[CLARITY_GRENADES];
static Vector3 grenadeVisualTarget[CLARITY_GRENADES];
static ClarityEffect clarityEffects[CLARITY_EFFECTS];
static int grenadeInventory[MAX_COMBATANTS];
static float grenadeCooldown[MAX_COMBATANTS],killConfirm[MAX_COMBATANTS],enemyAudioTimer[MAX_COMBATANTS];
static bool explosionDamageContext;
static uint32_t clarityNextId=1;
static float clarityClock,swarmHintTimer,swarmPulseTimer;
static Vector3 swarmHints[MAX_COMBATANTS];
static bool swarmHintActive[MAX_COMBATANTS];
static void clarity_reset(void);
static void clarity_tick(float dt);
static void clarity_draw(void);
static void clarity_hud(int slot,int w,int h);
static void clarity_death(int slot,int attacker);
static void clarity_explode(Vector3 pos,int owner);
static void enemy_blast(EnemyCharge c);
static bool clarity_throw(int slot);
static bool clarity_evade(int slot,float dt);
static bool clarity_shield_voxel(int v);
static bool clarity_shield_melee(int v,int attacker);
static bool clarity_debris_hit(int v,int target);
static void play_sfx_at(SfxId id,Vector3 pos,int source,int priority,float pitch);
static void clarity_net_snapshot(int peer);
static void clarity_net_reset(void);
static void clarity_net_event(ClarityEffect e);
static bool clarity_net_receive(NetMessageType type,NetReader *r,uint32_t tick);
static void clarity_supply_layout(void);
static void clarity_supply_update(void);
static void clarity_menu(void);
static void clarity_menu_shutdown(void);
static int couch_axis(bool horizontal);
#endif
