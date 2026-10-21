/* Gameplay equipment is independent of renderer enum values and bot loadouts. */
enum { EQ_DEFAULT,EQ_VOXEL,EQ_SHOTGUN,EQ_CUTTER,EQ_RAIL,EQ_ROCKET,EQ_GUN_COUNT,EQ_EMPTY=255 };
enum { EQ_SWORD=1,EQ_ADHESIVE=2,EQ_BUILD=4 };
typedef struct {const char *name;WeaponKind visual;float cost,cooldown,range,damage,zoom;} GunDefinition;
static const GunDefinition eqGuns[EQ_GUN_COUNT]={
    {"Default",WEAPON_STANDARD,MATTER_SHOT_COST,BULLET_COOLDOWN_SECONDS,100,34,1.5f},{"Voxel",WEAPON_DYNAMIC_LAUNCHER,MATTER_SHOT_COST,BULLET_COOLDOWN_SECONDS,100,34,1.25f},
    {"Shotgun",WEAPON_IMPULSE_SHOTGUN,6,.8f,12,5,1.2f},{"Cutter",WEAPON_CUTTER,8,0,12,20,2},
    {"Railgun",WEAPON_RAILGUN,15,1.8f,100,50,4},{"Rocket",WEAPON_ROCKET_LAUNCHER,25,1.5f,110,100,1.75f}};
enum { EQ_PELLETS=8, EQ_PANEL_CUBES=3, EQ_SWORD_CUTS=4 };
static const float eqSwitchSeconds=.25f,eqRailCharge=.65f,eqCutterUnitSeconds=.25f;
static const float eqShotgunCone=12,eqRocketSpeed=22,eqRocketLife=5,eqRocketInner=1,eqRocketOuter=4;
static const float eqShotgunImpulse=48,eqShotgunActivationBelief=1;
static const float eqSwordReach=2.7f,eqSwordSweep=70,eqSwordCycle=.8f,eqSwordDamage=40;
typedef struct {int gun[2],selected,upgrades;float cooldown[EQ_GUN_COUNT],switchTime,charge,cutTime,buildTime;bool initialized,held;Vector3 beamEnd;} Equipment;
static Equipment equipment[MAX_PLAYERS];
static bool equipmentLab;
static bool equipmentEmission;
typedef struct {bool active;int owner;Vector3 pos,vel;float life;} EquipmentRocket;
static EquipmentRocket eqRockets[32];
typedef struct {Vector3 start,end;float life;int gun;} EquipmentTrace;
static EquipmentTrace eqTraces[32];
static uint32_t eqReceiveTick,eqPartTick[3];
static void eq_reset(int i);
static void eq_reset_all(void);
static int eq_gun(int i);
static bool eq_upgrade(int i,int flag);
static bool eq_fire(int i);
static void eq_command(int i,const NetInputCommand *c,float dt);
static void eq_tick(float dt);
static void eq_local(float dt);
static bool eq_pickup(int i,int type,bool replace);
static int eq_pickup_kind(int type);
static void eq_layout(void);
static void eq_world(void);
static void eq_hud(int i,int w,int h);
static void eq_sword(int i);
static bool eq_build(int i);
static void eq_snapshot(int peer);
static bool eq_receive(NetMessageType type,NetReader *r,uint32_t tick);
static void eq_lab(void);

/* Local blend is never overwritten by network snapshots. */
typedef struct {bool active;float blend,tapTime;bool tapDown,tapArmed,revive,toggleArmed;uint32_t pending;} AimState;
static AimState aimStates[MAX_PLAYERS];
static uint32_t netLastAimInput[MAX_PLAYERS];
static void aim_cancel(int i);
static void aim_command(int i,const NetInputCommand *c);
static void aim_sample(int seat,int slot,float dt,NetInputCommand *c);
static void aim_interact(int i,bool down,float dt,NetInputCommand *c);
static void aim_tick(float dt);
static void aim_focus_tick(void);
static float aim_zoom(int i);
static float aim_fov(int i);
static float aim_scale(int i);
static float aim_player_scale(const Player *p);
static float aim_seat_scale(int seat);
static void aim_overlay(int i,int w,int h);
static void aim_shutdown(void);
