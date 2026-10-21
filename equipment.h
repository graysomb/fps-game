/* Gameplay equipment is independent of renderer enum values and bot loadouts. */
enum { EQ_DEFAULT,EQ_VOXEL,EQ_SHOTGUN,EQ_CUTTER,EQ_RAIL,EQ_ROCKET,EQ_BLADE,EQ_GUN_COUNT,EQ_EMPTY=255 };
enum { EQ_SWORD=1,EQ_ADHESIVE=2,EQ_BUILD=4 };
typedef struct {const char *name;WeaponKind visual;float cost,cooldown,range,damage,zoom;} GunDefinition;
static const GunDefinition eqGuns[EQ_GUN_COUNT]={
    {"Default",WEAPON_STANDARD,MATTER_SHOT_COST,BULLET_COOLDOWN_SECONDS,100,34,1.5f},{"Voxel",WEAPON_DYNAMIC_LAUNCHER,MATTER_SHOT_COST,BULLET_COOLDOWN_SECONDS,100,34,1.25f},
    {"Shotgun",WEAPON_IMPULSE_SHOTGUN,6,.8f,12,6.25f,1.2f},{"Cutter",WEAPON_CUTTER,4,0,24,20,2},
    {"Railgun",WEAPON_RAILGUN,15,1.8f,100,50,4},{"Rocket",WEAPON_ROCKET_LAUNCHER,25,1.5f,110,100,1.75f},{"Sword",WEAPON_SWORD,25,.8f,6,10000,1}};
enum { EQ_PELLETS=8, EQ_PANEL_SIDE=5, EQ_PANEL_CUBES=EQ_PANEL_SIDE*EQ_PANEL_SIDE*EQ_PANEL_SIDE, EQ_SWORD_CUTS=96 };
static const float eqSwitchSeconds=.25f,eqRailCharge=.65f,eqCutterUnitSeconds=1.f/24;
static const float eqShotgunCone=12,eqRocketSpeed=22,eqRocketLife=5,eqRocketInner=1,eqRocketOuter=4;
static const float eqShotgunImpulse=TETHER_THROW_IMPULSE,eqShotgunActivationBelief=1;
static const float eqGoldGrabCost=25;
static const float eqSwordReach=SWORD_REACH,eqSwordSweep=SWORD_SWEEP_DEGREES,eqSwordCycle=.8f,eqSwordDamage=10000;
typedef struct {int gun[2],selected,upgrades;float cooldown[EQ_GUN_COUNT],switchTime,charge,cutTime,buildTime;bool initialized,held,buildHeld,wasEmpty;float depleted,depletionPulse;Vector3 beamEnd;} Equipment;
static Equipment equipment[MAX_PLAYERS];
static bool equipmentLab;
static bool equipmentEmission;
typedef struct {bool active;int owner;Vector3 pos,vel;float life;} EquipmentRocket;
static EquipmentRocket eqRockets[32];
typedef struct {Vector3 start,end;float life;int gun,owner;} EquipmentTrace;
static EquipmentTrace eqTraces[32];
enum { SWORD_LANES=41 };
typedef struct {
    bool active,emitted,rewarded;
    uint32_t sequence;
    int direction,cuts;
    float age,stop[SWORD_LANES];
    Vector3 origin,forward;
    bool hit[MAX_COMBATANTS];
    uint64_t cutIds[EQ_SWORD_CUTS];
    bool terrainChanged;int minx,maxx,miny,maxy,minz,maxz;
} SwordWave;
static SwordWave swordWaves[MAX_PLAYERS];
static uint32_t swordSequences[MAX_PLAYERS];
static int swordNextDirection[MAX_PLAYERS],swordDirection[MAX_PLAYERS];
static void eq_sword_tick(float dt);
static uint32_t eqReceiveTick,eqPartTick[4];
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
static int eq_pickup_gun(int type);
static void eq_layout(void);
static void eq_world(void);
static void eq_hud(int i,int w,int h);
static bool eq_sword(int i);
static bool eq_build(int i);
static void eq_build_command(int i,const NetInputCommand *c);
static void eq_build_sample(int seat,int slot,NetInputCommand *c);
static void eq_snapshot(int peer);
static bool eq_receive(NetMessageType type,NetReader *r,uint32_t tick);
static void eq_lab(void);

/* Local blend is never overwritten by network snapshots. */
typedef struct {bool active;float blend,tapTime;bool tapDown,tapArmed,revive,toggleArmed,holdUsed;uint32_t pending;} AimState;
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
