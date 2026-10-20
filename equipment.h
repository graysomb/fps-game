/* Gameplay equipment is independent of renderer enum values and bot loadouts. */
enum { EQ_DEFAULT,EQ_VOXEL,EQ_SHOTGUN,EQ_CUTTER,EQ_RAIL,EQ_ROCKET,EQ_GUN_COUNT,EQ_EMPTY=255 };
enum { EQ_SWORD=1,EQ_ADHESIVE=2,EQ_BUILD=4 };
typedef struct {const char *name;WeaponKind visual;float cost,cooldown,range,damage;} GunDefinition;
static const GunDefinition eqGuns[EQ_GUN_COUNT]={
    {"Default",WEAPON_STANDARD,MATTER_SHOT_COST,BULLET_COOLDOWN_SECONDS,100,34},{"Voxel",WEAPON_DYNAMIC_LAUNCHER,MATTER_SHOT_COST,BULLET_COOLDOWN_SECONDS,100,34},
    {"Shotgun",WEAPON_IMPULSE_SHOTGUN,6,.8f,12,5},{"Cutter",WEAPON_CUTTER,8,0,12,20},
    {"Railgun",WEAPON_RAILGUN,15,1.8f,100,50},{"Rocket",WEAPON_ROCKET_LAUNCHER,25,1.5f,110,100}};
enum { EQ_PELLETS=8, EQ_PANEL_CUBES=3, EQ_SWORD_CUTS=4 };
static const float eqSwitchSeconds=.25f,eqRailCharge=.65f,eqCutterUnitSeconds=.25f;
static const float eqShotgunCone=12,eqRocketSpeed=22,eqRocketLife=5,eqRocketInner=1,eqRocketOuter=4;
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
