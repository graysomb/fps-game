/* Per-player feedback; simulation and network state are sampled outside viewports. */
typedef struct {Vector3 origin; float until; bool known;} DamageCue;
static DamageCue damageCues[MAX_PLAYERS][8];
static float reviveComplete[MAX_PLAYERS];
static bool damageOriginKnown;
static Vector3 damageOrigin;
static Shader reviveShader;
static bool reviveShaderAttempted;
static uint32_t combatUiReceiveTick,chargeReceiveTick[4];
static void combat_ui_reset(void);
static void combat_ui_damage(int victim,int attacker);
static void combat_ui_world(void);
static void combat_ui_hud(int viewer,Camera3D camera,int w,int h);
static void combat_ui_abilities(int viewer,int w,int h);
static void combat_ui_shutdown(void);
static void combat_ui_snapshot(int peer);
static bool combat_ui_receive(NetMessageType type,NetReader *r,uint32_t tick);
