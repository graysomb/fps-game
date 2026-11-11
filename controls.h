/* Local device profiles. Stable action names are the persistence interface. */
#define CT_PROFILES 16
#define CT_NAME 48
#define CT_COMMON 1
#define CT_GAME 2
#define CT_EDIT 4
#define CT_LAB 8
#define CT_PREP 16
#define CT_TEST 32
#define CT_ACTIONS(X) \
 X(FORWARD,"Move forward",CT_COMMON,KEY_W,KEY_I,0,-1) \
 X(BACK,"Move backward",CT_COMMON,KEY_S,KEY_K,0,-1) \
 X(LEFT,"Move left",CT_COMMON,KEY_A,KEY_J,0,-1) \
 X(RIGHT,"Move right",CT_COMMON,KEY_D,KEY_L,0,-1) \
 X(LOOK_LEFT,"Look left",CT_COMMON,KEY_F,KEY_LEFT,0,-1) \
 X(LOOK_RIGHT,"Look right",CT_COMMON,KEY_H,KEY_RIGHT,0,-1) \
 X(LOOK_UP,"Look up",CT_COMMON,KEY_T,KEY_UP,0,-1) \
 X(LOOK_DOWN,"Look down",CT_COMMON,KEY_G,KEY_DOWN,0,-1) \
 X(FIRE,"Fire / charged slash",CT_GAME,KEY_LEFT_CONTROL,KEY_RIGHT_CONTROL,GAMEPAD_BUTTON_RIGHT_TRIGGER_2,MOUSE_BUTTON_LEFT) \
 X(JUMP,"Jump",CT_GAME,KEY_SPACE,KEY_RIGHT_SHIFT,GAMEPAD_BUTTON_RIGHT_FACE_RIGHT,-1) \
 X(MELEE,"Melee / harvest",CT_GAME,KEY_Z,KEY_M,GAMEPAD_BUTTON_RIGHT_FACE_DOWN,MOUSE_BUTTON_MIDDLE) \
 X(TETHER,"Tether",CT_GAME,KEY_R,KEY_P,GAMEPAD_BUTTON_LEFT_TRIGGER_2,MOUSE_BUTTON_RIGHT) \
 X(BUILD,"Build (release to place)",CT_GAME,KEY_E,KEY_O,GAMEPAD_BUTTON_RIGHT_TRIGGER_1,-1) \
 X(AIM,"Aim / zoom",CT_GAME,KEY_LEFT_SHIFT,KEY_ENTER,GAMEPAD_BUTTON_RIGHT_THUMB,-1) \
 X(SWITCH,"Switch weapon",CT_GAME,KEY_Q,KEY_U,GAMEPAD_BUTTON_LEFT_FACE_LEFT,-1) \
 X(EQUIP,"Pick up / replace weapon",CT_GAME,KEY_B,KEY_Y,0,-1) \
 X(REVIVE,"Revive teammate",CT_GAME,KEY_X,KEY_COMMA,0,-1) \
 X(INTERACT,"Tap switch / hold interact",CT_GAME,0,0,GAMEPAD_BUTTON_RIGHT_FACE_LEFT,-1) \
 X(GRENADE,"Throw grenade",CT_GAME,KEY_C,KEY_PERIOD,GAMEPAD_BUTTON_LEFT_FACE_UP,-1) \
 X(CYCLE,"Cycle grenade",CT_GAME,KEY_TAB,KEY_SEMICOLON,GAMEPAD_BUTTON_LEFT_FACE_DOWN,-1) \
 X(HOLOGRAM,"Hologram",CT_GAME,KEY_V,KEY_SLASH,GAMEPAD_BUTTON_LEFT_TRIGGER_1,-1) \
 X(PREP_LEFT,"Previous upgrade",CT_PREP,KEY_F,KEY_LEFT,GAMEPAD_BUTTON_LEFT_FACE_LEFT,-1) \
 X(PREP_RIGHT,"Next upgrade",CT_PREP,KEY_H,KEY_RIGHT,GAMEPAD_BUTTON_LEFT_FACE_RIGHT,-1) \
 X(READY,"Ready / skip loot",CT_PREP,KEY_SPACE,KEY_RIGHT_SHIFT,GAMEPAD_BUTTON_RIGHT_FACE_DOWN,-1) \
 X(PREP_MENU,"Upgrade menu",CT_PREP,KEY_V,KEY_N,GAMEPAD_BUTTON_RIGHT_FACE_UP,-1) \
 X(UP,"Fly up",CT_EDIT,KEY_E,KEY_O,GAMEPAD_BUTTON_RIGHT_FACE_DOWN,-1) \
 X(DOWN,"Fly down",CT_EDIT,KEY_Q,KEY_U,GAMEPAD_BUTTON_RIGHT_FACE_RIGHT,-1) \
 X(FAST,"Fly faster",CT_EDIT,KEY_LEFT_SHIFT,KEY_RIGHT_SHIFT,GAMEPAD_BUTTON_LEFT_THUMB,-1) \
 X(PLACE,"Apply / place",CT_EDIT,KEY_LEFT_CONTROL,KEY_KP_0,GAMEPAD_BUTTON_RIGHT_TRIGGER_2,MOUSE_BUTTON_LEFT) \
 X(ERASE,"Erase",CT_EDIT,KEY_LEFT_ALT,KEY_KP_DECIMAL,GAMEPAD_BUTTON_LEFT_TRIGGER_2,MOUSE_BUTTON_RIGHT) \
 X(PALETTE,"Palette",CT_EDIT,KEY_TAB,KEY_KP_MULTIPLY,GAMEPAD_BUTTON_RIGHT_FACE_UP,-1) \
 X(PROPERTIES,"Tool properties",CT_EDIT,KEY_R,KEY_KP_ENTER,GAMEPAD_BUTTON_RIGHT_FACE_LEFT,-1) \
 X(SAMPLE,"Sample / select object",CT_EDIT,0,0,GAMEPAD_BUTTON_RIGHT_THUMB,MOUSE_BUTTON_MIDDLE) \
 X(UNDO,"Undo",CT_EDIT,KEY_Z,KEY_KP_DIVIDE,GAMEPAD_BUTTON_LEFT_TRIGGER_1,-1) \
 X(REDO,"Redo",CT_EDIT,KEY_Y,KEY_KP_DIVIDE,GAMEPAD_BUTTON_RIGHT_TRIGGER_1,-1) \
 X(SIZE_UP,"Increase brush",CT_EDIT,KEY_RIGHT_BRACKET,KEY_KP_ADD,GAMEPAD_BUTTON_LEFT_FACE_UP,-1) \
 X(SIZE_DOWN,"Decrease brush",CT_EDIT,KEY_LEFT_BRACKET,KEY_KP_SUBTRACT,GAMEPAD_BUTTON_LEFT_FACE_DOWN,-1) \
 X(DIST_UP,"Increase distance",CT_EDIT,0,KEY_KP_9,GAMEPAD_BUTTON_LEFT_FACE_RIGHT,-1) \
 X(DIST_DOWN,"Decrease distance",CT_EDIT,0,KEY_KP_7,GAMEPAD_BUTTON_LEFT_FACE_LEFT,-1) \
 X(LAB_NEXT,"Next enemy / golem action",CT_LAB,KEY_F3,KEY_F3,GAMEPAD_BUTTON_LEFT_FACE_RIGHT,-1) \
 X(LAB_RESET,"Reset encounter",CT_LAB,KEY_F4,KEY_F4,GAMEPAD_BUTTON_RIGHT_FACE_UP,-1) \
 X(LAB_PAUSE,"Pause simulation",CT_LAB,KEY_F5,KEY_F5,GAMEPAD_BUTTON_RIGHT_FACE_LEFT,-1) \
 X(LAB_STEP,"Step simulation",CT_LAB,KEY_F6,KEY_F6,GAMEPAD_BUTTON_RIGHT_TRIGGER_1,-1) \
 X(LAB_GOD,"Invulnerability",CT_LAB,KEY_F7,KEY_F7,GAMEPAD_BUTTON_MIDDLE_LEFT,-1) \
 X(LAB_FIXTURE,"Next fixture / regeneration",CT_LAB,KEY_F8,KEY_F8,GAMEPAD_BUTTON_LEFT_FACE_DOWN,-1) \
 X(LAB_OVERLAY,"Diagnostic overlay",CT_LAB,KEY_F9,KEY_F9,0,-1) \
 X(LAB_ARM,"Remove arm / Matter",CT_LAB,KEY_F10,KEY_F10,GAMEPAD_BUTTON_LEFT_TRIGGER_1,-1) \
 X(LAB_LEG,"Remove leg",CT_LAB,KEY_F11,KEY_F11,0,-1) \
 X(LAB_CHEST,"Remove chest",CT_LAB,KEY_F12,KEY_F12,0,-1) \
 X(JOIN_EDITOR,"Join keyboard editor",CT_EDIT,KEY_F9,KEY_F9,0,-1) \
 X(TEST_RETURN,"Return to editor",CT_TEST,KEY_F10,KEY_F10,0,-1) \
 X(TEST_TOOLS,"Sandbox tools",CT_TEST,KEY_F5,KEY_F5,0,-1) \
 X(CONSOLE,"Debug console",CT_COMMON,KEY_GRAVE,KEY_GRAVE,0,-1) \
 X(SETTINGS,"Controls & Settings",CT_COMMON,0,0,0,-1) \
 X(PAUSE,"Pause / session menu",CT_COMMON,KEY_PAUSE,KEY_PAUSE,GAMEPAD_BUTTON_MIDDLE_RIGHT,-1)
#define CT_ENUM(id,label,ctx,k1,k2,p,m) CT_##id,
typedef enum {CT_ACTIONS(CT_ENUM) CT_COUNT} CtAction;
#undef CT_ENUM
/* Device: 0 keyboard, 1 mouse, 2 controller. Kind: none/key/mouse/wheel/button. */
typedef struct {int kind,code,mods;} CtBinding;
typedef struct {float x,y,zoom,accel,brake,inner,outer,curve,trigger,assist,smooth,moveInner,moveOuter;int invertX,invertY,direct,swap;} CtTuning;
typedef struct {char name[CT_NAME];CtBinding bindings[3][CT_COUNT][2];CtTuning tuning[3];int aimToggle[3],tetherToggle[3],combined[3],layout;} CtProfile;
typedef struct {int version,count,assigned[4];CtProfile profiles[CT_PROFILES];} CtLibrary;
typedef struct {bool down,pressed,released;} CtActionState;
static CtLibrary ctLibrary;
static bool ctReady,ctMenuOpen,ctPauseOpened;
static CtActionState ctState[4][CT_COUNT];
static void ct_load(void);
static void ct_frame(float dt);
static bool ct_menu_frame(void);
static void ct_open(void);
static bool ct_down(int seat,int action);
static bool ct_pressed(int seat,int action);
static bool ct_actor_pressed(int actor,int action);
static const char *ct_hint(int actor,int action);
static int ct_seat(int actor);
static NetInputCommand ct_command(int seat,int slot,float dt);
static void ct_move(int slot,float dt);
static void ct_local_actions(float dt);
static void ct_edit_look(int seat,float dt);
static Vector2 ct_move_axes(int seat);
static Vector2 ct_look_axes(int seat);
static float ct_assist(int slot);
static void ct_cancel_all(void);

static void ct_mouse_reset(void);
