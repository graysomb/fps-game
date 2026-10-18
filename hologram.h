/* Target references [0,MAX_COMBATANTS) are actors; the next four are decoys.
   Decoys never occupy actor slots or participate in spawn/score accounting. */
typedef struct {
    Player body;
    float health,life,cooldown;
    uint32_t generation;
    Color color;
    WeaponKind weapon;
    bool armed;
    Vector3 direction,visualTarget;
    CharacterAnimation animation;
    CharacterPose pose;
} Hologram;
static Hologram holograms[MAX_PLAYERS];
static uint32_t hologramGeneration,hologramReceiveTick;
static bool target_is_hologram(int target){return target>=MAX_COMBATANTS && target<MAX_COMBATANTS+MAX_PLAYERS;}
static Player *ai_target_body(int target){return target_is_hologram(target)?&holograms[target-MAX_COMBATANTS].body:&players[target];}
static bool hologram_deploy(int owner);
static void hologram_tick(float dt);
static void hologram_draw(bool emission);
static void hologram_hud(int owner,int w,int h);
static void hologram_damage(int target,float damage);
static void hologram_projectile(int voxel,Vector3 start,Vector3 end,float damage);
static void hologram_melee(int owner,Vector3 start,Vector3 end);
static void hologram_blast(Vector3 pos,float radius,float damage);
static void hologram_net_snapshot(int peer);
static bool hologram_net_receive(NetMessageType type,NetReader *r,uint32_t tick);
