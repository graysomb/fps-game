#define SWARMER_CONTACT_DAMAGE 25.0f
#define SWARMER_RALLY_RADIUS 18.0f

/* Shared host-side perception and gun discipline; never serialized. */
typedef enum { AI_ALERT_NONE,AI_ALERT_SEEN,AI_ALERT_HEARD,AI_ALERT_NEAR,AI_ALERT_HIT } AiAlert;
typedef struct {
    int target,pendingTarget,intent,burstShots;
    Vector3 known,pendingPosition,lastPosition,lastTargetVelocity;
    AiAlert source,pendingSource;
    bool visible,wasVisible,initialized;
    double clock,lastFire;
    float rallyCooldown;
    float memory,commit,overrideDelay,reaction,recovery,tracking,intentTime;
    float aimYaw,aimPitch,goalYaw,goalPitch,aimTimer,error,stuck,jumpDelay;
    uint32_t random;
} AiCombat;
static AiCombat aiCombat[MAX_COMBATANTS];
static bool aiProjectileDamage;
static Vector3 aiProjectileOrigin;
static void ai_reset(int i);
static void ai_alert(int victim,int attacker,Vector3 position,AiAlert source);
static void ai_gunfire(int shooter,Vector3 origin);
static void ai_near_miss(int projectile,Vector3 start,Vector3 end);
static bool ai_shot_allowed(int i);
static void ai_shot_fired(int i);
static void ai_tick(int i,float dt);
static bool ai_ready_to_fire(int i);
static void ai_aim(int i,Vector3 target,float dt);
