#ifndef FPS_WEAPON_RENDERER_H
#define FPS_WEAPON_RENDERER_H

#include "raylib.h"
#include <stdint.h>

#define WEAPON_SHOT_SECONDS 0.12f
#define WEAPON_MELEE_ACTIVE_START 0.28f
#define WEAPON_MELEE_ACTIVE_END 0.62f
#define WEAPON_MELEE_PEAK 0.52f
#define SWORD_ANIM_SECONDS .6f
#define SWORD_EMIT_SECONDS .15f
#define SWORD_SLASH_END .35f
#define SWORD_WAVE_END .5f
#define SWORD_REACH 6.f
#define SWORD_SWEEP_DEGREES 100.f
typedef struct SwordSwingSample {float sweep,fade;Vector3 direction,edge_direction;} SwordSwingSample;
SwordSwingSample weapon_sword_sample(float age,int direction);
float weapon_sword_radius(float age,float lane);

typedef enum WeaponKind {
    WEAPON_STANDARD,
    WEAPON_DYNAMIC_LAUNCHER,
    WEAPON_IMPULSE_SHOTGUN,
    WEAPON_CUTTER,
    WEAPON_ADHESIVE_GRENADE,
    WEAPON_RAILGUN,
    WEAPON_ROCKET_LAUNCHER,
    WEAPON_SWORD,
    WEAPON_REINFORCED_BUILDER,
    WEAPON_GLOB_GUN,
    WEAPON_BLACK_HOLE_BOMB,
    WEAPON_GLOCK,
    WEAPON_IGLOO_GRENADE,
    WEAPON_DMR,
    WEAPON_KIND_COUNT
} WeaponKind;

typedef enum WeaponRenderPass {
    WEAPON_COLOR,
    WEAPON_EMISSION,
    WEAPON_GLASS
} WeaponRenderPass;

#define WEAPON_LAUNCHER_CUBES 3

/* Presentation only: callers own gameplay, aim, and projectile spawning. */
typedef struct WeaponVisual {
    uint32_t shot_sequence;
    float shot_time;
    float claw_open;
    /* Cosmetic state supplied by gameplay: charge and active effect in [0,1]. */
    float charge, activity, depleted;
    /* Normalized remaining cooldown, rigid tool deployment, attachment, and
     * successful build pulse. Updated by the owner once, never by a viewport. */
    float cooldown, deployment, attached, build_pulse;
} WeaponVisual;

typedef struct WeaponPose {
    Matrix transform;
    Vector3 muzzle;
    Vector3 grip_heel;
    WeaponKind kind;
} WeaponPose;

/* The caller supplies the existing hit sweep's endpoint, in world space.
 * No target or hit-result state: mining, attacking and missing share a swing. */
typedef struct WeaponMeleePose {
    float progress;
    Vector3 strike_endpoint;
    int sword_direction; /* 0: upper-right to lower-left; 1: upper-left to lower-right. */
} WeaponMeleePose;

typedef struct WeaponBloom {
    RenderTexture2D emission;
    RenderTexture2D blur[2];
    int width, height;
} WeaponBloom;

/* Local model-space anchors, independent of inventory and gameplay. */
typedef struct WeaponAnchors { Vector3 grip, support, muzzle; } WeaponAnchors;
WeaponAnchors weapon_anchors(WeaponKind kind);
const char *weapon_visual_name(WeaponKind kind);
/* Effects are stateless samples; caller supplies the same age/time to every view.
 * Attachment transform is target-world * contact-local composed by the caller.
 * No particles, forces, or game entities are allocated by these calls. */
typedef struct WeaponEffect {
    Matrix transform;
    float age, lifetime, progress;
    uint32_t seed;
    bool attached, cutaway;
} WeaponEffect;
void weapon_draw_glob(const WeaponEffect *effect,Vector3 eye,float time,WeaponRenderPass pass);
void weapon_draw_gravity(const WeaponEffect *effect,Vector3 eye,float time,WeaponRenderPass pass);
void weapon_draw_construction(const WeaponEffect *effect,Vector3 eye,float time,WeaponRenderPass pass);
void weapon_draw_dmr_impact(const WeaponEffect *effect,Vector3 eye,float time,WeaponRenderPass pass);
void weapon_draw_dmr_reticle(int width,int height,float opacity);

bool weapon_renderer_init(void);
void weapon_renderer_shutdown(void);
void weapon_visual_reset(WeaponVisual *visual);
void weapon_visual_shot(WeaponVisual *visual, float now);
void weapon_visual_receive(WeaponVisual *visual, uint32_t sequence, uint8_t age_ms, float now);
void weapon_visual_update(WeaponVisual *visual, bool tether, float dt);
WeaponPose weapon_pose(WeaponKind kind, Vector3 position, Vector3 forward, Vector3 right,
                       Vector3 up, bool first_person, float aspect,
                       const WeaponMeleePose *melee, float recoil);
float weapon_shot_amount(const WeaponVisual *visual, float now);
void weapon_draw(WeaponPose pose, const WeaponVisual *visual, Vector3 eye,
                 float now, bool gold, WeaponRenderPass pass);
/* World pickup presentation, sharing cached weapon meshes and materials. */
void weapon_draw_rocket(Vector3 position,Vector3 direction,Vector3 eye,float now,WeaponRenderPass pass);
void weapon_draw_pickup(bool launcher, Vector3 position, Vector3 eye, float time,
                        WeaponRenderPass pass);
void weapon_draw_drain(Vector3 start, Vector3 target, float time);
void weapon_draw_drain_state(Vector3 start,Vector3 target,float time,bool active);
bool weapon_drain_shader_ready(void);
/* Local visual particle transform; never creates gameplay particles. */
Matrix weapon_cube_transform(int cube, float now);
bool weapon_bloom_resize(WeaponBloom *bloom, int width, int height);
void weapon_bloom_unload(WeaponBloom *bloom);
/* Call outside BeginTextureMode. The mask retains the world's depth buffer. */
void weapon_bloom_begin(WeaponBloom *bloom, RenderTexture2D world);
void weapon_clear_depth(void);
/* Ends the mask, blurs it, then leaves the requested scene target active. */
void weapon_bloom_composite(WeaponBloom *bloom, RenderTexture2D scene);

void weapon_view_hide(float amount);

#endif
