#ifndef FPS_WEAPON_RENDERER_H
#define FPS_WEAPON_RENDERER_H

#include "raylib.h"
#include <stdint.h>

#define WEAPON_SHOT_SECONDS 0.12f

/* Presentation only: callers own gameplay, aim, and projectile spawning. */
typedef struct WeaponVisual {
    uint32_t shot_sequence;
    float shot_time;
    float claw_open;
} WeaponVisual;

typedef struct WeaponPose {
    Matrix transform;
    Vector3 muzzle;
} WeaponPose;

typedef struct WeaponBloom {
    RenderTexture2D emission;
    RenderTexture2D blur[2];
    int width, height;
} WeaponBloom;

bool weapon_renderer_init(void);
void weapon_renderer_shutdown(void);
void weapon_visual_reset(WeaponVisual *visual);
void weapon_visual_shot(WeaponVisual *visual, float now);
void weapon_visual_receive(WeaponVisual *visual, uint32_t sequence, uint8_t age_ms, float now);
void weapon_visual_update(WeaponVisual *visual, bool tether, float dt);
WeaponPose weapon_pose(Vector3 position, Vector3 forward, Vector3 right,
                       Vector3 up, bool first_person, float aspect, float melee, float recoil);
float weapon_shot_amount(const WeaponVisual *visual, float now);
void weapon_draw(WeaponPose pose, const WeaponVisual *visual, Vector3 eye,
                 float now, bool gold, bool powered_shot, bool emission_only);
bool weapon_bloom_resize(WeaponBloom *bloom, int width, int height);
void weapon_bloom_unload(WeaponBloom *bloom);
/* Call outside BeginTextureMode. The mask retains the world's depth buffer. */
void weapon_bloom_begin(WeaponBloom *bloom, RenderTexture2D world);
void weapon_clear_depth(void);
/* Ends the mask, blurs it, then leaves the requested scene target active. */
void weapon_bloom_composite(WeaponBloom *bloom, RenderTexture2D scene);

#endif
