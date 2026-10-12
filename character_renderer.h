#ifndef CHARACTER_RENDERER_H
#define CHARACTER_RENDERER_H
#include "weapon_renderer.h"

#define CHARACTER_MODEL_SCALE (1.0f/.537f)
#define CHARACTER_HALF_WIDTH .39f
#define CHARACTER_BODY_TOP .12f
/* Position is the eye/camera anchor; soles are one unit below it. */
typedef struct CharacterPose {
    Vector3 position;
    float yaw, pitch, stride, breath, time;
    bool airborne, downed;
    Color gear, effect;
    float effectAmount;
    Vector3 soles[2];
    bool feetResolved;

} CharacterPose;
typedef struct CharacterSupport {
    bool valid;
    float height;
    unsigned long long identity;
    Vector3 origin;
} CharacterSupport;
typedef CharacterSupport (*CharacterGroundQuery)(Vector3 sole,void *context);
typedef struct CharacterAnimation {
    bool initialized, grounded, swinging[2];
    Vector3 previous, soles[2], swingStart[2], landingStart[2];
    CharacterSupport support[2];
    float phase, yaw, landing;
} CharacterAnimation;
void character_animate(CharacterAnimation *state,CharacterPose *pose,float dt,
                       CharacterGroundQuery query,void *context);
Matrix character_body_transform(CharacterPose pose);
bool character_renderer_init(void);
void character_renderer_shutdown(void);
CharacterPose character_pose(Vector3 position, float yaw, float pitch, float speed,
                             bool grounded, bool downed, Color gear, float time);
void character_draw(CharacterPose pose, bool emission);
/* Draw after opaque characters, back-to-front, with depth testing enabled. */
void character_draw_visor(CharacterPose pose);
/* Unarmed presentation pose; gameplay supplies its own weapon/melee arm pass. */
void character_draw_idle_arms(CharacterPose pose);
void character_draw_paws(WeaponPose weapon, Color gear, bool emission);
void character_draw_arm(Vector3 shoulder, Vector3 hand, Color gear);
#endif
