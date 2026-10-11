#ifndef CHARACTER_RENDERER_H
#define CHARACTER_RENDERER_H
#include "weapon_renderer.h"

/* Presentation only. Position is the existing collision-box center. */
typedef struct CharacterPose {
    Vector3 position;
    float yaw, pitch, stride, breath, time;
    bool airborne, downed;
    Color gear, effect;
    float effectAmount;
} CharacterPose;
bool character_renderer_init(void);
void character_renderer_shutdown(void);
CharacterPose character_pose(Vector3 position, float yaw, float pitch, float speed,
                             bool grounded, bool downed, Color gear, float time);
void character_draw(CharacterPose pose, bool emission);
/* Draw after opaque characters, back-to-front, with depth testing enabled. */
void character_draw_visor(CharacterPose pose);
void character_draw_paws(WeaponPose weapon, Color gear, bool emission);
void character_draw_arm(Vector3 shoulder, Vector3 hand, Color gear);
#endif
