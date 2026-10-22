#ifndef GOLEM_RENDERER_H
#define GOLEM_RENDERER_H
#include "weapon_renderer.h"

typedef enum GolemJoint {
    GOLEM_ROOT, GOLEM_PELVIS, GOLEM_TORSO, GOLEM_HEAD,
    GOLEM_ARM_L, GOLEM_FOREARM_L, GOLEM_HAND_L,
    GOLEM_ARM_R, GOLEM_FOREARM_R, GOLEM_HAND_R,
    GOLEM_THIGH_L, GOLEM_SHIN_L, GOLEM_FOOT_L,
    GOLEM_THIGH_R, GOLEM_SHIN_R, GOLEM_FOOT_R, GOLEM_JOINT_COUNT
} GolemJoint;
typedef enum GolemAction { GOLEM_IDLE, GOLEM_WALK, GOLEM_SLAM, GOLEM_SWEEP, GOLEM_THROW } GolemAction;
typedef enum GolemSection { GOLEM_LEFT_ARM, GOLEM_RIGHT_ARM, GOLEM_LEFT_LEG, GOLEM_RIGHT_LEG, GOLEM_CHEST, GOLEM_SECTION_COUNT } GolemSection;
/* Event data is caller-owned. Capture before removing a section; sampling never
 * mutates it. A throw uses forearmOnly; ordinary arm loss removes at the shoulder. */
typedef struct GolemDamage {
    bool removed, captured, forearmOnly;
    float age, regeneration, recovery;
    Matrix detached[3];
    Vector3 origin;
    float yaw, floorY;
} GolemDamage;
typedef struct GolemInput {
    Vector3 position;
    float yaw, scale, time, travel, speed, actionProgress;
    GolemAction action;
    int side;
    GolemDamage damage[GOLEM_SECTION_COUNT];
    int hitJoint; /* -1 disables the localized hit flash. */
    float hitAge;
    uint32_t seed;
    /* Optional caller-resolved contacts enable planted feet on arbitrary paths.
     * Contacts are ankle positions, in world units, not sole positions. */
    bool feetResolved;
    Vector3 feet[2];
    float footYaw[2];
} GolemInput;
typedef struct GolemPose {
    GolemInput input;
    Matrix joint[GOLEM_JOINT_COUNT], chest, core, eye;
    float visibility[GOLEM_JOINT_COUNT], chestVisibility;
    float kneel, impact, sweep;
    Vector3 contact, hand[2], knee[2];
} GolemPose;
typedef struct GolemStats { int cubes, triangles, models; bool shader; } GolemStats;

GolemInput golem_input(void);
float golem_action_seconds(GolemAction action);
GolemPose golem_sample(const GolemInput *input);
void golem_capture_damage(const GolemPose *pose,GolemSection section,bool forearmOnly,GolemDamage *event);
bool golem_renderer_init(const char *shaderDirectory);
void golem_renderer_shutdown(void);
GolemStats golem_renderer_stats(void);
/* Sample once and reuse exactly the same pose in every viewport/pass. */
void golem_draw(const GolemPose *pose,Vector3 eye,WeaponRenderPass pass,bool effects);
void golem_draw_skeleton(const GolemPose *pose);
#endif
