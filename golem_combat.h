/* Local Enemy Lab reactor proxy. Never serialized as a network enemy. */
#ifndef GOLEM_COMBAT_H
#define GOLEM_COMBAT_H
static bool goCombatEnabled,goCombatDead;
static bool go_combat_actor(int slot){return goCombatEnabled && slot==MAX_PLAYERS;}
static bool go_combat_motor_grace(int part);
static void go_combat_damage(float amount,bool shieldOnly);
static void go_combat_pose(void);
static void go_combat_reset(void);
static void go_combat_draw(void);
static void go_combat_launch(Vector3 origin,Vector3 *target,Vector3 *velocity);
#endif
