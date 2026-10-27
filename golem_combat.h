/* Local Enemy Lab reactor proxy. Never serialized as a network enemy. */
#ifndef GOLEM_COMBAT_H
#define GOLEM_COMBAT_H
static bool goCombatEnabled,goCombatDead;
static int drawViewPlayerIndex;
static bool go_combat_actor(int slot){return goCombatEnabled && slot==MAX_PLAYERS;}
static bool go_combat_motor_grace(int part);
static int go_shield_part(int voxel);
static bool go_shield_hit(int voxel,float damage);
static bool go_shielded(int voxel);
static void go_shield_blast(Vector3 origin,float radius,float damage);
static bool go_drain_start(int player,int voxel);
static void go_drain_tick(int player,float dt);
static void go_drain_cancel(int player);
static Vector3 go_tether_anchor(void);
static Vector3 go_tether_aim(void);
static bool spawn_default_projectile(int owner,Vector3 start,Vector3 dir,bool dynamic);
static void acquire_tether_hit(int idx,int hit_id,Vector3 normal);
static void release_tether(int idx);
static void release_tether_velocity(int idx,const Vector3 *launch);
static void go_combat_damage(float amount,bool shieldOnly);
static void go_combat_pose(void);
static void go_combat_reset(void);
static void go_combat_draw(void);
static void go_combat_launch(Vector3 origin,Vector3 *target,Vector3 *velocity);
#endif
