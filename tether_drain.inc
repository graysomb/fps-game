/* Actor draining is retired. Keep clearing transient locks for reset/network compatibility. */
static void update_tether_drains(float dt) {
    (void)dt;
    for(int i=0;i<activePlayers;i++) {
        players[i].drainTarget=0;
        go_drain_cancel(i);
    }
}
