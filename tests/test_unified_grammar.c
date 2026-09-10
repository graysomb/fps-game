#include "../unified_sanctum_grammar.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

static void print_sanctum_ascii_blueprint(const SanctumCitadelPlan *plan) {
    printf("\n=======================================================\n");
    printf(" PRECURSOR CITADEL TOP-DOWN BLUEPRINT (Nodes: %d, Seed: %u, Steps: %d)\n",
           plan->node_count, plan->seed, plan->growth_steps);
    printf(" Census: Stone=%d, Marble=%d, Titanium=%d, HardLight=%d, Fluid=%d, Voids=%d\n",
           plan->count_stone, plan->count_marble, plan->count_titanium,
           plan->count_hardlight, plan->count_fluid, plan->count_void);
    printf(" Key: [S] Megalith, [M] Marble Colonnade, [T] Titanium Pylon,\n");
    printf("      [#] Hard-Light Bridge, [~] PBF Water, [!] Spire, [ ] Abyss Void\n");
    printf("=======================================================\n");

    int grid_w = 40;
    int grid_h = 40;
    char grid[41][41];
    for (int y = 0; y <= grid_h; ++y) {
        for (int x = 0; x <= grid_w; ++x) {
            grid[y][x] = '.';
        }
    }

    int span = 50; // Map [-50, 50] to [0, 40]
    for (int i = 0; i < plan->node_count; ++i) {
        const SanctumNode *n = &plan->nodes[i];
        int gx = (int)((float)(n->x + span) / (2.0f * span) * grid_w);
        int gz = (int)((float)(n->z + span) / (2.0f * span) * grid_h);
        if (gx < 0 || gx > grid_w || gz < 0 || gz > grid_h) continue;

        char glyph = '.';
        if (n->is_void) glyph = ' ';
        else if (n->is_fluid) glyph = '~';
        else if (n->prim_type == SANCTUM_PRIM_HARDLIGHT_BRIDGE) glyph = '#';
        else if (n->prim_type == SANCTUM_PRIM_STONE_ORTHOSTAT || n->prim_type == SANCTUM_PRIM_STONE_LINTEL) glyph = 'S';
        else if (n->prim_type == SANCTUM_PRIM_OBELISK || n->prim_type == SANCTUM_PRIM_ALTAR) glyph = 'O';
        else if (n->prim_type == SANCTUM_PRIM_MARBLE_COLUMN || n->prim_type == SANCTUM_PRIM_MARBLE_STYLOBATE ||
                 n->prim_type == SANCTUM_PRIM_MARBLE_ARCHITRAVE || n->prim_type == SANCTUM_PRIM_MARBLE_PEDIMENT) glyph = 'M';
        else if (n->prim_type == SANCTUM_PRIM_TITANIUM_PYLON || n->prim_type == SANCTUM_PRIM_TITANIUM_BUTTRESS) glyph = 'T';
        else if (n->prim_type == SANCTUM_PRIM_COFFERED_CEILING) glyph = 'C';
        else if (n->prim_type == SANCTUM_PRIM_BALUSTRADE) glyph = 'B';
        else if (n->prim_type == SANCTUM_PRIM_BRAZIER || n->prim_type == SANCTUM_PRIM_LIGHT_CHANNEL) glyph = '*';
        else if (n->prim_type == SANCTUM_PRIM_GRAVITY_CORE) glyph = '!';

        grid[gz][gx] = glyph;
    }

    for (int y = grid_h; y >= 0; y -= 2) {
        int world_z = (int)((float)y / grid_h * (2.0f * span) - span);
        printf("%3d | ", world_z);
        for (int x = 0; x <= grid_w; x += 2) {
            printf("%c ", grid[y][x]);
        }
        printf("\n");
    }
    printf("      ----------------------------------------\n");
    printf("      -50     -25      0       25      50     \n\n");
}

int main(void) {
    printf("=======================================================\n");
    printf(" RUNNING UNIFIED FIREFIGHT SYNCRETIC CITADEL TEST SUITE\n");
    printf("=======================================================\n\n");

    // [1/5] Test 4 Distinct Macro-Topologies
    printf("[1/5] Testing 4 Seed-Driven Macro-Topologies & Tactical Initializations...\n");
    const char *topo_names[] = {
        "SANCTUM_TOPO_STRONGHOLD",
        "SANCTUM_TOPO_ABYSSAL_RIFT",
        "SANCTUM_TOPO_SUNKEN_CRUCIBLE",
        "SANCTUM_TOPO_ASYMMETRIC_OUTPOST"
    };

    for (int t = 0; t < SANCTUM_TOPO_COUNT; ++t) {
        SanctumCitadelPlan plan_t;
        sanctum_plan_init(&plan_t, (uint32_t)t);
        assert(plan_t.topology == (SanctumTopology)t);
        assert(plan_t.node_count > 0);
        assert(plan_t.count_stone > 0);
        assert(plan_t.count_marble > 0);
        assert(plan_t.count_titanium > 0);
        assert(plan_t.spawn_point_count >= 3);
        assert(plan_t.supply_count >= 2);
        assert(verify_sanctum_invariants(&plan_t));

        printf("  [%d] %-32s -> Nodes: %3d, Sockets: %2d, Spawns: %d (Player: %d), Supplies: %d, Pads: %d\n",
               t, topo_names[t], plan_t.node_count, plan_t.socket_count,
               plan_t.spawn_point_count, plan_t.spawn_points[0].is_player ? 1 : 0,
               plan_t.supply_count, plan_t.jump_pad_count);
    }
    printf("-> All 4 distinct macro-topologies initialized with tactical combat elements!\n\n");

    // [2/5] Test Coalgebra Frontier Discovery across Topologies
    printf("[2/5] Testing Frontier Sockets & Affordance Inspection...\n");
    for (int t = 0; t < SANCTUM_TOPO_COUNT; ++t) {
        SanctumCitadelPlan plan_t;
        sanctum_plan_init(&plan_t, (uint32_t)t + 100);
        SanctumOpportunity opps[MAX_SANCTUM_OPPS];
        int opp_count = sanctum_coalgebra_frontier(&plan_t, opps, MAX_SANCTUM_OPPS);
        assert(opp_count > 0);
        for (int i = 0; i < opp_count; ++i) {
            assert(opps[i].socket_idx >= 0 && opps[i].socket_idx < plan_t.socket_count);
            assert(opps[i].score > 0.0f);
        }
    }
    printf("-> Frontier Coalgebra discovered valid continuation opportunities across all topologies!\n\n");

    // [3/5] Test Step-by-Step Open-Ended Growth for each topology
    printf("[3/5] Testing Open-Ended Expansion Loop (Growth Steps = 12)...\n");
    for (int t = 0; t < SANCTUM_TOPO_COUNT; ++t) {
        SanctumCitadelPlan grown = generate_unified_sanctum((uint32_t)(t * 100 + 42), 12);
        assert(grown.node_count > 0);
        assert(verify_sanctum_invariants(&grown));
        printf("  [%d] Grown %s: %d nodes, %d supplies, %d spawns\n",
               t, topo_names[t], grown.node_count, grown.supply_count, grown.spawn_point_count);
        if (t == 1) {
            print_sanctum_ascii_blueprint(&grown);
        }
    }
    printf("-> Growth expansion loops validated across all 4 map topologies!\n\n");

    // [4/5] Monte-Carlo Sweep over 50 randomized seeds across various growth steps
    printf("[4/5] Running Monte-Carlo sweep over 50 randomized seeds and growth budgets (Steps 5 to 25)...\n");
    for (uint32_t seed = 1; seed <= 50; ++seed) {
        int steps = 5 + (seed % 20);
        SanctumCitadelPlan test_plan = generate_unified_sanctum(seed * 7919 + 31, steps);
        assert(verify_sanctum_invariants(&test_plan));
        assert(test_plan.node_count > 0);
        assert(test_plan.socket_count > 0);
        assert(test_plan.spawn_point_count >= 3);
        assert(test_plan.supply_count >= 2);
    }
    printf("-> 50/50 randomized seeds across dynamic growth steps passed all invariants with 100%% stability!\n\n");

    printf("ALL UNIFIED SYNCRETIC CITADEL FIREFIGHT TESTS PASSED SUCCESSFULLY.\n");
    return 0;
}
