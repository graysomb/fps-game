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
    printf(" RUNNING OPEN-ENDED UNIFIED SYNCRETIC CITADEL TEST SUITE\n");
    printf("=======================================================\n\n");

    // [1/4] Test Plan Initialization & Root Nexus
    printf("[1/4] Testing Root Nexus Initialization...\n");
    SanctumCitadelPlan plan0;
    sanctum_plan_init(&plan0, 1337);
    assert(plan0.node_count > 0);
    assert(plan0.socket_count >= 4);
    assert(plan0.count_stone > 0);
    assert(plan0.count_marble > 0);
    assert(plan0.count_titanium > 0);
    assert(plan0.count_fluid > 0);
    assert(verify_sanctum_invariants(&plan0));
    printf("-> Root Nexus passed all initialization and multi-style census invariants!\n\n");

    // [2/4] Test Coalgebra Frontier Discovery
    printf("[2/4] Testing Frontier Sockets & Affordance Inspection...\n");
    SanctumOpportunity opps[MAX_SANCTUM_OPPS];
    int opp_count = sanctum_coalgebra_frontier(&plan0, opps, MAX_SANCTUM_OPPS);
    assert(opp_count > 0);
    for (int i = 0; i < opp_count; ++i) {
        assert(opps[i].socket_idx >= 0 && opps[i].socket_idx < plan0.socket_count);
        assert(opps[i].score > 0.0f);
    }
    printf("-> Frontier Coalgebra discovered %d valid continuation opportunities!\n\n", opp_count);

    // [3/4] Test Step-by-Step Open-Ended Growth
    printf("[3/4] Testing Open-Ended Expansion Loop (Growth Steps = 15)...\n");
    SanctumCitadelPlan plan15 = generate_unified_sanctum(2552, 15);
    assert(plan15.node_count > plan0.node_count);
    assert(verify_sanctum_invariants(&plan15));
    print_sanctum_ascii_blueprint(&plan15);
    printf("-> 15-step expansion successfully grew citadel from %d to %d nodes!\n\n",
           plan0.node_count, plan15.node_count);

    // [4/4] Monte-Carlo Sweep over 50 randomized seeds across various growth steps
    printf("[4/4] Running Monte-Carlo sweep over 50 randomized seeds and growth budgets (Steps 5 to 30)...\n");
    for (uint32_t seed = 1; seed <= 50; ++seed) {
        int steps = 5 + (seed % 25);
        SanctumCitadelPlan test_plan = generate_unified_sanctum(seed * 7919 + 31, steps);
        assert(verify_sanctum_invariants(&test_plan));
        assert(test_plan.node_count > 0);
        assert(test_plan.socket_count > 0);
    }
    printf("-> 50/50 randomized seeds across dynamic growth steps passed all invariants with 100%% stability!\n\n");

    printf("ALL UNIFIED SYNCRETIC CITADEL TESTS PASSED SUCCESSFULLY.\n");
    return 0;
}
