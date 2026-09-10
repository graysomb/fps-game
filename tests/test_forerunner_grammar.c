#include "../forerunner_grammar.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// Helper to verify bilateral reflection symmetry across X = 0
static bool verify_forerunner_symmetry(const ForerunnerPlan *plan) {
    for (int i = 0; i < plan->node_count; ++i) {
        const ForerunnerNode *n1 = &plan->nodes[i];
        if (n1->is_void) continue;

        int cm1_x = (n1->box.min_x + n1->box.max_x) / 2;
        if (abs(cm1_x) <= 1) continue; // Centered on axis, intrinsically symmetric

        // Find symmetric partner
        bool found_partner = false;
        for (int j = 0; j < plan->node_count; ++j) {
            const ForerunnerNode *n2 = &plan->nodes[j];
            if (n2->type != n1->type) continue;
            int cm2_x = (n2->box.min_x + n2->box.max_x) / 2;
            if (abs(cm1_x + cm2_x) <= 2 &&
                n1->box.min_y == n2->box.min_y && n1->box.max_y == n2->box.max_y &&
                abs(n1->box.min_z - n2->box.min_z) <= 2) {
                found_partner = true;
                break;
            }
        }

        if (!found_partner) {
            printf("Error: Node %d (type %d) at X=%d has no symmetric reflection mate!\n",
                   i, n1->type, cm1_x);
            return false;
        }
    }
    return true;
}

// Helper to verify mass/void polarity balance
static bool verify_mass_void_polarity(const ForerunnerPlan *plan) {
    bool has_void = false;
    bool has_mass = false;
    for (int i = 0; i < plan->node_count; ++i) {
        if (plan->nodes[i].is_void) has_void = true;
        else has_mass = true;
    }
    return has_void && has_mass;
}

// Helper to verify skybridge alignment across the chasm
static bool verify_skybridge_alignment(const ForerunnerPlan *plan) {
    if (!plan->has_bridge) return true;
    for (int i = 0; i < plan->node_count; ++i) {
        const ForerunnerNode *n = &plan->nodes[i];
        if (n->type == FORERUNNER_PRIM_HARDLIGHT_BRIDGE) {
            int cm_x = (n->box.min_x + n->box.max_x) / 2;
            if (abs(cm_x) > 1) {
                printf("Error: Hard-light bridge deviated from axis! X=%d\n", cm_x);
                return false;
            }
        }
    }
    return true;
}

int main(void) {
    printf("=======================================================\n");
    printf(" RUNNING FORERUNNER MASS-VOID TRI-GRAMMAR TEST SUITE   \n");
    printf("=======================================================\n\n");

    // Test 1: Cartographer Archetype
    printf("[1/4] Testing Archetype: Cartographer Megastructure...\n");
    ForerunnerPlan p0 = generate_forerunner_structure(117, FORERUNNER_ARCHETYPE_CARTOGRAPHER, FORERUNNER_STAGE_APEX);
    assert(p0.has_chasm);
    assert(p0.has_gateway);
    assert(p0.has_bridge);
    assert(p0.has_core);
    assert(verify_mass_void_polarity(&p0));
    assert(verify_forerunner_symmetry(&p0));
    assert(verify_skybridge_alignment(&p0));
    printf("-> Cartographer Archetype passed all invariants!\n\n");

    // Test 2: Crossroads Archetype (Cruciform Abyss & Bi-Level Bridges)
    printf("[2/4] Testing Archetype: Crossroads (Bi-Level Bridges)...\n");
    ForerunnerPlan p1 = generate_forerunner_structure(343, FORERUNNER_ARCHETYPE_CROSSROADS, FORERUNNER_STAGE_APEX);
    assert(p1.has_chasm);
    assert(p1.has_bridge);
    assert(p1.has_core);
    assert(verify_mass_void_polarity(&p1));
    assert(verify_forerunner_symmetry(&p1));
    assert(verify_skybridge_alignment(&p1));
    printf("-> Crossroads Archetype passed bi-level crossing & symmetry invariants!\n\n");

    // Test 3: Crucible Archetype (Radial Pylon Ring & Gravity Pit)
    printf("[3/4] Testing Archetype: Crucible (Radial Pylon Ring)...\n");
    ForerunnerPlan p2 = generate_forerunner_structure(777, FORERUNNER_ARCHETYPE_CRUCIBLE, FORERUNNER_STAGE_APEX);
    assert(p2.has_chasm);
    assert(p2.has_core);
    assert(verify_mass_void_polarity(&p2));
    assert(verify_forerunner_symmetry(&p2));
    printf("-> Crucible Archetype passed radial ring & pit invariants!\n\n");

    // Test 4: Spire Citadel Archetype (Pinnacle & Flying Buttresses)
    printf("[4/4] Testing Archetype: Spire Citadel (Flying Buttresses)...\n");
    ForerunnerPlan p3 = generate_forerunner_structure(2552, FORERUNNER_ARCHETYPE_SPIRE, FORERUNNER_STAGE_APEX);
    assert(p3.has_chasm);
    assert(p3.has_core);
    assert(verify_mass_void_polarity(&p3));
    assert(verify_forerunner_symmetry(&p3));
    print_forerunner_blueprint(&p3);
    printf("-> Spire Citadel Archetype passed flying buttress & pinnacle invariants!\n\n");

    // Monte Carlo 50-Seed Invariant Sweep across all 4 archetypes and 5 stages
    printf("Running Monte-Carlo sweep over 50 random seeds and all 4 archetypes...\n");
    for (uint32_t seed = 1; seed <= 50; ++seed) {
        ForerunnerArchetype arch = (ForerunnerArchetype)(seed % FORERUNNER_ARCHETYPE_COUNT);
        ForerunnerStage stage = (ForerunnerStage)(seed % 5);
        ForerunnerPlan plan = generate_forerunner_structure(seed, arch, stage);
        assert(verify_mass_void_polarity(&plan));
        assert(verify_forerunner_symmetry(&plan));
        if (arch != FORERUNNER_ARCHETYPE_CRUCIBLE) {
            assert(verify_skybridge_alignment(&plan));
        }
    }
    printf("-> 50/50 randomized seeds across 4 archetypes passed all invariants with 100%% stability!\n\n");

    printf("ALL FORERUNNER TRI-GRAMMAR TESTS PASSED SUCCESSFULLY.\n");
    return 0;
}
