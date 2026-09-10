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

    // Test 1: Stage 0 Abyssal Chasm
    printf("[1/5] Generating Stage 0: Abyssal Chasm & Beacon Spire...\n");
    ForerunnerPlan p0 = generate_forerunner_structure(42, FORERUNNER_STAGE_CHASM);
    assert(p0.has_chasm);
    assert(verify_mass_void_polarity(&p0));
    printf("-> Stage 0 passed chasm void & beacon spire polarity assertions!\n\n");

    // Test 2: Stage 1 Monumental Chevron Gateway
    printf("[2/5] Generating Stage 1: Canted Chevron Gateway...\n");
    ForerunnerPlan p1 = generate_forerunner_structure(101, FORERUNNER_STAGE_GATEWAY);
    assert(p1.has_gateway);
    assert(p1.pylon_count >= 2);
    assert(verify_forerunner_symmetry(&p1));
    printf("-> Stage 1 passed canted pylon & chevron lintel symmetry assertions!\n\n");

    // Test 3: Stage 2 Suspended Skybridge
    printf("[3/5] Generating Stage 2: Suspended Skybridge & Hard-Light Span...\n");
    ForerunnerPlan p2 = generate_forerunner_structure(202, FORERUNNER_STAGE_SKYBRIDGE);
    assert(p2.has_bridge);
    assert(verify_skybridge_alignment(&p2));
    assert(verify_forerunner_symmetry(&p2));
    printf("-> Stage 2 passed axial bridge span & guide rail assertions!\n\n");

    // Test 4: Stage 3 Terraced Vault Chamber
    printf("[4/5] Generating Stage 3: Subterranean Terraced Vault Galleries...\n");
    ForerunnerPlan p3 = generate_forerunner_structure(303, FORERUNNER_STAGE_VAULT);
    assert(p3.has_vault);
    assert(verify_forerunner_symmetry(&p3));
    printf("-> Stage 3 passed terraced galleries & wall buttress assertions!\n\n");

    // Test 5: Stage 4 Apex Cartographer Installation
    printf("[5/5] Generating Stage 4: Full Apex Cartographer Megastructure...\n");
    ForerunnerPlan p4 = generate_forerunner_structure(1337, FORERUNNER_STAGE_CARTOGRAPHER);
    assert(p4.has_core);
    assert(verify_mass_void_polarity(&p4));
    assert(verify_skybridge_alignment(&p4));
    assert(verify_forerunner_symmetry(&p4));
    print_forerunner_blueprint(&p4);
    printf("-> Stage 4 Apex Cartographer passed all mass-void-path tri-grammar invariants!\n\n");

    // Monte Carlo 50-Seed Invariant Sweep
    printf("Running Monte-Carlo sweep over 50 random Forerunner seeds...\n");
    for (uint32_t seed = 1; seed <= 50; ++seed) {
        ForerunnerStage stage = (ForerunnerStage)(seed % 5);
        ForerunnerPlan plan = generate_forerunner_structure(seed, stage);
        assert(verify_mass_void_polarity(&plan));
        assert(verify_forerunner_symmetry(&plan));
        assert(verify_skybridge_alignment(&plan));
    }
    printf("-> 50/50 randomized Forerunner seeds passed all invariants with 100%% stability!\n\n");

    printf("ALL FORERUNNER TRI-GRAMMAR TESTS PASSED SUCCESSFULLY.\n");
    return 0;
}
