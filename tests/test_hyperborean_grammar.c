#include "../hyperborean_grammar.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// Helper to check load-bearing support polygon invariants
static bool verify_hyper_supports(const HyperPlan *plan) {
    for (int i = 0; i < plan->node_count; ++i) {
        const HyperNode *n = &plan->nodes[i];
        if (n->type == HYPER_PRIM_ROUGH_LINTEL || n->type == HYPER_PRIM_ARCHITRAVE ||
            n->type == HYPER_PRIM_CYCLOPEAN_ARCHITRAVE || n->type == HYPER_PRIM_PEDIMENT) {
            if (n->parent_count < 2) {
                printf("Error: Elevated lintel/pediment %d has < 2 supporting parents!\n", i);
                return false;
            }
            int min_sup_x = 9999, max_sup_x = -9999;
            int min_sup_z = 9999, max_sup_z = -9999;
            for (int p = 0; p < n->parent_count; ++p) {
                int pid = n->parent_ids[p];
                if (pid < 0 || pid >= plan->node_count) {
                    printf("Error: Invalid parent id %d on node %d\n", pid, i);
                    return false;
                }
                const HyperNode *parent = &plan->nodes[pid];
                if (parent->box.min_x < min_sup_x) min_sup_x = parent->box.min_x;
                if (parent->box.max_x > max_sup_x) max_sup_x = parent->box.max_x;
                if (parent->box.min_z < min_sup_z) min_sup_z = parent->box.min_z;
                if (parent->box.max_z > max_sup_z) max_sup_z = parent->box.max_z;
            }

            int cm_x = (n->box.min_x + n->box.max_x) / 2;
            int cm_z = (n->box.min_z + n->box.max_z) / 2;
            if (cm_x < min_sup_x - 2 || cm_x > max_sup_x + 2 ||
                cm_z < min_sup_z - 2 || cm_z > max_sup_z + 2) {
                printf("Error: Center of mass (%d, %d) outside support polygon [%d..%d, %d..%d]\n",
                       cm_x, cm_z, min_sup_x, max_sup_x, min_sup_z, max_sup_z);
                return false;
            }
        }
    }
    return true;
}

// Helper to verify solstice axis alignment
static bool verify_solstice_alignment(const HyperPlan *plan) {
    if (!plan->has_avenue) return true;
    for (int i = 0; i < plan->node_count; ++i) {
        const HyperNode *n = &plan->nodes[i];
        if (n->type == HYPER_PRIM_HEEL_STONE) {
            int cm_x = (n->box.min_x + n->box.max_x) / 2;
            if (abs(cm_x) > 1) {
                printf("Error: Heel stone deviated from solstice axis! X=%d\n", cm_x);
                return false;
            }
        }
    }
    return true;
}

// Helper to verify concentric ring geometry
static bool verify_concentric_rings(const HyperPlan *plan) {
    for (int i = 0; i < plan->node_count; ++i) {
        const HyperNode *n = &plan->nodes[i];
        int cm_x = (n->box.min_x + n->box.max_x) / 2;
        int cm_z = (n->box.min_z + n->box.max_z) / 2;
        float dist = sqrtf((float)(cm_x * cm_x + cm_z * cm_z));

        // Skip avenue nodes at +Z >= 16
        if (cm_z >= 16) continue;

        if (n->type == HYPER_PRIM_MENHIR && dist > 10.0f) {
            // Outer henge ring
            if (fabsf(dist - (float)plan->r_outer_henge) > 3.0f) {
                printf("Error: Outer menhir %d at dist=%.1f violates r_outer_henge=%d\n",
                       i, dist, plan->r_outer_henge);
                return false;
            }
        } else if (n->type == HYPER_PRIM_FLUTED_COLUMN) {
            // Mid peristyle ring
            if (fabsf(dist - (float)plan->r_mid_peristyle) > 3.0f) {
                printf("Error: Marble column %d at dist=%.1f violates r_mid_peristyle=%d\n",
                       i, dist, plan->r_mid_peristyle);
                return false;
            }
        } else if (n->type == HYPER_PRIM_THOLOS_PODIUM || n->type == HYPER_PRIM_FLUID_BASIN) {
            // Core epicenter
            if (dist > 4.5f) {
                printf("Error: Tholos node %d at dist=%.1f outside core radius!\n", i, dist);
                return false;
            }
        }
    }
    return true;
}

int main(void) {
    printf("=======================================================\n");
    printf(" RUNNING HYPERBOREAN SUN-HENGE BI-ALGEBRA TEST SUITE   \n");
    printf("=======================================================\n\n");

    // Test 1: Stage 0 Solstice Avenue
    printf("[1/5] Generating Stage 0: Solstice Avenue & Heel Stone...\n");
    HyperPlan p0 = generate_hyperborean_structure(42, HYPER_STAGE_AVENUE);
    assert(p0.has_avenue);
    assert(verify_solstice_alignment(&p0));
    printf("-> Stage 0 passed solstice axis & Heel Stone assertions!\n\n");

    // Test 2: Stage 1 Outer Cromlech Henge
    printf("[2/5] Generating Stage 1: Outer Sarsen Cromlech Henge...\n");
    HyperPlan p1 = generate_hyperborean_structure(101, HYPER_STAGE_OUTER_HENGE);
    assert(p1.outer_stones_placed >= 14);
    assert(verify_hyper_supports(&p1));
    assert(verify_concentric_rings(&p1));
    printf("-> Stage 1 passed outer cromlech lintel & radial bounds assertions!\n\n");

    // Test 3: Stage 2 Classical Marble Peristyle
    printf("[3/5] Generating Stage 2: Concentric Classical Marble Peristyle...\n");
    HyperPlan p2 = generate_hyperborean_structure(202, HYPER_STAGE_CLASSICAL_PERISTYLE);
    assert(p2.peristyle_cols_placed >= 10);
    assert(verify_hyper_supports(&p2));
    assert(verify_concentric_rings(&p2));
    printf("-> Stage 2 passed marble fluted peristyle & architrave support assertions!\n\n");

    // Test 4: Stage 3 Great Hybrid Trilithons with Doric Pediment
    printf("[4/5] Generating Stage 3: Great Hybrid Trilithons with Doric Pediment...\n");
    HyperPlan p3 = generate_hyperborean_structure(303, HYPER_STAGE_GREAT_TRILITHONS);
    assert(p3.trilithons_placed >= 3);
    assert(verify_hyper_supports(&p3));
    printf("-> Stage 3 passed monumental trilithon & classical pediment load-bearing assertions!\n\n");

    // Test 5: Stage 4 Full Hyperborean Sanctum
    printf("[5/5] Generating Stage 4: Full Hyperborean Sanctum Complex...\n");
    HyperPlan p4 = generate_hyperborean_structure(1337, HYPER_STAGE_FULL_SANCTUM);
    assert(p4.has_tholos);
    assert(p4.has_moat);
    assert(p4.has_brazier);
    assert(verify_hyper_supports(&p4));
    assert(verify_solstice_alignment(&p4));
    assert(verify_concentric_rings(&p4));
    print_hyperborean_blueprint(&p4);
    printf("-> Stage 4 Full Sanctum passed all syncretic bi-algebra invariants!\n\n");

    // Monte-Carlo 50-Seed Invariant Sweep
    printf("Running Monte-Carlo sweep over 50 random Hyperborean seeds...\n");
    for (uint32_t seed = 1; seed <= 50; ++seed) {
        HyperStage stage = (HyperStage)(seed % 5);
        HyperPlan plan = generate_hyperborean_structure(seed, stage);
        assert(verify_hyper_supports(&plan));
        assert(verify_solstice_alignment(&plan));
        assert(verify_concentric_rings(&plan));
    }
    printf("-> 50/50 randomized Hyperborean seeds passed all invariants with 100%% stability!\n\n");

    printf("ALL HYPERBOREAN SUN-HENGE TESTS PASSED SUCCESSFULLY.\n");
    return 0;
}
