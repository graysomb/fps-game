#include "../megalith_grammar.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>

// Helper to verify structural support invariant: every capstone must have supporting uprights
static bool verify_capstone_supports(const MegalithPlan *plan) {
    for (int i = 0; i < plan->node_count; ++i) {
        const MegalithNode *n = &plan->nodes[i];
        if (n->type == MEGALITH_PRIMITIVE_CAPSTONE || n->type == MEGALITH_PRIMITIVE_CORBEL) {
            if (n->support_count < 1) {
                printf("Error: Unsupported capstone node %d!\n", i);
                return false;
            }
            // Verify that support parents exist and touch/underlie the capstone
            for (int s = 0; s < n->support_count; ++s) {
                int parent_id = n->support_parents[s];
                if (parent_id < 0 || parent_id >= plan->node_count) {
                    printf("Error: Invalid parent ID %d for capstone %d\n", parent_id, i);
                    return false;
                }
                const MegalithNode *parent = &plan->nodes[parent_id];
                if (parent->box.max_y < n->box.min_y - 1) {
                    printf("Error: Parent %d does not reach capstone %d (parent max_y=%d, capstone min_y=%d)\n",
                           parent_id, i, parent->box.max_y, n->box.min_y);
                    return false;
                }
            }
        }
    }
    return true;
}

// Helper to verify solstice axis alignment of the passage
static bool verify_passage_alignment(const MegalithPlan *plan) {
    if (!plan->has_passage) return true;
    for (int i = 0; i < plan->node_count; ++i) {
        const MegalithNode *n = &plan->nodes[i];
        if (n->box.min_z > plan->chamber_cz + 3 && n->type == MEGALITH_PRIMITIVE_WALL_SLAB) {
            // Passage walls should be centered near X=0 along the sunrise axis
            if (abs(n->box.min_x) > 5 || abs(n->box.max_x) > 5) {
                printf("Error: Passage slab %d deviated from solstice axis! [%d..%d]\n",
                       i, n->box.min_x, n->box.max_x);
                return false;
            }
        }
    }
    return true;
}

// Helper to verify stone circle radial geometry
static bool verify_circle_radius(const MegalithPlan *plan) {
    if (plan->circle_stones_placed < 4) return true;
    int r = plan->circle_radius;
    int tolerance = 3;
    for (int i = 0; i < plan->node_count; ++i) {
        const MegalithNode *n = &plan->nodes[i];
        if (n->type == MEGALITH_PRIMITIVE_ORTHOSTAT && n->box.min_y == 0 && n->box.min_z < plan->circle_radius + 1) {
            int mx = (n->box.min_x + n->box.max_x) / 2;
            int mz = (n->box.min_z + n->box.max_z) / 2;
            float dist = sqrtf((float)(mx * mx + mz * mz));
            // Either inner trilithon (dist < 8) or perimeter circle (dist ~ r)
            if (dist > 7.0f && fabsf(dist - (float)r) > (float)tolerance) {
                printf("Error: Stone %d at (%d, %d) violates circle radius %d (dist=%.1f)\n",
                       i, mx, mz, r, dist);
                return false;
            }
        }
    }
    return true;
}

int main(void) {
    printf("=======================================================\n");
    printf(" RUNNING PREHISTORIC MEGALITH BI-ALGEBRA TEST SUITE   \n");
    printf("=======================================================\n\n");

    // Test 1: Dolmen (Portal Tomb)
    printf("[1/4] Generating Dolmen (Portal Tomb)...\n");
    MegalithPlan dolmen = generate_megalith_structure(42, MEGALITH_ARCHETYPE_PASSAGE_GRAVE, MEGALITH_STAGE_DOLMEN);
    assert(dolmen.has_capstone);
    assert(verify_capstone_supports(&dolmen));
    print_megalith_blueprint(&dolmen);
    printf("-> Dolmen passed gravity support & capstone bridging assertions!\n\n");

    // Test 2: Passage Grave & Earthen Tumulus Barrow
    printf("[2/4] Generating Passage Grave & Tumulus Barrow...\n");
    MegalithPlan barrow = generate_megalith_structure(1337, MEGALITH_ARCHETYPE_PASSAGE_GRAVE, MEGALITH_STAGE_TUMULUS);
    assert(barrow.has_capstone);
    assert(barrow.has_passage);
    assert(barrow.has_mound);
    assert(barrow.has_hearth);
    assert(verify_capstone_supports(&barrow));
    assert(verify_passage_alignment(&barrow));
    print_megalith_blueprint(&barrow);
    printf("-> Tumulus Barrow passed corbel vault, solstice alignment & mound enclosure assertions!\n\n");

    // Test 3: Megalithic Stone Circle / Cromlech
    printf("[3/4] Generating Megalithic Stone Circle (Cromlech)...\n");
    MegalithPlan circle = generate_megalith_structure(777, MEGALITH_ARCHETYPE_STONE_CIRCLE, MEGALITH_STAGE_DOLMEN);
    assert(circle.circle_stones_placed >= 8);
    assert(verify_circle_radius(&circle));
    print_megalith_blueprint(&circle);
    printf("-> Stone Circle passed radial boundary expansion assertions!\n\n");

    // Test 4: Monumental Henge with Trilithons & Avenue
    printf("[4/4] Generating Monumental Henge with Trilithons & Avenue...\n");
    MegalithPlan henge = generate_megalith_structure(999, MEGALITH_ARCHETYPE_STONE_CIRCLE, MEGALITH_STAGE_HENGE);
    assert(henge.has_corbel); // Inner trilithons
    assert(henge.has_avenue);
    assert(henge.has_hearth);
    assert(verify_capstone_supports(&henge));
    assert(verify_circle_radius(&henge));
    print_megalith_blueprint(&henge);
    printf("-> Monumental Henge passed inner trilithon bridging & processional avenue assertions!\n\n");

    // Test 5: Monte-Carlo 50-Seed Stability & Invariant Sweep
    printf("Running Monte-Carlo sweep over 50 random prehistoric seeds...\n");
    for (uint32_t seed = 1; seed <= 50; ++seed) {
        MegalithArchetype arch = (seed % 2 == 0) ? MEGALITH_ARCHETYPE_PASSAGE_GRAVE : MEGALITH_ARCHETYPE_STONE_CIRCLE;
        MegalithStage stage = (MegalithStage)(seed % 6);
        MegalithPlan p = generate_megalith_structure(seed, arch, stage);
        assert(verify_capstone_supports(&p));
        assert(verify_passage_alignment(&p));
    }
    printf("-> 50/50 randomized megalithic seeds passed structural gravity support & alignment checks!\n\n");

    printf("ALL MEGALITH BI-ALGEBRA TESTS PASSED SUCCESSFULLY.\n");
    return 0;
}
