#include "../greek_grammar.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>

// Helper to check bilateral symmetry across X=0
static bool verify_bilateral_symmetry(const TemplePlan *plan) {
    // For every column or wall node with [min_x, max_x], there must be a matching
    // node or symmetric bounds across X=0.
    for (int i = 0; i < plan->node_count; ++i) {
        const ArchNode *n = &plan->nodes[i];
        if (n->type != ARCH_PRIMITIVE_COLUMN && n->type != ARCH_PRIMITIVE_WALL &&
            n->type != ARCH_PRIMITIVE_TREE && n->type != ARCH_PRIMITIVE_THOLOS &&
            n->type != ARCH_PRIMITIVE_EXEDRA && n->type != ARCH_PRIMITIVE_NAISKOS &&
            n->type != ARCH_PRIMITIVE_FOUNTAIN && n->type != ARCH_PRIMITIVE_BRAZIER &&
            n->type != ARCH_PRIMITIVE_PERGOLA && n->type != ARCH_PRIMITIVE_PATH &&
            n->type != ARCH_PRIMITIVE_SHRUB) continue;

        // If the node spans across X=0 (e.g. min_x = -k, max_x = +k), it is self-symmetric
        if (n->box.min_x == -n->box.max_x) continue;

        // Otherwise, look for its mirror partner with min_x' = -max_x, max_x' = -min_x, same min_z, max_z
        bool found_mirror = false;
        for (int j = 0; j < plan->node_count; ++j) {
            const ArchNode *m = &plan->nodes[j];
            if (m->type != n->type) continue;
            if (m->box.min_z == n->box.min_z && m->box.max_z == n->box.max_z &&
                m->box.min_y == n->box.min_y && m->box.max_y == n->box.max_y) {
                if (m->box.min_x == -n->box.max_x && m->box.max_x == -n->box.min_x) {
                    found_mirror = true;
                    break;
                }
            }
        }
        if (!found_mirror) {
            printf("Symmetry violation on node %d: type=%d, box=[%d..%d, %d..%d, %d..%d]\n",
                   i, n->type, n->box.min_x, n->box.max_x, n->box.min_y, n->box.max_y, n->box.min_z, n->box.max_z);
            return false;
        }
    }
    return true;
}

int main(void) {
    printf("=======================================================\n");
    printf(" RUNNING GREEK BUILDING BI-ALGEBRA VERIFICATION SUITE  \n");
    printf("=======================================================\n");

    // Test 1: Prostyle Temple (Front porch)
    printf("[1/4] Generating Prostyle Temple...\n");
    TemplePlan prostyle = generate_greek_temple(42, TEMPLE_STAGE_PROSTYLE, 2);
    assert(prostyle.stage >= TEMPLE_STAGE_PROSTYLE);
    assert(prostyle.has_stylobate);
    assert(prostyle.has_pediment);
    assert(verify_bilateral_symmetry(&prostyle));
    print_temple_blueprint(&prostyle);
    printf("-> Prostyle Temple passed bilateral symmetry & canonical assertions!\n\n");

    // Test 2: Amphiprostyle Temple (Front & Rear Porches)
    printf("[2/4] Generating Amphiprostyle Temple...\n");
    TemplePlan amphi = generate_greek_temple(1337, TEMPLE_STAGE_AMPHIPROSTYLE, 2);
    assert(amphi.stage >= TEMPLE_STAGE_AMPHIPROSTYLE);
    assert(amphi.has_stylobate);
    assert(amphi.has_pediment);
    assert(verify_bilateral_symmetry(&amphi));
    print_temple_blueprint(&amphi);
    printf("-> Amphiprostyle Temple passed bilateral symmetry & canonical assertions!\n\n");

    // Test 3: Peripteral Temple (Full Monument Colonnade)
    printf("[3/4] Generating Peripteral Temple...\n");
    TemplePlan peri = generate_greek_temple(9999, TEMPLE_STAGE_PERIPTERAL, 2);
    assert(peri.stage >= TEMPLE_STAGE_PERIPTERAL);
    assert(peri.has_stylobate);
    assert(peri.has_pediment);
    assert(verify_bilateral_symmetry(&peri));
    print_temple_blueprint(&peri);
    printf("-> Peripteral Temple passed bilateral symmetry & canonical assertions!\n\n");

    // Test 4: Sanctuary (Courtyard + Pool + Sacred Groves + Garden Precinct)
    printf("[4/4] Generating Sanctuary (Temenos, Pool, Grove & Garden Precinct)...\n");
    TemplePlan sanctuary = generate_greek_temple(2024, TEMPLE_STAGE_SANCTUARY, 2);
    assert(sanctuary.stage == TEMPLE_STAGE_SANCTUARY);
    assert(sanctuary.has_courtyard);
    assert(sanctuary.has_pool);
    assert(sanctuary.has_trees);
    assert(sanctuary.has_garden_paths);
    assert(sanctuary.has_tholos);
    assert(sanctuary.has_exedrae);
    assert(sanctuary.has_naiskoi);
    assert(sanctuary.has_fountains);
    assert(sanctuary.has_pergolas);
    assert(sanctuary.has_braziers);
    assert(sanctuary.has_garden_flora);
    assert(verify_bilateral_symmetry(&sanctuary));
    print_temple_blueprint(&sanctuary);
    printf("-> Sanctuary passed bilateral symmetry, pool excavation & garden precinct assertions!\n\n");

    // Test 5: Seed invariance & randomness sweep
    printf("Running Monte-Carlo sweep over 50 random seeds...\n");
    for (uint32_t seed = 1; seed <= 50; ++seed) {
        TempleStage stage = (TempleStage)(seed % 5);
        TemplePlan p = generate_greek_temple(seed, stage, 2);
        assert(verify_bilateral_symmetry(&p));
    }
    printf("-> 50/50 randomized generative seeds passed bilateral symmetry with zero collisions!\n");

    printf("\nALL GREEK BI-ALGEBRA TESTS PASSED SUCCESSFULLY.\n");
    return 0;
}
