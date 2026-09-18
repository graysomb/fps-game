#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
mkdir -p .build/buildings
cc -std=c11 -O2 -Wall -Wextra -Werror -I. tests/test_building_generation.c building_generation.c -o .build/buildings/test_generation
cc -std=c11 -O2 -Wall -Wextra -Werror -I. -I.build/raylib-macos/src \
  tests/test_building_diversity.c greek_grammar.c megalith_grammar.c hyperborean_grammar.c \
  forerunner_grammar.c unified_sanctum_grammar.c building_generation.c -lm \
  -o .build/buildings/test_diversity
.build/buildings/test_generation
.build/buildings/test_diversity
for name in greek megalith hyperborean forerunner; do
  cc -std=c11 -O2 -Wall -Wextra -Werror -I. -I.build/raylib-macos/src \
    "tests/test_${name}_grammar.c" "${name}_grammar.c" building_generation.c -lm \
    -o ".build/buildings/test_${name}"
  ".build/buildings/test_${name}" >/dev/null
done
cc -std=c11 -O2 -Wall -Wextra -Werror -I. -I.build/raylib-macos/src \
  tests/test_unified_grammar.c unified_sanctum_grammar.c building_generation.c -lm \
  -o .build/buildings/test_unified
.build/buildings/test_unified >/dev/null
