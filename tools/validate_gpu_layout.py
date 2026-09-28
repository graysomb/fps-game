#!/usr/bin/env python3
"""Fail when the C, GLSL, and Metal pipeline modes/buffer slots drift apart."""

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]


def constants(path: pathlib.Path, pattern: str) -> dict[str, int]:
    text = path.read_text(encoding="utf-8")
    return {name: int(value) for name, value in re.findall(pattern, text)}


layout = constants(ROOT / "physics_gpu_layout.h", r"GPU_MODE_([A-Z_]+)\s*=\s*(\d+)")
glsl = constants(ROOT / "shaders/pbd/pbd_pipeline.comp", r"const int MODE_([A-Z_]+)\s*=\s*(\d+)")
metal = constants(ROOT / "shaders/pbd/pbd_pipeline.metal", r"constant int MODE_([A-Z_]+)\s*=\s*(\d+)")
control = constants(ROOT / "physics_gpu_layout.h", r"FPS_GPU_CONTROL_([A-Z_]+)\s*=\s*(\d+)")

errors: list[str] = []
for shader_name, shader in (("GLSL", glsl), ("Metal", metal)):
    if shader != layout:
        missing = sorted(set(layout) - set(shader))
        extra = sorted(set(shader) - set(layout))
        wrong = sorted(name for name in set(layout) & set(shader) if layout[name] != shader[name])
        errors.append(f"{shader_name} modes differ: missing={missing} extra={extra} wrong={wrong}")

for shader_name, path, pattern in (
    ("GLSL", ROOT / "shaders/pbd/pbd_pipeline.comp", r"const int CONTROL_FLAG_([A-Z_]+)\s*=\s*(\d+)"),
    ("Metal", ROOT / "shaders/pbd/pbd_pipeline.metal", r"constant int CONTROL_FLAG_([A-Z_]+)\s*=\s*(\d+)"),
):
    shader_flags = constants(path, pattern)
    if shader_flags != control:
        errors.append(f"{shader_name} control flags differ: C={control} shader={shader_flags}")

glsl_text = (ROOT / "shaders/pbd/pbd_pipeline.comp").read_text(encoding="utf-8")
metal_text = (ROOT / "shaders/pbd/pbd_pipeline.metal").read_text(encoding="utf-8")
glsl_slots = {int(value) for value in re.findall(r"binding\s*=\s*(\d+)", glsl_text)}
metal_slots = {int(value) for value in re.findall(r"\[\[buffer\((\d+)\)\]\]", metal_text)}
expected_slots = set(range(22))
if glsl_slots != expected_slots:
    errors.append(f"GLSL buffer slots differ: {sorted(glsl_slots)}")
if not expected_slots.issubset(metal_slots):
    errors.append(f"Metal buffer slots differ: {sorted(metal_slots)}")

amr_c = constants(ROOT / "physics_gpu_common.inc", r"GPU_AMR_META_([A-Z_]+)\s*=\s*(\d+)")
amr_glsl = constants(ROOT / "shaders/pbd/pbd_pipeline.comp", r"const uint AMR_([A-Z_]+)\s*=\s*(\d+)u")
amr_metal = constants(ROOT / "shaders/pbd/pbd_pipeline.metal", r"constant uint AMR_([A-Z_]+)\s*=\s*(\d+)u")
for name, c_name in (("OCCUPANCY_SHIFT", "OCCUPANCY_SHIFT"),
                     ("REAL_MASK_SHIFT", "REAL_MASK_SHIFT")):
    expected = amr_c.get(c_name)
    if expected is None or amr_glsl.get(name) != expected or amr_metal.get(name) != expected:
        errors.append(f"AMR metadata {name} differs: C={expected} GLSL={amr_glsl.get(name)} Metal={amr_metal.get(name)}")

for mode in sorted(set(layout) - {"SCENE_COLLISIONS"}):
    if glsl_text.count(f"MODE_{mode}") < 2:
        errors.append(f"GLSL does not dispatch MODE_{mode}")
    if metal_text.count(f"MODE_{mode}") < 2:
        errors.append(f"Metal does not dispatch MODE_{mode}")

common_text = (ROOT / "physics_gpu_common.inc").read_text(encoding="utf-8")
if "gpu_dispatch_indirect_offset(GPU_MODE_VGS, 52)" not in common_text:
    errors.append("resident VGS indirect offset changed from byte 52")
if "dispatchArgs[13]" not in glsl_text or "dispatchArgs[13]" not in metal_text:
    errors.append("shader resident VGS indirect slot changed from dispatchArgs[13]")

if errors:
    print("\n".join(errors), file=sys.stderr)
    raise SystemExit(1)
print("GPU layout validation passed")
