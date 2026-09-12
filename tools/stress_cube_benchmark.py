#!/usr/bin/env python3
"""
Scalable Active Cube Stress Test: GPU (Metal) vs Multithreaded CPU.
Runs cube sizes from 4x4x4 (64 voxels) up to 32x32x32 (32,768 voxels)
and captures GPU exec time, CPU wait time, and CPU step time.
"""

import json
import os
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
GPU_BIN = ROOT / ".build" / "bin" / "fps_ray_gpu"
CPU_BIN = ROOT / ".build" / "bin" / "fps_ray_cpu"

SIZES = [4, 8, 12, 16, 20, 24, 28, 32]
STEPS = 60

def run_test(binary, size, steps=STEPS, extra_env=None):
    env = os.environ.copy()
    env["FPS_CUBE_SIZE"] = str(size)
    if extra_env:
        env.update(extra_env)

    cmd = [
        str(binary),
        "--debug-scenario=cube-stress",
        f"--debug-steps={steps}",
        "--debug-capture-steps=60",
        "--gpu-transfer-stats"
    ]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True, env=env)
    wall_ms = (time.perf_counter() - t0) * 1000.0

    output = proc.stdout + "\n" + proc.stderr
    stats = {
        "passed": proc.returncode == 0,
        "wall_ms": wall_ms,
        "gpuExecMs": 0.0,
        "waitMs": 0.0,
        "packMs": 0.0,
        "uploadMs": 0.0,
        "commitMs": 0.0,
        "step_ms": 0.0,
        "voxels": size ** 3,
        "particles": 0,
    }

    for line in output.splitlines():
        if "gpu-transfer" in line:
            for token in line.split():
                if "=" in token:
                    k, v = token.split("=", 1)
                    if k in stats:
                        try:
                            stats[k] = float(v)
                        except ValueError:
                            pass
        if "report=" in line:
            for token in line.split():
                if token.startswith("report="):
                    report_path = token.split("=", 1)[1]
                    if os.path.exists(report_path):
                        try:
                            with open(report_path, "r") as f:
                                rep = json.load(f)
                                samples = rep.get("samples", [])
                                if samples:
                                    last = samples[-1]
                                    stats["step_ms"] = last.get("physicsStepMs", 0.0)
                                    stats["particles"] = last.get("simParticles", 0)
                        except Exception:
                            pass

    return stats

def main():
    print(f"Running Scalable Active Cube Stress Test ({STEPS} steps per size)...")
    print(f"{'Size':<8} | {'Voxels':<8} | {'Particles':<10} | {'CPU Step':<10} | {'GPU Exec':<10} | {'GPU Step':<10} | {'GPU Speedup':<12} | {'GPU Wall':<10}")
    print("-" * 88)

    results = []
    for s in SIZES:
        gpu = run_test(GPU_BIN, s)
        cpu = run_test(CPU_BIN, s)

        gpu_exec_step = gpu["gpuExecMs"] / STEPS if STEPS else 0.0
        gpu_step = gpu["step_ms"]
        cpu_step = cpu["step_ms"]
        speedup = (cpu_step / gpu_exec_step) if gpu_exec_step > 0 else 0.0

        print(f"{s}x{s}x{s:<4} | {s**3:<8} | {gpu['particles']:<10} | {cpu_step:6.2f} ms  | {gpu_exec_step:6.2f} ms  | {gpu_step:6.2f} ms  | {speedup:6.2f}x      | {gpu['wall_ms']:6.1f} ms")
        results.append({
            "size": s,
            "voxels": s**3,
            "particles": gpu["particles"],
            "cpu_step_ms": cpu_step,
            "gpu_exec_step_ms": gpu_exec_step,
            "gpu_step_ms": gpu_step,
            "speedup": speedup,
            "gpu_wall_ms": gpu["wall_ms"],
            "gpu_wait_ms": gpu["waitMs"] / STEPS
        })

    with open(ROOT / "tools" / "stress_cube_results.json", "w") as f:
        json.dump(results, f, indent=2)

if __name__ == "__main__":
    main()
