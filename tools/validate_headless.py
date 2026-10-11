#!/usr/bin/env python3
"""Run rules, packet-format, and shader-layout checks without raylib or a display."""

import argparse
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]


def run(command: list[str]) -> None:
    print("+ " + " ".join(map(str, command)), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cc", default="gcc", help="C compiler executable (default: gcc)")
    parser.add_argument("--sanitize", action="store_true",
                        help="enable AddressSanitizer and UndefinedBehaviorSanitizer")
    args = parser.parse_args()
    flags = ["-std=c11", "-O1", "-g", "-Wall", "-Wextra", "-Werror", "-I."]
    if args.sanitize:
        flags += ["-fsanitize=address,undefined", "-fno-omit-frame-pointer",
                  "-fno-sanitize-recover=all"]
    build = ROOT / ".build"
    build.mkdir(exist_ok=True)
    try:
        with tempfile.TemporaryDirectory(prefix="headless-", dir=build) as directory:
            for name, sources in (
                ("firefight_rules", ["tests/firefight_rules_test.c"]),
                ("net_protocol", ["tests/net_protocol_test.c", "net_protocol.c"]),
            ):
                executable = str(pathlib.Path(directory) / (name + ".exe"))
                run([args.cc, *flags, *sources, "-o", executable])
                run([executable])
            run([sys.executable, "tools/validate_gpu_layout.py"])
    except (OSError, subprocess.CalledProcessError) as error:
        print(f"Headless validation failed: {error}", file=sys.stderr)
        return 1
    print("All headless checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
