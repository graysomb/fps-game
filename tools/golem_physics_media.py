"""Turn deterministic physical-golem captures into 20 fps GIFs (requires Pillow)."""
from pathlib import Path
from PIL import Image

root = Path(__file__).resolve().parents[1] / "artifacts/golem/physics"
names = {"idle": "idle", "motion-1": "walking", "motion-2": "slam",
         "motion-3": "sweep", "motion-4": "throw", "leg-loss": "kneel",
         "regeneration": "regeneration", "leg-regrowth": "leg-regrowth"}
for folder, name in names.items():
    paths = sorted((root / folder).glob("*.png"))
    if folder == "leg-regrowth" and not paths:
        continue  # Optional -Repair -Capture run.
    if len(paths) != 120:
        raise RuntimeError(f"{folder}: expected 120 frames, got {len(paths)}; run -Capture")
    frames = []
    for path in paths:
        with Image.open(path) as source:
            frames.append(source.convert("RGB").quantize(colors=128))
    output = root / f"physical-{name}.gif"
    frames[0].save(output, save_all=True, append_images=frames[1:], duration=50,
                   loop=0, disposal=2, optimize=False)
    print(output)
