"""Encode captures from validate_golem_combat.ps1 -Preview (Pillow required)."""
from pathlib import Path
from PIL import Image

root = Path(__file__).resolve().parents[1] / 'artifacts/golem/combat'
for pattern, name, duration in [
    ('preview-0-*.png', 'slam.gif', 100),
    ('preview-1-*.png', 'sweep.gif', 100),
    ('preview-2-*.png', 'hand-launch.gif', 100),
    ('shield-[0-9]*.png', 'shield-regeneration.gif', 200),
]:
    paths = sorted(root.glob(pattern))
    if not paths:
        raise SystemExit(f'No frames for {pattern}; run the preview harness first.')
    frames = []
    for path in paths:
        with Image.open(path) as image:
            frames.append(image.convert('RGB').resize((720, 540)).quantize(colors=160))
    frames[0].save(root / name, save_all=True, append_images=frames[1:],
                   duration=duration, loop=0, disposal=2)
    print(root / name)
images = ['shield-full.png', 'shield-broken.png', 'shield-regenerating.png', 'defeated.png']
sheet = Image.new('RGB', (960, 720))
for i, name in enumerate(images):
    with Image.open(root / name) as image:
        sheet.paste(image.convert('RGB').resize((480, 360)), ((i % 2) * 480, (i // 2) * 360))
sheet.save(root / 'reactor-comparison.jpg', quality=92)
