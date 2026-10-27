"""Encode captures from validate_golem_combat.ps1 -Preview (Pillow required)."""
from pathlib import Path
from PIL import Image

root = Path(__file__).resolve().parents[1] / 'artifacts/golem/combat'
for name in ['part-shields', 'shoulder-guns', 'head-rocket', 'structure-throw']:
    paths = sorted(root.glob(f'{name}-*.png'))
    if not paths:
        raise SystemExit(f'No frames for {name}; run the preview harness first.')
    frames = []
    for path in paths:
        with Image.open(path) as image:
            frames.append(image.convert('RGB').resize((720, 540)).quantize(colors=160))
    frames[0].save(root / f'{name}.gif', save_all=True, append_images=frames[1:],
                   duration=100, loop=0, disposal=2)
    print(root / f'{name}.gif')
