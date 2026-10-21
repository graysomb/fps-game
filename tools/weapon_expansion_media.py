"""Package deterministic renderer PNGs. Requires Pillow; never regenerates assets."""
from pathlib import Path
import argparse
from PIL import Image, ImageDraw

p = argparse.ArgumentParser()
p.add_argument('directory', type=Path, nargs='?', default=Path('artifacts/weapon-expansion/rendered'))
args = p.parse_args()
root = args.directory
names = ['glob', 'black-hole', 'glock', 'igloo', 'dmr']
concepts = ['glob-gun', 'black-hole-bomb', 'glock-pistol', 'igloo-grenade-chimney', 'dmr']
for name in names:
    files = sorted(root.glob(f'{name}-[0-9][0-9][0-9].png'))
    if len(files) != 160:
        raise RuntimeError(f'{name}: expected 160 deterministic frames, found {len(files)}')
    # One shared palette avoids per-frame palette flicker.
    palette = Image.open(files[40]).convert('RGB').quantize(colors=160)
    frames = []
    for path in files:
        with Image.open(path) as im:
            frames.append(im.convert('RGB').quantize(palette=palette, dither=Image.Dither.NONE))
    frames[0].save(root / f'{name}.gif', save_all=True, append_images=frames[1:], duration=50, loop=0, optimize=False)
    for im in frames:
        im.close()
sheet = Image.new('RGB', (1200, 5*350), '#101b2b')
d = ImageDraw.Draw(sheet)
for i, (name, concept) in enumerate(zip(names, concepts)):
    with Image.open(Path('artifacts/weapon-concepts/expansion-01') / f'{concept}.png') as im:
        im.thumbnail((340, 320)); sheet.paste(im, (0, i*350+24))
    for j, filename in enumerate([f'{name}-view-0.png', f'{name}-first-person-full.png']):
        with Image.open(root / filename) as im:
            im.thumbnail((420, 320)); sheet.paste(im, (350+j*425, i*350+24))
    d.text((10, i*350+5), f'{name.upper()}   |   CONCEPT / PROCEDURAL MODEL / FIRST PERSON', fill='white')
sheet.save(root / 'concept-comparison.png')
print('Wrote five 8-second GIFs and concept-comparison.png')
