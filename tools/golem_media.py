"""Package deterministic Foundry presentation captures. Requires Pillow."""
from pathlib import Path
import argparse
from PIL import Image, ImageDraw

p = argparse.ArgumentParser()
p.add_argument('directory', type=Path, nargs='?', default=Path('artifacts/golem/rendered'))
root = p.parse_args().directory
names = ['idle', 'walk', 'slam', 'sweep', 'throw', 'arm-loss', 'leg-loss', 'chest-loss', 'both-legs', 'hit']
for name in names:
    files = sorted(root.glob(f'{name}-[0-9][0-9][0-9].png'))
    if len(files) != 120:
        raise RuntimeError(f'{name}: expected 120 frames, found {len(files)}')
    # A palette sampled across the clip preserves gold and regeneration cyan.
    atlas = Image.new('RGB', (450*4, 340))
    for j, index in enumerate([0, 30, 65, 90]):
        with Image.open(files[index]) as im:
            atlas.paste(im.resize((450, 340)), (450*j, 0))
    palette = atlas.quantize(colors=192)
    frames = []
    for path in files:
        with Image.open(path) as im:
            frames.append(im.convert('RGB').resize((720, 544)).quantize(palette=palette, dither=Image.Dither.NONE))
    frames[0].save(root / f'{name}.gif', save_all=True, append_images=frames[1:], duration=50, loop=0, optimize=False)
    for frame in frames:
        frame.close()

sheet = Image.new('RGB', (1500, 1130), '#101822')
draw = ImageDraw.Draw(sheet)
images = [('docs/concepts/foundry-golem.png', 'APPROVED CONCEPT'),
          (root / 'view-0.png', 'PROCEDURAL MODEL / RAT SCALE'),
          (root / 'arm-loss-065.png', 'ARM REGENERATION'),
          (root / 'leg-loss-040.png', 'LOWER LEG LOSS / BRACE'),
          (root / 'chest-loss-040.png', 'EXPOSED CORE'),
          (root / 'skeleton.png', '16 JOINT RIG')]
for index, (path, label) in enumerate(images):
    x, y = index % 2 * 750, index // 2 * 375
    draw.text((x+15, y+10), label, fill='white')
    with Image.open(path) as im:
        im.thumbnail((730, 340))
        sheet.paste(im, (x+(750-im.width)//2, y+30))
sheet.save(root / 'comparison.png')
print('Wrote ten GIFs and comparison.png')
