"""Package rendered mallard frames into a loop; requires Pillow."""
from pathlib import Path
import argparse
from PIL import Image

parser = argparse.ArgumentParser()
parser.add_argument('directory', nargs='?', type=Path, default=Path('artifacts/duck'))
root = parser.parse_args().directory
files = sorted(root.glob('swim-[0-9][0-9][0-9].png'))
if len(files) != 160:
    raise RuntimeError(f'Expected 160 rendered frames, found {len(files)}')
with Image.open(files[20]) as source:
    palette = source.convert('RGB').quantize(colors=160)
frames = []
for path in files:
    with Image.open(path) as source:
        frames.append(source.convert('RGB').quantize(palette=palette, dither=Image.Dither.NONE))
frames[0].save(root / 'mallard-swimming.gif', save_all=True, append_images=frames[1:],
               duration=50, loop=0, optimize=False)
for frame in frames:
    frame.close()
palette.close()
print(root / 'mallard-swimming.gif')
