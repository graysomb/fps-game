"""Encode the shared scene manifest emitted by render_victory_town."""
import json
from pathlib import Path
from PIL import Image, ImageDraw
root=Path(__file__).resolve().parents[1]/'artifacts/victory-town'
scenes=json.loads((root/'scenes.json').read_text())
def gif(paths,out,width=768,seconds=None):
    frames=[]
    for path in paths:
        with Image.open(path) as im:
            im=im.convert('RGB');im.thumbnail((width,width*9//16),Image.Resampling.LANCZOS)
            frames.append(im.quantize(colors=160,method=Image.Quantize.FASTOCTREE))
    # GIF delays use hundredths of a second; distribute rounding over frames.
    seconds=seconds if seconds is not None else len(frames)/20
    ticks=[round(i*seconds*100/len(frames)) for i in range(len(frames)+1)]
    delays=[(ticks[i+1]-ticks[i])*10 for i in range(len(frames))]
    frames[0].save(out,save_all=True,append_images=frames[1:],duration=delays,loop=0,optimize=False,disposal=2)
offset=0
for scene in scenes:
    count=scene['frames']
    gif([root/'scenes'/f'{j:04}.png' for j in range(offset,offset+count)],root/f"{scene['id']}.gif",seconds=scene.get("playback_seconds",count/20))
    offset+=count
gif([root/'frames'/f'{j:04}.png' for j in range(offset)],root/'reclaimed-foundry-eight-scenes.gif',seconds=sum(s.get('playback_seconds',s['frames']/20) for s in scenes))
canvas=Image.new('RGB',(1280,385*((len(scenes)+1)//2)),(14,23,31));draw=ImageDraw.Draw(canvas)
for i,scene in enumerate(scenes):
    with Image.open(root/f'scene-{i}.png') as im:
        x=(i%2)*640;y=(i//2)*385
        canvas.paste(im.resize((640,360),Image.Resampling.LANCZOS),(x,y))
        draw.text((x+12,y+366),scene['title'],fill=(244,216,166))
canvas.save(root/'eight-scene-overview.png')
print('Saved',len(scenes),'scene GIFs and soundtrack-length montage.')
# A compact overview focused on this addition; full overview remains available.
new=Image.new('RGB',(1280,770),(14,23,31));labels=ImageDraw.Draw(new)
for slot,index in enumerate((0,2,4,6)):
    with Image.open(root/f'scene-{index}.png') as im:
        x=(slot%2)*640;y=(slot//2)*385
        new.paste(im.resize((640,360),Image.Resampling.LANCZOS),(x,y))
        labels.text((x+12,y+366),scenes[index]['title'],fill=(244,216,166))
new.save(root/'new-scenes-overview.png')
