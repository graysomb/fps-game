"""Encode deterministic 20 fps town captures. Requires Pillow; no game assets generated."""
from pathlib import Path
from PIL import Image, ImageDraw
root=Path(__file__).resolve().parents[1]/'artifacts/victory-town'
names=['supper','homes','garden-workshop','town-under-stars']
def gif(paths,out,width=768):
    frames=[]
    for path in paths:
        with Image.open(path) as im:
            im=im.convert('RGB');im.thumbnail((width,width*9//16),Image.Resampling.LANCZOS)
            frames.append(im.quantize(colors=160,method=Image.Quantize.FASTOCTREE))
    frames[0].save(out,save_all=True,append_images=frames[1:],duration=50,loop=0,optimize=False,disposal=2)
for i,name in enumerate(names):
    gif([root/'scenes'/f'{j:04}.png' for j in range(i*160,(i+1)*160)],root/f'{name}.gif')
gif(sorted((root/'frames').glob('*.png')),root/'reclaimed-foundry-montage.gif')
canvas=Image.new('RGB',(1280,770),(14,23,31));draw=ImageDraw.Draw(canvas)
for i,name in enumerate(names):
    with Image.open(root/f'scene-{i}.png') as im:
        im=im.resize((640,360),Image.Resampling.LANCZOS);x=(i%2)*640;y=(i//2)*385;canvas.paste(im,(x,y));draw.text((x+12,y+366),name.replace('-',' ').title(),fill=(244,216,166))
canvas.save(root/'overview.png')
concept=root.parent/'concepts/firefight-victory-mouse-homestead-v1.png'
if concept.exists():
    comparison=Image.new('RGB',(1280,1490),(14,23,31))
    with Image.open(concept) as im:
        im.thumbnail((1280,700),Image.Resampling.LANCZOS);comparison.paste(im,((1280-im.width)//2,0))
    comparison.paste(canvas,(0,720));comparison.save(root/'concept-comparison.png')
print('Saved four 8-second GIFs, 32-second montage, overview and concept comparison.')

previous=root/'before-detail/overview.png'
if previous.exists():
    with Image.open(previous) as old:
        sheet=Image.new('RGB',(2560,800),(14,23,31));sheet.paste(old,(0,30));sheet.paste(canvas,(1280,30));d=ImageDraw.Draw(sheet);d.text((20,10),'Previous preview',fill='white');d.text((1300,10),'Detailed town',fill='white');sheet.save(root/'before-after.png')
