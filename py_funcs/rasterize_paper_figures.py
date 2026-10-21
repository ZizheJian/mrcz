import math
import subprocess
from pathlib import Path

from PIL import Image


project_root=Path(__file__).resolve().parents[1]
paper_root=project_root/"论文"
output_root=paper_root/"figs"/"raster300"
column_width_pt=252
points_per_inch=72.27


def pixels(width_pt):
    return math.floor(width_pt/points_per_inch*300)


figures=[
    ("cross-element.png","figs/5/cross-element.pdf",pixels(column_width_pt)),
]


def rasterize_pdf(source,target,width):
    output_stem=target.with_suffix("")
    subprocess.run(["pdftoppm","-png","-singlefile","-f","1","-l","1","-scale-to-x",str(width),"-scale-to-y","-1",str(source),str(output_stem)],check=True)
    with Image.open(target) as image:
        if image.width!=width:
            height=round(image.height*width/image.width)
            image=image.resize((width,height),Image.Resampling.LANCZOS)
        image.convert("RGB").save(target,dpi=(300,300),optimize=True)


def rasterize_image(source,target,width):
    with Image.open(source) as image:
        image=image.convert("RGBA")
        height=round(image.height*width/image.width)
        image=image.resize((width,height),Image.Resampling.LANCZOS)
        background=Image.new("RGB",image.size,"white")
        background.paste(image,mask=image.getchannel("A"))
        background.save(target,dpi=(300,300),optimize=True)


def main():
    output_root.mkdir(parents=True,exist_ok=True)
    for output_name,source_name,width in figures:
        source=paper_root/source_name
        target=output_root/output_name
        if source.suffix.lower()==".pdf":rasterize_pdf(source,target,width)
        else:rasterize_image(source,target,width)
        with Image.open(target) as image:
            effective_dpi=image.width/(width/300)
            print(f"{target.relative_to(paper_root)}: {image.width}x{image.height},effective width DPI={effective_dpi:.2f}")


if __name__=="__main__":main()
