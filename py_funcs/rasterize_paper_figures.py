import math
import subprocess
from pathlib import Path

from PIL import Image


project_root=Path(__file__).resolve().parents[1]
paper_root=project_root/"论文"
output_root=paper_root/"figs"/"raster300"
column_width_pt=252
text_width_pt=516
points_per_inch=72.27


def pixels(width_pt):
    return math.floor(width_pt/points_per_inch*300)


figures=[
    ("d3-three-elements.png","figs/3/three_elements.pdf",pixels(column_width_pt)),
    ("d3-end-to-end.png","figs/3/End to End.pdf",pixels(column_width_pt)),
    ("fdwz-diagram.png","figs/5/fdwz diagram.pdf",pixels(column_width_pt)),
    ("cubic-interpolation.png","figs/5/cubic interpolation.pdf",pixels(0.8*column_width_pt)),
    ("fdwz-tree-ablation.png","figs/5/fdwz-tree-ablation.pdf",pixels(column_width_pt)),
    ("cross-element.png","figs/5/cross-element.pdf",pixels(column_width_pt)),
    ("mrcz-diagram.png","figs/6/mrcz diagram.pdf",pixels(column_width_pt)),
    ("transform-comparison-stacked.png","figs/6/transform-comparison-stacked.pdf",pixels(column_width_pt)),
    ("transform-mode-statistics.png","figs/6/transform-mode-statistics.pdf",pixels(column_width_pt)),
    ("mrcz-reordering-ablation.png","figs/6/mrcz-reordering-ablation.pdf",pixels(column_width_pt)),
    ("small-rd.png","figs/7/small-rd.pdf",pixels(text_width_pt)),
    ("large-rd.png","figs/7/large-rd.pdf",pixels(0.8333*text_width_pt)),
    ("combined-throughput.png","figs/7/combined-throughput.pdf",pixels(text_width_pt)),
    ("e2e-time.png","figs/7/e2e_time.pdf",pixels(column_width_pt)),
    ("visualization-cr33.png","figs/7/pb0-T-cr33-visualization.pdf",pixels(text_width_pt)),
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
