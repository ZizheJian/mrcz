import os,sys
import matplotlib.pyplot as plt

project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
paperfig_source=os.path.join(project_root,"third_party/paperfig/src")
if paperfig_source not in sys.path:sys.path.insert(0,paperfig_source)

from paperfig import DEFAULT_LEGEND_GAP_PT,DEFAULT_LEGEND_ROW_GAP_PT,DEFAULT_PANEL_ASPECT_RATIO,MPL_POINTS_PER_INCH,PAPER_COLUMN_WIDTH_PT,PAPER_FONT_SIZE,PAPER_TEXT_WIDTH_IN,PAPER_TEXT_WIDTH_PT,RASTER_DPI,TEX_POINTS_PER_INCH,PaperFigureSpec,add_panel_captions,add_shared_xlabel,align_xlabels,apply_paper_style as _apply_paper_style,layout_panel_grid,layout_panel_row,legend_above,panel_aspect_ratios,paper_subplots,save_paper_png,save_paper_png_vertical_crop,set_legend_row_gap,style_plot_axis,validate_figure_bounds

METHOD_COLORS={
    "MRCZ":"#111827",
    "FDWZ":"#e15759",
    "SZ3-4D":"#4e79a7",
    "SZ2-3D-fused":"#76b7b2",
    "QoZ-4D":"#59a14f",
    "ZFP-3D-fused":"#f28e2b",
    "SPERR-2D":"#b07aa1",
    "MGARD-1D":"#9c755f"}
METHOD_LABELS={
    "MRCZ":"MRCZ",
    "FDWZ":"FDWZ",
    "SZ3-4D":"SZ3",
    "SZ2-3D-fused":"SZ2",
    "QoZ-4D":"QoZ",
    "ZFP-3D-fused":"ZFP",
    "SPERR-2D":"SPERR",
    "MGARD-1D":"MGARD-X"}
METHOD_MARKERS={"MRCZ":"*","FDWZ":"o","SZ3-4D":"^","QoZ-4D":"D","ZFP-3D-fused":"v","SPERR-2D":"P"}
TIMING_METHODS=("MRCZ","FDWZ","SZ3-4D","QoZ-4D","ZFP-3D-fused","SPERR-2D")
RD_METHODS=("MRCZ","FDWZ","SZ3-4D","SZ2-3D-fused","QoZ-4D","ZFP-3D-fused","SPERR-2D","MGARD-1D")
TRANSFER_METHODS=("MRCZ","FDWZ","SZ3-4D","QoZ-4D","ZFP-3D-fused","SPERR-2D","MGARD-1D")

def apply_paper_style(svg_hashsalt=None,font_size=PAPER_FONT_SIZE,tick_size=None,legend_size=None,serif_fonts=None):
    if serif_fonts is None:serif_fonts=["Nimbus Roman"]
    return _apply_paper_style(svg_hashsalt,font_size,tick_size,legend_size,serif_fonts)

def save_svg(figure,path,**kwargs):
    figure.savefig(path,format="svg",**kwargs)
    with open(path,encoding="utf-8") as file_obj:content=file_obj.read()
    with open(path,"w",encoding="utf-8",newline="\n") as file_obj:file_obj.write("\n".join(line.rstrip() for line in content.splitlines())+"\n")
