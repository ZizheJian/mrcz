import csv,math,os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,LogLocator,MaxNLocator,NullLocator
from plot_style import METHOD_COLORS,MPL_POINTS_PER_INCH,PANEL_GAP_WITH_YTICKS_PT,PaperFigureSpec,add_shared_xlabel,apply_paper_style,layout_panel_grid,legend_spacing_kwargs,panel_aspect_ratios,save_paper_png_vertical_crop,save_svg,style_plot_axis

ROW_CONTENT_GAP_PT=6.0
PANEL_CAPTION_GAP_PT=1.0
SHARED_XLABEL_GAP_PT=2.0
SHARED_YLABEL_GAP_PT=1.0
LEGEND_GAP_PT=1.0
OUTER_EDGE_GAP_PT=3.0
BOTTOM_PLOT_PT=30.0
VERTICAL_EDGE_GAP_PT=3.0

array_layouts=["1D","2D","3D-fused","3D-per-element","4D"]
layouts=array_layouts+["4D+GLL"]
methods=["MRCZ","FDWZ"]+[f"{compressor}-{layout}" for compressor in ["SZ3","SZ2","QoZ","ZFP"] for layout in array_layouts]+[f"MGARD-{layout}" for layout in layouts]+["SPERR-2D","SPERR-3D-fused","SPERR-3D-per-element"]
family_colors={"MRCZ":METHOD_COLORS["MRCZ"],"FDWZ":METHOD_COLORS["FDWZ"],"SZ3":METHOD_COLORS["SZ3-4D"],"SZ2":METHOD_COLORS["SZ2-3D-fused"],"QoZ":METHOD_COLORS["QoZ-4D"],"ZFP":METHOD_COLORS["ZFP-3D-fused"],"MGARD":METHOD_COLORS["MGARD-1D"],"SPERR":METHOD_COLORS["SPERR-2D"]}
layout_styles={"1D":(0,(5,1,1,1)),"2D":"--","3D-fused":"-.","3D-per-element":":","4D":"-","4D+GLL":(0,(3,1,1,1))}
fields=[
    ("1UX","mixlay0/mixlay0.f00001_UX.csv"),("1UY","mixlay0/mixlay0.f00001_UY.csv"),("1P","mixlay0/mixlay0.f00001_P.csv"),("1T","mixlay0/mixlay0.f00001_T.csv"),("1S00","mixlay0/mixlay0.f00001_S00.csv"),("1S01","mixlay0/mixlay0.f00001_S01.csv"),
    ("2UZ","test0/test0.f00001_UZ.csv"),("2T","test0/test0.f00001_T.csv"),("2S00","test0/test0.f00001_S00.csv"),("2S01","test0/test0.f00001_S01.csv"),("3UX","tc_channel0/tc_channel0.f00001_UX.csv"),("3UY","tc_channel0/tc_channel0.f00001_UY.csv"),
    ("3UZ","tc_channel0/tc_channel0.f00001_UZ.csv"),("Dataset 4","mfem_ex0/sol.gf_no_name_field.csv"),("Dataset 5","mfem_ex1/sol.gf_no_name_field.csv"),("6P","nektar_cylinder/Cylinder-Base.vtu_p.csv"),("6U","nektar_cylinder/Cylinder-Base.vtu_u.csv"),("6V","nektar_cylinder/Cylinder-Base.vtu_v.csv")]

def read_rows(path):
    with open(path,newline="") as f:
        return list(csv.DictReader(f))

def read_method_rows(result_dir,relative_path,method,default_rows):
    if not method.startswith("MGARD-"):
        return default_rows,method
    coordinate_type="GLL" if method.endswith("+GLL") else "uniform"
    stored_method="MGARD-4D" if method.endswith("+GLL") else method
    path=os.path.join(result_dir,"mgard_layout_selection",os.path.dirname(relative_path),coordinate_type,os.path.basename(relative_path))
    return read_rows(path),stored_method

def prepare_axes():
    spec=PaperFigureSpec.text(1)
    figure,axes=plt.subplots(3,6,figsize=spec.size_in,dpi=spec.dpi)
    return figure,axes,spec

def method_style(method):
    if method in ["MRCZ","FDWZ"]:
        return family_colors[method],"-",0.4
    family,layout=method.split("-",1)
    return family_colors[family],layout_styles[layout],0.4

def finite_pair(first,second):
    try:
        pair=float(first),float(second)
    except (TypeError,ValueError):
        return None
    return pair if all(math.isfinite(value) for value in pair) else None

def percentile(values,fraction):
    values=sorted(values)
    position=(len(values)-1)*fraction
    lower=int(position)
    upper=min(lower+1,len(values)-1)
    weight=position-lower
    return values[lower]*(1-weight)+values[upper]*weight

def focus_rd_axes(axis,method_points):
    all_points=[point for points in method_points.values() for point in points]
    primary_points=method_points["MRCZ"]+method_points["FDWZ"]
    br_left=min(point[0] for point in all_points)*0.9
    br_right=max(percentile([point[0] for point in all_points],0.8),max(point[0] for point in primary_points))*1.05
    psnr_values=[point[1] for point in all_points]
    psnr_top=max(percentile(psnr_values,0.95),max(point[1] for point in primary_points))
    psnr_bottom=min(psnr_values)
    psnr_padding=max(1,(psnr_top-psnr_bottom)*0.03)
    axis.set_xscale("log")
    axis.set_xlim(br_left,br_right)
    axis.set_ylim(psnr_bottom-psnr_padding,psnr_top+psnr_padding)

def visible_bboxes(artists,renderer):
    return [artist.get_window_extent(renderer) for artist in artists if artist.get_visible() and artist.get_text()]

def place_shared_ylabel(figure,axes,label):
    figure.canvas.draw()
    renderer=figure.canvas.get_renderer()
    tick_left=min(bbox.x0 for axis in axes[:,0] for bbox in visible_bboxes(axis.get_yticklabels(which="both"),renderer))
    grid_bottom=min(axis.get_position().y0 for axis in axes.flat)
    grid_top=max(axis.get_position().y1 for axis in axes.flat)
    text=figure.text(0,(grid_bottom+grid_top)/2,label,rotation=90,ha="center",va="center")
    figure.canvas.draw()
    bbox=text.get_window_extent(figure.canvas.get_renderer())
    center_x=tick_left-SHARED_YLABEL_GAP_PT*figure.dpi/MPL_POINTS_PER_INCH-bbox.width/2
    text.set_x(center_x/figure.bbox.width)
    figure.canvas.draw()
    return text

def finish_figure(figure,axes,spec,xlabel,ylabel,output_path,paper_output):
    family_handles=[]
    for family in ["MRCZ","FDWZ","SZ3","SZ2","QoZ","ZFP","SPERR","MGARD"]:
        family_handles.append(plt.Line2D([0],[0],color=family_colors[family],linestyle="-",linewidth=0.8,label="MGARD-X" if family=="MGARD" else family))
    layout_handles=[plt.Line2D([0],[0],color="#4b5563",linestyle=layout_styles[layout],linewidth=0.8,label=layout) for layout in layouts]
    ylabel_artist=figure.text(0,0,ylabel,rotation=90,ha="center",va="center")
    figure.canvas.draw()
    ylabel_width_pt=ylabel_artist.get_window_extent(figure.canvas.get_renderer()).width*MPL_POINTS_PER_INCH/figure.dpi
    ylabel_artist.remove()
    layout_panel_grid(figure,axes,panel_aspect_ratios(1,1,1,1,1,1),bottom_pt=BOTTOM_PLOT_PT,column_gap_pt=PANEL_GAP_WITH_YTICKS_PT,row_content_gap_pt=ROW_CONTENT_GAP_PT,left_edge_gap_pt=OUTER_EDGE_GAP_PT+ylabel_width_pt+SHARED_YLABEL_GAP_PT,right_edge_gap_pt=OUTER_EDGE_GAP_PT)
    shared_xlabel=add_shared_xlabel(figure,axes[-1,:],xlabel,gap_pt=SHARED_XLABEL_GAP_PT,source_artists=[axis.xaxis.label for axis in axes[-1,:]])
    shared_ylabel=place_shared_ylabel(figure,axes,ylabel)
    figure_height_pt=figure.get_figheight()*MPL_POINTS_PER_INCH
    top_plot_pt=max(axis.get_position().y1 for axis in axes[0,:])*figure_height_pt
    layout_legend=figure.legend(handles=layout_handles,loc="lower center",ncol=6,frameon=False,bbox_to_anchor=(0.5,(top_plot_pt+LEGEND_GAP_PT)/figure_height_pt),borderaxespad=0,borderpad=0,**legend_spacing_kwargs(9.6))
    figure.canvas.draw()
    layout_top_pt=layout_legend.get_window_extent(figure.canvas.get_renderer()).y1*MPL_POINTS_PER_INCH/figure.dpi
    family_legend=figure.legend(handles=family_handles,loc="lower center",ncol=8,frameon=False,bbox_to_anchor=(0.5,(layout_top_pt+LEGEND_GAP_PT)/figure_height_pt),borderaxespad=0,borderpad=0,**legend_spacing_kwargs(8))
    figure.canvas.draw()
    save_svg(figure,output_path,metadata={"Date":None},bbox_inches="tight",pad_inches=0.02)
    os.makedirs(os.path.dirname(paper_output),exist_ok=True)
    save_paper_png_vertical_crop(figure,paper_output,spec,edge_gap_pt=VERTICAL_EDGE_GAP_PT)
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:
        figure.savefig(os.path.join(preview_dir,"small_rd.png"),format="png",dpi=160,bbox_inches="tight")
    plt.close(figure)

def plot_rd(result_dir):
    apply_paper_style("mrcz-small-rd",serif_fonts=["STIXGeneral"])
    figure,axes,spec=prepare_axes()
    for panel_index,(axis,(title,relative_path)) in enumerate(zip(axes.flat,fields)):
        rows=read_rows(os.path.join(result_dir,relative_path))
        method_points={}
        for method in methods:
            method_rows,stored_method=read_method_rows(result_dir,relative_path,method,rows)
            points=sorted(point for row in method_rows if (point:=finite_pair(row.get(f"{stored_method}_BR"),row.get(f"{stored_method}_PSNR"))) is not None)
            method_points[method]=points
            color,linestyle,linewidth=method_style(method)
            axis.plot([point[0] for point in points],[point[1] for point in points],color=color,linestyle=linestyle,linewidth=linewidth,zorder=1)
        focus_rd_axes(axis,method_points)
        axis.set_xlabel(f"({chr(ord('a')+panel_index)}) {title}")
        axis.grid(True,which="major",linestyle="--",linewidth=0.35,alpha=0.35)
        log_locator=LogLocator(base=10,numticks=4)
        x_left,x_right=axis.get_xlim()
        axis.xaxis.set_major_locator(FixedLocator([tick for tick in log_locator.tick_values(x_left,x_right) if x_left<=tick<=x_right]))
        axis.xaxis.set_minor_locator(NullLocator())
        axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
        style_plot_axis(axis,x_label_pad=PANEL_CAPTION_GAP_PT)
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    finish_figure(figure,axes,spec,"Bit rate (log scale)","PSNR (dB)",os.path.join(result_dir,"small_rd.svg"),os.path.join(project_root,"论文","figs","7","small-rd.png"))

def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","RD")
    plot_rd(result_dir)

if __name__=="__main__":
    main()
