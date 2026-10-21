import argparse,csv,math,os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,LogLocator,MaxNLocator,NullLocator
from plot_style import METHOD_COLORS,MPL_POINTS_PER_INCH,PaperFigureSpec,RD_METHODS,add_shared_xlabel,apply_paper_style,layout_panel_grid,legend_spacing_kwargs,panel_aspect_ratios,save_paper_png_vertical_crop,save_svg,style_plot_axis

PANEL_GAP_PT=4.0
ROW_CONTENT_GAP_PT=6.0
PANEL_CAPTION_GAP_PT=1.0
SHARED_XLABEL_GAP_PT=2.0
SHARED_YLABEL_GAP_PT=1.0
LEGEND_GAP_PT=1.0
OUTER_EDGE_GAP_PT=3.0
BOTTOM_PLOT_PT=30.0
VERTICAL_EDGE_GAP_PT=3.0

methods=RD_METHODS
colors=METHOD_COLORS
datasets={
    "pb0":("pb0",["UX","UY","UZ","P","T"],"pb0.f00001"),
    "jicf0":("jicf0",["UX","UY","UZ","P","S00"],"jicf0.f00000")}

# 读取一个field的rate-distortion CSV。
def read_rows(path):
    with open(path,newline="") as f:
        return list(csv.DictReader(f))

# MGARD的大数据RD来自可运行的8线程配置，BR和PSNR不依赖线程数。
def method_rows(project_root,dataset,path,method,default_rows):
    if method!="MGARD-1D":return default_rows
    mgard_path=os.path.join(project_root,"results","small_core_throughput",dataset,os.path.basename(path))
    return read_rows(mgard_path)

# 忽略缺失或非有限的CSV值。
def finite_pair(first,second):
    try:
        pair=float(first),float(second)
    except (TypeError,ValueError):
        return None
    return pair if all(math.isfinite(value) for value in pair) else None

# 返回一组数的线性插值分位数，用于截断明显偏离主要比较范围的曲线端点。
def percentile(values,fraction):
    values=sorted(values)
    position=(len(values)-1)*fraction
    lower=int(position)
    upper=min(lower+1,len(values)-1)
    weight=position-lower
    return values[lower]*(1-weight)+values[upper]*weight

# 用对数BR轴展开低码率区域，并聚焦MRCZ、FDWZ及大多数baseline共同覆盖的范围。
def focus_axes(axis,method_points):
    all_points=[point for points in method_points.values() for point in points]
    primary_points=method_points["MRCZ"]+method_points["FDWZ"]
    br_left=min(point[0] for point in all_points)*0.9
    br_right=max(percentile([point[0] for point in all_points],0.9),max(point[0] for point in primary_points))*1.05
    psnr_values=[point[1] for point in all_points]
    psnr_top=max(percentile(psnr_values,0.95),max(point[1] for point in primary_points))
    psnr_bottom=min(psnr_values)
    padding=max(1,(psnr_top-psnr_bottom)*0.03)
    axis.set_xscale("log")
    axis.set_xlim(br_left,br_right)
    axis.set_ylim(psnr_bottom-padding,psnr_top+padding)

# 为一个大数据集的五个field绘制rate-distortion曲线。
def plot_dataset(result_dir,dataset,project_root):
    directory,fields,prefix=datasets[dataset]
    figure,axes=plt.subplots(1,len(fields),figsize=(19,4.5))
    for axis,field in zip(axes,fields):
        rows=read_rows(os.path.join(result_dir,directory,f"{prefix}_{field}.csv"))
        method_points={}
        for method in methods:
            source_rows=method_rows(project_root,dataset,os.path.join(result_dir,directory,f"{prefix}_{field}.csv"),method,rows)
            points=sorted(point for row in source_rows if (point:=finite_pair(row.get(f"{method}_BR"),row.get(f"{method}_PSNR"))) is not None)
            method_points[method]=points
            axis.plot([point[0] for point in points],[point[1] for point in points],color=colors[method],linewidth=1.6,label=method)
        focus_axes(axis,method_points)
        axis.set_title(field)
        axis.grid(True,alpha=0.25,linewidth=0.5)
    handles=[plt.Line2D([0],[0],color=colors[method],linewidth=1.6,label="MGARD-X-1D" if method=="MGARD-1D" else method) for method in methods]
    figure.legend(handles=handles,loc="upper center",ncol=len(methods),frameon=False,bbox_to_anchor=(0.5,0.99),fontsize=8)
    figure.suptitle(f"{dataset} rate-distortion",fontsize=16,y=0.9)
    figure.supxlabel("Bit rate")
    figure.supylabel("PSNR (dB)")
    figure.tight_layout(rect=(0.015,0.04,1,0.82))
    output_path=os.path.join(result_dir,f"{dataset}_rd.svg")
    save_svg(figure,output_path,bbox_inches="tight",metadata={"Date":None})
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:
        figure.savefig(os.path.join(preview_dir,f"{dataset}_rd.png"),format="png",dpi=160,bbox_inches="tight")
    plt.close(figure)
    print(output_path)

# 将两个大数据集绘制为论文中的2x5双栏图。
def plot_combined(result_dir,project_root):
    spec=PaperFigureSpec.text(1,width_fraction=5/6)
    figure,axes=plt.subplots(2,5,figsize=spec.size_in,dpi=spec.dpi)
    panels=[]
    for dataset_index,dataset in ((7,"pb0"),(8,"jicf0")):
        directory,fields,prefix=datasets[dataset]
        panels.extend((f"{dataset_index}{field}",os.path.join(result_dir,directory,f"{prefix}_{field}.csv")) for field in fields)
    for panel_index,(axis,(title,path)) in enumerate(zip(axes.flat,panels)):
        rows=read_rows(path)
        method_points={}
        for method in methods:
            dataset="pb0" if title.startswith("7") else "jicf0"
            source_rows=method_rows(project_root,dataset,path,method,rows)
            points=sorted(point for row in source_rows if (point:=finite_pair(row.get(f"{method}_BR"),row.get(f"{method}_PSNR"))) is not None)
            method_points[method]=points
            axis.plot([point[0] for point in points],[point[1] for point in points],color=colors[method],linewidth=0.8,label=method,zorder=1)
        focus_axes(axis,method_points)
        axis.set_xlabel(f"({chr(ord('a')+panel_index)}) {title}")
        axis.grid(True,which="major",linestyle="--",linewidth=0.35,alpha=0.35)
        log_locator=LogLocator(base=10,numticks=4)
        x_left,x_right=axis.get_xlim()
        axis.xaxis.set_major_locator(FixedLocator([tick for tick in log_locator.tick_values(x_left,x_right) if x_left<=tick<=x_right]))
        axis.xaxis.set_minor_locator(NullLocator())
        axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
        style_plot_axis(axis,x_label_pad=PANEL_CAPTION_GAP_PT)
    handles=[plt.Line2D([0],[0],color=colors[method],linewidth=0.8,label="MGARD-X-1D" if method=="MGARD-1D" else method) for method in methods]
    provisional_ylabel=figure.text(0,0,"PSNR (dB)",rotation=90,ha="center",va="center")
    figure.canvas.draw()
    ylabel_width_pt=provisional_ylabel.get_window_extent(figure.canvas.get_renderer()).width*MPL_POINTS_PER_INCH/figure.dpi
    provisional_ylabel.remove()
    layout_panel_grid(figure,axes,panel_aspect_ratios(1,1,1,1,1),bottom_pt=BOTTOM_PLOT_PT,column_gap_pt=PANEL_GAP_PT,row_content_gap_pt=ROW_CONTENT_GAP_PT,left_edge_gap_pt=OUTER_EDGE_GAP_PT+ylabel_width_pt+SHARED_YLABEL_GAP_PT,right_edge_gap_pt=OUTER_EDGE_GAP_PT)
    add_shared_xlabel(figure,axes[-1,:],"Bit rate (log scale)",gap_pt=SHARED_XLABEL_GAP_PT,source_artists=[axis.xaxis.label for axis in axes[-1,:]])
    figure.canvas.draw()
    renderer=figure.canvas.get_renderer()
    tick_left=min(bbox.x0 for axis in axes[:,0] for bbox in [label.get_window_extent(renderer) for label in axis.get_yticklabels(which="both") if label.get_visible() and label.get_text()])
    grid_bottom=min(axis.get_position().y0 for axis in axes.flat)
    grid_top=max(axis.get_position().y1 for axis in axes.flat)
    shared_ylabel=figure.text(0,(grid_bottom+grid_top)/2,"PSNR (dB)",rotation=90,ha="center",va="center")
    figure.canvas.draw()
    ylabel_bbox=shared_ylabel.get_window_extent(figure.canvas.get_renderer())
    shared_ylabel.set_x((tick_left-SHARED_YLABEL_GAP_PT*figure.dpi/MPL_POINTS_PER_INCH-ylabel_bbox.width/2)/figure.bbox.width)
    figure_height_pt=figure.get_figheight()*MPL_POINTS_PER_INCH
    top_plot_pt=max(axis.get_position().y1 for axis in axes[0,:])*figure_height_pt
    figure.legend(handles=handles,loc="lower center",ncol=len(methods),frameon=False,bbox_to_anchor=(0.5,(top_plot_pt+LEGEND_GAP_PT)/figure_height_pt),borderaxespad=0,borderpad=0,**legend_spacing_kwargs(8))
    figure.canvas.draw()
    save_svg(figure,os.path.join(result_dir,"large_rd.svg"),metadata={"Date":None},bbox_inches="tight",pad_inches=0.02)
    paper_output=os.path.join(project_root,"论文","figs","7","large-rd.png")
    os.makedirs(os.path.dirname(paper_output),exist_ok=True)
    save_paper_png_vertical_crop(figure,paper_output,spec,edge_gap_pt=VERTICAL_EDGE_GAP_PT)
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:figure.savefig(os.path.join(preview_dir,"large_rd.png"),format="png",dpi=180)
    plt.close(figure)
    print(paper_output)

# 选择并绘制一个大数据集。
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("dataset",nargs="?",choices=["combined"]+sorted(datasets),default="combined")
    args=parser.parse_args()
    apply_paper_style("mrcz-large-rd",serif_fonts=["STIXGeneral"])
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","RD")
    if args.dataset=="combined":plot_combined(result_dir,project_root)
    else:plot_dataset(result_dir,args.dataset,project_root)

if __name__=="__main__":
    main()
