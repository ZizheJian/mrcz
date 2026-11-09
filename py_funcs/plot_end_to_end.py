import csv,os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from plot_style import MPL_POINTS_PER_INCH,PANEL_GAP_WITHOUT_Y_DECORATION_PT,PaperFigureSpec,TRANSFER_METHODS,add_panel_captions,apply_paper_style,layout_panel_row,legend_spacing_kwargs,save_paper_png_vertical_crop,save_svg,set_legend_row_gap,style_plot_axis

apply_paper_style("mrcz-end-to-end",serif_fonts=["STIXGeneral"])

methods=("Direct",)+TRANSFER_METHODS
components=[("Compression Time","Compression","#e3427a"),("Transfer Time Median","Transfer","#4bacc6"),("Decompression Time","Decompression","#f5aa36")]
titles={"pb0":"(a) Dataset 7","jicf0":"(b) Dataset 8"}
x_limits={"pb0":(0,78),"jicf0":(0,230)}
PLOT_ASPECT_RATIOS=(1.0,1.0)
LEFT_EDGE_GAP_PT=3.0
RIGHT_EDGE_GAP_PT=3.0
BOTTOM_PLOT_PT=34.0
VERTICAL_EDGE_GAP_PT=3.0

def read_rows(path):
    with open(path,newline="") as input_file:
        return list(csv.DictReader(input_file))

def draw_panel(axis,rows,dataset):
    rows_by_method={row["Method"]:row for row in rows}
    ordered_methods=[method for method in methods if method in rows_by_method]
    ordered=[rows_by_method[method] for method in ordered_methods]
    y=list(range(len(ordered_methods)))
    base=[0.0]*len(ordered_methods)
    for column,label,color in components:
        values=[float(row[column]) for row in ordered]
        bars=axis.barh(y,values,left=base,color=color,edgecolor="black",linewidth=0.45,height=0.70,label=label,zorder=1)
        if column!="Transfer Time Median" and "MGARD-1D" in ordered_methods:
            bars[ordered_methods.index("MGARD-1D")].set_hatch("////")
        base=[first+second for first,second in zip(base,values)]
    axis.grid(True,axis="x",linestyle="--",linewidth=0.35,alpha=0.35)
    axis.set_axisbelow(True)
    style_plot_axis(axis)
    axis.set_xlim(*x_limits[dataset])
    axis.set_yticks(y)
    axis.set_yticklabels(["MGARD-X-1D" if method=="MGARD-1D" else method for method in ordered_methods],fontsize=8)
    axis.set_xlabel("Time (s)")
    for row,y_position in zip(ordered,y):
        total=float(row["End-to-end Time Median"])
        axis.annotate(f"{total:.1f}",(total,y_position),xytext=(2,0),textcoords="offset points",ha="left",va="center",fontsize=8)

def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","end_to_end")
    rows=read_rows(os.path.join(result_dir,"end_to_end_summary.csv"))
    spec=PaperFigureSpec.column(1)
    figure,axes=plt.subplots(1,2,figsize=spec.size_in,dpi=spec.dpi,sharey=True)
    for axis,dataset in zip(axes,("pb0","jicf0")):
        draw_panel(axis,[row for row in rows if row["Dataset"]==dataset],dataset)
    axes[0].invert_yaxis()
    axes[1].tick_params(axis="y",labelleft=False)
    handles=[Patch(facecolor=color,edgecolor="black",linewidth=0.45,label=label) for _,label,color in components]
    layout_panel_row(figure,axes,PLOT_ASPECT_RATIOS,bottom_pt=BOTTOM_PLOT_PT,panel_gap_pt=PANEL_GAP_WITHOUT_Y_DECORATION_PT,left_edge_gap_pt=LEFT_EDGE_GAP_PT,right_edge_gap_pt=RIGHT_EDGE_GAP_PT)
    add_panel_captions(figure,axes,[titles[dataset] for dataset in ("pb0","jicf0")])
    figure_height_pt=figure.get_figheight()*MPL_POINTS_PER_INCH
    top_plot_pt=max(axis.get_position().y1 for axis in axes)*figure_height_pt
    legend=figure.legend(handles=handles,loc="lower center",ncol=3,frameon=False,bbox_to_anchor=(0.5,(top_plot_pt+1)/figure_height_pt),borderaxespad=0,borderpad=0,**legend_spacing_kwargs(8))
    set_legend_row_gap(legend)
    figure.canvas.draw()
    output_path=os.path.join(result_dir,"end_to_end_total.svg")
    png_output=os.path.join(result_dir,"end_to_end_total.png")
    paper_output=os.path.join(project_root,"overleaf","figs","7","e2e_time.png")
    save_svg(figure,output_path,metadata={"Date":None},bbox_inches="tight",pad_inches=0.02)
    save_paper_png_vertical_crop(figure,png_output,spec,edge_gap_pt=VERTICAL_EDGE_GAP_PT)
    save_paper_png_vertical_crop(figure,paper_output,spec,edge_gap_pt=VERTICAL_EDGE_GAP_PT)
    plt.close(figure)
    print(paper_output)

if __name__=="__main__":
    main()
