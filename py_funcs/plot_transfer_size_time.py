import csv,os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_style import MPL_POINTS_PER_INCH,PaperFigureSpec,METHOD_COLORS,METHOD_LABELS,TRANSFER_METHODS,apply_paper_style,layout_panel_row,save_paper_png_vertical_crop,save_svg,set_legend_row_gap,style_plot_axis

apply_paper_style("mrcz-transfer-size-time",serif_fonts=["STIXGeneral"])

methods=TRANSFER_METHODS
colors=METHOD_COLORS
labels=METHOD_LABELS
datasets=(("pb0",7),("jicf0",8))
PLOT_ASPECT_RATIO=3.0
LEFT_EDGE_GAP_PT=3.0
RIGHT_EDGE_GAP_PT=3.0
BOTTOM_PLOT_PT=24.0
VERTICAL_EDGE_GAP_PT=3.0


def read_rows(path):
    with open(path,newline="") as input_file:return list(csv.DictReader(input_file))


def draw_panel(axis,rows):
    rows_by_dataset_method={(row["Dataset"],row["Method"]):row for row in rows}
    sizes=[]
    medians=[]
    for dataset,dataset_index in datasets:
        marker="o" if dataset_index==7 else "^"
        for method in methods:
            row=rows_by_dataset_method[(dataset,method)]
            size=float(row["Compressed Bytes"])/(1024**3)
            trials=[float(row[f"Transfer Time Trial {index}"]) for index in range(1,4)]
            throughputs=[size/time for time in trials]
            median=size/float(row["Transfer Time Median"])
            lower=median-min(throughputs)
            upper=max(throughputs)-median
            sizes.append(size)
            medians.append(median)
            axis.errorbar(size,median,yerr=[[lower],[upper]],fmt=marker,markersize=4.0,color=colors[method],markeredgecolor="black",markeredgewidth=0.3,elinewidth=0.5,capsize=1.4,zorder=3)
    size_order=sorted(range(len(sizes)),key=sizes.__getitem__)
    throughput_order=sorted(range(len(medians)),key=medians.__getitem__)
    size_ranks=[0]*len(sizes)
    throughput_ranks=[0]*len(medians)
    for rank,index in enumerate(size_order,1):size_ranks[index]=rank
    for rank,index in enumerate(throughput_order,1):throughput_ranks[index]=rank
    sample_count=len(sizes)
    rank_correlation=1-6*sum((size_rank-throughput_rank)**2 for size_rank,throughput_rank in zip(size_ranks,throughput_ranks))/(sample_count*(sample_count**2-1))
    axis.set_xlim(0,42)
    axis.set_ylim(0,1.2)
    axis.set_xticks([0,20,40])
    axis.set_yticks([0,0.3,0.6,0.9,1.2])
    axis.text(0.97,0.08,f"Spearman $\\rho={rank_correlation:.3f}$",transform=axis.transAxes,ha="right",va="bottom",fontsize=8.0)
    axis.grid(True,which="major",linestyle="--",linewidth=0.35,alpha=0.35)
    style_plot_axis(axis)
    axis.set_xlabel("Transferred Bitstream Size (GiB)")


def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","end_to_end")
    rows=read_rows(os.path.join(result_dir,"end_to_end_summary.csv"))
    spec=PaperFigureSpec.column(1)
    figure,axis=plt.subplots(1,1,figsize=spec.size_in,dpi=spec.dpi)
    draw_panel(axis,rows)
    method_handles=[plt.Line2D([0],[0],marker="o",markersize=3.8,color=colors[method],markeredgecolor="black",markeredgewidth=0.3,linewidth=0,label=labels[method]) for method in methods]
    blank_handle=plt.Line2D([0],[0],linewidth=0,alpha=0,label="")
    dataset_handles=[plt.Line2D([0],[0],marker=marker,markersize=4.0,color="black",markerfacecolor="white",markeredgewidth=0.5,linewidth=0,label=f"Dataset {dataset_index}") for marker,dataset_index in (("o",7),("^",8))]
    axis.set_ylabel("Transfer Throughput (GiB/s)")
    layout_panel_row(figure,[axis],[PLOT_ASPECT_RATIO],bottom_pt=BOTTOM_PLOT_PT,panel_gap_pt=0,left_edge_gap_pt=LEFT_EDGE_GAP_PT,right_edge_gap_pt=RIGHT_EDGE_GAP_PT)
    figure_height_pt=figure.get_figheight()*MPL_POINTS_PER_INCH
    top_plot_pt=axis.get_position().y1*figure_height_pt
    legend=figure.legend(handles=method_handles+[blank_handle]+dataset_handles,loc="lower center",ncol=5,frameon=False,bbox_to_anchor=(0.5,(top_plot_pt+1)/figure_height_pt),handletextpad=0.12,columnspacing=0.20,borderaxespad=0,borderpad=0)
    set_legend_row_gap(legend)
    figure.canvas.draw()
    svg_output=os.path.join(result_dir,"transfer_size_time.svg")
    png_output=os.path.join(result_dir,"transfer_size_time.png")
    paper_output=os.path.join(project_root,"论文","figs","7","transfer-size-time.png")
    save_svg(figure,svg_output,metadata={"Date":None},bbox_inches="tight",pad_inches=0.02)
    save_paper_png_vertical_crop(figure,png_output,spec,edge_gap_pt=VERTICAL_EDGE_GAP_PT)
    save_paper_png_vertical_crop(figure,paper_output,spec,edge_gap_pt=VERTICAL_EDGE_GAP_PT)
    plt.close(figure)
    print(paper_output)


if __name__=="__main__":main()
