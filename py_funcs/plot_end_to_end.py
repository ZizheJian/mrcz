import csv,os
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from plot_style import TRANSFER_METHODS,apply_paper_style,style_plot_axis

apply_paper_style("mrcz-end-to-end")

methods=("Direct",)+TRANSFER_METHODS
components=[("Compression Time","Compression","#e3427a"),("Transfer Time Median","Transfer","#4bacc6"),("Decompression Time","Decompression","#f5aa36")]
titles={"pb0":"(a) Dataset 7","jicf0":"(b) Dataset 8"}
x_limits={"pb0":(0,78),"jicf0":(0,230)}

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
    style_plot_axis(axis,tick_length=2,tick_width=0.5,tick_pad=0.8,spine_width=0.6)
    axis.set_xlim(*x_limits[dataset])
    axis.set_yticks(y)
    axis.set_yticklabels(["MGARD-X-1D" if method=="MGARD-1D" else method for method in ordered_methods],fontsize=7)
    axis.set_xlabel("Time (s)",labelpad=1.0)
    axis.text(0.5,-0.22,titles[dataset],transform=axis.transAxes,ha="center",va="top",fontsize=8)
    for row,y_position in zip(ordered,y):
        total=float(row["End-to-end Time Median"])
        axis.annotate(f"{total:.1f}",(total,y_position),xytext=(2,0),textcoords="offset points",ha="left",va="center",fontsize=7)

def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","end_to_end")
    rows=read_rows(os.path.join(result_dir,"end_to_end_summary.csv"))
    figure,axes=plt.subplots(1,2,figsize=(3.5,2.35),sharey=True)
    for axis,dataset in zip(axes,("pb0","jicf0")):
        draw_panel(axis,[row for row in rows if row["Dataset"]==dataset],dataset)
    axes[0].invert_yaxis()
    axes[1].tick_params(axis="y",labelleft=False)
    handles=[Patch(facecolor=color,edgecolor="black",linewidth=0.45,label=label) for _,label,color in components]
    figure.legend(handles=handles,loc="upper center",ncol=3,frameon=False,bbox_to_anchor=(0.5,0.97),handlelength=1.0,columnspacing=0.7,fontsize=8)
    figure.subplots_adjust(left=0.21,right=0.98,bottom=0.22,top=0.84,wspace=0.15)
    output_path=os.path.join(result_dir,"end_to_end_total.svg")
    figure.savefig(output_path,format="svg",metadata={"Date":None},bbox_inches="tight",pad_inches=0.01)
    with open(output_path) as input_file:
        svg=input_file.read()
    with open(output_path,"w") as output_file:
        output_file.write("\n".join(line.rstrip() for line in svg.splitlines())+"\n")
    paper_output=os.path.join(project_root,"论文","figs","7","e2e_time.pdf")
    for axis in axes:axis.set_rasterization_zorder(1.5)
    figure.savefig(paper_output,format="pdf",dpi=300,bbox_inches="tight",pad_inches=0.01)
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:figure.savefig(os.path.join(preview_dir,"end_to_end_total.png"),format="png",dpi=180,bbox_inches="tight",pad_inches=0.01)
    plt.close(figure)
    print(paper_output)

if __name__=="__main__":
    main()
