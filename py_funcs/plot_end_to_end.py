import csv,os
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams["svg.hashsalt"]="mrcz-end-to-end"
plt.rcParams.update({"font.family":"serif","font.serif":["Nimbus Roman"],"mathtext.fontset":"stix","font.size":8.0,"axes.labelsize":8.0,"axes.titlesize":8.0,"xtick.labelsize":8.0,"ytick.labelsize":8.0,"legend.fontsize":8.0,"pdf.fonttype":42,"ps.fonttype":42})

methods=["Direct","MRCZ","FDWZ","SZ3-4D","QoZ-4D","ZFP-3D-fused","SPERR-2D"]
components=[("Compression Time","Compression","#e3427a"),("Transfer Time Median","Transfer","#4bacc6"),("Decompression Time","Decompression","#f5aa36")]
titles={"pb0":"(a) Dataset 7","jicf0":"(b) Dataset 8"}
y_limits={"pb0":(0,65),"jicf0":(0,195)}

def read_rows(path):
    with open(path,newline="") as input_file:
        return list(csv.DictReader(input_file))

def draw_panel(axis,rows,dataset):
    rows_by_method={row["Method"]:row for row in rows}
    ordered_methods=[method for method in methods if method in rows_by_method]
    ordered=[rows_by_method[method] for method in ordered_methods]
    x=list(range(len(ordered_methods)))
    base=[0.0]*len(ordered_methods)
    top_bars=None
    for column,label,color in components:
        values=[float(row[column]) for row in ordered]
        top_bars=axis.bar(x,values,bottom=base,color=color,edgecolor="black",linewidth=0.45,width=0.74,label=label,zorder=1)
        base=[first+second for first,second in zip(base,values)]
    axis.grid(True,axis="y",linestyle="--",linewidth=0.35,alpha=0.35)
    axis.set_axisbelow(True)
    axis.tick_params(direction="in",length=2,width=0.5,pad=0.8)
    for spine in axis.spines.values():spine.set_linewidth(0.6)
    axis.set_ylim(*y_limits[dataset])
    axis.set_xticks(x)
    axis.set_xticklabels(ordered_methods,rotation=55,ha="right",rotation_mode="anchor",fontsize=7)
    axis.text(0.5,-0.47,titles[dataset],transform=axis.transAxes,ha="center",va="top",fontsize=8)
    axis.bar_label(top_bars,labels=[f'{float(row["End-to-end Time Median"]):.1f}' for row in ordered],padding=0.5,fontsize=7)

def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","end_to_end")
    rows=read_rows(os.path.join(result_dir,"end_to_end_summary.csv"))
    figure,axes=plt.subplots(1,2,figsize=(3.5,2.25))
    for axis,dataset in zip(axes,("pb0","jicf0")):
        draw_panel(axis,[row for row in rows if row["Dataset"]==dataset],dataset)
    axes[0].set_ylabel("Time (s)",labelpad=2.0)
    handles=[Patch(facecolor=color,edgecolor="black",linewidth=0.45,label=label) for _,label,color in components]
    figure.legend(handles=handles,loc="upper center",ncol=3,frameon=False,bbox_to_anchor=(0.5,0.97),handlelength=1.0,columnspacing=0.7,fontsize=8)
    figure.subplots_adjust(left=0.105,right=0.995,bottom=0.36,top=0.86,wspace=0.20)
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
