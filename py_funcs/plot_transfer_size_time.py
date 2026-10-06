import csv,os

import matplotlib.pyplot as plt
from plot_style import METHOD_COLORS,METHOD_LABELS,TRANSFER_METHODS,apply_paper_style,style_plot_axis

apply_paper_style("mrcz-transfer-size-time")

methods=TRANSFER_METHODS
colors=METHOD_COLORS
labels=METHOD_LABELS
datasets=(("pb0",7),("jicf0",8))


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
    style_plot_axis(axis,spine_width=0.6)
    axis.set_xlabel("Transferred Data (GiB)",labelpad=1.0)


def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","end_to_end")
    rows=read_rows(os.path.join(result_dir,"end_to_end_summary.csv"))
    figure,axis=plt.subplots(1,1,figsize=(3.5,1.9))
    draw_panel(axis,rows)
    handles=[plt.Line2D([0],[0],marker="o",markersize=3.8,color=colors[method],markeredgecolor="black",markeredgewidth=0.3,linewidth=0,label=labels[method]) for method in methods]
    figure.legend(handles=handles,loc="upper center",ncol=4,frameon=False,bbox_to_anchor=(0.5,0.99),handletextpad=0.12,columnspacing=0.45,borderaxespad=0)
    dataset_handles=[plt.Line2D([0],[0],marker=marker,markersize=4.0,color="black",markerfacecolor="white",markeredgewidth=0.5,linewidth=0,label=f"Dataset {dataset_index}") for marker,dataset_index in (("o",7),("^",8))]
    figure.legend(handles=dataset_handles,loc="upper center",ncol=2,frameon=False,bbox_to_anchor=(0.5,0.82),handletextpad=0.15,columnspacing=0.65,borderaxespad=0)
    axis.set_ylabel("Transfer Throughput (GiB/s)",labelpad=1.5)
    figure.tight_layout(rect=(0,0,1,0.70),pad=0.2)
    svg_output=os.path.join(result_dir,"transfer_size_time.svg")
    pdf_output=os.path.join(project_root,"论文","figs","7","transfer-size-time.pdf")
    figure.savefig(svg_output,format="svg",metadata={"Date":None},bbox_inches="tight",pad_inches=0.02)
    with open(svg_output) as input_file:svg=input_file.read()
    with open(svg_output,"w") as output_file:output_file.write("\n".join(line.rstrip() for line in svg.splitlines())+"\n")
    figure.savefig(pdf_output,format="pdf",bbox_inches="tight",pad_inches=0.02)
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:figure.savefig(os.path.join(preview_dir,"transfer_size_time.png"),format="png",dpi=180,bbox_inches="tight",pad_inches=0.02)
    plt.close(figure)
    print(pdf_output)


if __name__=="__main__":main()
