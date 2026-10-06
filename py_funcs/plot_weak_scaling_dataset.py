import argparse,csv,os
import matplotlib.pyplot as plt
from plot_style import METHOD_COLORS,METHOD_MARKERS,TIMING_METHODS,apply_paper_style

apply_paper_style()

methods=TIMING_METHODS
colors=METHOD_COLORS
markers=METHOD_MARKERS
dataset_titles={"pb0":"pb0 (100 dB)","jicf":"jicf0 (150 dB)"}

# 读取单个数据集的weak-scaling汇总结果。
def read_rows(path):
    with open(path,newline="") as f:
        return list(csv.DictReader(f))

# 绘制压缩或解压的综合parallel efficiency。
def draw_panel(axis,rows,task):
    metric=f"{task} Efficiency"
    cpu_values=sorted({int(row["CPUs"]) for row in rows})
    positions={cpus:index for index,cpus in enumerate(cpu_values)}
    for method in methods:
        method_rows=sorted((row for row in rows if row["Method"]==method),key=lambda row:int(row["CPUs"]))
        if len(method_rows)<2:
            continue
        x=[positions[int(row["CPUs"])] for row in method_rows]
        y=[float(row[metric]) for row in method_rows]
        axis.plot(x,y,color=colors[method],marker=markers[method],markersize=7 if method=="MRCZ" else 4.5,linewidth=2.8 if method=="MRCZ" else 1.6,label=method)
    axis.axhline(1,color="#9ca3af",linestyle="--",linewidth=1)
    axis.set_xticks(range(len(cpu_values)),cpu_values)
    axis.set_ylim(bottom=0)
    axis.set_xlabel("CPU count")
    axis.set_ylabel("Parallel efficiency")
    axis.set_title(task)
    axis.grid(True,alpha=0.25,linewidth=0.5)

# 生成单个数据集的压缩和解压weak-scaling efficiency图。
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("dataset",choices=sorted(dataset_titles))
    args=parser.parse_args()
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","weak_scaling")
    rows=read_rows(os.path.join(result_dir,f"weak_{args.dataset}_summary.csv"))
    figure,axes=plt.subplots(1,2,figsize=(13,5.2))
    draw_panel(axes[0],rows,"Compression")
    draw_panel(axes[1],rows,"Decompression")
    handles=[plt.Line2D([0],[0],color=colors[method],marker=markers[method],linewidth=2.8 if method=="MRCZ" else 1.6,markersize=7 if method=="MRCZ" else 4.5,label=method) for method in methods]
    figure.legend(handles=handles,loc="upper center",ncol=len(methods),frameon=False,bbox_to_anchor=(0.5,0.99),fontsize=9)
    figure.suptitle(f"{dataset_titles[args.dataset]} weak scaling",fontsize=16,y=0.92)
    figure.tight_layout(rect=(0,0,1,0.84))
    output_path=os.path.join(result_dir,f"weak_{args.dataset}_efficiency.svg")
    figure.savefig(output_path,format="svg",bbox_inches="tight")
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:
        figure.savefig(os.path.join(preview_dir,f"weak_{args.dataset}_efficiency.png"),format="png",dpi=160,bbox_inches="tight")
    plt.close(figure)
    print(output_path)

if __name__=="__main__":
    main()
