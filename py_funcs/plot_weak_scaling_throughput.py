import csv,os
import matplotlib.pyplot as plt

methods=["MRCZ","FDWZ","SZ3-4D","QoZ-4D","ZFP-3D-fused","SPERR-2D"]
colors={"MRCZ":"#111827","FDWZ":"#e15759","SZ3-4D":"#4e79a7","QoZ-4D":"#59a14f","ZFP-3D-fused":"#f28e2b","SPERR-2D":"#b07aa1"}
markers={"MRCZ":"*","FDWZ":"o","SZ3-4D":"^","QoZ-4D":"D","ZFP-3D-fused":"v","SPERR-2D":"P"}

# 读取weak-scaling汇总结果。
def read_rows(path):
    with open(path,newline="") as f:
        return list(csv.DictReader(f))

# 用等间距CPU标签绘制一个数据集的一项总吞吐量。
def draw_panel(axis,rows,title,metric):
    cpu_values=sorted({int(row["CPUs"]) for row in rows})
    positions={cpus:index for index,cpus in enumerate(cpu_values)}
    for method in methods:
        method_rows=sorted((row for row in rows if row["Method"]==method),key=lambda row:int(row["CPUs"]))
        if len(method_rows)<2:
            continue
        axis.plot([positions[int(row["CPUs"])] for row in method_rows],[float(row[metric]) for row in method_rows],color=colors[method],marker=markers[method],markersize=5 if method=="MRCZ" else 3.5,linewidth=2.8 if method=="MRCZ" else 1.6,label=method)
    axis.set_title(title,fontsize=11)
    axis.set_xticks(range(len(cpu_values)),cpu_values)
    axis.set_xlabel("CPU count")
    axis.set_ylabel("Aggregate throughput (MB/s)")
    axis.grid(True,alpha=0.25,linewidth=0.5)
    axis.tick_params(labelsize=8)

# 生成pb0和jicf0压缩、解压总吞吐量的四面板图。
def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","weak_scaling")
    datasets={"pb0":read_rows(os.path.join(result_dir,"weak_pb0_summary.csv")),"jicf0":read_rows(os.path.join(result_dir,"weak_jicf_summary.csv"))}
    figure,axes=plt.subplots(2,2,figsize=(13,9))
    draw_panel(axes[0,0],datasets["pb0"],"pb0 compression","Compression Throughput MB/s")
    draw_panel(axes[0,1],datasets["jicf0"],"jicf0 compression","Compression Throughput MB/s")
    draw_panel(axes[1,0],datasets["pb0"],"pb0 decompression","Decompression Throughput MB/s")
    draw_panel(axes[1,1],datasets["jicf0"],"jicf0 decompression","Decompression Throughput MB/s")
    handles=[plt.Line2D([0],[0],color=colors[method],marker=markers[method],linewidth=2.8 if method=="MRCZ" else 1.6,markersize=5,label=method) for method in methods]
    figure.legend(handles=handles,loc="upper center",ncol=len(methods),frameon=False,bbox_to_anchor=(0.5,0.995),fontsize=9)
    figure.suptitle("Weak scaling aggregate throughput",fontsize=16,y=0.95)
    figure.tight_layout(rect=(0,0,1,0.91))
    output_path=os.path.join(result_dir,"weak_scaling_throughput.svg")
    figure.savefig(output_path,format="svg",bbox_inches="tight")
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:
        figure.savefig(os.path.join(preview_dir,"weak_scaling_throughput.png"),format="png",dpi=150,bbox_inches="tight")
    plt.close(figure)
    print(output_path)

if __name__=="__main__":
    main()
