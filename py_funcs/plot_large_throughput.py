import csv,math,os
import matplotlib.pyplot as plt

methods=["MRCZ","FDWZ","SZ3-4D","SZ2-3D-fused","QoZ-4D","ZFP-3D-fused","SPERR-2D"]
colors={"MRCZ":"#111827","FDWZ":"#e15759","SZ3-4D":"#4e79a7","SZ2-3D-fused":"#76b7b2","QoZ-4D":"#59a14f","ZFP-3D-fused":"#f28e2b","SPERR-2D":"#b07aa1"}
markers={"MRCZ":"*","FDWZ":"o","SZ3-4D":"^","SZ2-3D-fused":"s","QoZ-4D":"D","ZFP-3D-fused":"v","SPERR-2D":"P"}
fields=[
    ("pb0 / UX","pb0/pb0.f00001_UX.csv",64),("pb0 / UY","pb0/pb0.f00001_UY.csv",64),("pb0 / UZ","pb0/pb0.f00001_UZ.csv",64),("pb0 / P","pb0/pb0.f00001_P.csv",64),("pb0 / T","pb0/pb0.f00001_T.csv",64),
    ("jicf0 / UX","jicf0/jicf0.f00000_UX.csv",128),("jicf0 / UY","jicf0/jicf0.f00000_UY.csv",128),("jicf0 / UZ","jicf0/jicf0.f00000_UZ.csv",128),("jicf0 / P","jicf0/jicf0.f00000_P.csv",128),("jicf0 / S00","jicf0/jicf0.f00000_S00.csv",128)]

# 读取一个CSV中的全部记录。
def read_rows(path):
    with open(path,newline="") as f:
        return list(csv.DictReader(f))

# 将缺失值和非有限值过滤掉。
def finite_pair(first,second):
    try:
        pair=float(first),float(second)
    except (TypeError,ValueError):
        return None
    return pair if all(math.isfinite(value) and value>0 for value in pair) else None

# 绘制pb0和jicf0全部field在不同PSNR下的压缩或解压吞吐量。
def plot_throughput(result_dir,task):
    figure,axes=plt.subplots(2,5,figsize=(20,8))
    metric="Compression Throughput" if task=="compression" else "Decompression Throughput"
    for axis,(title,relative_path,max_cpus) in zip(axes.flat,fields):
        rd_rows=read_rows(os.path.join(result_dir,relative_path))
        throughput_path=os.path.join(result_dir,relative_path.removesuffix(".csv")+"_throughput.csv")
        throughput_rows={float(row["REL EB"]):row for row in read_rows(throughput_path) if int(row["Max CPUs"])==max_cpus}
        for method in methods:
            points=[]
            for row in rd_rows:
                throughput_row=throughput_rows.get(float(row["REL EB"]))
                if throughput_row is None:
                    continue
                point=finite_pair(row.get(f"{method}_PSNR"),throughput_row.get(f"{method} {metric}"))
                if point is not None:
                    points.append(point)
            points.sort()
            axis.plot([point[0] for point in points],[point[1] for point in points],color=colors[method],marker=markers[method],markersize=4 if method=="MRCZ" else 3,linewidth=2.8 if method=="MRCZ" else 1.5,label=method)
        axis.set_title(title,fontsize=10)
        axis.set_yscale("log")
        axis.grid(True,which="both",alpha=0.25,linewidth=0.5)
        axis.tick_params(labelsize=8)
    handles=[plt.Line2D([0],[0],color=colors[method],marker=markers[method],linewidth=2.8 if method=="MRCZ" else 1.5,markersize=5,label=method) for method in methods]
    figure.legend(handles=handles,loc="upper center",ncol=len(methods),frameon=False,bbox_to_anchor=(0.5,0.995),fontsize=9)
    figure.supxlabel("PSNR (dB)",fontsize=13)
    figure.supylabel(f"{task.capitalize()} throughput (MB/s per CPU, log scale)",fontsize=13)
    figure.tight_layout(rect=(0.025,0.035,1,0.93))
    output_path=os.path.join(result_dir,f"large_{task}_throughput.svg")
    figure.savefig(output_path,format="svg",bbox_inches="tight")
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:
        figure.savefig(os.path.join(preview_dir,f"large_{task}_throughput.png"),format="png",dpi=150,bbox_inches="tight")
    plt.close(figure)
    print(output_path)

# 生成压缩和解压两张大数据throughput图。
def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","RD")
    plot_throughput(result_dir,"compression")
    plot_throughput(result_dir,"decompression")

if __name__=="__main__":
    main()
