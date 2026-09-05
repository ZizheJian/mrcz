import argparse,csv,math,os
import matplotlib.pyplot as plt

plt.rcParams["svg.hashsalt"]="mrcz-large-rd"

methods=["MRCZ","FDWZ","SZ3-4D","SZ2-3D-fused","QoZ-4D","ZFP-3D-fused","SPERR-2D"]
colors={"MRCZ":"#111827","FDWZ":"#e15759","SZ3-4D":"#4e79a7","SZ2-3D-fused":"#76b7b2","QoZ-4D":"#59a14f","ZFP-3D-fused":"#f28e2b","SPERR-2D":"#b07aa1"}
markers={"MRCZ":"*","FDWZ":"o","SZ3-4D":"^","SZ2-3D-fused":"s","QoZ-4D":"D","ZFP-3D-fused":"v","SPERR-2D":"P"}
datasets={
    "pb0":("pb0",["UX","UY","UZ","P","T"],"pb0.f00001"),
    "jicf0":("jicf0",["UX","UY","UZ","P","S00"],"jicf0.f00000")}

# 读取一个field的rate-distortion CSV。
def read_rows(path):
    with open(path,newline="") as f:
        return list(csv.DictReader(f))

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

# 将BR左边界固定为0，并聚焦MRCZ、FDWZ及大多数baseline共同覆盖的范围。
def focus_axes(axis,method_points):
    all_points=[point for points in method_points.values() for point in points]
    primary_points=method_points["MRCZ"]+method_points["FDWZ"]
    br_right=max(percentile([point[0] for point in all_points],0.9),max(point[0] for point in primary_points))*1.05
    psnr_values=[point[1] for point in all_points]
    psnr_top=max(percentile(psnr_values,0.95),max(point[1] for point in primary_points))
    psnr_bottom=min(psnr_values)
    padding=max(1,(psnr_top-psnr_bottom)*0.03)
    axis.set_xlim(0,br_right)
    axis.set_ylim(psnr_bottom-padding,psnr_top+padding)

# 为一个大数据集的五个field绘制rate-distortion曲线。
def plot_dataset(result_dir,dataset):
    directory,fields,prefix=datasets[dataset]
    figure,axes=plt.subplots(1,len(fields),figsize=(19,4.5))
    for axis,field in zip(axes,fields):
        rows=read_rows(os.path.join(result_dir,directory,f"{prefix}_{field}.csv"))
        method_points={}
        for method in methods:
            points=sorted(point for row in rows if (point:=finite_pair(row.get(f"{method}_BR"),row.get(f"{method}_PSNR"))) is not None)
            method_points[method]=points
            axis.plot([point[0] for point in points],[point[1] for point in points],color=colors[method],marker=markers[method],markersize=7 if method=="MRCZ" else 4.5,linewidth=2.8 if method=="MRCZ" else 1.6,label=method)
        focus_axes(axis,method_points)
        axis.set_title(field)
        axis.grid(True,alpha=0.25,linewidth=0.5)
    handles=[plt.Line2D([0],[0],color=colors[method],marker=markers[method],linewidth=2.8 if method=="MRCZ" else 1.6,markersize=7 if method=="MRCZ" else 4.5,label=method) for method in methods]
    figure.legend(handles=handles,loc="upper center",ncol=len(methods),frameon=False,bbox_to_anchor=(0.5,0.99),fontsize=9)
    figure.suptitle(f"{dataset} rate-distortion",fontsize=16,y=0.9)
    figure.supxlabel("Bit rate (bits/value)")
    figure.supylabel("PSNR (dB)")
    figure.tight_layout(rect=(0.015,0.04,1,0.82))
    output_path=os.path.join(result_dir,f"{dataset}_rd.svg")
    figure.savefig(output_path,format="svg",bbox_inches="tight",metadata={"Date":None})
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:
        figure.savefig(os.path.join(preview_dir,f"{dataset}_rd.png"),format="png",dpi=160,bbox_inches="tight")
    plt.close(figure)
    print(output_path)

# 选择并绘制一个大数据集。
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("dataset",choices=sorted(datasets))
    args=parser.parse_args()
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    plot_dataset(os.path.join(project_root,"results","RD"),args.dataset)

if __name__=="__main__":
    main()
