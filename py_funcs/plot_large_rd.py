import argparse,csv,math,os
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

plt.rcParams["svg.hashsalt"]="mrcz-large-rd"

methods=["MRCZ","FDWZ","SZ3-4D","SZ2-3D-fused","QoZ-4D","ZFP-3D-fused","SPERR-2D"]
colors={"MRCZ":"#111827","FDWZ":"#e15759","SZ3-4D":"#4e79a7","SZ2-3D-fused":"#76b7b2","QoZ-4D":"#59a14f","ZFP-3D-fused":"#f28e2b","SPERR-2D":"#b07aa1"}
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
            axis.plot([point[0] for point in points],[point[1] for point in points],color=colors[method],linewidth=1.6,label=method)
        focus_axes(axis,method_points)
        axis.set_title(field)
        axis.grid(True,alpha=0.25,linewidth=0.5)
    handles=[plt.Line2D([0],[0],color=colors[method],linewidth=1.6,label=method) for method in methods]
    figure.legend(handles=handles,loc="upper center",ncol=len(methods),frameon=False,bbox_to_anchor=(0.5,0.99),fontsize=9)
    figure.suptitle(f"{dataset} rate-distortion",fontsize=16,y=0.9)
    figure.supxlabel("Bit rate")
    figure.supylabel("PSNR (dB)")
    figure.tight_layout(rect=(0.015,0.04,1,0.82))
    output_path=os.path.join(result_dir,f"{dataset}_rd.svg")
    figure.savefig(output_path,format="svg",bbox_inches="tight",metadata={"Date":None})
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:
        figure.savefig(os.path.join(preview_dir,f"{dataset}_rd.png"),format="png",dpi=160,bbox_inches="tight")
    plt.close(figure)
    print(output_path)

# 将两个大数据集绘制为论文中的2x5双栏图。
def plot_combined(result_dir,project_root):
    plt.rcParams.update({"font.family":"serif","font.serif":["Nimbus Roman"],"mathtext.fontset":"stix","font.size":8.0,"axes.labelsize":8.0,"axes.titlesize":8.0,"xtick.labelsize":8.0,"ytick.labelsize":8.0,"legend.fontsize":8.0,"pdf.fonttype":42,"ps.fonttype":42})
    figure,axes=plt.subplots(2,5,figsize=(7.16*5/6,2.56))
    panels=[]
    for dataset_index,dataset in ((7,"pb0"),(8,"jicf0")):
        directory,fields,prefix=datasets[dataset]
        panels.extend((f"{dataset_index}{field}",os.path.join(result_dir,directory,f"{prefix}_{field}.csv")) for field in fields)
    for panel_index,(axis,(title,path)) in enumerate(zip(axes.flat,panels)):
        rows=read_rows(path)
        method_points={}
        for method in methods:
            points=sorted(point for row in rows if (point:=finite_pair(row.get(f"{method}_BR"),row.get(f"{method}_PSNR"))) is not None)
            method_points[method]=points
            axis.plot([point[0] for point in points],[point[1] for point in points],color=colors[method],linewidth=0.8,label=method,zorder=1)
        focus_axes(axis,method_points)
        axis.set_xlabel(f"({chr(ord('a')+panel_index)}) {title}",labelpad=1.5)
        axis.grid(True,alpha=0.25,linewidth=0.35)
        axis.xaxis.set_major_locator(MaxNLocator(nbins=3))
        axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
        axis.tick_params(direction="in",length=1.8,width=0.45,pad=1)
        for spine in axis.spines.values():spine.set_linewidth(0.5)
    handles=[plt.Line2D([0],[0],color=colors[method],linewidth=0.8,label=method) for method in methods]
    figure.legend(handles=handles,loc="upper center",ncol=len(methods),frameon=False,bbox_to_anchor=(0.5,0.995),handlelength=1.0,handletextpad=0.25,columnspacing=0.65,borderaxespad=0)
    figure.supxlabel("Bit rate",x=0.53,y=0.025)
    figure.supylabel("PSNR (dB)",x=0.008,y=0.47)
    figure.tight_layout(rect=(0.015,0.05,1,0.92),pad=0.25,w_pad=0.3,h_pad=0.65)
    figure.savefig(os.path.join(result_dir,"large_rd.svg"),format="svg",metadata={"Date":None},bbox_inches="tight",pad_inches=0.02)
    paper_output=os.path.join(project_root,"论文","figs","7","large-rd.pdf")
    os.makedirs(os.path.dirname(paper_output),exist_ok=True)
    for axis in figure.axes:axis.set_rasterization_zorder(1.5)
    figure.savefig(paper_output,format="pdf",dpi=300,bbox_inches="tight",pad_inches=0.02)
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:figure.savefig(os.path.join(preview_dir,"large_rd.png"),format="png",dpi=180)
    plt.close(figure)
    print(paper_output)

# 选择并绘制一个大数据集。
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("dataset",nargs="?",choices=["combined"]+sorted(datasets),default="combined")
    args=parser.parse_args()
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","RD")
    if args.dataset=="combined":plot_combined(result_dir,project_root)
    else:plot_dataset(result_dir,args.dataset)

if __name__=="__main__":
    main()
