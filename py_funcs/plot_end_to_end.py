import csv,os
import matplotlib.pyplot as plt

plt.rcParams["svg.hashsalt"]="mrcz-end-to-end"

colors={"MRCZ":"#111827","FDWZ":"#e15759","SZ3-4D":"#4e79a7","SZ2-3D-fused":"#76b7b2","QoZ-4D":"#59a14f","ZFP-3D-fused":"#f28e2b","SPERR-2D":"#b07aa1"}
markers={"MRCZ":"*","FDWZ":"o","SZ3-4D":"^","SZ2-3D-fused":"s","QoZ-4D":"D","ZFP-3D-fused":"v","SPERR-2D":"P"}
titles={"pb0":"pb0 (REL EB = 1e-4, 64 CPUs)","jicf0":"jicf0 (REL EB = 1e-6, 128 CPUs)"}

def read_rows(path):
    with open(path,newline="") as input_file:
        return list(csv.DictReader(input_file))

def draw_panel(axis,rows,dataset):
    rows=sorted(rows,key=lambda row:float(row["End-to-end Time Median"]),reverse=True)
    for position,row in enumerate(rows):
        method=row["Method"]
        median=float(row["End-to-end Time Median"])
        minimum=float(row["End-to-end Time Min"])
        maximum=float(row["End-to-end Time Max"])
        axis.errorbar(median,position,xerr=[[median-minimum],[maximum-median]],fmt=markers[method],color=colors[method],markersize=10 if method=="MRCZ" else 7,elinewidth=2,capsize=3,zorder=3)
        axis.annotate(f"{median:.1f} s",(median,position),xytext=(7,0),textcoords="offset points",va="center",fontsize=9)
    axis.set_yticks(range(len(rows)),[row["Method"] for row in rows])
    axis.set_xscale("log")
    axis.set_xlabel("End-to-end time (s, log scale)")
    axis.set_title(titles[dataset])
    axis.grid(True,axis="x",which="both",alpha=0.25,linewidth=0.6)
    axis.set_axisbelow(True)
    axis.set_xlim(min(float(row["End-to-end Time Min"]) for row in rows)*0.78,max(float(row["End-to-end Time Max"]) for row in rows)*1.75)

def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","end_to_end")
    rows=read_rows(os.path.join(result_dir,"end_to_end_summary.csv"))
    figure,axes=plt.subplots(1,2,figsize=(14,5.4))
    for axis,dataset in zip(axes,["pb0","jicf0"]):
        draw_panel(axis,[row for row in rows if row["Dataset"]==dataset],dataset)
    figure.suptitle("Compression + Globus transfer + decompression",fontsize=16,y=0.98)
    figure.text(0.5,0.015,"Points use the median of three transfers; error bars show the transfer-derived min--max range. Disk I/O is excluded.",ha="center",fontsize=9)
    figure.tight_layout(rect=(0,0.05,1,0.93))
    output_path=os.path.join(result_dir,"end_to_end_total.svg")
    figure.savefig(output_path,format="svg",bbox_inches="tight",metadata={"Date":None})
    with open(output_path) as input_file:
        svg=input_file.read()
    with open(output_path,"w") as output_file:
        output_file.write("\n".join(line.rstrip() for line in svg.splitlines())+"\n")
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:
        figure.savefig(os.path.join(preview_dir,"end_to_end_total.png"),format="png",dpi=160,bbox_inches="tight")
    plt.close(figure)
    print(output_path)

if __name__=="__main__":
    main()
