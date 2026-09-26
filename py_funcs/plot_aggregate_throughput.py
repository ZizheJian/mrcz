import csv,glob,math,os
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator,NullLocator

plt.rcParams["svg.hashsalt"]="mrcz-aggregate-throughput"
plt.rcParams.update({"font.family":"serif","font.serif":["Nimbus Roman"],"mathtext.fontset":"stix","font.size":8.0,"axes.labelsize":8.0,"axes.titlesize":8.0,"xtick.labelsize":8.0,"ytick.labelsize":8.0,"legend.fontsize":8.0,"pdf.fonttype":42,"ps.fonttype":42})

methods=["MRCZ","FDWZ","SZ3-4D","QoZ-4D","ZFP-3D-fused","SPERR-2D"]
colors={"MRCZ":"#111827","FDWZ":"#e15759","SZ3-4D":"#4e79a7","QoZ-4D":"#59a14f","ZFP-3D-fused":"#f28e2b","SPERR-2D":"#b07aa1"}
datasets=[("pb0",7,64),("jicf0",8,128)]

def read_rows(path):
    with open(path,newline="") as input_file:return list(csv.DictReader(input_file))

def harmonic_mean(values):
    return len(values)/sum(1/value for value in values)

def aggregate_dataset(result_dir,dataset,cpus):
    tables=[read_rows(path) for path in sorted(glob.glob(os.path.join(result_dir,dataset,"*_throughput.csv")))]
    if len(tables)!=5:raise RuntimeError(f"{dataset} requires five field throughput files,found {len(tables)}")
    error_bounds=sorted({float(row["REL EB"]) for table in tables for row in table if int(row["Max CPUs"])==cpus},reverse=True)
    aggregate={task:{method:[] for method in methods} for task in ("Compression","Decompression")}
    for error_bound in error_bounds:
        for task in aggregate:
            for method in methods:
                field_values=[]
                for table in tables:
                    row=next(row for row in table if int(row["Max CPUs"])==cpus and math.isclose(float(row["REL EB"]),error_bound,rel_tol=1e-12))
                    field_values.append(float(row[f"{method} {task} Throughput"]))
                aggregate[task][method].append(harmonic_mean(field_values)*cpus/1024)
    return error_bounds,aggregate

def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","RD")
    figure,axes=plt.subplots(2,2,figsize=(3.5,2.8))
    panel_index=0
    for row,(dataset,dataset_index,cpus) in enumerate(datasets):
        error_bounds,aggregate=aggregate_dataset(result_dir,dataset,cpus)
        for column,task in enumerate(("Compression","Decompression")):
            axis=axes[row,column]
            for method in methods:
                axis.plot(error_bounds,aggregate[task][method],color=colors[method],linewidth=0.8,label=method)
            axis.set_xscale("log")
            axis.invert_xaxis()
            axis.xaxis.set_minor_locator(NullLocator())
            axis.set_ylim(0,18 if row==0 else 12)
            axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
            axis.grid(True,which="major",alpha=0.25,linewidth=0.35)
            axis.tick_params(direction="in",length=1.8,width=0.45,pad=1)
            for spine in axis.spines.values():spine.set_linewidth(0.5)
            panel=chr(ord("a")+panel_index)
            axis.set_xlabel(f"({panel}) Dataset {dataset_index}, {task}",labelpad=1.0)
            panel_index+=1
    handles=[plt.Line2D([0],[0],color=colors[method],linewidth=0.8,label=method) for method in methods]
    figure.legend(handles=handles,loc="upper center",ncol=3,frameon=False,bbox_to_anchor=(0.5,0.975),handlelength=0.9,handletextpad=0.2,columnspacing=0.55,borderaxespad=0,fontsize=8.0)
    figure.supxlabel("Relative Error Bound",x=0.54,y=0.028)
    figure.supylabel("Aggregate Throughput (GiB/s)",x=0.02,y=0.48)
    figure.tight_layout(rect=(0.035,0.04,1,0.85),pad=0.2,w_pad=0.35,h_pad=0.55)
    svg_output=os.path.join(result_dir,"aggregate_throughput.svg")
    pdf_output=os.path.join(project_root,"论文","figs","7","aggregate-throughput.pdf")
    figure.savefig(svg_output,format="svg",metadata={"Date":None},bbox_inches="tight",pad_inches=0.02)
    figure.savefig(pdf_output,format="pdf",bbox_inches="tight",pad_inches=0.02)
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:figure.savefig(os.path.join(preview_dir,"aggregate_throughput.png"),format="png",dpi=180)
    plt.close(figure)
    print(pdf_output)

if __name__=="__main__":main()
