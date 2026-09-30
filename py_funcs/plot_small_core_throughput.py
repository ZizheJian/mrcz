import csv,glob,math,os
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator,NullLocator

plt.rcParams["svg.hashsalt"]="mrcz-small-core-throughput"
plt.rcParams.update({"font.family":"serif","font.serif":["Nimbus Roman"],"mathtext.fontset":"stix","font.size":8.0,"axes.labelsize":8.0,"axes.titlesize":8.0,"xtick.labelsize":8.0,"ytick.labelsize":8.0,"legend.fontsize":8.0,"pdf.fonttype":42,"ps.fonttype":42})

methods=["MRCZ","FDWZ","SZ3-4D","QoZ-4D","ZFP-3D-fused","SPERR-2D"]
colors={"MRCZ":"#111827","FDWZ":"#e15759","SZ3-4D":"#4e79a7","QoZ-4D":"#59a14f","ZFP-3D-fused":"#f28e2b","SPERR-2D":"#b07aa1"}
datasets=[("pb0",7,8,64,[1e-2,1e-3,1e-4,1e-5,1e-6]),("jicf0",8,16,128,[1e-4,1e-5,1e-6,1e-7,1e-8])]
split_names={"FDWZ":"fdwz","SZ3-4D":"sz3","QoZ-4D":"qoz","ZFP-3D-fused":"zfp","SPERR-2D":"sperr"}

def read_rows(path):
    with open(path,newline="") as input_file:return list(csv.DictReader(input_file))

def harmonic_mean(values):
    return len(values)/sum(1/value for value in values)

def supplemental_dir(project_root,dataset,method):
    root=os.path.join(project_root,"results","small_core_throughput")
    if dataset=="pb0":return os.path.join(root,"pb0")
    if method=="MRCZ":return os.path.join(root,"jicf0")
    return os.path.join(root,"jicf0_split",split_names[method])

def aggregate_dataset(project_root,dataset,cpus,partition_num,error_bounds):
    middle_dir=os.path.join(project_root,"results","strong_scaling","fixed_rel_eb",dataset)
    field_files=[os.path.basename(path) for path in sorted(glob.glob(os.path.join(middle_dir,"*_throughput.csv")))]
    if len(field_files)!=5:raise RuntimeError(f"{dataset} requires five field throughput files,found {len(field_files)}")
    aggregate={task:{method:[] for method in methods} for task in ("Compression","Decompression")}
    for method in methods:
        source_dir=supplemental_dir(project_root,dataset,method)
        tables=[]
        for field_file in field_files:
            supplement_path=os.path.join(source_dir,field_file)
            middle_path=os.path.join(middle_dir,field_file)
            if not os.path.isfile(supplement_path):raise RuntimeError(f"missing {supplement_path}")
            tables.append(read_rows(supplement_path)+read_rows(middle_path))
        for error_bound in error_bounds:
            for task in aggregate:
                field_values=[]
                for table in tables:
                    column=f"{method} {task} Throughput"
                    row=next((row for row in table if int(row["Max CPUs"])==cpus and int(row["Partition Num"])==partition_num and math.isclose(float(row["REL EB"]),error_bound,rel_tol=1e-12) and row.get(column)),None)
                    if row is None:raise RuntimeError(f"missing {dataset}/{method}/{task}/{error_bound:g}")
                    field_values.append(float(row[f"{method} {task} Throughput"]))
                aggregate[task][method].append(harmonic_mean(field_values)*cpus/1024)
    return aggregate

def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    aggregates={}
    row_limits=[]
    for dataset,_,cpus,partition_num,error_bounds in datasets:
        aggregate=aggregate_dataset(project_root,dataset,cpus,partition_num,error_bounds)
        aggregates[dataset]=aggregate
        maximum=max(value for task in aggregate.values() for method_values in task.values() for value in method_values)
        row_limits.append(math.ceil(maximum)+0.5)

    figure,axes=plt.subplots(2,2,figsize=(3.5,2.8))
    panel_index=0
    for row,(dataset,dataset_index,cpus,_,error_bounds) in enumerate(datasets):
        aggregate=aggregates[dataset]
        for column,task in enumerate(("Compression","Decompression")):
            axis=axes[row,column]
            for method in methods:axis.plot(error_bounds,aggregate[task][method],color=colors[method],linewidth=0.8,label=method)
            axis.set_xscale("log")
            axis.invert_xaxis()
            axis.xaxis.set_minor_locator(NullLocator())
            axis.set_ylim(0,row_limits[row])
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
    result_dir=os.path.join(project_root,"results","small_core_throughput")
    svg_output=os.path.join(result_dir,"small_core_throughput.svg")
    png_output=os.path.join(result_dir,"small_core_throughput.png")
    pdf_output=os.path.join(project_root,"论文","figs","7","small-core-throughput.pdf")
    figure.savefig(svg_output,format="svg",metadata={"Date":None},bbox_inches="tight",pad_inches=0.02)
    figure.savefig(png_output,format="png",dpi=240,bbox_inches="tight",pad_inches=0.02)
    figure.savefig(pdf_output,format="pdf",bbox_inches="tight",pad_inches=0.02)
    plt.close(figure)
    for dataset,dataset_index,cpus,_,_ in datasets:
        aggregate=aggregates[dataset]
        for task in ("Compression","Decompression"):
            wins=sum(aggregate[task]["FDWZ"][index]==max(aggregate[task][method][index] for method in methods) for index in range(5))
            print(f"Dataset {dataset_index}, {cpus} cores, {task}: FDWZ wins {wins}/5")
    print(pdf_output)
    print(png_output)

if __name__=="__main__":main()
