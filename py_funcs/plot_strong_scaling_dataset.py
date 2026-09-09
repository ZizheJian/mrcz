import argparse,csv,os
import matplotlib.pyplot as plt

methods=["MRCZ","FDWZ","SZ3-4D","QoZ-4D","ZFP-3D-fused","SPERR-2D"]
colors={"MRCZ":"#111827","FDWZ":"#e15759","SZ3-4D":"#4e79a7","QoZ-4D":"#59a14f","ZFP-3D-fused":"#f28e2b","SPERR-2D":"#b07aa1"}
markers={"MRCZ":"*","FDWZ":"o","SZ3-4D":"^","QoZ-4D":"D","ZFP-3D-fused":"v","SPERR-2D":"P"}
dataset_titles={"pb0":"pb0 (100 dB)"}

def read_rows(path):
    with open(path,newline="") as file:
        return list(csv.DictReader(file))

def draw_panel(axis,rows,task):
    column=f"{task} Throughput MB/s"
    cpu_values=sorted({int(row["CPUs"]) for row in rows})
    base_cpus=min(cpu_values)
    positions={cpus:index for index,cpus in enumerate(cpu_values)}
    axis.plot(range(len(cpu_values)),[cpus/base_cpus for cpus in cpu_values],color="#9ca3af",linestyle="--",linewidth=1.2,label="Ideal")
    for method in methods:
        method_rows=sorted((row for row in rows if row["Method"]==method),key=lambda row:int(row["CPUs"]))
        if len(method_rows)<2:
            continue
        baseline=float(method_rows[0][column])
        x=[positions[int(row["CPUs"])] for row in method_rows]
        y=[float(row[column])/baseline for row in method_rows]
        axis.plot(x,y,color=colors[method],marker=markers[method],markersize=7 if method=="MRCZ" else 4.5,linewidth=2.8 if method=="MRCZ" else 1.6,label=method)
    axis.set_xticks(range(len(cpu_values)),cpu_values)
    axis.set_xlabel("CPU count")
    axis.set_ylabel(f"Speedup relative to {base_cpus} CPUs")
    axis.set_title(task)
    axis.grid(True,alpha=0.25,linewidth=0.5)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("dataset",choices=sorted(dataset_titles))
    args=parser.parse_args()
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","strong_scaling")
    rows=read_rows(os.path.join(result_dir,f"strong_{args.dataset}_summary.csv"))
    figure,axes=plt.subplots(1,2,figsize=(13,5.2))
    draw_panel(axes[0],rows,"Compression")
    draw_panel(axes[1],rows,"Decompression")
    handles=[plt.Line2D([0],[0],color=colors[method],marker=markers[method],linewidth=2.8 if method=="MRCZ" else 1.6,markersize=7 if method=="MRCZ" else 4.5,label=method) for method in methods]
    handles.append(plt.Line2D([0],[0],color="#9ca3af",linestyle="--",linewidth=1.2,label="Ideal"))
    figure.legend(handles=handles,loc="upper center",ncol=len(handles),frameon=False,bbox_to_anchor=(0.5,0.99),fontsize=9)
    figure.suptitle(f"{dataset_titles[args.dataset]} strong scaling",fontsize=16,y=0.92)
    figure.tight_layout(rect=(0,0,1,0.84))
    output_path=os.path.join(result_dir,f"strong_{args.dataset}_speedup.svg")
    figure.savefig(output_path,format="svg",bbox_inches="tight")
    with open(output_path) as input_file:
        svg=input_file.read()
    with open(output_path,"w") as output_file:
        output_file.write("\n".join(line.rstrip() for line in svg.splitlines())+"\n")
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:
        figure.savefig(os.path.join(preview_dir,f"strong_{args.dataset}_speedup.png"),format="png",dpi=160,bbox_inches="tight")
    plt.close(figure)
    print(output_path)

if __name__=="__main__":
    main()
