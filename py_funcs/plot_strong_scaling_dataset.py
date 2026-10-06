import argparse,csv,glob,json,math,os
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from plot_style import METHOD_COLORS,TIMING_METHODS,apply_paper_style,style_plot_axis

apply_paper_style("mrcz-strong-scaling")

methods=TIMING_METHODS
colors=METHOD_COLORS
datasets=[("pb0",7,1e-4,[8,16,32,64]),("jicf0",8,1e-6,[16,32,64,128])]

def read_rows(path):
    with open(path,newline="") as input_file:return list(csv.DictReader(input_file))

def harmonic_mean(values):
    return len(values)/sum(1/value for value in values)

def aggregate_speedups(result_dir,dataset,error_bound,cpu_values):
    tables=[read_rows(path) for path in sorted(glob.glob(os.path.join(result_dir,dataset,"*_throughput.csv")))]
    if len(tables)!=5:raise RuntimeError(f"{dataset} requires five field throughput files,found {len(tables)}")
    throughputs={task:{method:[] for method in methods} for task in ("Compression","Decompression")}
    for cpus in cpu_values:
        for task in throughputs:
            for method in methods:
                field_values=[]
                for table in tables:
                    row=next((row for row in table if int(row["Max CPUs"])==cpus and int(row["Partition Num"])==max(cpu_values) and math.isclose(float(row["REL EB"]),error_bound,rel_tol=1e-12)),None)
                    if row is None:raise RuntimeError(f"missing {dataset}/{method}/{task}/{cpus}-CPU result")
                    field_values.append(float(row[f"{method} {task} Throughput"]))
                throughputs[task][method].append(harmonic_mean(field_values)*cpus)
    speedups={task:{method:[value/throughputs[task][method][0] for value in throughputs[task][method]] for method in methods} for task in throughputs}
    return speedups

def read_placeholder(path):
    with open(path,encoding="utf-8") as input_file:return json.load(input_file)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--placeholder-json")
    parser.add_argument("--output")
    args=parser.parse_args()
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","strong_scaling","fixed_rel_eb")
    os.makedirs(result_dir,exist_ok=True)
    placeholder=read_placeholder(args.placeholder_json) if args.placeholder_json else None
    figure,axes=plt.subplots(2,2,figsize=(3.5,2.8),sharey=True)
    panel_index=0
    for row,(dataset,dataset_index,error_bound,cpu_values) in enumerate(datasets):
        speedups=placeholder[dataset] if placeholder else aggregate_speedups(result_dir,dataset,error_bound,cpu_values)
        positions=list(range(len(cpu_values)))
        ideal=[cpus/cpu_values[0] for cpus in cpu_values]
        for column,task in enumerate(("Compression","Decompression")):
            axis=axes[row,column]
            axis.plot(positions,ideal,color="#7f7f7f",linestyle="--",linewidth=0.8,label="Ideal")
            for method in methods:axis.plot(positions,speedups[task][method],color=colors[method],linewidth=0.8,label=method)
            axis.set_xticks(positions,cpu_values)
            axis.set_xlim(positions[0],positions[-1])
            axis.set_ylim(0,8.4)
            axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
            axis.grid(True,which="major",alpha=0.25,linewidth=0.35)
            style_plot_axis(axis)
            panel=chr(ord("a")+panel_index)
            axis.set_xlabel(f"({panel}) Dataset {dataset_index}, {task}",labelpad=1.0)
            panel_index+=1
    handles=[plt.Line2D([0],[0],color=colors[method],linewidth=0.8,label=method) for method in methods]
    handles.append(plt.Line2D([0],[0],color="#7f7f7f",linestyle="--",linewidth=0.8,label="Ideal"))
    figure.legend(handles=handles,loc="upper center",ncol=4,frameon=False,bbox_to_anchor=(0.5,0.975),handlelength=0.9,handletextpad=0.2,columnspacing=0.55,borderaxespad=0,fontsize=8.0)
    figure.supxlabel("CPU Cores",x=0.54,y=0.028)
    figure.supylabel("Speedup",x=0.02,y=0.48)
    if placeholder:figure.text(0.5,0.5,"Placeholder Data",ha="center",va="center",fontsize=17,color="#7f7f7f",alpha=0.18,rotation=25)
    figure.tight_layout(rect=(0.035,0.04,1,0.85),pad=0.2,w_pad=0.35,h_pad=0.55)
    svg_output=os.path.join(result_dir,"strong_scaling.svg")
    pdf_output=args.output or os.path.join(project_root,"论文","figs","7","strong-scaling.pdf")
    if placeholder and args.output is None:pdf_output=os.path.join(project_root,"论文","figs","7","strong-scaling-placeholder.pdf")
    os.makedirs(os.path.dirname(os.path.abspath(pdf_output)),exist_ok=True)
    figure.savefig(svg_output,format="svg",metadata={"Date":None},bbox_inches="tight",pad_inches=0.02)
    figure.savefig(pdf_output,format="pdf",bbox_inches="tight",pad_inches=0.02)
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:figure.savefig(os.path.join(preview_dir,"strong_scaling.png"),format="png",dpi=180)
    plt.close(figure)
    print(pdf_output)

if __name__=="__main__":main()
