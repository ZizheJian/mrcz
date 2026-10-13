import glob,math,os

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator,NullLocator

import plot_aggregate_throughput as large
import plot_small_core_throughput as small
from plot_style import METHOD_COLORS,TIMING_METHODS,apply_paper_style,style_plot_axis

apply_paper_style("mrcz-combined-throughput",tick_size=7.0)

methods=TIMING_METHODS
colors=METHOD_COLORS
tasks=("Compression","Decompression")


def group_limit(aggregate,extra_aggregate):
    maximum_value=max(value for task in tasks for values in aggregate[task].values() for value in values)
    maximum_value=max(maximum_value,max(value for task in tasks for value in extra_aggregate[task]))
    return math.ceil(maximum_value*1.05*2)/2


def mgard_8core(project_root,dataset,error_bounds):
    directory=os.path.join(project_root,"results","small_core_throughput",dataset)
    tables=[small.read_rows(path) for path in sorted(glob.glob(os.path.join(directory,"*_throughput.csv")))]
    if len(tables)!=5:raise RuntimeError(f"{dataset} requires five MGARD throughput files,found {len(tables)}")
    aggregate={task:[] for task in tasks}
    for error_bound in error_bounds:
        for task in tasks:
            field_values=[]
            for table in tables:
                row=next((row for row in table if int(row["Max CPUs"])==8 and math.isclose(float(row["REL EB"]),error_bound,rel_tol=1e-12)),None)
                if row is None or not row.get(f"MGARD-1D {task} Throughput"):raise RuntimeError(f"missing {dataset}/MGARD-1D/{task}/{error_bound:g}")
                field_values.append(float(row[f"MGARD-1D {task} Throughput"]))
            aggregate[task].append(small.harmonic_mean(field_values)*8/1024)
    return aggregate


def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    large_result_dir=os.path.join(project_root,"results","RD")
    data={}
    limits={}
    for dataset,dataset_index,small_cpus,partition_num,error_bounds in small.datasets:
        small_aggregate=small.aggregate_dataset(project_root,dataset,small_cpus,partition_num,error_bounds)
        large_cpus=partition_num
        large_error_bounds,large_aggregate=large.aggregate_dataset(large_result_dir,dataset,large_cpus)
        mgard_aggregate=mgard_8core(project_root,dataset,error_bounds)
        data[(0,dataset)]=(small_cpus,error_bounds,small_aggregate,mgard_aggregate)
        data[(1,dataset)]=(large_cpus,large_error_bounds,large_aggregate,mgard_aggregate)
        limits[(0,dataset)]=group_limit(small_aggregate,mgard_aggregate)
        projected_mgard={task:[value*large_cpus/8 for value in mgard_aggregate[task]] for task in tasks}
        limits[(1,dataset)]=group_limit(large_aggregate,projected_mgard)

    columns=[("pb0",7,"Compression"),("pb0",7,"Decompression"),("jicf0",8,"Compression"),("jicf0",8,"Decompression")]
    figure,axes=plt.subplots(2,4,figsize=(7.16,2.99))
    for row in range(2):
        for column,(dataset,dataset_index,task) in enumerate(columns):
            axis=axes[row,column]
            cpus,error_bounds,aggregate,mgard_aggregate=data[(row,dataset)]
            for method in methods:axis.plot(error_bounds,aggregate[task][method],color=colors[method],linewidth=0.8,label=method)
            mgard_values=[value*cpus/8 for value in mgard_aggregate[task]]
            axis.plot(error_bounds,mgard_values,color=colors["MGARD-1D"],linestyle="-" if cpus==8 else "--",linewidth=1.0,zorder=2)
            axis.set_xscale("log")
            axis.invert_xaxis()
            axis.xaxis.set_minor_locator(NullLocator())
            axis.set_ylim(0,limits[(row,dataset)])
            axis.yaxis.set_major_locator(MaxNLocator(nbins=4,integer=True))
            axis.grid(True,which="major",linestyle="--",linewidth=0.35,alpha=0.35)
            style_plot_axis(axis)
            panel=chr(ord("a")+row*4+column)
            axis.set_xlabel(f"({panel}) Dataset {dataset_index}, {task}\n{cpus} cores",labelpad=1.0)

    handles=[plt.Line2D([0],[0],color=colors[method],linewidth=0.8,label=method) for method in methods]
    handles.extend([plt.Line2D([0],[0],color=colors["MGARD-1D"],linewidth=1.0,label="MGARD-X-1D measured"),plt.Line2D([0],[0],color=colors["MGARD-1D"],linestyle="--",linewidth=1.0,label="MGARD-X-1D ideal")])
    figure.legend(handles=handles,loc="upper center",ncol=4,frameon=False,bbox_to_anchor=(0.5,0.995),handlelength=1.0,handletextpad=0.25,columnspacing=0.75,borderaxespad=0,fontsize=8.0)
    figure.supxlabel("Relative Error Bound",x=0.515,y=0.018)
    figure.supylabel("Aggregate Throughput (GiB/s)",x=0.007,y=0.47)
    figure.tight_layout(rect=(0.025,0.045,1,0.86),pad=0.25,w_pad=0.45,h_pad=0.75)

    result_dir=os.path.join(project_root,"results","RD")
    svg_output=os.path.join(result_dir,"combined_throughput.svg")
    png_output=os.path.join(result_dir,"combined_throughput.png")
    pdf_output=os.path.join(project_root,"论文","figs","7","combined-throughput.pdf")
    figure.savefig(svg_output,format="svg",metadata={"Date":None},bbox_inches="tight",pad_inches=0.02)
    figure.savefig(png_output,format="png",dpi=240,bbox_inches="tight",pad_inches=0.02)
    figure.savefig(pdf_output,format="pdf",bbox_inches="tight",pad_inches=0.02)
    plt.close(figure)
    print(pdf_output)
    print(png_output)


if __name__=="__main__":main()
