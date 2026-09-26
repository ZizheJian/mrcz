import math,os

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator,NullLocator

import plot_aggregate_throughput as large
import plot_small_core_throughput as small


plt.rcParams["svg.hashsalt"]="mrcz-combined-throughput"
plt.rcParams.update({"font.family":"serif","font.serif":["Nimbus Roman"],"mathtext.fontset":"stix","font.size":8.0,"axes.labelsize":8.0,"axes.titlesize":8.0,"xtick.labelsize":7.0,"ytick.labelsize":7.0,"legend.fontsize":8.0,"pdf.fonttype":42,"ps.fonttype":42})

methods=large.methods
colors=large.colors
tasks=("Compression","Decompression")


def maximum(aggregate):
    return max(value for task in aggregate.values() for method_values in task.values() for value in method_values)


def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    large_result_dir=os.path.join(project_root,"results","RD")
    data={}
    limits={}
    for dataset,dataset_index,small_cpus,partition_num,error_bounds in small.datasets:
        small_aggregate=small.aggregate_dataset(project_root,dataset,small_cpus,partition_num,error_bounds)
        large_cpus=partition_num
        large_error_bounds,large_aggregate=large.aggregate_dataset(large_result_dir,dataset,large_cpus)
        data[(0,dataset)]=(small_cpus,error_bounds,small_aggregate)
        data[(1,dataset)]=(large_cpus,large_error_bounds,large_aggregate)
        limits[(0,dataset)]=math.ceil(maximum(small_aggregate))+0.5
        limits[(1,dataset)]=18 if dataset_index==7 else 12

    columns=[("pb0",7,"Compression"),("pb0",7,"Decompression"),("jicf0",8,"Compression"),("jicf0",8,"Decompression")]
    figure,axes=plt.subplots(2,4,figsize=(7.16,3.25))
    for row in range(2):
        for column,(dataset,dataset_index,task) in enumerate(columns):
            axis=axes[row,column]
            cpus,error_bounds,aggregate=data[(row,dataset)]
            for method in methods:axis.plot(error_bounds,aggregate[task][method],color=colors[method],linewidth=0.8,label=method)
            axis.set_xscale("log")
            axis.invert_xaxis()
            axis.xaxis.set_minor_locator(NullLocator())
            axis.set_ylim(0,limits[(row,dataset)])
            axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
            axis.grid(True,which="major",alpha=0.25,linewidth=0.35)
            axis.tick_params(direction="in",length=1.8,width=0.45,pad=1)
            for spine in axis.spines.values():spine.set_linewidth(0.5)
            panel=chr(ord("a")+row*4+column)
            axis.set_xlabel(f"({panel}) Dataset {dataset_index}, {task}\n{cpus} cores",labelpad=1.0)

    handles=[plt.Line2D([0],[0],color=colors[method],linewidth=0.8,label=method) for method in methods]
    figure.legend(handles=handles,loc="upper center",ncol=6,frameon=False,bbox_to_anchor=(0.5,0.99),handlelength=1.0,handletextpad=0.25,columnspacing=0.75,borderaxespad=0,fontsize=8.0)
    figure.supxlabel("Relative Error Bound",x=0.515,y=0.018)
    figure.supylabel("Aggregate Throughput (GiB/s)",x=0.007,y=0.47)
    figure.tight_layout(rect=(0.025,0.045,1,0.91),pad=0.25,w_pad=0.45,h_pad=0.75)

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
