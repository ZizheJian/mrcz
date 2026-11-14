import csv,glob,math,os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator,NullLocator

from plot_style import METHOD_COLORS,TIMING_METHODS,apply_paper_style,save_svg
from paperfig import DEFAULT_PANEL_CAPTION_GAP_PT,PANEL_GAP_WITH_YTICKS_PT,PaperFigureSpec,add_shared_xlabel,layout_panel_grid,legend_above,legend_spacing_kwargs,save_paper_png_vertical_crop,style_plot_axis

apply_paper_style("mrcz-combined-throughput",serif_fonts=["STIXGeneral"])

methods=TIMING_METHODS
colors=METHOD_COLORS
tasks=("Compression","Decompression")
datasets=[("pb0",7,8,64,[1e-2,1e-3,1e-4,1e-5,1e-6]),("jicf0",8,16,128,[1e-4,1e-5,1e-6,1e-7,1e-8])]
PLOT_ASPECT_RATIOS=(2.0,2.0,2.0,2.0)
WIDTH_FRACTION=2*math.sqrt(2)/3
ROW_CONTENT_GAP_PT=6.0
OUTER_EDGE_GAP_PT=3.0
BOTTOM_PLOT_PT=42.0
VERTICAL_EDGE_GAP_PT=3.0


def read_rows(path):
    with open(path,newline="") as input_file:return list(csv.DictReader(input_file))


def harmonic_mean(values):
    return len(values)/sum(1/value for value in values)


def aggregate_small_dataset(project_root,dataset,cpus,partition_num,error_bounds):
    source_dir=os.path.join(project_root,"results","small_core_throughput",dataset)
    tables=[read_rows(path) for path in sorted(glob.glob(os.path.join(source_dir,"*_throughput.csv")))]
    if len(tables)!=5:raise RuntimeError(f"{dataset} requires five field throughput files,found {len(tables)}")
    aggregate={task:{method:[] for method in methods} for task in tasks}
    for error_bound in error_bounds:
        for task in tasks:
            for method in methods:
                field_values=[]
                for table in tables:
                    column=f"{method} {task} Throughput"
                    row=next((row for row in table if int(row["Max CPUs"])==cpus and int(row["Partition Num"])==partition_num and math.isclose(float(row["REL EB"]),error_bound,rel_tol=1e-12) and row.get(column)),None)
                    if row is None:raise RuntimeError(f"missing {dataset}/{method}/{task}/{error_bound:g}")
                    field_values.append(float(row[column]))
                aggregate[task][method].append(harmonic_mean(field_values)*cpus/1024)
    return aggregate


def aggregate_large_dataset(result_dir,dataset,cpus):
    tables=[read_rows(path) for path in sorted(glob.glob(os.path.join(result_dir,dataset,"*_throughput.csv")))]
    if len(tables)!=5:raise RuntimeError(f"{dataset} requires five field throughput files,found {len(tables)}")
    error_bounds=sorted({float(row["REL EB"]) for table in tables for row in table if int(row["Max CPUs"])==cpus},reverse=True)
    aggregate={task:{method:[] for method in methods} for task in tasks}
    for error_bound in error_bounds:
        for task in tasks:
            for method in methods:
                field_values=[]
                for table in tables:
                    row=next(row for row in table if int(row["Max CPUs"])==cpus and math.isclose(float(row["REL EB"]),error_bound,rel_tol=1e-12))
                    field_values.append(float(row[f"{method} {task} Throughput"]))
                aggregate[task][method].append(harmonic_mean(field_values)*cpus/1024)
    return error_bounds,aggregate


def group_limit(aggregate,extra_aggregate):
    maximum_value=max(value for task in tasks for values in aggregate[task].values() for value in values)
    maximum_value=max(maximum_value,max(value for task in tasks for value in extra_aggregate[task]))
    return math.ceil(maximum_value*1.05*2)/2


def mgard_8core(project_root,dataset,error_bounds):
    directory=os.path.join(project_root,"results","small_core_throughput",dataset)
    tables=[read_rows(path) for path in sorted(glob.glob(os.path.join(directory,"*_throughput.csv")))]
    if len(tables)!=5:raise RuntimeError(f"{dataset} requires five MGARD throughput files,found {len(tables)}")
    aggregate={task:[] for task in tasks}
    for error_bound in error_bounds:
        for task in tasks:
            field_values=[]
            for table in tables:
                row=next((row for row in table if int(row["Max CPUs"])==8 and math.isclose(float(row["REL EB"]),error_bound,rel_tol=1e-12)),None)
                if row is None or not row.get(f"MGARD-1D {task} Throughput"):raise RuntimeError(f"missing {dataset}/MGARD-1D/{task}/{error_bound:g}")
                field_values.append(float(row[f"MGARD-1D {task} Throughput"]))
            aggregate[task].append(harmonic_mean(field_values)*8/1024)
    return aggregate


def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    large_result_dir=os.path.join(project_root,"results","RD")
    data={}
    limits={}
    for dataset,dataset_index,small_cpus,partition_num,error_bounds in datasets:
        small_aggregate=aggregate_small_dataset(project_root,dataset,small_cpus,partition_num,error_bounds)
        large_cpus=partition_num
        large_error_bounds,large_aggregate=aggregate_large_dataset(large_result_dir,dataset,large_cpus)
        mgard_aggregate=mgard_8core(project_root,dataset,error_bounds)
        data[(0,dataset)]=(small_cpus,error_bounds,small_aggregate,mgard_aggregate)
        data[(1,dataset)]=(large_cpus,large_error_bounds,large_aggregate,mgard_aggregate)
        limits[(0,dataset)]=group_limit(small_aggregate,mgard_aggregate)
        projected_mgard={task:[value*large_cpus/8 for value in mgard_aggregate[task]] for task in tasks}
        limits[(1,dataset)]=group_limit(large_aggregate,projected_mgard)

    columns=[("pb0",7,"Compression"),("pb0",7,"Decompression"),("jicf0",8,"Compression"),("jicf0",8,"Decompression")]
    spec=PaperFigureSpec.text(1,width_fraction=WIDTH_FRACTION)
    figure,axes=plt.subplots(2,4,figsize=spec.size_in,dpi=spec.dpi)
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
            axis.set_xlabel(f"({panel}) Dataset {dataset_index}, {task}\n{cpus} cores")
            axis.xaxis.labelpad=DEFAULT_PANEL_CAPTION_GAP_PT

    handles=[plt.Line2D([0],[0],color=colors[method],linewidth=0.8,label=method) for method in methods]
    handles.extend([plt.Line2D([0],[0],color=colors["MGARD-1D"],linewidth=1.0,label="MGARD-X-1D measured"),plt.Line2D([0],[0],color=colors["MGARD-1D"],linestyle="--",linewidth=1.0,label="MGARD-X-1D ideal")])
    layout_panel_grid(figure,axes,PLOT_ASPECT_RATIOS,bottom_pt=BOTTOM_PLOT_PT,column_gap_pt=PANEL_GAP_WITH_YTICKS_PT,row_content_gap_pt=ROW_CONTENT_GAP_PT,left_edge_gap_pt=OUTER_EDGE_GAP_PT,right_edge_gap_pt=OUTER_EDGE_GAP_PT,shared_ylabel="Aggregate Throughput (GiB/s)")
    add_shared_xlabel(figure,axes[-1,:],"Relative Error Bound")
    legend_above(axes[0,:],handles=handles,center="figure",ncol=4,**legend_spacing_kwargs(8))

    result_dir=os.path.join(project_root,"results","RD")
    svg_output=os.path.join(result_dir,"combined_throughput.svg")
    png_output=os.path.join(result_dir,"combined_throughput.png")
    paper_output=os.path.join(project_root,"overleaf","figs","7","combined-throughput.png")
    save_svg(figure,svg_output,metadata={"Date":None},bbox_inches="tight",pad_inches=0.02)
    save_paper_png_vertical_crop(figure,png_output,spec,edge_gap_pt=VERTICAL_EDGE_GAP_PT)
    save_paper_png_vertical_crop(figure,paper_output,spec,edge_gap_pt=VERTICAL_EDGE_GAP_PT)
    plt.close(figure)
    print(paper_output)
    print(png_output)


if __name__=="__main__":main()
