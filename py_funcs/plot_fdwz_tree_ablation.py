import csv,os
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FixedFormatter

project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
max_tree_path=os.path.join(project_root,"results/RD/pb0_fdwz_max_tree/pb0.f00001_P.csv")
single_tree_path=os.path.join(project_root,"results/RD/pb0_fdwz_single_tree/pb0.f00001_P.csv")
tree_distribution_path=os.path.join(project_root,"results/RD/pb0_fdwz_tree_qb_distribution.csv")
prediction_error_distribution_path=os.path.join(project_root,"results/RD/pb0_fdwz_tree_prediction_error_distribution_1e-3.csv")
canonical_output_path=os.path.join(project_root,"论文/figs/5/fdwz-tree-ablation.pdf")
cdf_output_path=os.path.join(project_root,"论文/figs/5/fdwz-tree-ablation-cdf.pdf")
frequency_output_path=os.path.join(project_root,"论文/figs/5/fdwz-tree-ablation-frequency.pdf")
linear_full_output_path=os.path.join(project_root,"论文/figs/5/fdwz-tree-qb-frequency-full.pdf")
linear_prediction_output_path=os.path.join(project_root,"论文/figs/5/fdwz-tree-qb-frequency-prediction-focused.pdf")
prediction_error_full_output_path=os.path.join(project_root,"论文/figs/5/fdwz-tree-prediction-error-density-full.pdf")
prediction_error_focused_output_path=os.path.join(project_root,"论文/figs/5/fdwz-tree-prediction-error-density-focused.pdf")

tree_colors=["#000000","#0072B2","#D55E00","#009E73","#CC79A7"]
tree_styles=["-","--","-.",":",(0,(3,1,1,1))]

def read_curve(path):
    with open(path,encoding="utf-8") as csv_file:
        rows=list(csv.DictReader(csv_file))
    return sorted((float(row["FDWZ-ablation_BR"]),float(row["FDWZ-ablation_PSNR"])) for row in rows)

def read_tree_distributions(path):
    distributions=defaultdict(lambda:defaultdict(int))
    with open(path,encoding="utf-8") as csv_file:
        for row in csv.DictReader(csv_file):
            distributions[int(row["Tree ID"])][abs(int(row["QB"]))]+=int(row["Count"])
    return distributions

def read_signed_tree_distributions(path):
    distributions=defaultdict(lambda:defaultdict(int))
    with open(path,encoding="utf-8") as csv_file:
        for row in csv.DictReader(csv_file):
            distributions[int(row["Tree ID"])][int(row["QB"])]+=int(row["Count"])
    return distributions

def read_prediction_error_distributions(path):
    distributions=defaultdict(lambda:defaultdict(int))
    totals=defaultdict(int)
    bin_width=None
    with open(path,encoding="utf-8") as csv_file:
        for row in csv.DictReader(csv_file):
            tree_id=int(row["Tree ID"])
            count=int(row["Count"])
            totals[tree_id]+=count
            bin_width=float(row["Bin Width"])
            distributions[tree_id][int(row["Bin ID"])]+=count
    return distributions,totals,bin_width

def configure_common_axes(axes):
    axes.tick_params(axis="both",which="major",direction="in",length=2.5,width=0.6,pad=1)
    axes.tick_params(axis="both",which="minor",direction="in",length=1.5,width=0.5)
    for spine in axes.spines.values():
        spine.set_linewidth(0.6)

def draw_tree_distribution(axes,distributions,mode):
    for tree_id,counts in sorted(distributions.items()):
        total=sum(counts.values())
        cumulative=0
        points=[]
        for qb,count in sorted(counts.items()):
            cumulative+=count
            y=cumulative/total if mode=="cdf" else count/total
            points.append((qb,y))
        axes.plot([point[0] for point in points],[point[1] for point in points],label=f"Tree{tree_id}",color=tree_colors[tree_id],linestyle=tree_styles[tree_id],linewidth=0.9,drawstyle="steps-post")
    axes.set_xscale("symlog",linthresh=1)
    axes.set_xlim(0,4000)
    axes.xaxis.set_major_locator(FixedLocator([0,1,10,100,1000]))
    axes.xaxis.set_major_formatter(FixedFormatter(["0","1","10","100","1k"]))
    axes.set_xlabel(r"$|QB|$",labelpad=1)
    if mode=="cdf":
        axes.set_ylim(0,1.02)
        axes.yaxis.set_major_locator(FixedLocator([0,0.5,1]))
        axes.yaxis.set_major_formatter(FixedFormatter(["0",".5","1"]))
        axes.set_ylabel("CDF",labelpad=1)
    else:
        axes.set_yscale("log")
        axes.set_ylim(1e-8,1.2)
        axes.yaxis.set_major_locator(FixedLocator([1e-6,1e-3,1]))
        axes.set_ylabel("Frequency",labelpad=1)
    axes.legend(loc="lower center",bbox_to_anchor=(0.5,0.97),ncol=3,frameon=False,handlelength=1.1,handletextpad=0.25,columnspacing=0.6,borderaxespad=0,labelspacing=0.1)
    configure_common_axes(axes)

def draw_rate_distortion(axes):
    for path,label,color,marker in [(max_tree_path,"Multi-tree","#0072B2","o"),(single_tree_path,"Single-tree","#D55E00","s")]:
        points=read_curve(path)
        axes.plot([point[0] for point in points],[point[1] for point in points],label=label,color=color,marker=marker,markersize=2.8,linewidth=1,markeredgewidth=0.5)
    axes.set_xlim(0,4.5)
    axes.set_ylim(40,130)
    axes.xaxis.set_major_locator(FixedLocator([0,2,4]))
    axes.xaxis.set_major_formatter(FixedFormatter(["0","2","4"]))
    axes.yaxis.set_major_locator(FixedLocator([40,80,120]))
    axes.set_xlabel("Bit rate",labelpad=1)
    axes.set_ylabel("PSNR (dB)",labelpad=1)
    axes.legend(loc="lower center",bbox_to_anchor=(0.5,0.98),frameon=False,handlelength=1,handletextpad=0.25,borderaxespad=0,labelspacing=0.1)
    configure_common_axes(axes)

def draw_linear_qb_frequency(axes,distributions,x_limits,x_ticks):
    for tree_id,counts in sorted(distributions.items()):
        total=sum(counts.values())
        qbs=range(min(counts),max(counts)+1)
        axes.plot(qbs,[counts.get(qb,0)/total for qb in qbs],label=f"Tree{tree_id}",color=tree_colors[tree_id],linestyle=tree_styles[tree_id],linewidth=0.9)
    axes.set_xlim(*x_limits)
    axes.set_ylim(0,1.02)
    axes.xaxis.set_major_locator(FixedLocator(x_ticks))
    axes.yaxis.set_major_locator(FixedLocator([0,0.5,1]))
    axes.yaxis.set_major_formatter(FixedFormatter(["0",".5","1"]))
    axes.set_xlabel("QB",labelpad=1)
    axes.set_ylabel("Frequency",labelpad=1)
    axes.legend(loc="lower center",bbox_to_anchor=(0.5,0.97),ncol=3,frameon=False,handlelength=1.1,handletextpad=0.25,columnspacing=0.6,borderaxespad=0,labelspacing=0.1)
    configure_common_axes(axes)

def draw_prediction_error_density(axes,distributions,totals,bin_width,x_limits,x_ticks):
    first_visible_bin=int(x_limits[0]/bin_width)-1
    past_last_visible_bin=int(x_limits[1]/bin_width)+2
    maximum_density=0
    for tree_id,counts in sorted(distributions.items()):
        first_bin=max(min(counts),first_visible_bin)
        past_last_bin=min(max(counts)+1,past_last_visible_bin)
        bins=range(first_bin,past_last_bin)
        densities=[counts.get(bin_id,0)/(totals[tree_id]*bin_width) for bin_id in bins]
        maximum_density=max(maximum_density,max(densities))
        axes.plot([bin_id*bin_width for bin_id in bins],densities,label=f"Tree{tree_id}",color=tree_colors[tree_id],linestyle=tree_styles[tree_id],linewidth=0.9,drawstyle="steps-mid")
    axes.set_xlim(*x_limits)
    axes.set_ylim(0,maximum_density*1.08)
    axes.xaxis.set_major_locator(FixedLocator(x_ticks))
    axes.yaxis.set_major_locator(FixedLocator([0,5,10]))
    axes.set_xlabel("Prediction Error",labelpad=1)
    axes.set_ylabel("Probability Density",labelpad=1)
    handles,labels=axes.get_legend_handles_labels()
    order=[0,3,1,4,2]
    axes.legend([handles[i] for i in order],[labels[i] for i in order],loc="lower center",bbox_to_anchor=(0.5,0.98),ncol=3,frameon=False,handlelength=1,handletextpad=0.2,columnspacing=0.45,borderaxespad=0,labelspacing=0.05)
    configure_common_axes(axes)

def make_qb_figure(mode,output_paths):
    figure=plt.figure(figsize=(3.5,1.55))
    left_axes=figure.add_axes([0.12,0.34,0.49,0.48])
    right_axes=figure.add_axes([0.735,0.34,0.255,0.48])
    draw_tree_distribution(left_axes,read_tree_distributions(tree_distribution_path),mode)
    draw_rate_distortion(right_axes)
    figure.text(0.365,0.015,"(a) QB distributions",ha="center",va="bottom",fontsize=9)
    figure.text(0.8625,0.015,"(b) Rate-distortion",ha="center",va="bottom",fontsize=9)
    for output_path in output_paths:
        figure.savefig(output_path,format="pdf",dpi=300)
    plt.close(figure)

def make_ablation_figure(output_path):
    distributions,totals,bin_width=read_prediction_error_distributions(prediction_error_distribution_path)
    figure=plt.figure(figsize=(3.5,1.35))
    left_axes=figure.add_axes([0.105,0.31,0.505,0.51])
    right_axes=figure.add_axes([0.725,0.31,0.27,0.51])
    draw_prediction_error_density(left_axes,distributions,totals,bin_width,(-0.15,0.15),[-0.15,0,0.15])
    left_axes.xaxis.set_major_formatter(FixedFormatter(["-0.15","0","0.15"]))
    left_axes.get_xticklabels()[0].set_horizontalalignment("left")
    left_axes.get_xticklabels()[-1].set_horizontalalignment("right")
    draw_rate_distortion(right_axes)
    for axes in (left_axes,right_axes):
        axes.set_axisbelow(True)
        axes.grid(True,linestyle="--",linewidth=0.45,alpha=0.35)
    figure.text(0.3575,0.015,"(a) Prediction-Error Density",ha="center",va="bottom")
    figure.text(0.86,0.015,"(b) Rate-Distortion",ha="center",va="bottom")
    figure.savefig(output_path,format="pdf",dpi=300)
    plt.close(figure)

def make_linear_distribution_figure(x_limits,x_ticks,output_path):
    figure=plt.figure(figsize=(2.35,1.55))
    axes=figure.add_axes([0.18,0.25,0.80,0.57])
    draw_linear_qb_frequency(axes,read_signed_tree_distributions(tree_distribution_path),x_limits,x_ticks)
    figure.savefig(output_path,format="pdf",dpi=300)
    plt.close(figure)

def make_prediction_error_distribution_figure(x_limits,x_ticks,output_path):
    distributions,totals,bin_width=read_prediction_error_distributions(prediction_error_distribution_path)
    figure=plt.figure(figsize=(2.35,1.55))
    axes=figure.add_axes([0.22,0.25,0.76,0.57])
    draw_prediction_error_density(axes,distributions,totals,bin_width,x_limits,x_ticks)
    figure.savefig(output_path,format="pdf",dpi=300)
    plt.close(figure)

def make_focused_prediction_error_figure(output_path):
    distributions,totals,bin_width=read_prediction_error_distributions(prediction_error_distribution_path)
    figure=plt.figure(figsize=(2.35,1.55))
    axes=figure.add_axes([0.22,0.25,0.76,0.57])
    draw_prediction_error_density(axes,distributions,totals,bin_width,(-0.15,0.15),[-0.15,0,0.15])
    axes.xaxis.set_major_formatter(FixedFormatter([r"$-.15$","0",r"$.15$"]))
    axes.get_xticklabels()[0].set_horizontalalignment("left")
    axes.get_xticklabels()[-1].set_horizontalalignment("right")
    figure.savefig(output_path,format="pdf",dpi=300)
    plt.close(figure)

plt.rcParams.update({"font.family":"serif","font.serif":["Nimbus Roman","Times New Roman","Times"],"font.size":8,"axes.labelsize":8,"axes.titlesize":8,"xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":8,"pdf.fonttype":42,"ps.fonttype":42})
make_ablation_figure(canonical_output_path)
make_qb_figure("cdf",[cdf_output_path])
make_qb_figure("frequency",[frequency_output_path])
make_linear_distribution_figure((-3300,1900),[-3000,-1500,0,1500],linear_full_output_path)
make_linear_distribution_figure((-200,700),[-200,0,200,400,600],linear_prediction_output_path)
make_prediction_error_distribution_figure((-225,390),[-200,0,200],prediction_error_full_output_path)
make_focused_prediction_error_figure(prediction_error_focused_output_path)
