import csv,os
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FixedFormatter
from plot_style import PaperFigureSpec,add_panel_captions,align_xlabels,apply_paper_style,layout_panel_row,legend_above,panel_aspect_ratios,save_paper_png_vertical_crop,style_plot_axis

project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
tree_count_paths=[os.path.join(project_root,f"results/RD/pb0_fdwz_tree_count/tree{tree_num}/pb0.f00001_P.csv") for tree_num in range(1,6)]
prediction_error_distribution_path=os.path.join(project_root,"results/RD/pb0_fdwz_tree_prediction_error_distribution_1e-3.csv")
canonical_output_path=os.path.join(project_root,"论文/figs/5/fdwz-tree-ablation.png")

X_TICK_GAP_PT=1.5
Y_TICK_GAP_PT=1.0
X_LABEL_GAP_PT=1.0
Y_LABEL_GAP_PT=1.0
PLOT_ASPECT_RATIOS=panel_aspect_ratios(2,1)
PANEL_GAP_PT=9.0
LEFT_EDGE_GAP_PT=3.0
RIGHT_EDGE_GAP_PT=3.0
PANEL_CAPTION_GAP_PT=2.0
LEGEND_GAP_PT=0.0
VERTICAL_EDGE_GAP_PT=3.0

tree_colors=["#000000","#0072B2","#D55E00","#009E73","#CC79A7"]
tree_styles=["-","--","-.",":",(0,(3,1,1,1))]

def read_br_at_eb(path,target_eb):
    with open(path,encoding="utf-8") as csv_file:
        for row in csv.DictReader(csv_file):
            if float(row["REL EB"])==target_eb:
                return float(row["FDWZ-ablation_BR"])
    raise ValueError(f"{path} does not contain REL EB={target_eb}")

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
    style_plot_axis(axes,tick_length=2.5,tick_width=0.6,x_tick_pad=X_TICK_GAP_PT,y_tick_pad=Y_TICK_GAP_PT,x_label_pad=X_LABEL_GAP_PT,y_label_pad=Y_LABEL_GAP_PT,spine_width=0.6,minor_tick_length=1.5,minor_tick_width=0.5)

def draw_tree_count_bars(axes):
    tree_nums=list(range(1,6))
    bit_rates=[read_br_at_eb(path,1e-3) for path in tree_count_paths]
    axes.bar(tree_nums,bit_rates,width=0.68,color=["#9ecae1","#6baed6","#4292c6","#2171b5","#08519c"],edgecolor="black",linewidth=0.45)
    axes.set_xlim(0.45,5.55)
    axes.set_ylim(0,0.27)
    axes.xaxis.set_major_locator(FixedLocator(tree_nums))
    axes.yaxis.set_major_locator(FixedLocator([0,0.1,0.2]))
    axes.yaxis.set_major_formatter(FixedFormatter(["0","0.1","0.2"]))
    axes.set_xlabel("Number of Trees")
    axes.set_ylabel("BR")
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
    axes.set_xlabel("Prediction Error")
    axes.set_ylabel("Probability Density")
    handles,labels=axes.get_legend_handles_labels()
    order=[0,3,1,4,2]
    configure_common_axes(axes)
    return [handles[i] for i in order],[labels[i] for i in order]

def make_ablation_figure(output_path):
    distributions,totals,bin_width=read_prediction_error_distributions(prediction_error_distribution_path)
    spec=PaperFigureSpec.column(1)
    figure,axes=plt.subplots(1,2,figsize=spec.size_in,dpi=spec.dpi)
    left_axes,right_axes=axes
    legend_handles,legend_labels=draw_prediction_error_density(left_axes,distributions,totals,bin_width,(-0.15,0.15),[-0.15,0,0.15])
    left_axes.xaxis.set_major_formatter(FixedFormatter(["-0.15","0","0.15"]))
    left_axes.get_xticklabels()[0].set_horizontalalignment("left")
    left_axes.get_xticklabels()[-1].set_horizontalalignment("right")
    draw_tree_count_bars(right_axes)
    layout_panel_row(figure,axes,PLOT_ASPECT_RATIOS,bottom_pt=30,panel_gap_pt=PANEL_GAP_PT,left_edge_gap_pt=LEFT_EDGE_GAP_PT,right_edge_gap_pt=RIGHT_EDGE_GAP_PT)
    legend_above(left_axes,legend_handles,legend_labels,gap_pt=LEGEND_GAP_PT,ncol=3,frameon=False,borderpad=0,handlelength=1,handletextpad=0.2,columnspacing=0.45,labelspacing=0.05)
    for axis in axes:
        axis.set_axisbelow(True)
        axis.grid(True,linestyle="--",linewidth=0.45,alpha=0.35)
    align_xlabels(figure,axes,gap_pt=X_LABEL_GAP_PT)
    add_panel_captions(figure,axes,["(a) Prediction-Error Density","(b) Tree Count"],gap_pt=PANEL_CAPTION_GAP_PT)
    save_paper_png_vertical_crop(figure,output_path,spec,edge_gap_pt=VERTICAL_EDGE_GAP_PT)
    plt.close(figure)

apply_paper_style(serif_fonts=["STIXGeneral"])
make_ablation_figure(canonical_output_path)
