# 绘制MRCZ系数重排的rate--distortion和throughput消融结果。

import argparse,csv,os
import matplotlib
matplotlib.use("Agg")
import numpy as np
from matplotlib import pyplot as plt
from plot_style import PaperFigureSpec,add_panel_captions,align_xlabels,apply_paper_style,layout_panel_row,legend_above,panel_aspect_ratios,save_paper_png_vertical_crop,style_plot_axis

DEFAULT_ROOT="results/mrcz_reordering"
DEFAULT_OUTPUT="论文/figs/6/mrcz-reordering-ablation.png"
DEFAULT_SUMMARY="results/mrcz_reordering/summary.csv"

X_TICK_GAP_PT=1.5
Y_TICK_GAP_PT=1.0
X_LABEL_GAP_PT=1.0
Y_LABEL_GAP_PT=1.0
PLOT_ASPECT_RATIOS=panel_aspect_ratios(1,2)
PANEL_GAP_PT=9.0
LEFT_EDGE_GAP_PT=3.0
RIGHT_EDGE_GAP_PT=3.0
PANEL_CAPTION_GAP_PT=2.0
BOTTOM_PT=40.0
VERTICAL_EDGE_GAP_PT=3.0

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--root",default=DEFAULT_ROOT)
    parser.add_argument("--output",default=DEFAULT_OUTPUT)
    parser.add_argument("--summary",default=DEFAULT_SUMMARY)
    return parser.parse_args()

def read_rows(path):
    with open(path,newline="") as file_obj:
        return list(csv.DictReader(file_obj))

def rd_curve(root,layout):
    path=os.path.join(root,"rd",layout,"pb0.f00001_P.csv")
    rows=sorted(read_rows(path),key=lambda row:float(row["MRCZ-ablation_PSNR"]))
    br=np.asarray([float(row["MRCZ-ablation_BR"]) for row in rows])
    psnr=np.asarray([float(row["MRCZ-ablation_PSNR"]) for row in rows])
    eb=np.asarray([float(row["REL EB"]) for row in rows])
    return eb,br,psnr

def throughput(root,layout):
    path=os.path.join(root,"throughput",layout,"pb0.f00001_P_throughput.csv")
    rows=read_rows(path)
    compression=np.asarray([float(row["MRCZ-ablation Compression Throughput"]) for row in rows])
    decompression=np.asarray([float(row["MRCZ-ablation Decompression Throughput"]) for row in rows])
    return compression,decompression

def write_summary(path,element_br,mode_br,throughputs):
    os.makedirs(os.path.dirname(path) or ".",exist_ok=True)
    reductions=(element_br-mode_br)/element_br*100
    cr_improvements=(element_br/mode_br-1)*100
    with open(path,"w",newline="") as file_obj:
        writer=csv.writer(file_obj,lineterminator="\n")
        writer.writerow(["Metric","Value"])
        writer.writerow(["Minimum BR reduction (%)",reductions.min()])
        writer.writerow(["Maximum BR reduction (%)",reductions.max()])
        writer.writerow(["Mean BR reduction (%)",reductions.mean()])
        writer.writerow(["Minimum CR improvement (%)",cr_improvements.min()])
        writer.writerow(["Maximum CR improvement (%)",cr_improvements.max()])
        writer.writerow(["Mean CR improvement (%)",cr_improvements.mean()])
        for label,(compression,decompression) in throughputs.items():
            writer.writerow([f"{label} compression throughput mean (MB/s per CPU)",compression.mean()])
            writer.writerow([f"{label} compression throughput std (MB/s per CPU)",compression.std(ddof=1)])
            writer.writerow([f"{label} decompression throughput mean (MB/s per CPU)",decompression.mean()])
            writer.writerow([f"{label} decompression throughput std (MB/s per CPU)",decompression.std(ddof=1)])
        direct_decompression=throughputs["Direct"][1].mean()
        tiled_decompression=throughputs["Tiled"][1].mean()
        writer.writerow(["Tiled decompression improvement over direct (%)",(tiled_decompression/direct_decompression-1)*100])

def main():
    args=parse_args()
    _,element_br,element_psnr=rd_curve(args.root,"element_major")
    _,mode_br,mode_psnr=rd_curve(args.root,"mode_major")
    throughputs={
        "Element-major":throughput(args.root,"element_major"),
        "Direct":throughput(args.root,"direct_transpose"),
        "Tiled":throughput(args.root,"tiled_transpose"),
    }
    write_summary(args.summary,element_br,mode_br,throughputs)

    apply_paper_style(serif_fonts=["STIXGeneral"])
    spec=PaperFigureSpec.column(1)
    fig,axes=plt.subplots(1,2,figsize=spec.size_in,dpi=spec.dpi)
    rd_ax,throughput_ax=axes
    axes=(rd_ax,throughput_ax)
    colors={"Element-major":"#ff5a5a","Mode-major":"#3b7ded"}
    rd_ax.plot(element_br,element_psnr,color=colors["Element-major"],linewidth=0.9,label="Element-major")
    rd_ax.plot(mode_br,mode_psnr,color=colors["Mode-major"],linewidth=0.9,label="Mode-major")
    rd_ax.set_xlabel("Bit rate")
    rd_ax.set_ylabel("PSNR (dB)")
    rd_ax.set_xlim(0,3.0)
    rd_ax.set_ylim(60,135)
    rd_ax.set_xticks([0,1,2,3])
    rd_legend_handles,rd_legend_labels=rd_ax.get_legend_handles_labels()

    labels=list(throughputs)
    x=np.arange(len(labels))
    width=0.34
    compression_mean=np.asarray([throughputs[label][0].mean() for label in labels])
    decompression_mean=np.asarray([throughputs[label][1].mean() for label in labels])
    throughput_ax.bar(x-width/2,compression_mean,width,color="#e3912e",edgecolor="black",linewidth=0.4,label="Compression")
    throughput_ax.bar(x+width/2,decompression_mean,width,color="#4c91cf",edgecolor="black",linewidth=0.4,label="Decompression")
    throughput_ax.set_ylabel("Throughput\n(MB/s per CPU)")
    throughput_ax.set_xticks(x)
    throughput_ax.set_xticklabels(["Original","Direct\ntranspose","Tiled\ntranspose"])
    throughput_ax.set_ylim(0,70)
    throughput_legend_handles,throughput_legend_labels=throughput_ax.get_legend_handles_labels()

    for ax in axes:
        ax.grid(True,axis="y",linestyle="--",linewidth=0.4,alpha=0.35)
        style_plot_axis(ax,tick_length=2,tick_width=0.5,x_tick_pad=X_TICK_GAP_PT,y_tick_pad=Y_TICK_GAP_PT,x_label_pad=X_LABEL_GAP_PT,y_label_pad=Y_LABEL_GAP_PT,spine_width=0.6)
    layout_panel_row(fig,axes,PLOT_ASPECT_RATIOS,bottom_pt=BOTTOM_PT,panel_gap_pt=PANEL_GAP_PT,left_edge_gap_pt=LEFT_EDGE_GAP_PT,right_edge_gap_pt=RIGHT_EDGE_GAP_PT)
    legend_above(rd_ax,rd_legend_handles,rd_legend_labels,ncol=1,frameon=False,handlelength=1.0,handletextpad=0.25)
    legend_above(throughput_ax,throughput_legend_handles,throughput_legend_labels,ncol=2,frameon=False,handlelength=0.7,handletextpad=0.2,columnspacing=0.55)
    align_xlabels(fig,axes,gap_pt=X_LABEL_GAP_PT)
    add_panel_captions(fig,axes,["(a) Rate–Distortion","(b) Throughput"],gap_pt=PANEL_CAPTION_GAP_PT)
    os.makedirs(os.path.dirname(args.output) or ".",exist_ok=True)
    save_paper_png_vertical_crop(fig,args.output,spec,edge_gap_pt=VERTICAL_EDGE_GAP_PT)
    plt.close(fig)

if __name__=="__main__":
    main()
