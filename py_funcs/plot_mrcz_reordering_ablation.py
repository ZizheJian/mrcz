# 绘制MRCZ系数重排的rate--distortion和throughput消融结果。

import argparse,csv,os
import numpy as np
from matplotlib import pyplot as plt
from plot_style import apply_paper_style,style_plot_axis

DEFAULT_ROOT="results/mrcz_reordering"
DEFAULT_OUTPUT="论文/figs/6/mrcz-reordering-ablation.pdf"
DEFAULT_SUMMARY="results/mrcz_reordering/summary.csv"

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

    apply_paper_style(font_size=8.4,serif_fonts=["STIXGeneral"])
    fig=plt.figure(figsize=(3.5,1.65))
    rd_ax=fig.add_axes([0.12,0.28,0.25,0.50])
    throughput_ax=fig.add_axes([0.48,0.28,0.49,0.50])
    axes=(rd_ax,throughput_ax)
    colors={"Element-major":"#ff5a5a","Mode-major":"#3b7ded"}
    rd_ax.plot(element_br,element_psnr,color=colors["Element-major"],linewidth=0.9,label="Element-major")
    rd_ax.plot(mode_br,mode_psnr,color=colors["Mode-major"],linewidth=0.9,label="Mode-major")
    rd_ax.set_xlabel("Bit rate",labelpad=1)
    rd_ax.set_ylabel("PSNR (dB)",labelpad=1)
    rd_ax.set_xlim(0,3.0)
    rd_ax.set_ylim(60,135)
    rd_ax.set_xticks([0,1,2,3])
    rd_ax.legend(loc="lower center",bbox_to_anchor=(0.5,1.02),ncol=1,frameon=False,handlelength=1.0,handletextpad=0.25,labelspacing=0.05,borderaxespad=0)

    labels=list(throughputs)
    x=np.arange(len(labels))
    width=0.34
    compression_mean=np.asarray([throughputs[label][0].mean() for label in labels])
    decompression_mean=np.asarray([throughputs[label][1].mean() for label in labels])
    throughput_ax.bar(x-width/2,compression_mean,width,color="#e3912e",edgecolor="black",linewidth=0.4,label="Compression")
    throughput_ax.bar(x+width/2,decompression_mean,width,color="#4c91cf",edgecolor="black",linewidth=0.4,label="Decompression")
    throughput_ax.set_ylabel("Throughput\n(MB/s per CPU)",labelpad=1)
    throughput_ax.set_xticks(x)
    throughput_ax.set_xticklabels(["Original","Direct\ntranspose","Tiled\ntranspose"])
    throughput_ax.set_ylim(0,70)
    throughput_ax.legend(loc="lower center",bbox_to_anchor=(0.5,1.02),ncol=2,frameon=False,handlelength=0.7,handletextpad=0.2,columnspacing=0.55,borderaxespad=0)

    for ax in axes:
        ax.grid(True,axis="y",linestyle="--",linewidth=0.4,alpha=0.35)
        style_plot_axis(ax,tick_length=2,tick_width=0.5,tick_pad=1.5,spine_width=0.6)
    for label,ax in zip(("(a) Rate–Distortion","(b) Throughput"),axes):
        position=ax.get_position()
        fig.text((position.x0+position.x1)/2,0.04,label,ha="center",va="bottom",fontsize=8)
    os.makedirs(os.path.dirname(args.output) or ".",exist_ok=True)
    fig.savefig(args.output)

if __name__=="__main__":
    main()
