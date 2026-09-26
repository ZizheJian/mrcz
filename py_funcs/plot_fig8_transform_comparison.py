# 绘制固定relative error bound的alpha sweep与固定alpha的rate--distortion曲线。

import argparse,csv,glob,os
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

CONFIGS=[
    ("DLT_GLL","DLT (GLL)","#ff4d4d"),
    ("DChT_GLL","DChT (GLL)","#3b7ded"),
    ("DCT_GLL","DCT (GLL)","#2dbf63"),
    ("DLT_uniform","DLT (Uniform)","#bf83ff"),
    ("DChT_uniform","DChT (Uniform)","#f28c00"),
    ("DCT_uniform","DCT (Uniform)","#555555"),
    ("DST_uniform","DST (Uniform)","#12afc7"),
]

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--alpha-root",default="results/fig8/alpha")
    parser.add_argument("--alpha-csv",default="results/fig8/pb0_fig8_alpha_sweep.csv")
    parser.add_argument("--rd-root",default="results/fig8/rd")
    parser.add_argument("--rd-summary",default="results/fig8/pb0_fig8_rd.csv")
    parser.add_argument("--output",default="论文/figs/6/transform-comparison.pdf")
    parser.add_argument("--legend-layout",choices=("side","stacked"),default="stacked")
    return parser.parse_args()

def read_alpha_curves(root):
    curves={}
    rows_by_alpha={}
    for key,_,_ in CONFIGS:
        points=[]
        pattern=os.path.join(root,key,"alpha_*","pb0.f00001_P.csv")
        for path in glob.glob(pattern):
            alpha=float(os.path.basename(os.path.dirname(path))[len("alpha_"):])
            with open(path,newline="") as file_obj:
                rows=list(csv.DictReader(file_obj))
            if len(rows)!=1:
                raise ValueError(f"Expected one row in {path},found {len(rows)}")
            row=rows[0]
            br=float(row["MRCZ-ablation_BR"])
            psnr=float(row["MRCZ-ablation_PSNR"])
            points.append((alpha,br,psnr))
            rows_by_alpha.setdefault(alpha,{})[key]=(row["MRCZ-ablation_BR"],row["MRCZ-ablation_PSNR"])
        points.sort()
        curves[key]=([point[1] for point in points],[point[2] for point in points])
    return curves,rows_by_alpha

def write_alpha_summary(path,rows_by_alpha):
    os.makedirs(os.path.dirname(path) or ".",exist_ok=True)
    with open(path,"w",newline="") as file_obj:
        writer=csv.writer(file_obj,lineterminator="\n")
        header=["alpha"]
        for key,_,_ in CONFIGS:
            header.extend([f"{key}_BR",f"{key}_PSNR"])
        writer.writerow(header)
        for alpha in sorted(rows_by_alpha):
            row=[f"{alpha:g}"]
            for key,_,_ in CONFIGS:
                row.extend(rows_by_alpha[alpha].get(key,("","")))
            writer.writerow(row)

def read_rd_curves(root):
    curves={}
    rows_by_eb={}
    for key,_,_ in CONFIGS:
        path=os.path.join(root,key,"pb0.f00001_P.csv")
        with open(path,newline="") as file_obj:
            rows=sorted(csv.DictReader(file_obj),key=lambda row:float(row["REL EB"]),reverse=True)
        curves[key]=([float(row["MRCZ-ablation_BR"]) for row in rows],[float(row["MRCZ-ablation_PSNR"]) for row in rows])
        for row in rows:
            eb=float(row["REL EB"])
            rows_by_eb.setdefault(eb,{})[key]=(row["MRCZ-ablation_BR"],row["MRCZ-ablation_PSNR"])
    return curves,rows_by_eb

def write_rd_summary(path,rows_by_eb):
    os.makedirs(os.path.dirname(path) or ".",exist_ok=True)
    with open(path,"w",newline="") as file_obj:
        writer=csv.writer(file_obj,lineterminator="\n")
        header=["relative_eb"]
        for key,_,_ in CONFIGS:
            header.extend([f"{key}_BR",f"{key}_PSNR"])
        writer.writerow(header)
        for eb in sorted(rows_by_eb,reverse=True):
            row=[f"{eb:g}"]
            for key,_,_ in CONFIGS:
                row.extend(rows_by_eb[eb][key])
            writer.writerow(row)

def main():
    args=parse_args()
    alpha_curves,rows_by_alpha=read_alpha_curves(args.alpha_root)
    rd_curves,rows_by_eb=read_rd_curves(args.rd_root)
    write_alpha_summary(args.alpha_csv,rows_by_alpha)
    write_rd_summary(args.rd_summary,rows_by_eb)

    plt.rcParams.update({"font.family":"serif","font.serif":["Nimbus Roman"],"mathtext.fontset":"stix","font.size":8.4,"axes.labelsize":8.4,"xtick.labelsize":8.4,"ytick.labelsize":8.4,"legend.fontsize":8.4,"pdf.fonttype":42,"ps.fonttype":42})
    fig,ax=plt.subplots(figsize=(3.67,2.08) if args.legend_layout=="stacked" else (3.5,1.55))
    for key,_,color in CONFIGS:
        alpha_br,alpha_psnr=alpha_curves[key]
        rd_br,rd_psnr=rd_curves[key]
        linewidth=1.0 if key=="DChT_GLL" else 0.7
        zorder=4 if key=="DChT_GLL" else 2
        ax.plot(alpha_br,alpha_psnr,color=color,linestyle="-",linewidth=linewidth,zorder=zorder)
        ax.plot(rd_br,rd_psnr,color=color,linestyle="--",linewidth=linewidth,zorder=zorder)

    ax.set_xlim(0,3.5)
    ax.set_ylim(80,115)
    ax.set_xlabel("Bit rate",labelpad=1)
    ax.set_ylabel("PSNR (dB)",labelpad=1)
    ax.grid(True,linestyle="--",linewidth=0.45,alpha=0.35)
    ax.tick_params(direction="in",length=2.5,width=0.6)
    for spine in ax.spines.values():spine.set_linewidth(0.7)

    config_handles=[Line2D([0],[0],color=color,linewidth=0.9,label=label) for _,label,color in CONFIGS]
    if args.legend_layout=="stacked":
        config_by_key={key:Line2D([0],[0],color=color,linewidth=0.9,label=label) for key,label,color in CONFIGS}
        blank_handle=Line2D([0],[0],linestyle="none",label=" ")
        config_columns=[
            (0.045,[config_by_key["DLT_GLL"],config_by_key["DLT_uniform"]]),
            (0.285,[config_by_key["DChT_GLL"],config_by_key["DChT_uniform"]]),
            (0.53,[config_by_key["DCT_GLL"],config_by_key["DCT_uniform"]]),
            (0.762,[blank_handle,config_by_key["DST_uniform"]]),
        ]
        config_legends=[fig.legend(handles=handles,loc="upper left",bbox_to_anchor=(x,0.995),ncol=1,frameon=False,handlelength=0.35,handletextpad=0.15,borderaxespad=0.0) for x,handles in config_columns]
    else:
        config_legend=ax.legend(handles=config_handles,loc="lower left",bbox_to_anchor=(-0.08,1.01),ncol=2,frameon=False,columnspacing=0.35,handlelength=0.9,handletextpad=0.2,borderaxespad=0.0)
        ax.add_artist(config_legend)
    style_handles=[
        Line2D([0],[0],color="black",linestyle="-",linewidth=0.8,label=r"Fixed $\varepsilon_{\mathrm{rel}}$, vary $\alpha$"),
        Line2D([0],[0],color="black",linestyle="--",linewidth=0.8,label=r"Fixed $\alpha$, vary $\varepsilon_{\mathrm{rel}}$"),
    ]
    if args.legend_layout=="stacked":
        style_legends=[
            fig.legend(handles=[style_handles[0]],loc="upper left",bbox_to_anchor=(0.045,0.83),ncol=1,frameon=False,handlelength=1.2,handletextpad=0.3,borderaxespad=0.0),
            fig.legend(handles=[style_handles[1]],loc="upper left",bbox_to_anchor=(0.53,0.83),ncol=1,frameon=False,handlelength=1.2,handletextpad=0.3,borderaxespad=0.0),
        ]
    else:
        style_legend=ax.legend(handles=style_handles,loc="lower right",bbox_to_anchor=(1.05,1.01),ncol=1,frameon=False,handlelength=1.2,handletextpad=0.25,borderaxespad=0.0)

    os.makedirs(os.path.dirname(args.output) or ".",exist_ok=True)
    if args.legend_layout=="stacked":
        fig.subplots_adjust(left=0.12,right=0.98,bottom=0.24,top=0.72)
        fig.savefig(args.output,bbox_inches="tight",pad_inches=0.01)
        fig.savefig(os.path.splitext(args.output)[0]+".png",dpi=300,bbox_inches="tight",pad_inches=0.01)
    else:
        fig.tight_layout(pad=0.25)
        fig.savefig(args.output,bbox_inches="tight",bbox_extra_artists=(config_legend,style_legend))
        fig.savefig(os.path.splitext(args.output)[0]+".png",dpi=300,bbox_inches="tight",bbox_extra_artists=(config_legend,style_legend))

if __name__=="__main__":
    main()
