# 绘制固定relative error bound的alpha sweep与固定alpha的rate--distortion曲线。

import argparse,csv,glob,os
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from plot_style import MPL_POINTS_PER_INCH,PAPER_FONT_SIZE,PaperFigureSpec,apply_paper_style,layout_panel_row,save_paper_png_vertical_crop,style_plot_axis

CONFIGS=[
    ("DLT_GLL","DLT (GLL)","#ff4d4d"),
    ("DChT_GLL","DChT (GLL)","#3b7ded"),
    ("DCT_GLL","DCT (GLL)","#2dbf63"),
    ("DLT_uniform","DLT (Uniform)","#bf83ff"),
    ("DChT_uniform","DChT (Uniform)","#f28c00"),
    ("DCT_uniform","DCT (Uniform)","#555555"),
    ("DST_uniform","DST (Uniform)","#12afc7"),
]

X_TICK_GAP_PT=1.5
Y_TICK_GAP_PT=1.0
X_LABEL_GAP_PT=1.0
Y_LABEL_GAP_PT=1.0
PLOT_ASPECT_RATIO=3.0
LEFT_EDGE_GAP_PT=3.0
RIGHT_EDGE_GAP_PT=3.0
BOTTOM_PT=24.0
STYLE_LEGEND_GAP_PT=1.0
CONFIG_LEGEND_GAP_PT=1.0
CONFIG_LEGEND_ROW_GAP_PT=1.0
CONFIG_HANDLE_LENGTH_PT=2.8
CONFIG_HANDLE_TEXT_GAP_PT=1.2
CONFIG_COLUMN_GAP_PT=5.2
STYLE_HANDLE_LENGTH_PT=9.6
STYLE_HANDLE_TEXT_GAP_PT=2.4
STYLE_COLUMN_GAP_PT=16.0
VERTICAL_EDGE_GAP_PT=3.0

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--alpha-root",default="results/fig8/alpha")
    parser.add_argument("--alpha-csv",default="results/fig8/pb0_fig8_alpha_sweep.csv")
    parser.add_argument("--rd-root",default="results/fig8/rd")
    parser.add_argument("--rd-summary",default="results/fig8/pb0_fig8_rd.csv")
    parser.add_argument("--output",default="论文/figs/6/transform-comparison-stacked.png")
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

def add_legends(figure,axes):
    config_by_key={key:Line2D([0],[0],color=color,linewidth=0.9,label=label) for key,label,color in CONFIGS}
    blank_handle=Line2D([0],[0],linestyle="none",label=" ")
    config_handles=[config_by_key["DLT_GLL"],config_by_key["DLT_uniform"],config_by_key["DChT_GLL"],config_by_key["DChT_uniform"],config_by_key["DCT_GLL"],config_by_key["DCT_uniform"],blank_handle,config_by_key["DST_uniform"]]
    style_handles=[
        Line2D([0],[0],color="black",linestyle="-",linewidth=0.8,label=r"Fixed $\varepsilon_{\mathrm{rel}}$, vary $\alpha$"),
        Line2D([0],[0],color="black",linestyle="--",linewidth=0.8,label=r"Fixed $\alpha$, vary $\varepsilon_{\mathrm{rel}}$"),
    ]
    figure_height_pt=figure.get_figheight()*MPL_POINTS_PER_INCH
    axes_position=axes.get_position()
    style_y=axes_position.y1+STYLE_LEGEND_GAP_PT/figure_height_pt
    style_legends=[figure.legend(handles=[handle],loc="lower left",bbox_to_anchor=(0,style_y),bbox_transform=figure.transFigure,ncol=1,frameon=False,borderpad=0,handlelength=STYLE_HANDLE_LENGTH_PT/PAPER_FONT_SIZE,handletextpad=STYLE_HANDLE_TEXT_GAP_PT/PAPER_FONT_SIZE,columnspacing=STYLE_COLUMN_GAP_PT/PAPER_FONT_SIZE,borderaxespad=0,labelspacing=0) for handle in style_handles]
    config_legend=figure.legend(handles=config_handles,loc="lower center",bbox_to_anchor=(0.5,style_y),bbox_transform=figure.transFigure,ncol=4,frameon=False,borderpad=0,handlelength=CONFIG_HANDLE_LENGTH_PT/PAPER_FONT_SIZE,handletextpad=CONFIG_HANDLE_TEXT_GAP_PT/PAPER_FONT_SIZE,columnspacing=CONFIG_COLUMN_GAP_PT/PAPER_FONT_SIZE,borderaxespad=0,labelspacing=CONFIG_LEGEND_ROW_GAP_PT/PAPER_FONT_SIZE)
    figure.canvas.draw()
    renderer=figure.canvas.get_renderer()
    rendered_handles=getattr(config_legend,"legend_handles",None) or config_legend.legendHandles
    column_lefts=[figure.transFigure.inverted().transform((rendered_handles[index].get_window_extent(renderer).x0,0))[0] for index in [0,4]]
    for legend,left in zip(style_legends,column_lefts):legend.set_bbox_to_anchor((left,style_y),transform=figure.transFigure)
    figure.canvas.draw()
    renderer=figure.canvas.get_renderer()
    style_top=max(figure.transFigure.inverted().transform((0,legend.get_window_extent(renderer).y1))[1] for legend in style_legends)
    config_legend.set_bbox_to_anchor((0.5,style_top+CONFIG_LEGEND_GAP_PT/figure_height_pt),transform=figure.transFigure)

def make_figure(alpha_curves,rd_curves,output_path):
    spec=PaperFigureSpec.column(1)
    figure,axes=plt.subplots(figsize=spec.size_in,dpi=spec.dpi)
    for key,_,color in CONFIGS:
        alpha_br,alpha_psnr=alpha_curves[key]
        rd_br,rd_psnr=rd_curves[key]
        linewidth=1.0 if key=="DChT_GLL" else 0.7
        zorder=4 if key=="DChT_GLL" else 2
        axes.plot(alpha_br,alpha_psnr,color=color,linestyle="-",linewidth=linewidth,zorder=zorder)
        axes.plot(rd_br,rd_psnr,color=color,linestyle="--",linewidth=linewidth,zorder=zorder)
    axes.set_xlim(0,3.5)
    axes.set_ylim(80,115)
    axes.set_xlabel("Bit rate")
    axes.set_ylabel("PSNR (dB)")
    axes.grid(True,linestyle="--",linewidth=0.45,alpha=0.35)
    style_plot_axis(axes,tick_length=2.5,tick_width=0.6,x_tick_pad=X_TICK_GAP_PT,y_tick_pad=Y_TICK_GAP_PT,x_label_pad=X_LABEL_GAP_PT,y_label_pad=Y_LABEL_GAP_PT,spine_width=0.7)
    layout_panel_row(figure,[axes],[PLOT_ASPECT_RATIO],bottom_pt=BOTTOM_PT,panel_gap_pt=0,left_edge_gap_pt=LEFT_EDGE_GAP_PT,right_edge_gap_pt=RIGHT_EDGE_GAP_PT)
    add_legends(figure,axes)
    save_paper_png_vertical_crop(figure,output_path,spec,edge_gap_pt=VERTICAL_EDGE_GAP_PT)
    plt.close(figure)

def main():
    args=parse_args()
    alpha_curves,rows_by_alpha=read_alpha_curves(args.alpha_root)
    rd_curves,rows_by_eb=read_rd_curves(args.rd_root)
    write_alpha_summary(args.alpha_csv,rows_by_alpha)
    write_rd_summary(args.rd_summary,rows_by_eb)

    apply_paper_style(serif_fonts=["STIXGeneral"])
    os.makedirs(os.path.dirname(args.output) or ".",exist_ok=True)
    make_figure(alpha_curves,rd_curves,args.output)

if __name__=="__main__":
    main()
