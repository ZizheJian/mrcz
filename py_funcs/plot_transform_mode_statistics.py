# 比较DLT和DChT在不同tensor-product mode order上的系数幅值与量化符号熵。
# 默认从pb0的P field中均匀抽取多个连续element block，避免扫描完整数据文件。

import argparse
import csv
import os
import sys
import matplotlib
matplotlib.use("Agg")
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from plot_style import DEFAULT_LEGEND_GAP_PT,MPL_POINTS_PER_INCH,PANEL_GAP_WITH_YLABEL_PT,PaperFigureSpec,add_panel_captions,align_xlabels,apply_paper_style,layout_panel_row,legend_spacing_kwargs,save_paper_png_vertical_crop,style_plot_axis

sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from py_funcs.read_nekrs_header import read_nekrs_header

DEFAULT_DATA_PATH="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/element-data/hyuan/pb0.f00001"
DEFAULT_OUTPUT="overleaf/figs/6/transform-mode-statistics.png"
DEFAULT_CSV_OUTPUT="results/RD/transform_mode_statistics.csv"
DEFAULT_STREAM_ROOT="results/fig9/stream_breakdown"

OUTER_PANEL_ASPECT_RATIOS=(1/1.2,1/1.2,1/1.2)
STREAM_AXES_GAP_PT=2.0
LEFT_EDGE_GAP_PT=3.0
RIGHT_EDGE_GAP_PT=3.0
BOTTOM_PT=30.0
VERTICAL_EDGE_GAP_PT=3.0

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--data-path",default=DEFAULT_DATA_PATH)
    parser.add_argument("--output",default=DEFAULT_OUTPUT)
    parser.add_argument("--csv-output",default=DEFAULT_CSV_OUTPUT)
    parser.add_argument("--stream-root",default=DEFAULT_STREAM_ROOT)
    parser.add_argument("--block-count",type=int,default=64)
    parser.add_argument("--block-elements",type=int,default=128)
    parser.add_argument("--relative-eb",type=float,default=1e-4)
    parser.add_argument("--absolute-eb",type=float,default=0.060799235524628958)
    parser.add_argument("--mode-order",choices=("maximum","total"),default="total")
    return parser.parse_args()

def field_offset(header,field):
    base=136+header.file_element_num*4
    value_bytes=header.file_element_num*header.element_size*header.word_size
    if field=="P":
        if "X" in header.variables:
            base+=3*value_bytes
        if "U" in header.variables:
            base+=3*value_bytes
        return base
    raise ValueError(f"不支持的field: {field}")

def read_sample(path,block_count,block_elements):
    with open(path,"rb") as file_obj:
        header=read_nekrs_header(file_obj)
        if header.file_element_shape[0]!=header.file_element_shape[1] or header.file_element_shape[1]!=header.file_element_shape[2]:
            raise ValueError("当前统计脚本要求三个element维度长度相同")
        if block_elements>header.file_element_num:
            raise ValueError("block-elements不能超过element总数")
        starts=np.linspace(0,header.file_element_num-block_elements,block_count,dtype=np.int64)
        dtype=np.dtype(header.np_dtype).newbyteorder(header.endian_mode)
        offset=field_offset(header,"P")
        chunks=[]
        for start in starts:
            file_obj.seek(offset+int(start)*header.element_size*header.word_size)
            chunk=np.fromfile(file_obj,dtype=dtype,count=block_elements*header.element_size)
            if chunk.size!=block_elements*header.element_size:
                raise IOError(f"读取element block失败: start={start}")
            chunks.append(chunk.reshape(block_elements,*header.c_element_shape))
    return np.concatenate(chunks,axis=0).astype(np.float64),header

def build_forward_matrix(n,transform):
    if n==1:
        nodes=np.array([0.0])
    else:
        interior=np.polynomial.legendre.Legendre.basis(n-1).deriv().roots()
        nodes=np.concatenate(([-1.0],interior,[1.0]))
    if transform=="DLT":
        inverse_matrix=np.polynomial.legendre.legvander(nodes,n-1).T
    elif transform=="DChT":
        inverse_matrix=np.polynomial.chebyshev.chebvander(nodes,n-1).T
    else:
        raise ValueError(f"不支持的transform: {transform}")
    return np.linalg.inv(inverse_matrix)

def transform_elements(values,matrix):
    temp=np.einsum("eabc,ai->eibc",values,matrix,optimize=True)
    temp=np.einsum("eibc,bj->eijc",temp,matrix,optimize=True)
    return np.einsum("eijc,ck->eijk",temp,matrix,optimize=True)

def mode_statistics(coeff,value_range,coeff_eb,mode_order_kind):
    n=coeff.shape[1]
    indices=np.indices((n,n,n))
    if mode_order_kind=="maximum":
        mode_order=np.maximum.reduce(indices)
    else:
        mode_order=np.sum(indices,axis=0)
    medians=[]
    lower=[]
    upper=[]
    entropies=[]
    for order in range(int(mode_order.max())+1):
        values=np.abs(coeff[:,mode_order==order]).reshape(-1)/value_range
        medians.append(np.median(values))
        lower.append(np.quantile(values,0.25))
        upper.append(np.quantile(values,0.75))
        raw_values=coeff[:,mode_order==order].reshape(-1)
        quantized=np.rint(raw_values/(2*coeff_eb)).astype(np.int64)
        _,counts=np.unique(quantized,return_counts=True)
        probabilities=counts/counts.sum()
        entropy=-np.sum(probabilities*np.log2(probabilities))
        entropies.append(max(0.0,float(entropy)))
    return np.asarray(medians),np.asarray(lower),np.asarray(upper),np.asarray(entropies)

def read_stream_rows(root,transform):
    path=os.path.join(root,f"{transform}_GLL","pb0.f00001_P.csv")
    with open(path,newline="") as file_obj:
        return sorted(csv.DictReader(file_obj),key=lambda row:float(row["MRCZ-ablation_PSNR"]))

def main():
    args=parse_args()
    values,header=read_sample(args.data_path,args.block_count,args.block_elements)
    value_range=args.absolute_eb/args.relative_eb
    matrix_dlt=build_forward_matrix(header.c_element_shape[0],"DLT")
    matrix_dcht=build_forward_matrix(header.c_element_shape[0],"DChT")
    coeff_dlt=transform_elements(values,matrix_dlt)
    coeff_dcht=transform_elements(values,matrix_dcht)
    dlt=mode_statistics(coeff_dlt,value_range,args.absolute_eb*0.05,args.mode_order)
    dcht=mode_statistics(coeff_dcht,value_range,args.absolute_eb*0.05,args.mode_order)
    orders=np.arange(len(dlt[0]))
    order_label="Maximum mode order" if args.mode_order=="maximum" else "Total mode order"
    axis_order_label="Maximum Order" if args.mode_order=="maximum" else "Total Order"

    apply_paper_style(serif_fonts=["STIXGeneral"])
    spec=PaperFigureSpec.column(1)
    fig=plt.figure(figsize=spec.size_in,dpi=spec.dpi)
    coefficient_ax=fig.add_axes([0.1,0.1,0.2,0.2])
    residual_ax=fig.add_axes([0.1,0.1,0.2,0.2],sharex=coefficient_ax)
    magnitude_ax=fig.add_axes([0.4,0.1,0.2,0.2])
    entropy_ax=fig.add_axes([0.7,0.1,0.2,0.2])
    axes=(coefficient_ax,residual_ax,magnitude_ax,entropy_ax)
    colors={"DLT":"#ff5a5a","DChT":"#3b7ded"}
    mode_marker="o" if args.mode_order=="maximum" else None
    for label in ("DLT","DChT"):
        rows=read_stream_rows(args.stream_root,label)
        relative_eb=[float(row["REL EB"]) for row in rows]
        coefficient_br=[float(row["MRCZ-ablation_Coefficient_BR"]) for row in rows]
        residual_br=[float(row["MRCZ-ablation_Residual_BR"]) for row in rows]
        coefficient_ax.plot(relative_eb,coefficient_br,color=colors[label],linestyle="-",linewidth=0.8,zorder=1)
        residual_ax.plot(relative_eb,residual_br,color=colors[label],linestyle="--",linewidth=0.8,zorder=1)
    coefficient_ax.set_xscale("log")
    coefficient_ax.invert_xaxis()
    coefficient_ax.set_ylim(0,4)
    coefficient_ax.set_yticks([2,4])
    coefficient_ax.tick_params(labelbottom=False)
    residual_ax.set_ylim(0,0.045)
    residual_ax.set_yticks([0.02,0.04])
    residual_ax.set_xticks([1e-2,1e-4,1e-6])
    residual_ax.set_xlabel("Relative Error Bound")
    residual_ax.set_ylabel("BR")

    for label,stats in (("DLT",dlt),("DChT",dcht)):
        median,lower,upper,_=stats
        magnitude_ax.plot(orders,median,marker=mode_marker,linewidth=0.8,markersize=1.8,label=label,color=colors[label],zorder=1)
        magnitude_ax.fill_between(orders,lower,upper,color=colors[label],alpha=0.16,linewidth=0,zorder=1)
    magnitude_ax.set_yscale("log")
    magnitude_ax.set_xlim(0,orders[-1])
    magnitude_ax.set_xticks([0,10,20])
    magnitude_ax.set_ylim(1e-10,1)
    magnitude_ax.set_yticks([1e-9,1e-5,1e-1])
    magnitude_ax.set_xlabel(axis_order_label)
    magnitude_ax.set_ylabel("Normalized $|c|$")

    entropy_ax.plot(orders,dlt[3],marker=mode_marker,linewidth=0.8,markersize=1.8,label="DLT",color=colors["DLT"],zorder=1)
    entropy_ax.plot(orders,dcht[3],marker=mode_marker,linewidth=0.8,markersize=1.8,label="DChT",color=colors["DChT"],zorder=1)
    entropy_ax.set_xlim(0,orders[-1])
    entropy_ax.set_xticks([0,10,20])
    entropy_ax.set_xlabel(axis_order_label)
    entropy_ax.set_ylabel("Entropy")
    if args.mode_order=="total":
        entropy_ax.set_yscale("symlog",linthresh=1e-3,linscale=0.5)
        entropy_ax.set_ylim(-3e-4,20)
        entropy_ax.set_yticks([0,1e-3,1e-2,1e-1,1,10])
    else:
        entropy_ax.set_yscale("log")

    for ax in axes:
        ax.grid(True,which="major",linestyle="--",linewidth=0.4,alpha=0.35)
        style_plot_axis(ax)
        ax.tick_params(axis="x",which="minor",length=0)
    layout_panel_row(fig,[residual_ax,magnitude_ax,entropy_ax],OUTER_PANEL_ASPECT_RATIOS,bottom_pt=BOTTOM_PT,panel_gap_pt=PANEL_GAP_WITH_YLABEL_PT,left_edge_gap_pt=LEFT_EDGE_GAP_PT,right_edge_gap_pt=RIGHT_EDGE_GAP_PT)
    figure_height_pt=fig.get_figheight()*MPL_POINTS_PER_INCH
    stream_position=residual_ax.get_position()
    stream_gap=STREAM_AXES_GAP_PT/figure_height_pt
    stream_height=(stream_position.height-stream_gap)/2
    residual_ax.set_position([stream_position.x0,stream_position.y0,stream_position.width,stream_height])
    coefficient_ax.set_position([stream_position.x0,stream_position.y0+stream_height+stream_gap,stream_position.width,stream_height])
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    br_label=residual_ax.yaxis.label
    br_label_bbox=br_label.get_window_extent(renderer)
    br_label_center_x=fig.transFigure.inverted().transform(((br_label_bbox.x0+br_label_bbox.x1)/2,0))[0]
    br_label.set_visible(False)
    fig.text(br_label_center_x,(coefficient_ax.get_position().y1+residual_ax.get_position().y0)/2,"BR",ha="center",va="center",rotation=90)
    align_xlabels(fig,[residual_ax,magnitude_ax,entropy_ax])
    add_panel_captions(fig,[residual_ax,magnitude_ax,entropy_ax],["(a) BR Breakdown","(b) Magnitudes","(c) Entropy"])
    transform_handles=[Line2D([0],[0],color=colors[label],linewidth=0.8,label=label) for label in ("DLT","DChT")]
    stream_handles=[Line2D([0],[0],color="black",linestyle=linestyle,linewidth=0.8,label=label) for label,linestyle in (("Coefficients","-"),("Residuals","--"))]
    axes_top=max(coefficient_ax.get_position().y1,magnitude_ax.get_position().y1,entropy_ax.get_position().y1)
    fig.legend(handles=transform_handles+stream_handles,loc="lower center",bbox_to_anchor=(0.5,axes_top+DEFAULT_LEGEND_GAP_PT/figure_height_pt),bbox_transform=fig.transFigure,ncol=4,frameon=False,borderpad=0,borderaxespad=0,**legend_spacing_kwargs(9.6))
    os.makedirs(os.path.dirname(args.output) or ".",exist_ok=True)
    save_paper_png_vertical_crop(fig,args.output,spec,edge_gap_pt=VERTICAL_EDGE_GAP_PT)
    plt.close(fig)
    os.makedirs(os.path.dirname(args.csv_output) or ".",exist_ok=True)
    with open(args.csv_output,"w",newline="") as file_obj:
        writer=csv.writer(file_obj,lineterminator="\n")
        writer.writerow([order_label.title(),"DLT Median Normalized Magnitude","DChT Median Normalized Magnitude","DLT Q1","DLT Q3","DChT Q1","DChT Q3","DLT Entropy","DChT Entropy"])
        for order in orders:
            writer.writerow([order,dlt[0][order],dcht[0][order],dlt[1][order],dlt[2][order],dcht[1][order],dcht[2][order],dlt[3][order],dcht[3][order]])
    print(f"sampled_elements={values.shape[0]}")
    print(f"element_shape={header.c_element_shape}")
    for order in orders:
        print(f"order={order} dlt_median={dlt[0][order]:.8e} dcht_median={dcht[0][order]:.8e} dlt_entropy={dlt[3][order]:.6f} dcht_entropy={dcht[3][order]:.6f}")

if __name__=="__main__":
    main()
