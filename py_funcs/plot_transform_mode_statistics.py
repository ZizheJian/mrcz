# 比较DLT和DChT在不同tensor-product mode order上的系数幅值与量化符号熵。
# 默认从pb0的P field中均匀抽取多个连续element block，避免扫描完整数据文件。

import argparse
import csv
import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from py_funcs.read_nekrs_header import read_nekrs_header

DEFAULT_DATA_PATH="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001"
DEFAULT_OUTPUT="论文/figs/6/transform-mode-statistics.pdf"
DEFAULT_CSV_OUTPUT="results/RD/transform_mode_statistics.csv"
DEFAULT_STREAM_ROOT="results/fig9/stream_breakdown"

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

    plt.rcParams.update({"font.family":"serif","font.serif":["Nimbus Roman"],"mathtext.fontset":"stix","font.size":8.4,"axes.labelsize":8.4,"axes.titlesize":8.4,"xtick.labelsize":8.4,"ytick.labelsize":8.4,"legend.fontsize":8.4,"pdf.fonttype":42,"ps.fonttype":42})
    fig=plt.figure(figsize=(3.67,1.55))
    outer_grid=fig.add_gridspec(1,3,left=0.115,right=0.995,bottom=0.31,top=0.84,wspace=0.62)
    stream_grid=outer_grid[0,0].subgridspec(2,1,hspace=0.10)
    coefficient_ax=fig.add_subplot(stream_grid[0])
    residual_ax=fig.add_subplot(stream_grid[1],sharex=coefficient_ax)
    magnitude_ax=fig.add_subplot(outer_grid[0,1])
    entropy_ax=fig.add_subplot(outer_grid[0,2])
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
    residual_ax.set_xlabel("Relative Error Bound",labelpad=1)
    fig.text(0.025,(coefficient_ax.get_position().y1+residual_ax.get_position().y0)/2,"BR",ha="center",va="center",rotation=90,fontsize=8.4)

    for label,stats in (("DLT",dlt),("DChT",dcht)):
        median,lower,upper,_=stats
        magnitude_ax.plot(orders,median,marker=mode_marker,linewidth=0.8,markersize=1.8,label=label,color=colors[label],zorder=1)
        magnitude_ax.fill_between(orders,lower,upper,color=colors[label],alpha=0.16,linewidth=0,zorder=1)
    magnitude_ax.set_yscale("log")
    magnitude_ax.set_xlabel(axis_order_label,labelpad=1)
    magnitude_ax.set_ylabel("Normalized $|c|$",labelpad=1)

    entropy_ax.plot(orders,dlt[3],marker=mode_marker,linewidth=0.8,markersize=1.8,label="DLT",color=colors["DLT"],zorder=1)
    entropy_ax.plot(orders,dcht[3],marker=mode_marker,linewidth=0.8,markersize=1.8,label="DChT",color=colors["DChT"],zorder=1)
    entropy_ax.set_xlabel(axis_order_label,labelpad=1)
    entropy_ax.set_ylabel("Entropy (bits/value)",labelpad=1)
    if args.mode_order=="total":
        entropy_ax.set_yscale("symlog",linthresh=1e-3,linscale=0.5)
        entropy_ax.set_ylim(-3e-4,20)
        entropy_ax.set_yticks([0,1e-3,1e-2,1e-1,1,10])
    else:
        entropy_ax.set_yscale("log")

    for ax in axes:
        ax.grid(True,which="both",linestyle="--",linewidth=0.4,alpha=0.35)
        ax.tick_params(direction="in",length=2,width=0.5,pad=1)
        ax.tick_params(axis="x",pad=2.5)
        for spine in ax.spines.values():spine.set_linewidth(0.6)
    transform_handles=[Line2D([0],[0],color=colors[label],linewidth=0.8,label=label) for label in ("DLT","DChT")]
    stream_handles=[Line2D([0],[0],color="black",linestyle=linestyle,linewidth=0.8,label=label) for label,linestyle in (("Coefficients","-"),("Residuals","--"))]
    fig.legend(handles=transform_handles+stream_handles,loc="upper center",bbox_to_anchor=(0.5,0.99),ncol=4,frameon=False,handlelength=1.2,handletextpad=0.25,columnspacing=0.65,borderaxespad=0)
    panel_labels=[("(a) BR Breakdown",coefficient_ax,0.04),("(b) Magnitudes",magnitude_ax,0.04),("(c) Entropy",entropy_ax,0.04)]
    for label,ax,y in panel_labels:
        position=ax.get_position()
        fig.text((position.x0+position.x1)/2,y,label,ha="center",va="bottom",fontsize=8.4,linespacing=0.9)
    os.makedirs(os.path.dirname(args.output) or ".",exist_ok=True)
    for ax in axes:ax.set_rasterization_zorder(1.5)
    fig.savefig(args.output,dpi=300)
    fig.savefig(os.path.splitext(args.output)[0]+".png",dpi=300)
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
