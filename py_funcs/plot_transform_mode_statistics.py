# 比较DLT和DChT在不同tensor-product mode order上的系数幅值与量化符号熵。
# 默认从pb0的P field中均匀抽取多个连续element block，避免扫描完整数据文件。

import argparse
import csv
import os
import sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from py_funcs.read_nekrs_header import read_nekrs_header

DEFAULT_DATA_PATH="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001"
DEFAULT_OUTPUT="论文/figs/6/transform-mode-statistics.pdf"
DEFAULT_CSV_OUTPUT="results/RD/transform_mode_statistics.csv"

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--data-path",default=DEFAULT_DATA_PATH)
    parser.add_argument("--output",default=DEFAULT_OUTPUT)
    parser.add_argument("--csv-output",default=DEFAULT_CSV_OUTPUT)
    parser.add_argument("--block-count",type=int,default=64)
    parser.add_argument("--block-elements",type=int,default=128)
    parser.add_argument("--relative-eb",type=float,default=1e-4)
    parser.add_argument("--absolute-eb",type=float,default=0.060799235524628958)
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

def mode_statistics(coeff,value_range,coeff_eb):
    n=coeff.shape[1]
    indices=np.indices((n,n,n))
    mode_order=np.maximum.reduce(indices)
    medians=[]
    lower=[]
    upper=[]
    entropies=[]
    for order in range(n):
        values=np.abs(coeff[:,mode_order==order]).reshape(-1)/value_range
        medians.append(np.median(values))
        lower.append(np.quantile(values,0.25))
        upper.append(np.quantile(values,0.75))
        raw_values=coeff[:,mode_order==order].reshape(-1)
        quantized=np.rint(raw_values/(2*coeff_eb)).astype(np.int64)
        _,counts=np.unique(quantized,return_counts=True)
        probabilities=counts/counts.sum()
        entropies.append(-np.sum(probabilities*np.log2(probabilities)))
    return np.asarray(medians),np.asarray(lower),np.asarray(upper),np.asarray(entropies)

def main():
    args=parse_args()
    values,header=read_sample(args.data_path,args.block_count,args.block_elements)
    value_range=args.absolute_eb/args.relative_eb
    matrix_dlt=build_forward_matrix(header.c_element_shape[0],"DLT")
    matrix_dcht=build_forward_matrix(header.c_element_shape[0],"DChT")
    coeff_dlt=transform_elements(values,matrix_dlt)
    coeff_dcht=transform_elements(values,matrix_dcht)
    dlt=mode_statistics(coeff_dlt,value_range,args.absolute_eb*0.05)
    dcht=mode_statistics(coeff_dcht,value_range,args.absolute_eb*0.05)
    orders=np.arange(header.c_element_shape[0])

    plt.rcParams.update({"font.size":10,"axes.labelsize":10,"xtick.labelsize":9,"ytick.labelsize":9,"legend.fontsize":9})
    fig,axes=plt.subplots(2,1,figsize=(3.5,4.8))
    colors={"DLT":"#ff5a5a","DChT":"#3b7ded"}
    for label,stats in (("DLT",dlt),("DChT",dcht)):
        median,lower,upper,_=stats
        axes[0].plot(orders,median,marker="o",linewidth=1.8,markersize=4,label=label,color=colors[label])
        axes[0].fill_between(orders,lower,upper,color=colors[label],alpha=0.16,linewidth=0)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Maximum mode order")
    axes[0].set_ylabel("Normalized $|c|$")
    axes[0].set_title("(a) Coefficient Magnitude")
    axes[0].grid(True,which="both",linestyle="--",alpha=0.35)
    axes[0].legend(frameon=False)

    axes[1].plot(orders,dlt[3],marker="o",linewidth=1.8,markersize=4,label="DLT",color=colors["DLT"])
    axes[1].plot(orders,dcht[3],marker="o",linewidth=1.8,markersize=4,label="DChT",color=colors["DChT"])
    axes[1].set_xlabel("Maximum mode order")
    axes[1].set_ylabel("Entropy (bits/coefficient)")
    axes[1].set_yscale("log")
    axes[1].set_title("(b) Quantized Entropy")
    axes[1].grid(True,which="both",linestyle="--",alpha=0.35)
    axes[1].legend(frameon=False)

    fig.tight_layout(pad=0.5,h_pad=1.0)
    os.makedirs(os.path.dirname(args.output) or ".",exist_ok=True)
    fig.savefig(args.output,bbox_inches="tight")
    fig.savefig(os.path.splitext(args.output)[0]+".png",dpi=300,bbox_inches="tight")
    os.makedirs(os.path.dirname(args.csv_output) or ".",exist_ok=True)
    with open(args.csv_output,"w",newline="") as file_obj:
        writer=csv.writer(file_obj,lineterminator="\n")
        writer.writerow(["Maximum Mode Order","DLT Median Normalized Magnitude","DChT Median Normalized Magnitude","DLT Q1","DLT Q3","DChT Q1","DChT Q3","DLT Entropy","DChT Entropy"])
        for order in orders:
            writer.writerow([order,dlt[0][order],dcht[0][order],dlt[1][order],dlt[2][order],dcht[1][order],dcht[2][order],dlt[3][order],dcht[3][order]])
    print(f"sampled_elements={values.shape[0]}")
    print(f"element_shape={header.c_element_shape}")
    for order in orders:
        print(f"order={order} dlt_median={dlt[0][order]:.8e} dcht_median={dcht[0][order]:.8e} dlt_entropy={dlt[3][order]:.6f} dcht_entropy={dcht[3][order]:.6f}")

if __name__=="__main__":
    main()
