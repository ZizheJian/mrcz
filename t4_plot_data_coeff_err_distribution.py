# 功能：比较一个NekRS标量field、对应变换系数以及误差的概率分布。
# 程序将三组分布画在同一张1x3图片中，用于观察变换系数和误差是否比原始数据
# 更集中在零附近。
#
# 使用方法：
# 1. 将data_path设为目标NekRS field文件。
# 2. 预先准备data_path+"_coeff"和data_path+"_err"两个float64二进制文件；每个文件都应包含header.file_element_num*header.element_size个值，并与原始field逐点对应。
# 3. 当前代码通过跳过坐标X和速度U来读取后面的第一个标量field；更换数据时应先确认header.variables中的field顺序与此假设一致。
# 4. 根据目标数据调整三个add_hist调用中的横轴范围和bin_width。
# 5. 运行python3 t4_plot_data_coeff_err_distribution.py。程序会打印各bin的中心和概率，并生成fig/data_coeff_err_distribution.png。
#
# 本程序假定原始数据、_coeff和_err文件均已存在且内容完整。

import os
import numpy as np
from matplotlib import pyplot as plt
from py_funcs.read_nekrs_header import read_nekrs_header

data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"

def read_exact(file_obj,dtype,count,label):
    values=np.fromfile(file_obj,dtype=dtype,count=count)
    if values.size!=count:
        raise ValueError(f"{label}数据不完整：期望{count}个值，实际读取{values.size}个")
    return values

with open(data_path,"rb") as f:
    header=read_nekrs_header(f)
    value_count=header.file_element_num*header.element_size
    file_dtype=np.dtype(header.np_dtype).newbyteorder(header.endian_mode)
    f.seek(header.file_element_num*4,1) # element_map
    f.seek(value_count*3*header.word_size,1) # X
    f.seek(value_count*3*header.word_size,1) # U
    data=read_exact(f,file_dtype,value_count,"原始场")

with open(f"{data_path}_coeff","rb") as f:
    coeff=read_exact(f,np.float64,value_count,"系数")

with open(f"{data_path}_err","rb") as f:
    err=read_exact(f,np.float64,value_count,"误差")

def add_hist(ax,data,plot_min,plot_max,label,bin_width):
    bin_max=(int(np.ceil(plot_max/bin_width-0.5))+0.5)*bin_width
    bin_min=(int(np.floor(plot_min/bin_width+0.5))-0.5)*bin_width
    bins=np.arange(bin_min,bin_max+bin_width,bin_width)
    bins=np.concatenate(([-np.inf],bins,[np.inf]))
    counts,_=np.histogram(data,bins=bins)
    counts=counts[1:-1].astype(np.float64)
    counts/=data.size
    centers=(bins[1:-2]+bins[2:-1])/2
    for i in range(0,len(centers)):
        print(f"{centers[i]} {counts[i]}")
    ax.step(centers,counts,where="mid")
    ax.set_xlim(plot_min,plot_max)
    ax.set_ylim(bottom=0)
    ax.set_xticks(np.linspace(plot_min,plot_max,5))
    ax.ticklabel_format(axis="x",style="sci",useOffset=False)
    ax.set_title(label)

fig,axs=plt.subplots(1,3,figsize=(12,4),sharey=False)

add_hist(axs[0],data,-0.2,0.2,"data",4e-4)
add_hist(axs[1],coeff,-2e-8,2e-8,"coeff",4e-11)
add_hist(axs[2],err,-1e-6,1e-6,"err",2e-9)
os.makedirs("fig",exist_ok=True)
fig.tight_layout()
fig.savefig("fig/data_coeff_err_distribution.png",dpi=300)
plt.close(fig)
