import numpy as np
from matplotlib import pyplot as plt
from py_funcs.read_nekrs_header import read_nekrs_header

data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"

f=open(data_path,"rb")
header=read_nekrs_header(f)
f.seek(header.file_element_num*4,1) # element_map
f.seek(header.file_element_num*header.element_size*3*4,1) # X
f.seek(header.file_element_num*header.element_size*3*4,1) # U
data=np.fromfile(f,dtype=np.float32,count=header.file_element_num*header.element_size)
f.close()

f=open(f"{data_path}_coeff","rb")
coeff=np.fromfile(f,dtype=np.float64,count=header.file_element_num*header.element_size)
f.close()

f=open(f"{data_path}_err","rb")
err=np.fromfile(f,dtype=np.float64,count=header.file_element_num*header.element_size)
f.close()

def add_hist(ax,data,plot_min,plot_max,label,bin_width):
    bin_max=(int(np.ceil(plot_max/bin_width-0.5))+0.5)*bin_width
    bin_min=(int(np.floor(plot_min/bin_width+0.5))-0.5)*bin_width
    bins=np.arange(bin_min,bin_max+bin_width,bin_width)
    bins=np.concatenate(([-np.inf],bins,[np.inf]))
    counts,_=np.histogram(data,bins=bins)
    counts=counts[1:-1].astype(np.float64)
    counts/=(header.file_element_num*header.element_size)
    centers=(bins[1:-2]+bins[2:-1])/2
    for i in range(0,len(centers)):
        print(f"{centers[i]} {counts[i]}")
    input()
    # ax.step(centers[mask],counts[mask],where="mid")
    # ax.set_xlim(plot_min,plot_max)
    # ax.set_ylim(bottom=0)
    # ax.set_xticks(np.arange(plot_min,plot_max+(plot_max-plot_min)/4,(plot_max-plot_min)/4))
    # ax.ticklabel_format(axis="x",style="sci",useOffset=False)
    # ax.set_title(label)

fig,axs=plt.subplots(1,3,figsize=(12,4),sharey=False)

add_hist(axs[0],data,-0.2,0.2,"data",4e-4)
add_hist(axs[1],coeff,-2e-8,2e-8,"coeff",4e-11)
add_hist(axs[2],err,-1e-6,1e-6,"err",2e-9)
plt.savefig("fig/data_coeff_err_distribution.png",dpi=300)