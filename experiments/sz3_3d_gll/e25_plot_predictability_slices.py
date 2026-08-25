import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,Normalize
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch

# XXX 每个CSV中的position行已经对所有被处理的element完成平均。
# XXX 偶数边长取最靠近中心的两个切面并逐点平均；奇数边长直接取中心切面。
experiment_dir=os.path.dirname(os.path.abspath(__file__))
input_cases=[
    ("test0 · T","e25_test0_f00001_T_predictability.csv"),
    ("mixlay0 · T","e25_mixlay0_f00001_T_predictability.csv"),
    ("tc_channel0 · UZ","e25_tc_channel0_f00001_UZ_predictability.csv"),
    ("pb0 · UZ · elements [1335789,1352173)","e25_pb0_f00001_UZ_predictability.csv"),
    ("jicf0 · UZ · elements [578353,594737)","e25_jicf0_f00000_UZ_predictability.csv"),
]
input_cases=[(name,os.path.join(experiment_dir,path)) for name,path in input_cases]
output_path=os.path.join(experiment_dir,"e25_oracle_central_slices.png")
cross_element_output_path=os.path.join(experiment_dir,"e25_cross_element_central_slices.png")
best_direction_output_path=os.path.join(experiment_dir,"e25_best_prediction_direction_slices.png")
axis_names="XYZ"
direction_names=["element","X","Y","Z"]

def read_predictability(path):
    oracle_values={}
    element_win_ratios={}
    direction_values={name:{} for name in direction_names}
    with open(path,newline="") as f:
        for row in csv.DictReader(f):
            if row["scope"]!="position":
                continue
            indices=(int(row["i0"]),int(row["i1"]),int(row["i2"]))
            if row["predictor"]=="oracle_best":
                oracle_values[indices]=float(row["NMAE"])
                element_win_ratios[indices]=float(row["element_win_ratio"])
            if row["predictor"] in direction_values:
                direction_values[row["predictor"]][indices]=float(row["NMAE"])
    if not oracle_values:
        raise ValueError(f"{path}中没有position/oracle_best结果")
    shape=tuple(max(indices[dim] for indices in oracle_values)+1 for dim in range(3))
    oracle_data=np.empty(shape)
    element_win_data=np.empty(shape)
    direction_data=np.full((len(direction_names),*shape),np.inf)
    for indices,value in oracle_values.items():
        oracle_data[indices]=value
        element_win_data[indices]=element_win_ratios[indices]
    for direction_id,name in enumerate(direction_names):
        for indices,value in direction_values[name].items():
            direction_data[(direction_id,*indices)]=value
    return oracle_data,element_win_data,direction_data

def central_slice(data,normal_dim):
    n=data.shape[normal_dim]
    center_indices=[n//2] if n%2==1 else [n//2-1,n//2]
    result=np.mean(np.take(data,center_indices,axis=normal_dim),axis=normal_dim)
    remaining_dims=[dim for dim in range(3) if dim!=normal_dim]
    vertical_dim,horizontal_dim=remaining_dims
    return result,center_indices,vertical_dim,horizontal_dim

datasets=[(name,*read_predictability(path)) for name,path in input_cases]

def plot_slices(data_id,title,colorbar_label,path,log_scale):
    fig,axes=plt.subplots(len(datasets),3,figsize=(10,3*len(datasets)),constrained_layout=True)
    for row_id,dataset in enumerate(datasets):
        name=dataset[0]
        data=dataset[data_id]
        norm=LogNorm(vmin=data.min(),vmax=data.max()) if log_scale else Normalize(vmin=data.min(),vmax=data.max())
        last_image=None
        for normal_dim in range(3):
            ax=axes[row_id,normal_dim]
            plane,center_indices,vertical_dim,horizontal_dim=central_slice(data,normal_dim)
            aspect="auto" if min(plane.shape)==1 else "equal"
            last_image=ax.imshow(plane,origin="lower",cmap="viridis",norm=norm,interpolation="nearest",aspect=aspect)
            layer_text="+".join(str(index) for index in center_indices)
            ax.set_title(f"normal {axis_names[normal_dim]}, layer {layer_text}")
            ax.set_xlabel(f"{axis_names[horizontal_dim]} index")
            ax.set_ylabel(f"{axis_names[vertical_dim]} index")
            ax.set_xticks(range(plane.shape[1]))
            ax.set_yticks(range(plane.shape[0]))
            if normal_dim==0:
                ax.text(-0.38,0.5,f"{name}\n{data.shape[0]}×{data.shape[1]}×{data.shape[2]}",transform=ax.transAxes,ha="center",va="center",rotation=90,fontweight="bold")
        colorbar=fig.colorbar(last_image,ax=axes[row_id,:],location="right",shrink=0.78,pad=0.015)
        colorbar.set_label(colorbar_label)
    fig.suptitle(title)
    fig.savefig(path,dpi=300,bbox_inches="tight")
    plt.close(fig)
    print(os.path.abspath(path))

plot_slices(1,"Oracle-best normalized prediction error on central element slices","Oracle-best NMAE (independent log scale)",output_path,True)
plot_slices(2,"Where cross-element prediction strictly beats every spatial direction","Cross-element unique-win ratio",cross_element_output_path,False)

direction_colors=["#7b3294","#008837","#c51b7d","#2c7bb6"]
direction_cmap=ListedColormap(direction_colors)
direction_norm=BoundaryNorm(np.arange(-0.5,len(direction_names)+0.5),direction_cmap.N)
fig,axes=plt.subplots(len(datasets),3,figsize=(10,3*len(datasets)),constrained_layout=True)
for row_id,(name,oracle_data,element_win_data,direction_data) in enumerate(datasets):
    for normal_dim in range(3):
        ax=axes[row_id,normal_dim]
        n=oracle_data.shape[normal_dim]
        center_indices=[n//2] if n%2==1 else [n//2-1,n//2]
        direction_planes=np.mean(np.take(direction_data,center_indices,axis=normal_dim+1),axis=normal_dim+1)
        best_direction=np.argmin(direction_planes,axis=0)
        remaining_dims=[dim for dim in range(3) if dim!=normal_dim]
        vertical_dim,horizontal_dim=remaining_dims
        aspect="auto" if min(best_direction.shape)==1 else "equal"
        ax.imshow(best_direction,origin="lower",cmap=direction_cmap,norm=direction_norm,interpolation="nearest",aspect=aspect)
        layer_text="+".join(str(index) for index in center_indices)
        ax.set_title(f"normal {axis_names[normal_dim]}, layer {layer_text}")
        ax.set_xlabel(f"{axis_names[horizontal_dim]} index")
        ax.set_ylabel(f"{axis_names[vertical_dim]} index")
        ax.set_xticks(range(best_direction.shape[1]))
        ax.set_yticks(range(best_direction.shape[0]))
        if normal_dim==0:
            ax.text(-0.38,0.5,f"{name}\n{oracle_data.shape[0]}×{oracle_data.shape[1]}×{oracle_data.shape[2]}",transform=ax.transAxes,ha="center",va="center",rotation=90,fontweight="bold")
fig.suptitle("Direction with the lowest element-averaged prediction MAE",y=1.025)
fig.legend(handles=[Patch(facecolor=direction_colors[i],label=direction_names[i]) for i in range(len(direction_names))],loc="outside upper center",ncol=4)
fig.savefig(best_direction_output_path,dpi=300,bbox_inches="tight")
plt.close(fig)
print(os.path.abspath(best_direction_output_path))
