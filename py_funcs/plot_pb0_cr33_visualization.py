import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap,Normalize

plt.rcParams.update({"font.family":"serif","font.serif":["Nimbus Roman"],"mathtext.fontset":"stix","font.size":8.0,"axes.labelsize":8.0,"axes.titlesize":8.0,"xtick.labelsize":7.0,"ytick.labelsize":7.0,"pdf.fonttype":42,"ps.fonttype":42})

methods=["MRCZ","FDWZ","SZ3-4D","SZ2-3D-fused","QoZ-4D","ZFP-3D-fused","SPERR-2D"]
labels={"MRCZ":"MRCZ","FDWZ":"FDWZ","SZ3-4D":"SZ3","SZ2-3D-fused":"SZ2","QoZ-4D":"QoZ","ZFP-3D-fused":"ZFP","SPERR-2D":"SPERR"}
crop_fractions=(1/6,1/2,1/3,2/3)

def load_slice(path):
    data=np.load(path)
    points=data["points"]
    values=data["values"]
    offsets=data["offsets"]
    connectivity=data["connectivity"]
    if offsets[0]!=0 or offsets[-1]!=len(connectivity) or np.any(np.diff(offsets)!=3):raise RuntimeError(f"{path}不只包含三角形")
    triangles=connectivity.reshape(-1,3)
    return points,values,offsets,connectivity,triangles

def validate_geometry(reference_points,reference_offsets,reference_connectivity,points,offsets,connectivity,method):
    if reference_points.shape!=points.shape or reference_offsets.shape!=offsets.shape or reference_connectivity.shape!=connectivity.shape:raise RuntimeError(f"{method}的切片网格形状与原图不一致")
    if not np.array_equal(reference_offsets,offsets) or not np.array_equal(reference_connectivity,connectivity):raise RuntimeError(f"{method}的切片连接关系与原图不一致")
    if not np.allclose(reference_points,points,rtol=0,atol=1e-7):raise RuntimeError(f"{method}的切片坐标与原图不一致")

def style_axis(axis):
    axis.set_aspect("equal")
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():spine.set_linewidth(0.45)

def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input_dir=os.path.join(project_root,"results","visualization","cr33_slices")
    output_dir=os.path.join(project_root,"results","visualization")
    original_points,original_values,original_offsets,original_connectivity,triangles=load_slice(os.path.join(input_dir,"pb0_T_x0_original.npz"))
    triangulation=mtri.Triangulation(original_points[:,1],original_points[:,2],triangles)
    horizontal_min=float(original_points[:,1].min())
    horizontal_max=float(original_points[:,1].max())
    vertical_min=float(original_points[:,2].min())
    vertical_max=float(original_points[:,2].max())
    horizontal_range=horizontal_max-horizontal_min
    vertical_range=vertical_max-vertical_min
    crop=(horizontal_min+crop_fractions[0]*horizontal_range,horizontal_min+crop_fractions[1]*horizontal_range,vertical_min+crop_fractions[2]*vertical_range,vertical_min+crop_fractions[3]*vertical_range)
    errors={}
    for method in methods:
        points,reconstructed_values,offsets,connectivity,_=load_slice(os.path.join(input_dir,f"pb0_T_x0_{method}.npz"))
        validate_geometry(original_points,original_offsets,original_connectivity,points,offsets,connectivity,method)
        errors[method]=reconstructed_values-original_values
    y_coordinates=original_points[:,1]
    z_coordinates=original_points[:,2]
    crop_mask=(y_coordinates>=crop[0])&(y_coordinates<=crop[1])&(z_coordinates>=crop[2])&(z_coordinates<=crop[3])
    if not np.any(crop_mask):raise RuntimeError("误差裁剪区域不包含切片点")
    error_limit=max(float(np.abs(values[crop_mask]).max()) for values in errors.values())
    if error_limit<=0:raise RuntimeError("所有重建误差均为零")
    color_threshold=0.01
    if color_threshold>=error_limit:raise RuntimeError("颜色阈值必须小于共享误差范围")
    threshold_offset=color_threshold/(2*error_limit)
    error_colormap=LinearSegmentedColormap.from_list("neutral-gray-diverging",[(0.0,"#000000"),(0.5-threshold_offset,"#2166ac"),(0.5,"#d9d9d9"),(0.5+threshold_offset,"#b2182b"),(1.0,"#000000")],N=4097)
    error_norm=Normalize(vmin=-error_limit,vmax=error_limit,clip=True)

    figure_width=7.16
    figure_height=2.55
    small_side=0.98
    bottom=0.25
    row_gap=0.17
    large_side=2*small_side+row_gap
    original_left=0.06
    original_bottom=bottom
    error_left=2.52
    panel_gap=0.07
    figure=plt.figure(figsize=(figure_width,figure_height))
    original_axis=figure.add_axes([original_left/figure_width,original_bottom/figure_height,large_side/figure_width,large_side/figure_height])
    original_image=original_axis.tripcolor(triangulation,original_values,shading="gouraud",cmap="viridis",rasterized=True)
    original_axis.add_patch(Rectangle((crop[0],crop[2]),crop[1]-crop[0],crop[3]-crop[2],fill=False,edgecolor="#d62728",linewidth=0.8))
    original_axis.set_xlim(horizontal_min,horizontal_max)
    original_axis.set_ylim(vertical_min,vertical_max)
    style_axis(original_axis)
    original_axis.set_xlabel("(a) Original 7T",labelpad=1)
    original_colorbar_axis=figure.add_axes([2.25/figure_width,original_bottom/figure_height,0.055/figure_width,large_side/figure_height])
    original_colorbar=figure.colorbar(original_image,cax=original_colorbar_axis)
    original_colorbar.ax.tick_params(length=2,width=0.4,pad=1)

    error_axes=[]
    for index,method in enumerate(methods):
        row=0 if index<4 else 1
        column=index if index<4 else index-4
        left=error_left+column*(small_side+panel_gap)
        panel_bottom=bottom+small_side+row_gap if row==0 else bottom
        axis=figure.add_axes([left/figure_width,panel_bottom/figure_height,small_side/figure_width,small_side/figure_height])
        image=axis.tripcolor(triangulation,errors[method],shading="gouraud",cmap=error_colormap,norm=error_norm,rasterized=True)
        axis.set_xlim(crop[0],crop[1])
        axis.set_ylim(crop[2],crop[3])
        style_axis(axis)
        axis.set_xlabel(f"({chr(ord('b')+index)}) {labels[method]}",labelpad=1)
        error_axes.append(axis)

    error_colorbar_bottom=bottom
    error_colorbar_axis=figure.add_axes([6.78/figure_width,error_colorbar_bottom/figure_height,0.055/figure_width,large_side/figure_height])
    error_colorbar=figure.colorbar(image,cax=error_colorbar_axis)
    error_colorbar.set_ticks([-error_limit,-color_threshold,0,color_threshold,error_limit])
    error_colorbar.set_ticklabels([f"{-error_limit:.3f}",f"{-color_threshold:.2f}","0",f"{color_threshold:.2f}",f"{error_limit:.3f}"])
    error_colorbar.ax.axhline(-color_threshold,color="black",linewidth=0.35)
    error_colorbar.ax.axhline(color_threshold,color="black",linewidth=0.35)
    error_colorbar.ax.tick_params(length=2,width=0.4,pad=1)

    pdf_output=os.path.join(project_root,"论文","figs","7","pb0-T-cr33-visualization.pdf")
    png_output=os.path.join(output_dir,"pb0_T_cr33_visualization.png")
    figure.savefig(pdf_output,format="pdf",dpi=300,bbox_inches="tight",pad_inches=0.02)
    figure.savefig(png_output,format="png",dpi=300,bbox_inches="tight",pad_inches=0.02)
    plt.close(figure)
    print(f"Shared symmetric linear error range: [{-error_limit:.9g},{error_limit:.9g}]")
    print(pdf_output)
    print(png_output)

if __name__=="__main__":main()
