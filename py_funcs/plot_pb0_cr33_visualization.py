import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap,Normalize
from plot_style import DEFAULT_PANEL_CAPTION_GAP_PT,PaperFigureSpec,METHOD_LABELS,RD_METHODS,apply_paper_style,save_paper_png_vertical_crop,style_plot_axis

apply_paper_style()

methods=RD_METHODS[:-1]
labels=METHOD_LABELS
crop_fractions=(1/6,1/2,1/3,2/3)
SMALL_SIDE_PT=70.0
ROW_GAP_PT=12.0
LARGE_SIDE_PT=2*SMALL_SIDE_PT+ROW_GAP_PT
BOTTOM_PT=16.0
ORIGINAL_LEFT_PT=4.0
ORIGINAL_COLORBAR_LEFT_PT=162.0
ERROR_LEFT_PT=181.0
PANEL_GAP_PT=5.0
ERROR_COLORBAR_LEFT_PT=485.0
COLORBAR_WIDTH_PT=4.0
VERTICAL_EDGE_GAP_PT=3.0

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
    style_plot_axis(axis,x_label_pad=DEFAULT_PANEL_CAPTION_GAP_PT)

def add_axis_pt(figure,spec,left,bottom,width,height):
    return figure.add_axes([left/spec.width_pt,bottom/spec.height_pt,width/spec.width_pt,height/spec.height_pt])

def style_colorbar(colorbar):
    colorbar.ax.tick_params(axis="y",direction="in",length=1.8,width=0.45,pad=1.0)
    for spine in colorbar.ax.spines.values():spine.set_linewidth(0.5)

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

    spec=PaperFigureSpec.text(1)
    figure=plt.figure(figsize=spec.size_in,dpi=spec.dpi)
    original_axis=add_axis_pt(figure,spec,ORIGINAL_LEFT_PT,BOTTOM_PT,LARGE_SIDE_PT,LARGE_SIDE_PT)
    original_image=original_axis.tripcolor(triangulation,original_values,shading="gouraud",cmap="viridis",rasterized=True)
    original_axis.add_patch(Rectangle((crop[0],crop[2]),crop[1]-crop[0],crop[3]-crop[2],fill=False,edgecolor="#d62728",linewidth=0.8))
    original_axis.set_xlim(horizontal_min,horizontal_max)
    original_axis.set_ylim(vertical_min,vertical_max)
    style_axis(original_axis)
    original_axis.set_xlabel("(a) Original 7T")
    original_colorbar_axis=add_axis_pt(figure,spec,ORIGINAL_COLORBAR_LEFT_PT,BOTTOM_PT,COLORBAR_WIDTH_PT,LARGE_SIDE_PT)
    original_colorbar=figure.colorbar(original_image,cax=original_colorbar_axis)
    style_colorbar(original_colorbar)

    error_axes=[]
    for index,method in enumerate(methods):
        row=0 if index<4 else 1
        column=index if index<4 else index-4
        left=ERROR_LEFT_PT+column*(SMALL_SIDE_PT+PANEL_GAP_PT)
        panel_bottom=BOTTOM_PT+SMALL_SIDE_PT+ROW_GAP_PT if row==0 else BOTTOM_PT
        axis=add_axis_pt(figure,spec,left,panel_bottom,SMALL_SIDE_PT,SMALL_SIDE_PT)
        image=axis.tripcolor(triangulation,errors[method],shading="gouraud",cmap=error_colormap,norm=error_norm,rasterized=True)
        axis.set_xlim(crop[0],crop[1])
        axis.set_ylim(crop[2],crop[3])
        style_axis(axis)
        axis.set_xlabel(f"({chr(ord('b')+index)}) {labels[method]}")
        error_axes.append(axis)

    error_colorbar_axis=add_axis_pt(figure,spec,ERROR_COLORBAR_LEFT_PT,BOTTOM_PT,COLORBAR_WIDTH_PT,LARGE_SIDE_PT)
    error_colorbar=figure.colorbar(image,cax=error_colorbar_axis)
    error_colorbar.set_ticks([-error_limit,-color_threshold,0,color_threshold,error_limit])
    error_colorbar.set_ticklabels([f"{-error_limit:.3f}",f"{-color_threshold:.2f}","0",f"{color_threshold:.2f}",f"{error_limit:.3f}"])
    error_colorbar.ax.axhline(-color_threshold,color="black",linewidth=0.35)
    error_colorbar.ax.axhline(color_threshold,color="black",linewidth=0.35)
    style_colorbar(error_colorbar)

    result_output=os.path.join(output_dir,"pb0_T_cr33_visualization.png")
    paper_output=os.path.join(project_root,"论文","figs","7","visualization-cr33.png")
    save_paper_png_vertical_crop(figure,result_output,spec,edge_gap_pt=VERTICAL_EDGE_GAP_PT)
    save_paper_png_vertical_crop(figure,paper_output,spec,edge_gap_pt=VERTICAL_EDGE_GAP_PT)
    plt.close(figure)
    print(f"Shared symmetric linear error range: [{-error_limit:.9g},{error_limit:.9g}]")
    print(result_output)
    print(paper_output)

if __name__=="__main__":main()
