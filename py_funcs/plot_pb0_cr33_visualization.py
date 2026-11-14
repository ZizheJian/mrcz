import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.patches import Rectangle
from plot_style import METHOD_LABELS,RD_METHODS,apply_paper_style
from paperfig import DEFAULT_PANEL_CAPTION_GAP_PT,PaperFigureSpec,TEX_POINTS_PER_INCH,save_paper_png,style_plot_axis

apply_paper_style()

methods=RD_METHODS[:-1]
labels=METHOD_LABELS
crop_fractions=(1/6,1/2,1/3,2/3)
MASTER_DPI=1200
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
COLOR_THRESHOLD=0.01

def parse_args(project_root):
    parser=argparse.ArgumentParser()
    parser.add_argument("--input-dir",default=os.path.join(project_root,"results","visualization","cr33_slices"))
    parser.add_argument("--master-dir",default=os.path.join(project_root,"results","visualization","cr33_panel_masters"))
    parser.add_argument("--result-output",default=os.path.join(project_root,"results","visualization","pb0_T_cr33_visualization.png"))
    parser.add_argument("--paper-output",default=os.path.join(project_root,"overleaf","figs","7","visualization-cr33.png"))
    parser.add_argument("--master-dpi",type=int,default=MASTER_DPI)
    parser.add_argument("--compose-only",action="store_true")
    return parser.parse_args()

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

def make_error_colormap(error_limit,color_threshold):
    threshold_offset=color_threshold/(2*error_limit)
    return LinearSegmentedColormap.from_list("neutral-gray-diverging",[(0.0,"#000000"),(0.5-threshold_offset,"#2166ac"),(0.5,"#d9d9d9"),(0.5+threshold_offset,"#b2182b"),(1.0,"#000000")],N=4097)

def load_source_data(input_dir):
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
    if COLOR_THRESHOLD>=error_limit:raise RuntimeError("颜色阈值必须小于共享误差范围")
    metadata={"master_dpi":None,"methods":list(methods),"crop":list(crop),"original_range":[float(original_values.min()),float(original_values.max())],"error_range":[-error_limit,error_limit],"color_threshold":COLOR_THRESHOLD}
    return triangulation,original_values,errors,(horizontal_min,horizontal_max,vertical_min,vertical_max),metadata

def save_master_panel(path,side_pt,dpi,draw_panel):
    side_in=side_pt/TEX_POINTS_PER_INCH
    figure=plt.figure(figsize=(side_in,side_in),dpi=dpi)
    axis=figure.add_axes([0,0,1,1])
    draw_panel(axis)
    axis.set_axis_off()
    figure.savefig(path,format="png",dpi=dpi,bbox_inches=None,pad_inches=0,facecolor="white")
    plt.close(figure)

def render_panel_masters(master_dir,master_dpi,triangulation,original_values,errors,bounds,metadata):
    os.makedirs(master_dir,exist_ok=True)
    horizontal_min,horizontal_max,vertical_min,vertical_max=bounds
    crop=metadata["crop"]
    error_limit=metadata["error_range"][1]
    error_colormap=make_error_colormap(error_limit,metadata["color_threshold"])
    error_norm=Normalize(vmin=-error_limit,vmax=error_limit,clip=True)
    def draw_original(axis):
        axis.tripcolor(triangulation,original_values,shading="gouraud",cmap="viridis")
        axis.add_patch(Rectangle((crop[0],crop[2]),crop[1]-crop[0],crop[3]-crop[2],fill=False,edgecolor="#d62728",linewidth=0.8))
        axis.set_xlim(horizontal_min,horizontal_max)
        axis.set_ylim(vertical_min,vertical_max)
    save_master_panel(os.path.join(master_dir,"original.png"),LARGE_SIDE_PT,master_dpi,draw_original)
    for method in methods:
        def draw_error(axis,method=method):
            axis.tripcolor(triangulation,errors[method],shading="gouraud",cmap=error_colormap,norm=error_norm)
            axis.set_xlim(crop[0],crop[1])
            axis.set_ylim(crop[2],crop[3])
        save_master_panel(os.path.join(master_dir,f"error_{method}.png"),SMALL_SIDE_PT,master_dpi,draw_error)
    metadata["master_dpi"]=master_dpi
    with open(os.path.join(master_dir,"metadata.json"),"w",encoding="utf-8",newline="\n") as file_obj:json.dump(metadata,file_obj,indent=2);file_obj.write("\n")

def load_master_metadata(master_dir):
    with open(os.path.join(master_dir,"metadata.json"),encoding="utf-8") as file_obj:metadata=json.load(file_obj)
    if metadata["methods"]!=list(methods):raise RuntimeError("母版中的方法顺序与当前绘图配置不一致")
    return metadata

def style_axis(axis):
    axis.set_aspect("equal")
    axis.set_xticks([])
    axis.set_yticks([])
    style_plot_axis(axis,x_label_gap_pt=DEFAULT_PANEL_CAPTION_GAP_PT)

def add_axis_pt(figure,spec,left,bottom,width,height):
    return figure.add_axes([left/spec.width_pt,bottom/spec.height_pt,width/spec.width_pt,height/spec.height_pt])

def style_colorbar(colorbar):
    colorbar.ax.tick_params(axis="y",direction="in",length=1.8,width=0.45,pad=1.0)
    for spine in colorbar.ax.spines.values():spine.set_linewidth(0.5)

def show_master(axis,path):
    axis.imshow(plt.imread(path),interpolation="lanczos",resample=True)
    style_axis(axis)

def compose_figure(master_dir,metadata,result_output,paper_output):
    error_limit=metadata["error_range"][1]
    color_threshold=metadata["color_threshold"]
    error_colormap=make_error_colormap(error_limit,color_threshold)
    error_norm=Normalize(vmin=-error_limit,vmax=error_limit,clip=True)
    spec=PaperFigureSpec.text(1)
    figure=plt.figure(figsize=spec.size_in,dpi=spec.dpi)
    original_axis=add_axis_pt(figure,spec,ORIGINAL_LEFT_PT,BOTTOM_PT,LARGE_SIDE_PT,LARGE_SIDE_PT)
    show_master(original_axis,os.path.join(master_dir,"original.png"))
    original_axis.set_xlabel("(a) Original 7T")
    original_colorbar_axis=add_axis_pt(figure,spec,ORIGINAL_COLORBAR_LEFT_PT,BOTTOM_PT,COLORBAR_WIDTH_PT,LARGE_SIDE_PT)
    original_mappable=ScalarMappable(norm=Normalize(vmin=metadata["original_range"][0],vmax=metadata["original_range"][1]),cmap="viridis")
    original_colorbar=figure.colorbar(original_mappable,cax=original_colorbar_axis)
    style_colorbar(original_colorbar)
    for index,method in enumerate(methods):
        row=0 if index<4 else 1
        column=index if index<4 else index-4
        left=ERROR_LEFT_PT+column*(SMALL_SIDE_PT+PANEL_GAP_PT)
        panel_bottom=BOTTOM_PT+SMALL_SIDE_PT+ROW_GAP_PT if row==0 else BOTTOM_PT
        axis=add_axis_pt(figure,spec,left,panel_bottom,SMALL_SIDE_PT,SMALL_SIDE_PT)
        show_master(axis,os.path.join(master_dir,f"error_{method}.png"))
        axis.set_xlabel(f"({chr(ord('b')+index)}) {labels[method]}")
    error_colorbar_axis=add_axis_pt(figure,spec,ERROR_COLORBAR_LEFT_PT,BOTTOM_PT,COLORBAR_WIDTH_PT,LARGE_SIDE_PT)
    error_colorbar=figure.colorbar(ScalarMappable(norm=error_norm,cmap=error_colormap),cax=error_colorbar_axis)
    error_colorbar.set_ticks([-error_limit,-color_threshold,0,color_threshold,error_limit])
    error_colorbar.set_ticklabels([f"{-error_limit:.3f}",f"{-color_threshold:.2f}","0",f"{color_threshold:.2f}",f"{error_limit:.3f}"])
    error_colorbar.ax.axhline(-color_threshold,color="black",linewidth=0.35)
    error_colorbar.ax.axhline(color_threshold,color="black",linewidth=0.35)
    style_colorbar(error_colorbar)
    os.makedirs(os.path.dirname(result_output) or ".",exist_ok=True)
    os.makedirs(os.path.dirname(paper_output) or ".",exist_ok=True)
    save_paper_png(figure,result_output,spec)
    save_paper_png(figure,paper_output,spec)
    plt.close(figure)

def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    args=parse_args(project_root)
    if args.master_dpi<300:raise ValueError("master-dpi不能低于论文图的300 dpi")
    if not args.compose_only:
        triangulation,original_values,errors,bounds,metadata=load_source_data(args.input_dir)
        render_panel_masters(args.master_dir,args.master_dpi,triangulation,original_values,errors,bounds,metadata)
    metadata=load_master_metadata(args.master_dir)
    compose_figure(args.master_dir,metadata,args.result_output,args.paper_output)
    print(f"Panel masters: {args.master_dir} ({metadata['master_dpi']} dpi)")
    print(f"Shared symmetric linear error range: [{metadata['error_range'][0]:.9g},{metadata['error_range'][1]:.9g}]")
    print(args.result_output)
    print(args.paper_output)

if __name__=="__main__":main()
