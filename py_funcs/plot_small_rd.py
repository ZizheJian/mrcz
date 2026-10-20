import csv,math,os
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator,MaxNLocator
from plot_style import METHOD_COLORS,PAPER_TEXT_WIDTH_IN,apply_paper_style,save_svg,style_plot_axis

apply_paper_style("mrcz-small-rd",serif_fonts=["STIXGeneral"])

array_layouts=["1D","2D","3D-fused","3D-per-element","4D"]
layouts=array_layouts+["4D+GLL"]
methods=["MRCZ","FDWZ"]+[f"{compressor}-{layout}" for compressor in ["SZ3","SZ2","QoZ","ZFP"] for layout in array_layouts]+[f"MGARD-{layout}" for layout in layouts]+["SPERR-2D","SPERR-3D-fused","SPERR-3D-per-element"]
family_colors={"MRCZ":METHOD_COLORS["MRCZ"],"FDWZ":METHOD_COLORS["FDWZ"],"SZ3":METHOD_COLORS["SZ3-4D"],"SZ2":METHOD_COLORS["SZ2-3D-fused"],"QoZ":METHOD_COLORS["QoZ-4D"],"ZFP":METHOD_COLORS["ZFP-3D-fused"],"MGARD":METHOD_COLORS["MGARD-1D"],"SPERR":METHOD_COLORS["SPERR-2D"]}
layout_styles={"1D":(0,(5,1,1,1)),"2D":"--","3D-fused":"-.","3D-per-element":":","4D":"-","4D+GLL":(0,(3,1,1,1))}
fields=[
    ("1UX","mixlay0/mixlay0.f00001_UX.csv"),("1UY","mixlay0/mixlay0.f00001_UY.csv"),("1P","mixlay0/mixlay0.f00001_P.csv"),("1T","mixlay0/mixlay0.f00001_T.csv"),("1S00","mixlay0/mixlay0.f00001_S00.csv"),("1S01","mixlay0/mixlay0.f00001_S01.csv"),
    ("2UZ","test0/test0.f00001_UZ.csv"),("2T","test0/test0.f00001_T.csv"),("2S00","test0/test0.f00001_S00.csv"),("2S01","test0/test0.f00001_S01.csv"),("3UX","tc_channel0/tc_channel0.f00001_UX.csv"),("3UY","tc_channel0/tc_channel0.f00001_UY.csv"),
    ("3UZ","tc_channel0/tc_channel0.f00001_UZ.csv"),("Dataset 4","mfem_ex0/sol.gf_no_name_field.csv"),("Dataset 5","mfem_ex1/sol.gf_no_name_field.csv"),("6P","nektar_cylinder/Cylinder-Base.vtu_p.csv"),("6U","nektar_cylinder/Cylinder-Base.vtu_u.csv"),("6V","nektar_cylinder/Cylinder-Base.vtu_v.csv")]

def read_rows(path):
    with open(path,newline="") as f:
        return list(csv.DictReader(f))

def read_method_rows(result_dir,relative_path,method,default_rows):
    if not method.startswith("MGARD-"):
        return default_rows,method
    coordinate_type="GLL" if method.endswith("+GLL") else "uniform"
    stored_method="MGARD-4D" if method.endswith("+GLL") else method
    path=os.path.join(result_dir,"mgard_layout_selection",os.path.dirname(relative_path),coordinate_type,os.path.basename(relative_path))
    return read_rows(path),stored_method

def prepare_axes():
    figure,axes=plt.subplots(3,6,figsize=(PAPER_TEXT_WIDTH_IN,3.85))
    return figure,axes.flat

def method_style(method):
    if method in ["MRCZ","FDWZ"]:
        return family_colors[method],"-",0.4
    family,layout=method.split("-",1)
    return family_colors[family],layout_styles[layout],0.4

def finite_pair(first,second):
    try:
        pair=float(first),float(second)
    except (TypeError,ValueError):
        return None
    return pair if all(math.isfinite(value) for value in pair) else None

def percentile(values,fraction):
    values=sorted(values)
    position=(len(values)-1)*fraction
    lower=int(position)
    upper=min(lower+1,len(values)-1)
    weight=position-lower
    return values[lower]*(1-weight)+values[upper]*weight

def focus_rd_axes(axis,method_points):
    all_points=[point for points in method_points.values() for point in points]
    primary_points=method_points["MRCZ"]+method_points["FDWZ"]
    br_left=min(point[0] for point in all_points)*0.9
    br_right=max(percentile([point[0] for point in all_points],0.8),max(point[0] for point in primary_points))*1.05
    psnr_values=[point[1] for point in all_points]
    psnr_top=max(percentile(psnr_values,0.95),max(point[1] for point in primary_points))
    psnr_bottom=min(psnr_values)
    psnr_padding=max(1,(psnr_top-psnr_bottom)*0.03)
    axis.set_xscale("log")
    axis.set_xlim(br_left,br_right)
    axis.set_ylim(psnr_bottom-psnr_padding,psnr_top+psnr_padding)

def finish_figure(figure,xlabel,ylabel,output_path,paper_output):
    family_handles=[]
    for family in ["MRCZ","FDWZ","SZ3","SZ2","QoZ","ZFP","SPERR","MGARD"]:
        family_handles.append(plt.Line2D([0],[0],color=family_colors[family],linestyle="-",linewidth=0.8,label="MGARD-X" if family=="MGARD" else family))
    layout_handles=[plt.Line2D([0],[0],color="#4b5563",linestyle=layout_styles[layout],linewidth=0.8,label=layout) for layout in layouts]
    figure.legend(handles=family_handles,loc="upper center",ncol=8,frameon=False,bbox_to_anchor=(0.5,0.995),handlelength=1.0,handletextpad=0.25,columnspacing=0.6,borderaxespad=0)
    figure.legend(handles=layout_handles,loc="upper center",ncol=6,frameon=False,bbox_to_anchor=(0.5,0.945),handlelength=1.2,handletextpad=0.25,columnspacing=0.65,borderaxespad=0)
    figure.supxlabel(xlabel,x=0.53,y=0.025)
    figure.supylabel(ylabel,x=0.008,y=0.48)
    figure.tight_layout(rect=(0.015,0.045,1,0.88),pad=0.25,w_pad=0.2,h_pad=0.55)
    save_svg(figure,output_path,metadata={"Date":None},bbox_inches="tight",pad_inches=0.02)
    os.makedirs(os.path.dirname(paper_output),exist_ok=True)
    figure.savefig(paper_output,format="png",dpi=300)
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:
        figure.savefig(os.path.join(preview_dir,"small_rd.png"),format="png",dpi=160,bbox_inches="tight")
    plt.close(figure)

def plot_rd(result_dir):
    figure,axes=prepare_axes()
    for panel_index,(axis,(title,relative_path)) in enumerate(zip(axes,fields)):
        rows=read_rows(os.path.join(result_dir,relative_path))
        method_points={}
        for method in methods:
            method_rows,stored_method=read_method_rows(result_dir,relative_path,method,rows)
            points=sorted(point for row in method_rows if (point:=finite_pair(row.get(f"{stored_method}_BR"),row.get(f"{stored_method}_PSNR"))) is not None)
            method_points[method]=points
            color,linestyle,linewidth=method_style(method)
            axis.plot([point[0] for point in points],[point[1] for point in points],color=color,linestyle=linestyle,linewidth=linewidth,zorder=1)
        focus_rd_axes(axis,method_points)
        axis.set_xlabel(f"({chr(ord('a')+panel_index)}) {title}",labelpad=1.5)
        axis.grid(True,which="major",linestyle="--",linewidth=0.35,alpha=0.35)
        axis.xaxis.set_major_locator(LogLocator(base=10,numticks=4))
        axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
        style_plot_axis(axis)
        axis.tick_params(axis="x",which="minor",length=0)
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    finish_figure(figure,"Bit rate (log scale)","PSNR (dB)",os.path.join(result_dir,"small_rd.svg"),os.path.join(project_root,"论文","figs","7","small-rd.png"))

def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","RD")
    plot_rd(result_dir)

if __name__=="__main__":
    main()
