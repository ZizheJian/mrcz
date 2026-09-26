import csv,math,os
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

plt.rcParams["svg.hashsalt"]="mrcz-small-rd"

layouts=["1D","2D","3D-fused","3D-per-element","4D"]
methods=["MRCZ","FDWZ"]+[f"{compressor}-{layout}" for compressor in ["SZ3","SZ2","QoZ","ZFP"] for layout in layouts]+["SPERR-2D","SPERR-3D-fused","SPERR-3D-per-element"]
family_colors={"MRCZ":"#111827","FDWZ":"#e15759","SZ3":"#4e79a7","SZ2":"#76b7b2","QoZ":"#59a14f","ZFP":"#f28e2b","SPERR":"#b07aa1"}
layout_styles={"1D":(0,(5,1,1,1)),"2D":"--","3D-fused":"-.","3D-per-element":":","4D":"-"}
fields=[
    ("1UX","mixlay0/mixlay0.f00001_UX.csv"),("1UY","mixlay0/mixlay0.f00001_UY.csv"),("1P","mixlay0/mixlay0.f00001_P.csv"),("1T","mixlay0/mixlay0.f00001_T.csv"),("1S00","mixlay0/mixlay0.f00001_S00.csv"),("1S01","mixlay0/mixlay0.f00001_S01.csv"),
    ("2UZ","test0/test0.f00001_UZ.csv"),("2T","test0/test0.f00001_T.csv"),("2S00","test0/test0.f00001_S00.csv"),("2S01","test0/test0.f00001_S01.csv"),("3UX","tc_channel0/tc_channel0.f00001_UX.csv"),("3UY","tc_channel0/tc_channel0.f00001_UY.csv"),
    ("3UZ","tc_channel0/tc_channel0.f00001_UZ.csv"),("Dataset 4","mfem_ex0/sol.gf_no_name_field.csv"),("Dataset 5","mfem_ex1/sol.gf_no_name_field.csv"),("6P","nektar_cylinder/Cylinder-Base.vtu_p.csv"),("6U","nektar_cylinder/Cylinder-Base.vtu_u.csv"),("6V","nektar_cylinder/Cylinder-Base.vtu_v.csv")]

def read_rows(path):
    with open(path,newline="") as f:
        return list(csv.DictReader(f))

def prepare_axes():
    figure,axes=plt.subplots(3,6,figsize=(7.16,3.85))
    return figure,axes.flat

def method_style(method):
    if method in ["MRCZ","FDWZ"]:
        return family_colors[method],"-",0.8
    family,layout=method.split("-",1)
    return family_colors[family],layout_styles[layout],0.8

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
    br_right=max(percentile([point[0] for point in all_points],0.9),max(point[0] for point in primary_points))*1.05
    psnr_values=[point[1] for point in all_points]
    psnr_top=max(percentile(psnr_values,0.95),max(point[1] for point in primary_points))
    psnr_bottom=min(psnr_values)
    psnr_padding=max(1,(psnr_top-psnr_bottom)*0.03)
    axis.set_xlim(0,br_right)
    axis.set_ylim(psnr_bottom-psnr_padding,psnr_top+psnr_padding)

def finish_figure(figure,xlabel,ylabel,output_path,paper_output):
    family_handles=[]
    for family in ["MRCZ","FDWZ","SZ3","SZ2","QoZ","ZFP","SPERR"]:
        family_handles.append(plt.Line2D([0],[0],color=family_colors[family],linestyle="-",linewidth=0.8,label=family))
    layout_handles=[plt.Line2D([0],[0],color="#4b5563",linestyle=layout_styles[layout],linewidth=0.8,label=layout) for layout in layouts]
    figure.legend(handles=family_handles,loc="upper center",ncol=7,frameon=False,bbox_to_anchor=(0.5,0.995),handlelength=1.0,handletextpad=0.25,columnspacing=0.75,borderaxespad=0)
    figure.legend(handles=layout_handles,loc="upper center",ncol=5,frameon=False,bbox_to_anchor=(0.5,0.945),handlelength=1.2,handletextpad=0.25,columnspacing=0.8,borderaxespad=0)
    figure.supxlabel(xlabel,x=0.53,y=0.025)
    figure.supylabel(ylabel,x=0.008,y=0.48)
    figure.tight_layout(rect=(0.015,0.045,1,0.88),pad=0.25,w_pad=0.2,h_pad=0.55)
    figure.savefig(output_path,format="svg",metadata={"Date":None},bbox_inches="tight",pad_inches=0.02)
    os.makedirs(os.path.dirname(paper_output),exist_ok=True)
    for axis in figure.axes:axis.set_rasterization_zorder(1.5)
    figure.savefig(paper_output,format="pdf",dpi=300,bbox_inches="tight",pad_inches=0.02)
    preview_dir=os.environ.get("MRCZ_PLOT_PREVIEW_DIR")
    if preview_dir:
        figure.savefig(os.path.join(preview_dir,"small_rd.png"),format="png",dpi=160,bbox_inches="tight")
    plt.close(figure)

def plot_rd(result_dir):
    plt.rcParams.update({"font.family":"serif","font.serif":["Nimbus Roman"],"mathtext.fontset":"stix","font.size":8.0,"axes.labelsize":8.0,"axes.titlesize":8.0,"xtick.labelsize":8.0,"ytick.labelsize":8.0,"legend.fontsize":8.0,"pdf.fonttype":42,"ps.fonttype":42})
    figure,axes=prepare_axes()
    for panel_index,(axis,(title,relative_path)) in enumerate(zip(axes,fields)):
        rows=read_rows(os.path.join(result_dir,relative_path))
        method_points={}
        for method in methods:
            points=sorted(point for row in rows if (point:=finite_pair(row.get(f"{method}_BR"),row.get(f"{method}_PSNR"))) is not None)
            method_points[method]=points
            color,linestyle,linewidth=method_style(method)
            axis.plot([point[0] for point in points],[point[1] for point in points],color=color,linestyle=linestyle,linewidth=linewidth,zorder=1)
        focus_rd_axes(axis,method_points)
        axis.set_xlabel(f"({chr(ord('a')+panel_index)}) {title}",labelpad=1.5)
        axis.grid(True,which="both",alpha=0.25,linewidth=0.35)
        axis.xaxis.set_major_locator(MaxNLocator(nbins=3))
        axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
        axis.tick_params(direction="in",length=1.8,width=0.45,pad=1)
        for spine in axis.spines.values():spine.set_linewidth(0.5)
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    finish_figure(figure,"Bit rate","PSNR (dB)",os.path.join(result_dir,"small_rd.svg"),os.path.join(project_root,"论文","figs","7","small-rd.pdf"))

def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","RD")
    plot_rd(result_dir)

if __name__=="__main__":
    main()
