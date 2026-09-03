import csv,math,os
import matplotlib.pyplot as plt

layouts=["1D","2D","3D-fused","3D-per-element","4D"]
methods=["MRCZ","FDWZ"]+[f"{compressor}-{layout}" for compressor in ["SZ3","SZ2","QoZ","ZFP"] for layout in layouts]+["SPERR-2D","SPERR-3D-fused","SPERR-3D-per-element"]
family_colors={"MRCZ":"#111827","FDWZ":"#e15759","SZ3":"#4e79a7","SZ2":"#76b7b2","QoZ":"#59a14f","ZFP":"#f28e2b","SPERR":"#b07aa1"}
layout_styles={"1D":"-","2D":"--","3D-fused":"-.","3D-per-element":":","4D":(0,(5,1))}
layout_markers={"1D":"o","2D":"s","3D-fused":"^","3D-per-element":"x","4D":"D"}
fields=[
    ("mixlay0 / UX","mixlay0/mixlay0.f00001_UX.csv"),("mixlay0 / UY","mixlay0/mixlay0.f00001_UY.csv"),("mixlay0 / P","mixlay0/mixlay0.f00001_P.csv"),("mixlay0 / T","mixlay0/mixlay0.f00001_T.csv"),("mixlay0 / S00","mixlay0/mixlay0.f00001_S00.csv"),("mixlay0 / S01","mixlay0/mixlay0.f00001_S01.csv"),
    ("test0 / UZ","test0/test0.f00001_UZ.csv"),("test0 / T","test0/test0.f00001_T.csv"),("test0 / S00","test0/test0.f00001_S00.csv"),("test0 / S01","test0/test0.f00001_S01.csv"),("tc_channel0 / UX","tc_channel0/tc_channel0.f00001_UX.csv"),("tc_channel0 / UY","tc_channel0/tc_channel0.f00001_UY.csv"),
    ("tc_channel0 / UZ","tc_channel0/tc_channel0.f00001_UZ.csv"),("MFEM ex0 / no_name_field","mfem_ex0/sol.gf_no_name_field.csv"),("MFEM ex1 / no_name_field","mfem_ex1/sol.gf_no_name_field.csv"),("Nektar++ / p","nektar_cylinder/Cylinder-Base.vtu_p.csv"),("Nektar++ / u","nektar_cylinder/Cylinder-Base.vtu_u.csv"),("Nektar++ / v","nektar_cylinder/Cylinder-Base.vtu_v.csv")]

def read_rows(path):
    with open(path,newline="") as f:
        return list(csv.DictReader(f))

def prepare_axes():
    figure,axes=plt.subplots(3,6,figsize=(22,12))
    return figure,axes.flat

def method_style(method):
    if method in ["MRCZ","FDWZ"]:
        return family_colors[method],"-","*" if method=="MRCZ" else "o",2.8
    family,layout=method.split("-",1)
    return family_colors[family],layout_styles[layout],layout_markers[layout],1.3

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

def finish_figure(figure,xlabel,ylabel,output_path):
    handles=[]
    for method in methods:
        color,linestyle,marker,linewidth=method_style(method)
        handles.append(plt.Line2D([0],[0],color=color,linestyle=linestyle,marker=marker,linewidth=linewidth,markersize=4,label=method))
    figure.legend(handles=handles,loc="upper center",ncol=5,frameon=False,bbox_to_anchor=(0.5,0.995),fontsize=8)
    figure.supxlabel(xlabel,fontsize=13)
    figure.supylabel(ylabel,fontsize=13)
    figure.tight_layout(rect=(0.025,0.035,1,0.9))
    figure.savefig(output_path,format="svg",bbox_inches="tight")
    plt.close(figure)

def plot_rd(result_dir):
    figure,axes=prepare_axes()
    for axis,(title,relative_path) in zip(axes,fields):
        rows=read_rows(os.path.join(result_dir,relative_path))
        method_points={}
        for method in methods:
            points=sorted(point for row in rows if (point:=finite_pair(row.get(f"{method}_BR"),row.get(f"{method}_PSNR"))) is not None)
            method_points[method]=points
            color,linestyle,marker,linewidth=method_style(method)
            axis.plot([point[0] for point in points],[point[1] for point in points],color=color,linestyle=linestyle,marker=marker,markersize=2.5,linewidth=linewidth)
        focus_rd_axes(axis,method_points)
        axis.set_title(title,fontsize=9)
        axis.grid(True,which="both",alpha=0.25,linewidth=0.5)
        axis.tick_params(labelsize=7)
    finish_figure(figure,"Bit rate (bits/value)","PSNR (dB)",os.path.join(result_dir,"small_rd.svg"))

def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","RD")
    plot_rd(result_dir)

if __name__=="__main__":
    main()
