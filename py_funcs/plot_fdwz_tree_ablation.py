import csv,os
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FixedFormatter

project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
max_tree_path=os.path.join(project_root,"results/RD/pb0_fdwz_max_tree/pb0.f00001_P.csv")
single_tree_path=os.path.join(project_root,"results/RD/pb0_fdwz_single_tree/pb0.f00001_P.csv")
output_path=os.path.join(project_root,"论文/figs/5/fdwz-tree-ablation.pdf")

def read_curve(path):
    with open(path,encoding="utf-8") as csv_file:
        rows=list(csv.DictReader(csv_file))
    return sorted((float(row["FDWZ-ablation_BR"]),float(row["FDWZ-ablation_PSNR"])) for row in rows)

plt.rcParams.update({"font.family":"serif","font.serif":["Nimbus Roman","Times New Roman","Times"],"font.size":9,"axes.labelsize":9,"axes.titlesize":9,"xtick.labelsize":9,"ytick.labelsize":9,"legend.fontsize":9,"pdf.fonttype":42,"ps.fonttype":42})
figure=plt.figure(figsize=(3.5,1.55))
left_axes=figure.add_axes([0.0,0.0,0.64,1.0])
left_axes.axis("off")
axes=figure.add_axes([0.735,0.34,0.255,0.48])

for path,label,color,marker in [(max_tree_path,"Multi-tree","#0072B2","o"),(single_tree_path,"Single-tree","#D55E00","s")]:
    points=read_curve(path)
    axes.plot([point[0] for point in points],[point[1] for point in points],label=label,color=color,marker=marker,markersize=3.5,linewidth=1.2,markeredgewidth=0.6)

axes.set_xlim(0,4.5)
axes.set_ylim(40,130)
axes.xaxis.set_major_locator(FixedLocator([0,2,4]))
axes.xaxis.set_major_formatter(FixedFormatter(["0","2","4"]))
axes.yaxis.set_major_locator(FixedLocator([40,80,120]))
axes.set_xlabel("Bit rate",labelpad=1)
axes.set_ylabel("PSNR (dB)",labelpad=1)
axes.tick_params(axis="both",which="major",direction="in",length=2.5,width=0.6,pad=1)
axes.tick_params(axis="both",which="minor",direction="in",length=1.5,width=0.5)
axes.legend(loc="lower center",bbox_to_anchor=(0.5,0.97),frameon=False,handlelength=1.1,handletextpad=0.3,borderaxespad=0,labelspacing=0.15)
for spine in axes.spines.values():
    spine.set_linewidth(0.6)
figure.text(0.8625,0.015,"(b) Rate-distortion",ha="center",va="bottom",fontsize=9)

figure.savefig(output_path,format="pdf",dpi=300)
