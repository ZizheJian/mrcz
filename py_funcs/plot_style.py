import matplotlib.pyplot as plt

PAPER_FONT_SIZE=8.0
METHOD_COLORS={
    "MRCZ":"#111827",
    "FDWZ":"#e15759",
    "SZ3-4D":"#4e79a7",
    "SZ2-3D-fused":"#76b7b2",
    "QoZ-4D":"#59a14f",
    "ZFP-3D-fused":"#f28e2b",
    "SPERR-2D":"#b07aa1",
    "MGARD-1D":"#9c755f"}
METHOD_LABELS={
    "MRCZ":"MRCZ",
    "FDWZ":"FDWZ",
    "SZ3-4D":"SZ3",
    "SZ2-3D-fused":"SZ2",
    "QoZ-4D":"QoZ",
    "ZFP-3D-fused":"ZFP",
    "SPERR-2D":"SPERR",
    "MGARD-1D":"MGARD"}
METHOD_MARKERS={"MRCZ":"*","FDWZ":"o","SZ3-4D":"^","QoZ-4D":"D","ZFP-3D-fused":"v","SPERR-2D":"P"}
TIMING_METHODS=("MRCZ","FDWZ","SZ3-4D","QoZ-4D","ZFP-3D-fused","SPERR-2D")
RD_METHODS=("MRCZ","FDWZ","SZ3-4D","SZ2-3D-fused","QoZ-4D","ZFP-3D-fused","SPERR-2D","MGARD-1D")
TRANSFER_METHODS=("MRCZ","FDWZ","SZ3-4D","QoZ-4D","ZFP-3D-fused","SPERR-2D","MGARD-1D")

def apply_paper_style(svg_hashsalt=None,font_size=PAPER_FONT_SIZE,tick_size=None,legend_size=None,serif_fonts=None):
    if tick_size is None:tick_size=font_size
    if legend_size is None:legend_size=font_size
    if serif_fonts is None:serif_fonts=["Nimbus Roman"]
    if svg_hashsalt is not None:plt.rcParams["svg.hashsalt"]=svg_hashsalt
    plt.rcParams.update({"font.family":"serif","font.serif":serif_fonts,"mathtext.fontset":"stix","font.size":font_size,"axes.labelsize":font_size,"axes.titlesize":font_size,"xtick.labelsize":tick_size,"ytick.labelsize":tick_size,"legend.fontsize":legend_size,"pdf.fonttype":42,"ps.fonttype":42})

def style_plot_axis(axis,tick_length=1.8,tick_width=0.45,tick_pad=1,spine_width=0.5,minor_tick_length=None,minor_tick_width=None):
    axis.tick_params(direction="in",length=tick_length,width=tick_width,pad=tick_pad)
    if minor_tick_length is not None:
        axis.tick_params(axis="both",which="minor",direction="in",length=minor_tick_length,width=minor_tick_width if minor_tick_width is not None else tick_width)
    for spine in axis.spines.values():spine.set_linewidth(spine_width)
