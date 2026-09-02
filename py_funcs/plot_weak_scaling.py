import csv,html,os

width=1200
height=820
methods=["MRCZ","FDWZ","SZ3-4D","QoZ-4D","ZFP-3D-fused","SPERR-2D"]
colors={"MRCZ":"#111827","FDWZ":"#e15759","SZ3-4D":"#4e79a7","QoZ-4D":"#59a14f","ZFP-3D-fused":"#f28e2b","SPERR-2D":"#b07aa1"}

def load_rows(path):
    with open(path,newline="") as f:
        return list(csv.DictReader(f))

def svg_text(x,y,value,size=15,anchor="middle",weight="normal",rotate=None):
    transform=f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ""
    return f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}" font-weight="{weight}"{transform}>{html.escape(str(value))}</text>'

def draw_panel(parts,rows,title,metric,x,y,w,h):
    left=x+64
    right=x+w-18
    top=y+36
    bottom=y+h-48
    cpu_values=sorted({int(row["CPUs"]) for row in rows})
    x_index={cpus:index for index,cpus in enumerate(cpu_values)}
    x_scale=lambda cpus:left+(right-left)*x_index[cpus]/(len(cpu_values)-1)
    error_metric=metric+" Std Dev"
    y_max=max(1.05,max(float(row[metric])+float(row.get(error_metric,0) or 0) for row in rows)*1.05)
    y_scale=lambda value:bottom-(bottom-top)*value/y_max
    parts.append(svg_text((left+right)/2,y+20,title,18,weight="bold"))
    for index in range(5):
        value=y_max*index/4
        py=y_scale(value)
        parts.append(f'<line x1="{left}" y1="{py}" x2="{right}" y2="{py}" stroke="#d1d5db" stroke-width="1"/>')
        parts.append(svg_text(left-10,py+5,f"{value:.2f}",13,anchor="end"))
    for cpus in cpu_values:
        px=x_scale(cpus)
        parts.append(f'<line x1="{px}" y1="{bottom}" x2="{px}" y2="{bottom+5}" stroke="#374151"/>')
        parts.append(svg_text(px,bottom+23,cpus,13))
    parts.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#374151" stroke-width="1.5"/>')
    parts.append(f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#374151" stroke-width="1.5"/>')
    parts.append(svg_text((left+right)/2,bottom+43,"CPU count",14))
    parts.append(svg_text(x+15,(top+bottom)/2,"Efficiency",14,rotate=-90))
    for method in methods:
        method_rows=sorted((row for row in rows if row["Method"]==method),key=lambda row:int(row["CPUs"]))
        if len(method_rows)<2:
            continue
        points=" ".join(f'{x_scale(int(row["CPUs"])):.2f},{y_scale(float(row[metric])):.2f}' for row in method_rows)
        stroke_width=4 if method=="MRCZ" else 2.2
        parts.append(f'<polyline points="{points}" fill="none" stroke="{colors[method]}" stroke-width="{stroke_width}" stroke-linejoin="round" stroke-linecap="round"/>')
        for row in method_rows:
            px=x_scale(int(row["CPUs"]))
            value=float(row[metric])
            error=float(row.get(error_metric,0) or 0)
            py=y_scale(value)
            if error>0:
                error_top=y_scale(min(value+error,y_max))
                error_bottom=y_scale(max(value-error,0))
                parts.append(f'<line x1="{px}" y1="{error_top}" x2="{px}" y2="{error_bottom}" stroke="{colors[method]}" stroke-width="1.4"/>')
                parts.append(f'<line x1="{px-4}" y1="{error_top}" x2="{px+4}" y2="{error_top}" stroke="{colors[method]}" stroke-width="1.4"/>')
                parts.append(f'<line x1="{px-4}" y1="{error_bottom}" x2="{px+4}" y2="{error_bottom}" stroke="{colors[method]}" stroke-width="1.4"/>')
            radius=4.5 if method=="MRCZ" else 3.2
            parts.append(f'<circle cx="{px}" cy="{py}" r="{radius}" fill="{colors[method]}"/>')

def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","weak_scaling")
    datasets={"pb0":load_rows(os.path.join(result_dir,"weak_pb0_summary.csv")),"jicf0":load_rows(os.path.join(result_dir,"weak_jicf_summary.csv"))}
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">','<rect width="100%" height="100%" fill="white"/>','<g font-family="Arial,Helvetica,sans-serif" fill="#111827">',svg_text(width/2,38,"Weak scaling efficiency",25,weight="bold")]
    legend_width=1050
    legend_left=(width-legend_width)/2
    for index,method in enumerate(methods):
        x=legend_left+index*legend_width/(len(methods)-1)
        parts.append(f'<line x1="{x-27}" y1="76" x2="{x-5}" y2="76" stroke="{colors[method]}" stroke-width="{4 if method=="MRCZ" else 2.2}"/>')
        parts.append(f'<circle cx="{x-16}" cy="76" r="{4.5 if method=="MRCZ" else 3.2}" fill="{colors[method]}"/>')
        parts.append(svg_text(x,81,method,14,anchor="start",weight="bold" if method=="MRCZ" else "normal"))
    panel_width=545
    panel_height=315
    draw_panel(parts,datasets["pb0"],"pb0 compression","Compression Efficiency",35,105,panel_width,panel_height)
    draw_panel(parts,datasets["jicf0"],"jicf0 compression","Compression Efficiency",620,105,panel_width,panel_height)
    draw_panel(parts,datasets["pb0"],"pb0 decompression","Decompression Efficiency",35,470,panel_width,panel_height)
    draw_panel(parts,datasets["jicf0"],"jicf0 decompression","Decompression Efficiency",620,470,panel_width,panel_height)
    parts.extend(["</g>","</svg>"])
    output_name=os.environ.get("MRCZ_WEAK_SCALING_OUTPUT","weak_scaling_efficiency.svg")
    output_path=os.path.join(result_dir,output_name)
    with open(output_path,"w") as f:
        f.write("\n".join(parts)+"\n")
    print(output_path)

if __name__=="__main__":
    main()
