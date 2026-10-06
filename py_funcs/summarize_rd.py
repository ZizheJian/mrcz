import csv,math,os,statistics

small_fields=[
    ("mixlay0/UX","mixlay0/mixlay0.f00001_UX.csv"),("mixlay0/UY","mixlay0/mixlay0.f00001_UY.csv"),("mixlay0/P","mixlay0/mixlay0.f00001_P.csv"),("mixlay0/T","mixlay0/mixlay0.f00001_T.csv"),("mixlay0/S00","mixlay0/mixlay0.f00001_S00.csv"),("mixlay0/S01","mixlay0/mixlay0.f00001_S01.csv"),
    ("test0/UZ","test0/test0.f00001_UZ.csv"),("test0/T","test0/test0.f00001_T.csv"),("test0/S00","test0/test0.f00001_S00.csv"),("test0/S01","test0/test0.f00001_S01.csv"),("tc_channel0/UX","tc_channel0/tc_channel0.f00001_UX.csv"),("tc_channel0/UY","tc_channel0/tc_channel0.f00001_UY.csv"),
    ("tc_channel0/UZ","tc_channel0/tc_channel0.f00001_UZ.csv"),("MFEM ex0/no_name_field","mfem_ex0/sol.gf_no_name_field.csv"),("MFEM ex1/no_name_field","mfem_ex1/sol.gf_no_name_field.csv"),("Nektar++/p","nektar_cylinder/Cylinder-Base.vtu_p.csv"),("Nektar++/u","nektar_cylinder/Cylinder-Base.vtu_u.csv"),("Nektar++/v","nektar_cylinder/Cylinder-Base.vtu_v.csv")]
large_fields=[
    ("pb0/UX","pb0/pb0.f00001_UX.csv"),("pb0/UY","pb0/pb0.f00001_UY.csv"),("pb0/UZ","pb0/pb0.f00001_UZ.csv"),("pb0/P","pb0/pb0.f00001_P.csv"),("pb0/T","pb0/pb0.f00001_T.csv"),
    ("jicf0/UX","jicf0/jicf0.f00000_UX.csv"),("jicf0/UY","jicf0/jicf0.f00000_UY.csv"),("jicf0/UZ","jicf0/jicf0.f00000_UZ.csv"),("jicf0/P","jicf0/jicf0.f00000_P.csv"),("jicf0/S00","jicf0/jicf0.f00000_S00.csv")]
required_baselines=["SZ3-4D","SZ2-3D-fused","QoZ-4D","ZFP-3D-fused","SPERR-2D"]
optional_baselines=["MGARD-1D"]
baselines=required_baselines+optional_baselines
methods=["MRCZ","FDWZ"]+baselines

def read_rows(path):
    with open(path,newline="") as f:
        return list(csv.DictReader(f))

def method_rows(result_dir,scope,relative_path,method,default_rows):
    if method!="MGARD-1D":return default_rows
    if scope=="large":
        dataset=relative_path.split(os.sep,1)[0]
        path=os.path.join(os.path.dirname(result_dir),"small_core_throughput",dataset,os.path.basename(relative_path))
    else:
        path=os.path.join(result_dir,"mgard_layout_selection",os.path.dirname(relative_path),"uniform",os.path.basename(relative_path))
    return read_rows(path)

def prepare_curve(rows,method):
    by_psnr={}
    for row in rows:
        try:
            psnr=float(row[f"{method}_PSNR"])
            br=float(row[f"{method}_BR"])
        except (KeyError,TypeError,ValueError):
            continue
        if math.isfinite(psnr) and math.isfinite(br) and br>0:
            by_psnr[psnr]=min(br,by_psnr.get(psnr,math.inf))
    return sorted(by_psnr.items())

def interpolate(curve,target):
    if target==curve[0][0]:
        return curve[0][1]
    for index in range(1,len(curve)):
        left_psnr,left_br=curve[index-1]
        right_psnr,right_br=curve[index]
        if target<=right_psnr:
            weight=(target-left_psnr)/(right_psnr-left_psnr)
            return math.exp(math.log(left_br)+(math.log(right_br)-math.log(left_br))*weight)
    return curve[-1][1]

def select_target(psnr_low,psnr_high):
    return (psnr_low+psnr_high)/2

def summarize_field(result_dir,scope,title,relative_path):
    rows=read_rows(os.path.join(result_dir,relative_path))
    curves={method:prepare_curve(method_rows(result_dir,scope,relative_path,method,rows),method) for method in methods}
    required_methods=["MRCZ","FDWZ"]+required_baselines
    if any(len(curves[method])<2 for method in required_methods):
        raise RuntimeError(f"{title}缺少可插值的正式方法曲线")
    common_methods=required_methods+[method for method in optional_baselines if len(curves[method])>=2]
    psnr_low=max(curves[method][0][0] for method in common_methods)
    psnr_high=min(curves[method][-1][0] for method in common_methods)
    if psnr_high<=psnr_low:
        raise RuntimeError(f"{title}没有共同PSNR范围")
    target=select_target(psnr_low,psnr_high)
    comparable_methods=[method for method,curve in curves.items() if len(curve)>=2 and curve[0][0]<=target<=curve[-1][0]]
    brs={method:interpolate(curves[method],target) for method in comparable_methods}
    comparable_baselines=[method for method in baselines if method in brs]
    best_baseline=min(comparable_baselines,key=lambda method:brs[method])
    return {
        "Data":scope,"Field":title,"Common PSNR Min":psnr_low,"Common PSNR Max":psnr_high,"Target PSNR":target,
        "MRCZ vs FDWZ CR Growth":brs["FDWZ"]/brs["MRCZ"]-1,"FDWZ vs Best Official CR Growth":brs[best_baseline]/brs["FDWZ"]-1,"Best Official Baseline":best_baseline}

def summarize_scope(rows,scope):
    selected=rows if scope=="all" else [row for row in rows if row["Data"]==scope]
    mrcz_values=[row["MRCZ vs FDWZ CR Growth"] for row in selected]
    fdwz_values=[row["FDWZ vs Best Official CR Growth"] for row in selected]
    mrcz_max=max(selected,key=lambda row:row["MRCZ vs FDWZ CR Growth"])
    fdwz_max=max(selected,key=lambda row:row["FDWZ vs Best Official CR Growth"])
    return {
        "field_num":len(selected),"mrcz_mean":statistics.geometric_mean(1+value for value in mrcz_values)-1,"mrcz_max":mrcz_max["MRCZ vs FDWZ CR Growth"],"mrcz_max_field":mrcz_max["Field"],
        "fdwz_mean":statistics.geometric_mean(1+value for value in fdwz_values)-1,"fdwz_max":fdwz_max["FDWZ vs Best Official CR Growth"],"fdwz_max_field":fdwz_max["Field"]}

def main():
    project_root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    result_dir=os.path.join(project_root,"results","RD")
    rows=[]
    for scope,fields in [("large",large_fields),("small",small_fields)]:
        for title,relative_path in fields:
            if title=="test0/UZ":
                continue
            rows.append(summarize_field(result_dir,scope,title,relative_path))
    output_path=os.path.join(result_dir,"midpoint_cr_summary.csv")
    with open(output_path,"w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader()
        writer.writerows({name:f"{value:.10g}" if isinstance(value,float) else value for name,value in row.items()} for row in rows)
    for scope in ["large","small","all"]:
        summary=summarize_scope(rows,scope)
        print(f"{scope}: {summary['field_num']} fields")
        print(f"  MRCZ relative to FDWZ CR growth: geometric mean={summary['mrcz_mean']*100:.2f}%, max={summary['mrcz_max']*100:.2f}% ({summary['mrcz_max_field']})")
        print(f"  FDWZ relative to best official CR growth: geometric mean={summary['fdwz_mean']*100:.2f}%, max={summary['fdwz_max']*100:.2f}% ({summary['fdwz_max_field']})")
    print(output_path)

if __name__=="__main__":
    main()
