import csv,math

def prepare_br_psnr_curve(points):
    by_psnr={}
    for psnr,br in points:
        psnr=float(psnr)
        br=float(br)
        if math.isfinite(psnr) and math.isfinite(br) and br>0:
            by_psnr[psnr]=min(br,by_psnr.get(psnr,math.inf))
    return sorted(by_psnr.items())

def interpolate_br_from_curve(curve,target_psnr):
    target_psnr=float(target_psnr)
    if not math.isfinite(target_psnr):
        raise ValueError("target_psnr必须是有限数")
    if not curve:
        raise ValueError("BR--PSNR曲线没有有效数据点")
    if target_psnr<curve[0][0] or target_psnr>curve[-1][0]:
        raise ValueError(f"目标PSNR {target_psnr} dB超出曲线范围[{curve[0][0]},{curve[-1][0]}] dB")
    if target_psnr==curve[0][0]:
        return curve[0][1]
    for index in range(1,len(curve)):
        left_psnr,left_br=curve[index-1]
        right_psnr,right_br=curve[index]
        if target_psnr==right_psnr:
            return right_br
        if target_psnr<right_psnr:
            weight=(target_psnr-left_psnr)/(right_psnr-left_psnr)
            return math.exp(math.log(left_br)+(math.log(right_br)-math.log(left_br))*weight)
    raise RuntimeError("未能在已验证的PSNR范围内完成插值")

def interpolate_br_at_psnr(csv_path,method,target_psnr):
    br_column=f"{method}_BR"
    psnr_column=f"{method}_PSNR"
    points=[]
    with open(csv_path,newline="") as f:
        reader=csv.DictReader(f)
        fieldnames=reader.fieldnames or []
        missing=[column for column in [br_column,psnr_column] if column not in fieldnames]
        if missing:
            raise KeyError(f"CSV缺少列: {missing}")
        for row in reader:
            br=row.get(br_column,"")
            psnr=row.get(psnr_column,"")
            if br!="" and psnr!="":
                points.append((psnr,br))
    return interpolate_br_from_curve(prepare_br_psnr_curve(points),target_psnr)
