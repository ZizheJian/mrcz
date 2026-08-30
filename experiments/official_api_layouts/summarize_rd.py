# XXX 在各方法共同覆盖的PSNR范围内插值BR，以同PSNR的BR regret选择每个官方压缩器最好的布局。

import csv,math,os,statistics,sys

experiment_dir=os.path.dirname(os.path.abspath(__file__))
project_root=os.path.abspath(os.path.join(experiment_dir,"../.."))
sys.path.insert(0,project_root)
from run_small_rd import datasets,result_method_name
summary_path=os.path.join(experiment_dir,"layout_summary.csv")
case_path=os.path.join(experiment_dir,"case_summary.csv")
family_methods={
    "SZ3":["SZ3-1D-uniform","SZ3-2D-uniform","SZ3-3D-uniform-fused","SZ3-3D-uniform-per-element","SZ3-4D-uniform"],
    "SZ2":["SZ2-1D-uniform","SZ2-2D-uniform","SZ2-3D-uniform-fused","SZ2-3D-uniform-per-element","SZ2-4D-uniform"],
    "QoZ":["QoZ-1D-uniform","QoZ-2D-uniform","QoZ-3D-uniform-fused","QoZ-3D-uniform-per-element","QoZ-4D-uniform"],
    "ZFP":["ZFP-1D-uniform","ZFP-2D-uniform","ZFP-3D-uniform-fused","ZFP-3D-uniform-per-element","ZFP-4D-uniform"],
    "SPERR":["SPERR-2D-uniform","SPERR-3D-uniform-fused","SPERR-3D-uniform-per-element"]}
# XXX test0/UZ的RD曲线已知接近PSNR轴且不能区分大多数方法，保留原始结果但不让它支配布局选择。
excluded_cases={("test0","UZ")}

def percentile(values,p):
    values=sorted(values)
    if not values:
        return math.nan
    position=(len(values)-1)*p
    left=int(position)
    right=min(left+1,len(values)-1)
    weight=position-left
    return values[left]*(1-weight)+values[right]*weight

def prepare_curve(points):
    by_psnr={}
    for psnr,br in points:
        if math.isfinite(psnr) and math.isfinite(br) and br>0:
            by_psnr[psnr]=min(br,by_psnr.get(psnr,math.inf))
    return sorted(by_psnr.items())

def interpolate(curve,target):
    for index in range(1,len(curve)):
        left_psnr,left_br=curve[index-1]
        right_psnr,right_br=curve[index]
        if target<=right_psnr:
            if right_psnr==left_psnr:
                return min(left_br,right_br)
            weight=(target-left_psnr)/(right_psnr-left_psnr)
            return left_br+(right_br-left_br)*weight
    return curve[-1][1]

def main():
    curves={}
    methods={method for family in family_methods.values() for method in family}
    for dataset in datasets:
        for field in dataset["fields"]:
            path=f"{dataset['data']}_{field}.csv"
            if not os.path.isfile(path):
                continue
            with open(path,newline="") as f:
                for row in csv.DictReader(f):
                    for method in methods:
                        csv_method=result_method_name(method)
                        br=row.get(f"{csv_method}_BR","")
                        psnr=row.get(f"{csv_method}_PSNR","")
                        if br!="" and psnr!="":
                            curves.setdefault((dataset["name"],field,method),[]).append((float(psnr),float(br)))
    regrets={}
    wins={}
    cases=[]
    for family,methods in family_methods.items():
        dataset_fields=sorted({(dataset,field) for dataset,field,method in curves if method in methods and (dataset,field) not in excluded_cases})
        for dataset,field in dataset_fields:
            method_curves={method:prepare_curve(curves.get((dataset,field,method),[])) for method in methods}
            if any(len(curve)<2 for curve in method_curves.values()):
                continue
            low=max(curve[0][0] for curve in method_curves.values())
            high=min(curve[-1][0] for curve in method_curves.values())
            if high<=low:
                continue
            targets=[low+(high-low)*index/20 for index in range(21)]
            local={method:[] for method in methods}
            for target in targets:
                brs={method:interpolate(curve,target) for method,curve in method_curves.items()}
                best=min(brs.values())
                for method,br in brs.items():
                    regret=br/best-1
                    regrets.setdefault((family,method),[]).append(regret)
                    local[method].append(regret)
                winner=min(brs,key=brs.get)
                wins[(family,winner)]=wins.get((family,winner),0)+1
            ordered=sorted(methods,key=lambda method:statistics.mean(local[method]))
            cases.append({"Family":family,"Dataset":dataset,"Field":field,"PSNR Min":low,"PSNR Max":high,"Best Layout":result_method_name(ordered[0]),"Second Layout":result_method_name(ordered[1]),"Best Mean BR Regret":statistics.mean(local[ordered[0]]),"Second Mean BR Regret":statistics.mean(local[ordered[1]])})
    with open(summary_path,"w",newline="") as f:
        names=["Family","Method","Comparison Points","Win Rate","Median BR Regret","Mean BR Regret","P90 BR Regret"]
        writer=csv.DictWriter(f,fieldnames=names)
        writer.writeheader()
        for family,methods in family_methods.items():
            family_total=sum(len(regrets.get((family,method),[])) for method in methods)//len(methods)
            ordered=sorted(methods,key=lambda method:statistics.mean(regrets.get((family,method),[math.inf])))
            print(family)
            for method in ordered:
                values=regrets.get((family,method),[])
                if not values:
                    continue
                row={"Family":family,"Method":result_method_name(method),"Comparison Points":len(values),"Win Rate":wins.get((family,method),0)/family_total,"Median BR Regret":statistics.median(values),"Mean BR Regret":statistics.mean(values),"P90 BR Regret":percentile(values,0.9)}
                writer.writerow(row)
                print(f"  {result_method_name(method)}: mean={row['Mean BR Regret']*100:.2f}%,median={row['Median BR Regret']*100:.2f}%,p90={row['P90 BR Regret']*100:.2f}%,win={row['Win Rate']*100:.1f}%")
    with open(case_path,"w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(cases[0]) if cases else ["Family","Dataset","Field"])
        writer.writeheader()
        writer.writerows(cases)

if __name__=="__main__":
    main()
