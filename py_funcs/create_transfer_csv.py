import argparse,csv,os

methods=["MRCZ","FDWZ","SZ3-4D","SZ2-3D-fused","QoZ-4D","ZFP-3D-fused","SPERR-2D","MGARD-1D"]

parser=argparse.ArgumentParser()
parser.add_argument("data_path")
parser.add_argument("--rel-eb",required=True,type=float)
parser.add_argument("--fields",required=True)
args=parser.parse_args()
fields=[field for field in args.fields.split(",") if field]
if not fields:
    raise ValueError("--fields不能为空")

output_path=f"{args.data_path}_transfer_rel_eb_{args.rel_eb:g}.csv"
with open(output_path,"w",newline="") as output_file:
    writer=csv.DictWriter(output_file,fieldnames=["Method","REL EB","Trial","Files","Bytes Transferred","Effective Bytes Per Second"],lineterminator="\n")
    writer.writeheader()
    for method in methods:
        paths=[f"{args.data_path}_{field}.{method}" for field in fields]
        missing=[path for path in paths if not os.path.isfile(path)]
        if missing:
            raise FileNotFoundError(f"{method}缺少压缩文件: {missing}")
        transferred_bytes=sum(os.path.getsize(path) for path in paths)
        for trial in range(1,4):
            writer.writerow({"Method":method,"REL EB":f"{args.rel_eb:g}","Trial":trial,"Files":len(paths),"Bytes Transferred":transferred_bytes,"Effective Bytes Per Second":""})

print(output_path)
