import math,os
import pandas as pd

def read_rate_distortion_from_conclusion(experiment_id:int):
    br_dict={}
    cr_dict={}
    mse_dict={}
    psnr_dict={}
    with open(f"history/conclusion{experiment_id:04d}.txt","r") as f_conclusion:
        lines=f_conclusion.readlines()[2:]
        for line in lines:
            if line.startswith("abs_eb="):
                parts=line.strip().split("\t")
                abs_eb=float(parts[0].split("=")[1].strip())
                br=float(parts[1].split("=")[1].strip())
                cr=float(parts[2].split("=")[1].strip())
                mse=float(parts[3].split("=")[1].strip())
                psnr=float(parts[4].split("=")[1].strip())
                br_dict[abs_eb]=br
                cr_dict[abs_eb]=cr
                mse_dict[abs_eb]=mse
                psnr_dict[abs_eb]=psnr
    return br_dict,cr_dict,mse_dict,psnr_dict

def read_rate_distortion_from_history(experiment_id:int,block_list:dict,global_param):
    br_dict={}
    cr_dict={}
    mse_dict={}
    psnr_dict={}
    if not os.path.exists(f"history/history{experiment_id:04d}.csv"):
        return br_dict,cr_dict,mse_dict,psnr_dict
    total_element_num={}
    with open(f"history/history{experiment_id:04d}.csv","r") as f_history:
        df=pd.read_csv(f_history)
        for _,row in df.iterrows():
            block_id=int(row["block_id"])
            abs_eb=float(row["abs_eb"])
            if (block_list is not None) and (block_id not in block_list[float(row["abs_eb"])]):
                continue
            br=float(row["br"])
            mse=float(row["mse"])
            num=global_param.blocks[block_id]["element_num"]*global_param.header.element_size
            if abs_eb not in br_dict:
                br_dict[abs_eb]=br*num
                mse_dict[abs_eb]=mse*num
                total_element_num[abs_eb]=num
            else:
                br_dict[abs_eb]+=br*num
                mse_dict[abs_eb]+=mse*num
                total_element_num[abs_eb]+=num
    for abs_eb in br_dict:
        br_dict[abs_eb]/=total_element_num[abs_eb]
        cr_dict[abs_eb]=8*global_param.header.word_size/br_dict[abs_eb]
        mse_dict[abs_eb]/=total_element_num[abs_eb]
        psnr_dict[abs_eb]=10*math.log10((global_param.data_range["P"][1]-global_param.data_range["P"][0])**2/mse_dict[abs_eb])
    br_dict=dict(sorted(br_dict.items(),reverse=True))
    cr_dict=dict(sorted(cr_dict.items(),reverse=True))
    mse_dict=dict(sorted(mse_dict.items(),reverse=True))
    psnr_dict=dict(sorted(psnr_dict.items(),reverse=True))
    return br_dict,cr_dict,mse_dict,psnr_dict