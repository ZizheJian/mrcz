import os
import numpy as np
import pandas as pd
from py_funcs.global_param_class import global_param_class
from py_funcs.read_nekrs_header import nekrs_header_c
from py_funcs.run_1block_1eb import run_1block_1eb
from py_funcs.run_multiblock_1eb import run_multiblock_1eb

def get_global_param(global_param:global_param_class):
    # 如果有imcomplete的记录，则读取该记录
    # 否则，如果有best（包括baseline）的记录，则读取该记录，并修正current_state和experiment_id
    print("Reading latest conclusion...")
    global_param.initialize_param()
    if os.path.exists("history/conclusion0000.txt"):
        index=1
        best_index=0
        incomplete_index=None
        while os.path.exists(f"history/conclusion{index:04d}.txt"):
            with open(f"history/conclusion{index:04d}.txt","r") as f_conclusion:
                lines=f_conclusion.readlines()
                if "best" in lines[0]:
                    best_index=index
                if "incomplete" in lines[0]:
                    incomplete_index=index
            index+=1
        next_index=index
        if incomplete_index is not None:
            index=incomplete_index
        else:
            index=best_index
        global_param.read(index)
        if global_param.conclusion=="baseline":
            global_param.current_state=0
            global_param.experiment_id=0
        if global_param.conclusion=="best":
            global_param.current_state=0
            global_param.experiment_id=next_index

def generate_baseline_history(global_param:global_param_class):
    if not os.path.exists("history"):
        os.makedirs("history", exist_ok=True)
    if not os.path.exists("history/history0000.csv"):
        with open(f"history/history0000.csv","w") as f_history:
            df=pd.DataFrame(columns=["block_id","abs_eb","br","cr","mse","psnr"])
            df.to_csv(f_history,index=False)
    for abs_eb in global_param.abs_eb_list:
        df=pd.read_csv(f"history/history0000.csv")
        tested_blocks=[row["block_id"] for _,row in df.iterrows() if row["abs_eb"]==abs_eb]
        remaining_blocks=[block_id for block_id in global_param.block_list if block_id not in tested_blocks]
        if len(remaining_blocks)==0:
            continue
        if global_param.experiment_id!=0:
            raise ValueError("生成baseline history时使用的不是初始参数")
        print(f"Running blocks {remaining_blocks} with eb {abs_eb}")
        br_list,cr_list,mse_list,psnr_list=run_multiblock_1eb(global_param,remaining_blocks,abs_eb)
        df=append_history(df,remaining_blocks,abs_eb,br_list,cr_list,mse_list,psnr_list)
        df.to_csv(f"history/history0000.csv",index=False)
    global_param.write("baseline")

def append_history(df:pd.DataFrame,block_list:list,abs_eb:float,br_list:list,cr_list:list,mse_list:list,psnr_list:list):
    for i,block_id in enumerate(block_list):
        new_row=pd.DataFrame([[block_id,abs_eb,br_list[i],cr_list[i],mse_list[i],psnr_list[i]]],columns=["block_id","abs_eb","br","cr","mse","psnr"])
        if df.empty:
            df=new_row
        else:
            df=pd.concat([df,new_row],ignore_index=True)
    df=df.sort_values(by=["block_id","abs_eb"],ascending=[True,False])
    return df