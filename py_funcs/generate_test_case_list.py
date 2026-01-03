import os,random
import numpy as np
import pandas as pd
from py_funcs.global_param_class import global_param_class

def generate_test_case_list(global_param:global_param_class,case_num:int):
    if os.path.exists(f"history/history{global_param.experiment_id:04d}.csv"):
        df=pd.read_csv(f"history/history{global_param.experiment_id:04d}.csv")
        existing_cases=list((row["block_id"],row["abs_eb"]) for _,row in df.iterrows())
    else:
        existing_cases=[]
    case_num_per_eb=(case_num+len(global_param.abs_eb_list)-1)//len(global_param.abs_eb_list)
    case_num=case_num_per_eb*len(global_param.abs_eb_list)
    if len(existing_cases)>=case_num:
        return []
    ret=[]
    for abs_eb in global_param.abs_eb_list:
        tested_blocks=[block_id for block_id,eb in existing_cases if eb==abs_eb]
        if len(tested_blocks)>=case_num_per_eb:
            continue
        remaining_blocks=[i for i in range(global_param.block_num) if i not in tested_blocks]
        selected_blocks=random.sample(remaining_blocks,case_num_per_eb-len(tested_blocks))
        for block_id in selected_blocks:
            ret.append((block_id,abs_eb))
    return ret