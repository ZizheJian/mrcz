import os,random
import numpy as np
import pandas as pd
from py_funcs.global_param_class import global_param_class

def generate_test_case_dict(global_param:global_param_class,case_num_per_eb:int):
    existing_cases={abs_eb:[] for abs_eb in global_param.abs_eb_list}
    if os.path.exists(f"history/history{global_param.experiment_id:04d}.csv"):
        df=pd.read_csv(f"history/history{global_param.experiment_id:04d}.csv")
        for _,row in df.iterrows():
            existing_cases[row["abs_eb"]].append(int(row["block_id"]))
    ret={abs_eb:[] for abs_eb in global_param.abs_eb_list}
    for abs_eb in global_param.abs_eb_list:
        tested_blocks=existing_cases[abs_eb]
        if len(tested_blocks)>=case_num_per_eb:
            continue
        remaining_blocks=[i for i in range(global_param.block_num) if i not in tested_blocks]
        selected_blocks=random.sample(remaining_blocks,case_num_per_eb-len(tested_blocks))
        for block_id in selected_blocks:
            ret[abs_eb].append(block_id)
    return ret