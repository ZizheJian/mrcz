# 有五种状态：
#   state=0: 要随机一组新参数，并进入state=1
#   state=1: 对于每个eb，在1个block上测试
#   state=2: 对于每个eb，在少数block上测试
#   state=3: 对于每个eb，在多数block上测试
#   state=4: 测试所有case
#   state=5: 结束该实验，返回state=0

import os,yaml,math
import numpy as np
import pandas as pd
from py_funcs.global_param_class import global_param_class,generate_new_global_param
from py_funcs.read_nekrs_header import nekrs_header_c
from py_funcs.history import select_global_param,generate_baseline_history,append_history
from py_funcs.generate_test_case_dict import generate_test_case_dict
from py_funcs.run_multiblock_1eb import run_multiblock_1eb
from py_funcs.compare_with_best import compare_with_best
np.set_printoptions(linewidth=1000)

global_param=global_param_class()
global_param.original_data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"

with open(f"{global_param.original_data_path}.header.yaml","r") as f_header:
    header_and_block_info=yaml.load(f_header,Loader=yaml.UnsafeLoader)
    global_param.header=nekrs_header_c()
    global_param.header.__dict__.update(header_and_block_info["header"])
    global_param.data_range=header_and_block_info["data_range"]
    global_param.block_num=header_and_block_info["block_num"]
    global_param.blocks=header_and_block_info["blocks"]

global_param.abs_eb_list=[3.162e-4,1e-4,3.162e-5,1e-5,3.162e-6,1e-6,3.162e-7,1e-7,3.162e-8,1e-8]
global_param.block_list=[i for i in range(global_param.block_num)]

def run_and_judge(global_param:global_param_class,case_num_per_eb:int):
    if not os.path.exists(f"history/history{global_param.experiment_id:04d}.csv"):
        df=pd.DataFrame(columns=["block_id","abs_eb","br","cr","mse","psnr"])
        df.to_csv(f"history/history{global_param.experiment_id:04d}.csv",index=False)
    test_case_dict=generate_test_case_dict(global_param,case_num_per_eb)
    for abs_eb in global_param.abs_eb_list:
        block_list=test_case_dict[abs_eb]
        print(f"Running blocks {block_list} with eb {abs_eb}")
        br_list,cr_list,mse_list,psnr_list=run_multiblock_1eb(global_param,block_list,abs_eb)
        df=pd.read_csv(f"history/history{global_param.experiment_id:04d}.csv")
        df=append_history(df,block_list,abs_eb,br_list,cr_list,mse_list,psnr_list)
        df.to_csv(f"history/history{global_param.experiment_id:04d}.csv",index=False)
    global_param.write("incomplete")
    res=compare_with_best(global_param)
    print(f"Result= {res}")
    if res=="better" or res=="uncertain":
        global_param.current_state+=1
        global_param.write("incomplete")
    else:
        global_param.write("worse")
        select_global_param(global_param)
    return global_param

select_global_param(global_param)
if global_param.experiment_id==0:
    generate_baseline_history(global_param)
    global_param.experiment_id+=1
    global_param.write("incomplete")
while True:
    print(f"Experiment ID: {global_param.experiment_id}, State: {global_param.current_state}")
    if global_param.current_state==0:
        global_param=generate_new_global_param(global_param)
        global_param.write("incomplete")
    elif global_param.current_state==1:
        global_param=run_and_judge(global_param,1)
    elif global_param.current_state==2:
        global_param=run_and_judge(global_param,int(math.ceil((global_param.block_num)**(1/3))))
    elif global_param.current_state==3:
        global_param=run_and_judge(global_param,int(math.ceil((global_param.block_num)**(2/3))))
    elif global_param.current_state==4:
        global_param=run_and_judge(global_param,global_param.block_num)
    elif global_param.current_state==5:
        global_param.write("best")
        # global_param.current_state=0
        # global_param.experiment_id+=1
        select_global_param(global_param)