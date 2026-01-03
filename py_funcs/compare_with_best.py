import os,math
import pandas as pd
from py_funcs.global_param_class import global_param_class
from py_funcs.plot import plot_curves

def distance_point_line(cr:float,psnr:float,cr1:float,psnr1:float,cr2:float,psnr2:float)->float:
    dx=cr2-cr1
    dy=psnr2-psnr1
    if dx==0 and dy==0:
        return math.sqrt((cr-cr1)**2+(psnr-psnr1)**2)
    # 投影参数 t=((P-A)·(B-A))/|B-A|^2
    t=((cr-cr1)*dx+(psnr-psnr1)*dy)/(dx**2+dy**2)
    # 夹紧 t 到 [0,1]
    if t<0:
        # 最近到 A
        return math.sqrt((cr-cr1)**2+(psnr-psnr1)**2)
    elif t>1:
        # 最近到 B
        return math.sqrt((cr-cr2)**2+(psnr-psnr2)**2)
    else:
        # 投影在线段内，用直线距离
        proj_cr=cr1+t*dx
        proj_psnr=psnr1+t*dy
        return math.sqrt((cr-proj_cr)**2+(psnr-proj_psnr)**2)

def judge_point_curve(cr:float,psnr:float,reference_cr_list:list[float],reference_psnr_list:list[float])->str:
    dist_list=[distance_point_line(cr,psnr,reference_cr_list[i],reference_psnr_list[i],reference_cr_list[i+1],reference_psnr_list[i+1]) for i in range(len(reference_cr_list)-1)]
    min_line_index=dist_list.index(min(dist_list))
    A=(reference_cr_list[min_line_index],reference_psnr_list[min_line_index])
    B=(reference_cr_list[min_line_index+1],reference_psnr_list[min_line_index+1])
    P=(cr,psnr)
    if A[0]>B[0] or A[1]>B[1] or (A[0]==B[0] and A[1]==B[1]):
        return "uncertain"
    cross=(B[0]-A[0])*(P[1]-A[1])-(B[1]-A[1])*(P[0]-A[0])
    if cross>0:
        return "better"
    elif cross<0:
        return "worse"
    else:
        return "uncertain"

def compare_with_best(global_param:global_param_class):
    index=0
    best_index=0
    while os.path.exists(f"history/conclusion{index:04d}.txt"):
        with open(f"history/conclusion{index:04d}.txt","r") as f:
            lines=f.readlines()
            if "best" in lines[0]:
                best_index=index
        index+=1
    block_list=[]
    br_list=[]
    best_br_list=[]
    psnr_list=[]
    best_psnr_list=[]
    ######## 计算新实验的cr、psnr并记录每个abs_eb对应的blocks ########
    for abs_eb in global_param.abs_eb_list:
        df=pd.read_csv(f"history/history{global_param.experiment_id:04d}.csv")
        rows=[row for _,row in df.iterrows() if (row["abs_eb"]==abs_eb)]
        br=0
        mse=0
        total_num=0
        block_list.append([])
        for row in rows:
            block_list[-1].append(int(row["block_id"]))
            num=global_param.blocks[int(row["block_id"])]["element_num"]*global_param.header.element_size    
            br+=row["br"]*num
            mse+=row["mse"]*num
            total_num+=num
        br/=total_num
        mse/=total_num
        br_list.append(br)
        psnr_list.append(10*math.log10((global_param.data_range["P"][1]-global_param.data_range["P"][0])**2/mse))
    ######## 计算史上最佳实验在记录的blocks上的cr、psnr ########
    for abs_eb in global_param.abs_eb_list:
        df=pd.read_csv(f"history/history{best_index:04d}.csv")
        rows=[row for _,row in df.iterrows() if (row["abs_eb"]==abs_eb) and (row["block_id"] in block_list[global_param.abs_eb_list.index(abs_eb)])]
        br=0
        mse=0
        total_num=0
        for row in rows:
            num=global_param.blocks[int(row["block_id"])]["element_num"]*global_param.header.element_size    
            br+=row["br"]*num
            mse+=row["mse"]*num
            total_num+=num
        br/=total_num
        mse/=total_num
        best_br_list.append(br)
        best_psnr_list.append(10*math.log10((global_param.data_range["P"][1]-global_param.data_range["P"][0])**2/mse))
    
    plot_curves(global_param.experiment_id,[(f"{global_param.experiment_id:04d}",br_list,psnr_list),(f"Best{best_index:04d}",best_br_list,best_psnr_list)])

    if br_list[0]>=br_list[-1] or psnr_list[0]>=psnr_list[-1]:
        return "uncertain"
    score=0
    for i in range(len(br_list)):
        res=judge_point_curve(br_list[i],psnr_list[i],best_br_list,best_psnr_list)
        print(res)
        if res=="better":
            score+=1
        elif res=="worse":
            score-=1
    for i in range(len(best_br_list)):
        res=judge_point_curve(best_br_list[i],best_psnr_list[i],br_list,psnr_list)
        print(res)
        if res=="better":
            score-=1
        elif res=="worse":
            score+=1
    if score>2*len(global_param.abs_eb_list)/3:
        return "better"
    elif score<-2*len(global_param.abs_eb_list)/3:
        return "worse"
    else:
        return "uncertain"
        