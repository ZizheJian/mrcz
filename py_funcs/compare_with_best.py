import os,math
import pandas as pd
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

def judge_point_curve(abs_eb_list,cr:float,psnr:float,reference_cr_dict:dict,reference_psnr_dict:dict)->str:
    dist_list=[distance_point_line(cr,psnr,
                                   reference_cr_dict[abs_eb_list[i]],reference_psnr_dict[abs_eb_list[i]],
                                   reference_cr_dict[abs_eb_list[i+1]],reference_psnr_dict[abs_eb_list[i+1]]) for i in range(len(abs_eb_list)-1)]
    min_line_index=dist_list.index(min(dist_list))
    A=(reference_cr_dict[abs_eb_list[min_line_index]],reference_psnr_dict[abs_eb_list[min_line_index]])
    B=(reference_cr_dict[abs_eb_list[min_line_index+1]],reference_psnr_dict[abs_eb_list[min_line_index+1]])
    P=(cr,psnr)
    cross=(B[0]-A[0])*(P[1]-A[1])-(B[1]-A[1])*(P[0]-A[0])
    if A[0]==B[0] and A[1]==B[1]:
        if P[0]==A[0] and P[1]==A[1]:
            return "uncertain"
        if P[0]<=A[0] and P[1]>=A[1]:
            return "better"
        if P[0]>=A[0] and P[1]<=A[1]:
            return "worse"
        return "uncertain"
    if A[0]<=B[0] and A[1]<=B[1]:# AB朝向右上
        if cross>0:
            return "better"
        elif cross<0:
            return "worse"
        else:
            return "uncertain"
    if A[0]>=B[0] and A[1]<=B[1]:# AB朝向左上
        if P[1]>=A[1] and (cross>=0 or P[0]<=B[0]):
            return "better"
        if P[1]<=B[1] and (cross<=0 or P[0]>=A[0]):
            return "worse"
        return "uncertain"
    if A[0]<=B[0] and A[1]>=B[1]:# AB朝向右下
        if P[0]<=B[0] and (cross>=0 or P[1]>=A[1]):
            return "better"
        if P[0]>=A[0] and (cross<=0 or P[1]<=B[1]):
            return "worse"
        return "uncertain"
    if P[0]<=B[0] and P[1]>=A[1]:
        return "better"
    if P[0]>=A[0] and P[1]<=B[1]:
        return "worse"
    return "uncertain"

def compare_with_best(global_param):
    index=0
    best_index=0
    while os.path.exists(f"history/conclusion{index:04d}.txt"):
        with open(f"history/conclusion{index:04d}.txt","r") as f:
            lines=f.readlines()
            if "best" in lines[0]:
                best_index=index
        index+=1
    block_dict={}
    br_dict={}
    best_br_dict={}
    psnr_dict={}
    best_psnr_dict={}
    ######## 计算新实验的cr、psnr并记录每个abs_eb对应的blocks ########
    for abs_eb in global_param.abs_eb_list:
        df=pd.read_csv(f"history/history{global_param.experiment_id:04d}.csv")
        rows=[row for _,row in df.iterrows() if (row["abs_eb"]==abs_eb)]
        br=0
        mse=0
        total_num=0
        block_dict[abs_eb]=[]
        for row in rows:
            block_dict[abs_eb].append(int(row["block_id"]))
            num=global_param.blocks[int(row["block_id"])]["element_num"]*global_param.header.element_size    
            br+=row["br"]*num
            mse+=row["mse"]*num
            total_num+=num
        br/=total_num
        mse/=total_num
        br_dict[abs_eb]=br
        psnr_dict[abs_eb]=10*math.log10((global_param.data_range["P"][1]-global_param.data_range["P"][0])**2/mse)
    ######## 计算史上最佳实验在记录的blocks上的cr、psnr ########
    for abs_eb in global_param.abs_eb_list:
        df=pd.read_csv(f"history/history{best_index:04d}.csv")
        rows=[row for _,row in df.iterrows() if (row["abs_eb"]==abs_eb) and (row["block_id"] in block_dict[abs_eb])]
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
        best_br_dict[abs_eb]=br
        best_psnr_dict[abs_eb]=10*math.log10((global_param.data_range["P"][1]-global_param.data_range["P"][0])**2/mse)
    
    plot_curves(global_param.experiment_id,[(f"{global_param.experiment_id:04d}",br_dict,psnr_dict),(f"Best{best_index:04d}",best_br_dict,best_psnr_dict)])

    if br_dict[global_param.abs_eb_list[0]]>=br_dict[global_param.abs_eb_list[-1]] or psnr_dict[global_param.abs_eb_list[0]]>=psnr_dict[global_param.abs_eb_list[-1]]:
        return "uncertain"
    score=0
    for abs_eb in global_param.abs_eb_list:
        res=judge_point_curve(global_param.abs_eb_list,br_dict[abs_eb],psnr_dict[abs_eb],best_br_dict,best_psnr_dict)
        if res=="better":
            score+=1
        elif res=="worse":
            score-=1
    for abs_eb in global_param.abs_eb_list:
        res=judge_point_curve(global_param.abs_eb_list,best_br_dict[abs_eb],best_psnr_dict[abs_eb],br_dict,psnr_dict)
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
        