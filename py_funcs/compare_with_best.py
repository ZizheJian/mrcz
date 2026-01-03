import os,math
import pandas as pd
from py_funcs.global_param_class import global_param_class

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
    
    
    


def compare_with_best(global_param:global_param_class):
    index=0
    best_index=0
    while os.path.exists(f"history/conclusion{index:04d}.txt"):
        with open(f"history/conclusion{index:04d}.txt","r") as f:
            lines=f.readlines()
            if "best" in lines[0]:
                best_index=index
        index+=1
    cr_list=[]
    best_cr_list=[]
    psnr_list=[]
    best_psnr_list=[]
    for abs_eb in global_param.abs_eb_list:
        df=pd.read_csv(f"history/history{best_index:04d}.csv")
        rows=[row for _,row in df.iterrows() if (row["abs_eb"]==abs_eb) and (row["block_id"] in global_param.block_list)]
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
        best_cr_list.append(global_param.header.word_size/br)
        best_psnr_list.append(10*math.log10((global_param.data_range["P"][1]-global_param.data_range["P"][0])**2/mse))
    with open(f"history/conclusion{global_param.experiment_id:04d}.txt","r") as f:
        lines=f.readlines()[2:]
        for line in lines:
            if not line.startswith("abs_eb="):
                break
            cr_list.append(float(line.split()[5].strip()))
            psnr_list.append(float(line.split()[9].strip()))
    print(best_cr_list)
    print(cr_list)
    print(best_psnr_list)
    print(psnr_list)
    exit()
    if cr_list[0]>=cr_list[-1] or psnr_list[0]>=psnr_list[-1]:
        return "uncertain"
    score=0
    for i in range(len(cr_list)):
        judge_point_curve(cr_list[i],psnr_list[i],best_cr_list,best_psnr_list)
        