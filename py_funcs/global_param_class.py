import os,copy,random
import numpy as np
import pandas as pd
from typing import Any
from py_funcs.read_nekrs_header import nekrs_header_c

class global_param_class:
    def __init__(self):
        self.original_data_path:str=""
        if os.environ.get("SLURM_NTASKS") is not None:
            self.thread_num=int(os.environ.get("SLURM_NTASKS"))
        else:
            self.thread_num=1
        np.set_printoptions(precision=15,suppress=True)
        self.header:nekrs_header_c=None
        self.data_range:dict[str,list[list[float]]]=None
        self.block_num:int=None
        self.blocks:list[dict[str,Any]]=None
        self.block_list:list[int]=None
        self.abs_eb_list:list[float]=None

    def initialize_param(self):
        self.current_state:int=0
        self.experiment_id:int=0
        self.param:dict[str,Any]={}
        self.param["param_eb_scale"]=0.05 # 正数
        self.param["zero_freq_bound"]=2 # 正数。量化coefficients时量化值为0的阈值
        self.x=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(self.header.element_shape[0]-1).deriv().roots(),[1])) # 不变
        self.param["L"]=np.zeros((self.header.element_shape[0],self.header.element_shape[0]),dtype=np.float64) # 实数可逆矩阵
        self.param["L"][0]=1
        self.param["L"][1]=self.x
        for i in range(2,self.header.element_shape[0]):
            self.param["L"][i]=2*self.x*self.param["L"][i-1]-self.param["L"][i-2]
        L_norm=(self.param["L"]**2).sum()**0.5
        self.param["param_eb_scale"]*=(L_norm**3)
        self.param["L"]/=L_norm
        self.param["neg_L"]=np.linalg.inv(self.param["L"])
        self.param["element_base_weight"]=np.ones((8,*self.header.element_shape),dtype=np.float64) # 不变
        self.param["element_base_weight"][0]=np.einsum("a,b,c->abc",1-self.x,1-self.x,1-self.x)/8
        self.param["element_base_weight"][1]=np.einsum("a,b,c->abc",1-self.x,1-self.x,self.x+1)/8
        self.param["element_base_weight"][2]=np.einsum("a,b,c->abc",1-self.x,self.x+1,1-self.x)/8
        self.param["element_base_weight"][3]=np.einsum("a,b,c->abc",1-self.x,self.x+1,self.x+1)/8
        self.param["element_base_weight"][4]=np.einsum("a,b,c->abc",self.x+1,1-self.x,1-self.x)/8
        self.param["element_base_weight"][5]=np.einsum("a,b,c->abc",self.x+1,1-self.x,self.x+1)/8
        self.param["element_base_weight"][6]=np.einsum("a,b,c->abc",self.x+1,self.x+1,1-self.x)/8
        self.param["element_base_weight"][7]=np.einsum("a,b,c->abc",self.x+1,self.x+1,self.x+1)/8
        
        self.param["f_qb_postprocessing_method"]="None" # "None"/"bitmap"
    
    def read(self,file_index:int):
        with open(f"history/conclusion{file_index:04d}.txt","r") as f_conclusion:
            lines=f_conclusion.readlines()
            self.conclusion=lines.pop(0).split("=")[1].strip()
            self.current_state=int(lines.pop(0).split("=")[1].strip())
            self.experiment_id=file_index
            while lines[0].startswith("abs_eb="):
                lines.pop(0)
            self.param["param_eb_scale"]=float(lines.pop(0).split("=")[1].strip())
            self.param["zero_freq_bound"]=float(lines.pop(0).split("=")[1].strip())
            for i in range(self.header.element_shape[0]):
                if i==0:
                    line=lines.pop(0).split("=")[1]
                else:
                    line=lines.pop(0)
                self.param["L"][i]=list(map(float,line.replace("[","").replace("]","").strip().split()))
            self.param["f_qb_postprocessing_method"]=lines.pop(0).split("=")[1].strip()

    def write(self,conclusion):
        os.makedirs("history",exist_ok=True)
        br_dict={}
        cr_dict={}
        mse_dict={}
        psnr_dict={}
        if os.path.exists(f"history/history{self.experiment_id:04d}.csv"):
            total_element_num={}
            with open(f"history/history{self.experiment_id:04d}.csv","r") as f_history:
                df=pd.read_csv(f_history)
                for _,row in df.iterrows():
                    block_id=int(row["block_id"])
                    abs_eb=float(row["abs_eb"])
                    br=float(row["br"])
                    mse=float(row["mse"])
                    if abs_eb not in br_dict:
                        br_dict[abs_eb]=br*self.blocks[block_id]["element_num"]*self.header.element_size
                        mse_dict[abs_eb]=mse*self.blocks[block_id]["element_num"]*self.header.element_size
                        total_element_num[abs_eb]=self.blocks[block_id]["element_num"]*self.header.element_size
                    else:
                        br_dict[abs_eb]+=br*self.blocks[block_id]["element_num"]*self.header.element_size
                        mse_dict[abs_eb]+=mse*self.blocks[block_id]["element_num"]*self.header.element_size
                        total_element_num[abs_eb]+=self.blocks[block_id]["element_num"]*self.header.element_size
            for abs_eb in br_dict:
                br_dict[abs_eb]/=total_element_num[abs_eb]
                cr_dict[abs_eb]=self.header.word_size/br_dict[abs_eb]
                mse_dict[abs_eb]/=total_element_num[abs_eb]
                psnr_dict[abs_eb]=10*np.log10((self.data_range["P"][1]-self.data_range["P"][0])**2/mse_dict[abs_eb])
        with open(f"history/conclusion{self.experiment_id:04d}.txt","w") as f_conclusion:
            content=f"conclusion= {conclusion}\n"
            content+=f"state= {self.current_state}\n"
            for abs_eb in self.abs_eb_list:
                if abs_eb in br_dict:
                    content+=f"abs_eb= {abs_eb} \t br= {br_dict[abs_eb]} \t cr= {cr_dict[abs_eb]} \t mse= {mse_dict[abs_eb]} \t psnr= {psnr_dict[abs_eb]}\n"
            content+=f"param_eb_scale= {self.param['param_eb_scale']}\n"
            content+=f"zero_freq_bound= {self.param['zero_freq_bound']}\n"
            content+=f"L= {self.param['L']}\n"
            content+=f"f_qb_postprocessing_method= {self.param['f_qb_postprocessing_method']}\n"
            f_conclusion.write(content)

def generate_new_global_param(global_param:global_param_class)->global_param_class:
    index=0
    while os.path.exists(f"history/conclusion{index:04d}.txt"):
        index+=1
    conclusion_num=index
    if conclusion_num==1:
        new_global_param=copy.deepcopy(global_param)
        new_global_param.current_state=1
        new_global_param.experiment_id=conclusion_num
        new_global_param.param["param_eb_scale"]*=(1+random.uniform(-1e-3,1e-3))
        new_global_param.param["zero_freq_bound"]*=(1+random.uniform(-1e-3,1e-3))
        new_global_param.param["L"]+=np.random.uniform(-1e-3,1e-3,new_global_param.param["L"].shape)
        if random.random()<0.1:
            new_global_param.param["f_qb_postprocessing_method"]="bitmap" if new_global_param.param["f_qb_postprocessing_method"]=="None" else "None"
        return new_global_param
    else:
        raise NotImplementedError("Only support generating new parameters based on the initial parameters.")