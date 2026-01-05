import os,copy,random,itertools
import numpy as np
from typing import Any
from py_funcs.read_nekrs_header import nekrs_header_c
from py_funcs.read_rate_distortion import read_rate_distortion_from_history
from py_funcs.compare_with_best import compare_with_best

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
        br_dict,cr_dict,mse_dict,psnr_dict=read_rate_distortion_from_history(self.experiment_id,None,self)
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
    if global_param.experiment_id<=20:
        new_global_param=copy.deepcopy(global_param)
        new_global_param.current_state=1
        if global_param.experiment_id==1:
            new_global_param.param["param_eb_scale"]=0.05
            new_global_param.param["L"]=np.zeros((new_global_param.header.element_shape[0],new_global_param.header.element_shape[0]),dtype=np.float64)
            new_global_param.param["L"][0]=1
            new_global_param.param["L"][1]=new_global_param.x
            for i in range(2,new_global_param.header.element_shape[0]):
                new_global_param.param["L"][i]=((2*i-1)*new_global_param.x*new_global_param.param["L"][i-1]-(i-1)*new_global_param.param["L"][i-2])/i
            L_norm=(new_global_param.param["L"]**2).sum()**0.5
            new_global_param.param["param_eb_scale"]*=(L_norm**3)
            new_global_param.param["L"]/=L_norm
            new_global_param.param["neg_L"]=np.linalg.inv(new_global_param.param["L"])
        else:
            new_global_param.param["param_eb_scale"]*=(1+random.uniform(-1e-3,1e-3))
            new_global_param.param["zero_freq_bound"]*=(1+random.uniform(-1e-3,1e-3))
            new_global_param.param["L"]+=np.random.uniform(-1e-3,1e-3,new_global_param.param["L"].shape)
            L_norm=(new_global_param.param["L"]**2).sum()**0.5
            print(f"L_norm= {L_norm}")
            new_global_param.param["L"]/=L_norm
            new_global_param.param["neg_L"]=np.linalg.inv(new_global_param.param["L"])
        if random.random()<0.1:
            new_global_param.param["f_qb_postprocessing_method"]="bitmap" if new_global_param.param["f_qb_postprocessing_method"]=="None" else "None"
        return new_global_param
    else:
        index=0
        best_param_list=[]
        while os.path.exists(f"history/conclusion{index:04d}.txt"):
            temp_global_param=copy.deepcopy(global_param)
            temp_global_param.read(index)
            if temp_global_param.conclusion=="best":
                best_param_list.append(temp_global_param)
            index+=1
        best_param_list.reverse()
        monotonic=copy.deepcopy(global_param)
        max_monotonic=0
        for param_name in ["param_eb_scale","zero_freq_bound"]:
            monotonic.param[param_name]=0
            for param_index in range(len(best_param_list)-2):
                param0=best_param_list[param_index].param[param_name]
                param1=best_param_list[param_index+1].param[param_name]
                param2=best_param_list[param_index+2].param[param_name]
                if param0>=param1>=param2 or param0<=param1<=param2:
                    monotonic.param[param_name]+=1
                    continue
                break
            if monotonic.param[param_name]>max_monotonic:
                max_monotonic=monotonic.param[param_name]
        for pos in itertools.product(range(global_param.header.element_shape[0]),range(global_param.header.element_shape[0])):
            monotonic.param["L"][pos]=0
            for param_index in range(len(best_param_list)-2):
                param0=best_param_list[param_index].param["L"][pos]
                param1=best_param_list[param_index+1].param["L"][pos]
                param2=best_param_list[param_index+2].param["L"][pos]
                if param0>=param1>=param2 or param0<=param1<=param2:
                    monotonic.param["L"][pos]+=1
                    continue
                break
            if monotonic.param["L"][pos]>max_monotonic:
                max_monotonic=monotonic.param["L"][pos]
        lr_base=0.001*(0.9**global_param.experiment_id)*(1.1**(len(best_param_list)-1))
        print(f"global_param.experiment_id= {global_param.experiment_id}, len(best_param_list)= {len(best_param_list)}, lr_base= {lr_base}")
        print(f"max_monotonic= {max_monotonic}")
        print(f"monotonic.param[\"param_eb_scale\"]= {monotonic.param['param_eb_scale']}, monotonic.param[\"zero_freq_bound\"]= {monotonic.param['zero_freq_bound']}")
        print(f"monotonic.param[\"L\"]= {monotonic.param['L']}")
        monotonic.param["param_eb_scale"]=lr_base*0.5**(max_monotonic-monotonic.param["param_eb_scale"])
        monotonic.param["zero_freq_bound"]=lr_base*0.5**(max_monotonic-monotonic.param["zero_freq_bound"])
        for pos in zip(range(global_param.header.element_shape[0]),range(global_param.header.element_shape[0])):
            monotonic.param["L"][pos]=lr_base*0.5**(max_monotonic-monotonic.param["L"][pos])
        new_global_param=copy.deepcopy(global_param)
        new_global_param.current_state=1
        new_global_param.param["param_eb_scale"]*=(1+random.uniform(-monotonic.param["param_eb_scale"],monotonic.param["param_eb_scale"]))
        new_global_param.param["zero_freq_bound"]*=(1+random.uniform(-monotonic.param["zero_freq_bound"],monotonic.param["zero_freq_bound"]))
        for pos in zip(range(global_param.header.element_shape[0]),range(global_param.header.element_shape[0])):
            new_global_param.param["L"]+=random.uniform(-monotonic.param["L"][pos],monotonic.param["L"][pos])/global_param.header.element_shape[0]
        L_norm=(new_global_param.param["L"]**2).sum()**0.5
        print(f"L_norm= {L_norm}")
        new_global_param.param["L"]/=L_norm
        new_global_param.param["neg_L"]=np.linalg.inv(new_global_param.param["L"])
        if random.random()<0.1:
            new_global_param.param["f_qb_postprocessing_method"]="bitmap" if new_global_param.param["f_qb_postprocessing_method"]=="None" else "None"
        return new_global_param