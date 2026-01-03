import numpy as np
from typing import Tuple

def quantize_ndarray1(pred:np.ndarray,eb:float)->Tuple[np.ndarray,np.ndarray]:
    qb=np.round(pred/(2*eb)).astype(np.int64)
    return qb,2*eb*qb

def quantize_ndarray2(pred:np.ndarray,truth:np.ndarray,eb:float)->Tuple[np.ndarray,np.ndarray]:
    error=pred-truth
    qb,new_error=quantize_ndarray1(error,eb)
    return qb,pred-new_error

def quantize_ndarray3(pred:np.ndarray,eb:float,zero_freq_bound:float)->Tuple[np.ndarray,np.ndarray]:
    threshold=zero_freq_bound*eb
    qb=np.round(pred/(2*eb)).astype(np.int64)
    qb[np.abs(pred)<=threshold]=0
    return qb,2*eb*qb

def quantize_float1(pred:float,eb:float)->Tuple[int,float]:
    qb=int(round(pred/(2*eb)))
    return qb,2*eb*qb

def quantize_float2(pred:float,truth:float,eb:float)->Tuple[int,float]:
    error=pred-truth
    qb,new_error=quantize_float1(error,eb)
    return qb,pred-new_error