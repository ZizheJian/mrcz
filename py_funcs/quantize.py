import numpy as np

def quantize_ndarray1(pred:np.ndarray,eb:float)->np.ndarray:
    qb=np.round(pred/(2*eb)).astype(np.int64)
    return qb,2*eb*qb

def quantize_ndarray2(pred:np.ndarray,truth:np.ndarray,eb:float)->np.ndarray:
    error=pred-truth
    qb,new_error=quantize_ndarray1(error,eb)
    return qb,pred-new_error

def quantize_float1(pred:float,eb:float):
    qb=int(round(pred/(2*eb)))
    return qb,2*eb*qb

def quantize_float2(pred:float,truth:float,eb:float):
    error=pred-truth
    qb,new_error=quantize_float1(error,eb)
    return qb,pred-new_error