import numpy as np

def quantize_ndarray1(pred:np.ndarray,eb:float)->np.ndarray:
    mask_positive=(pred>=0)
    mask_negative=(~mask_positive)
    qb=np.zeros_like(pred,dtype=np.int32)
    qb[mask_positive]=np.ceil((pred[mask_positive]-eb)/(2*eb))
    qb[mask_negative]=np.floor((pred[mask_negative]+eb)/(2*eb))
    return qb,2*eb*qb

def quantize_ndarray2(pred:np.ndarray,truth:np.ndarray,eb:float)->np.ndarray:
    error=pred-truth
    qb,new_error=quantize_ndarray1(error,eb)
    return qb,pred-new_error