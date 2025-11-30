import numpy as np
from py_funcs.quantize_ndarray import quantize_ndarray1,quantize_ndarray2

def interpolation_compress_element(element,element_shape,x,abs_eb):
    new_element=np.zeros_like(element,dtype=np.float64)
    qb=np.zeros_like(element,dtype=np.int32)
    if element_shape[0]%2==0:
        left=element_shape[0]//2-1
        right=element_shape[0]//2
        mask=np.zeros_like(element,dtype=bool)
        mask[left:right+1,left:right+1,left:right+1]=True
        center=element[mask].mean()
        new_element[mask]=center
        qb[mask],new_element[mask]=quantize_ndarray2(new_element[mask],element[mask],abs_eb)
        left-=1
        right+=1
    else:
        center=element[element_shape[0]//2,element_shape[1]//2,element_shape[2]//2]
        new_element[element_shape[0]//2,element_shape[1]//2,element_shape[2]//2]=center
        left=element_shape[0]//2-1
        right=element_shape[0]//2+1
    while left>=0:
        ########surface########
        new_element[left,left+1:right,left+1:right]=new_element[left+1,left+1:right,left+1:right]
        new_element[right,left+1:right,left+1:right]=new_element[right-1,left+1:right,left+1:right]
        new_element[left:left+1,left,left+1:right]=new_element[left:left+1,left+1,left+1:right]
        new_element[left:left+1,right,left+1:right]=new_element[left:left+1,right-1,left+1:right]
        new_element[left:left+1,left:left+1,left]=new_element[left:left+1,left:left+1,left+1]
        new_element[left:left+1,left:left+1,right]=new_element[left:left+1,left:left+1,right-1]
        ########edge########
        new_element[left,left,left+1:right]=new_element[left+1,left+1,left+1:right]
        new_element[right,left,left+1:right]=new_element[right-1,left+1,left+1:right]
        new_element[left,right,left+1:right]=new_element[left+1,right-1,left+1:right]
        new_element[right,right,left+1:right]=new_element[right-1,right-1,left+1:right]
        new_element[left,left+1:right,left]=new_element[left+1,left+1:right,left+1]
        new_element[right,left+1:right,left]=new_element[right-1,left+1:right,left+1]
        new_element[left,left+1:right,right]=new_element[left+1,left+1:right,right-1]
        new_element[right,left+1:right,right]=new_element[right-1,left+1:right,right-1]
        new_element[left:left+1,left,left]=new_element[left:left+1,left+1,left+1]
        new_element[left:left+1,right,left]=new_element[left:left+1,right-1,left+1]
        new_element[left:left+1,left,right]=new_element[left:left+1,left+1,right-1]
        new_element[left:left+1,right,right]=new_element[left:left+1,right-1,right-1]
        ########corner########
        new_element[left,left,left]=new_element[left+1,left+1,left+1]
        new_element[right,left,left]=new_element[right-1,left+1,left+1]
        new_element[left,right,left]=new_element[left+1,right-1,left+1]
        new_element[right,right,left]=new_element[right-1,right-1,left+1]
        new_element[left,left,right]=new_element[left+1,left+1,right-1]
        new_element[right,left,right]=new_element[right-1,left+1,right-1]
        new_element[left,right,right]=new_element[left+1,right-1,right-1]
        new_element[right,right,right]=new_element[right-1,right-1,right-1]
        ########quantize########
        mask=np.zeros_like(element,dtype=bool)
        mask[left,left:right+1,left:right+1]=True
        mask[right,left:right+1,left:right+1]=True
        mask[left:right+1,left,left:right+1]=True
        mask[left:right+1,right,left:right+1]=True
        mask[left:right+1,left:right+1,left]=True
        mask[left:right+1,left:right+1,right]=True
        qb[mask],new_element[mask]=quantize_ndarray2(new_element[mask],element[mask],abs_eb)
        left-=1
        right+=1
    mse=np.sum((new_element-element)**2)
    max_err=np.abs(new_element-element).max()
    return qb,max_err,mse

    
    