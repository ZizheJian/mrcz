import numpy as np

def predict_and_overwrite_linear_dim2(element,x,eb,previous_indices,new_indices,next_indices):
    qb=np.array([],dtype=np.int32)
    max_err=0
    mse=0
    for i in range(len(new_indices)):
        if new_indices[i] is None:
            continue
        x_left=x[previous_indices[i]]
        x_mid=x[new_indices[i]]
        x_right=x[previous_indices[i+1]]
        y_left=element[np.ix_(np.array(previous_indices),np.array(previous_indices),np.array([previous_indices[i]]))]
        y_mid=element[np.ix_(np.array(previous_indices),np.array(previous_indices),np.array([new_indices[i]]))]
        y_right=element[np.ix_(np.array(previous_indices),np.array(previous_indices),np.array([previous_indices[i+1]]))]
        h=(y_left*(x_right-x_mid)+y_right*(x_mid-x_left))/(x_right-x_left)
        err=h-y_mid
        local_qb=np.round(err/(2*eb)).astype(int)
        qb=np.concatenate((qb,local_qb.flatten()))
        element[np.ix_(np.array(previous_indices),np.array(previous_indices),np.array([new_indices[i]]))]=h-2*eb*local_qb
        max_err=max(max_err,np.abs(err-2*eb*local_qb).max())
        mse+=((err-2*eb*local_qb)**2).sum()
    return qb,max_err,mse

def predict_and_overwrite_linear_dim1(element,x,eb,previous_indices,new_indices,next_indices):
    qb=np.array([],dtype=np.int32)
    max_err=0
    mse=0
    for i in range(len(new_indices)):
        if new_indices[i] is None:
            continue
        x_left=x[previous_indices[i]]
        x_mid=x[new_indices[i]]
        x_right=x[previous_indices[i+1]]
        y_left=element[np.ix_(np.array(previous_indices),np.array([previous_indices[i]]),np.array(next_indices))]
        y_mid=element[np.ix_(np.array(previous_indices),np.array([new_indices[i]]),np.array(next_indices))]
        y_right=element[np.ix_(np.array(previous_indices),np.array([previous_indices[i+1]]),np.array(next_indices))]
        h=(y_left*(x_right-x_mid)+y_right*(x_mid-x_left))/(x_right-x_left)
        err=h-y_mid
        local_qb=np.round(err/(2*eb)).astype(int)
        qb=np.concatenate((qb,local_qb.flatten()))
        element[np.ix_(np.array(previous_indices),np.array([new_indices[i]]),np.array(next_indices))]=h-2*eb*local_qb
        max_err=max(max_err,np.abs(err-2*eb*local_qb).max())
        mse+=((err-2*eb*local_qb)**2).sum()
    return qb,max_err,mse


def predict_and_overwrite_linear_dim0(element,x,eb,previous_indices,new_indices,next_indices):
    qb=np.array([],dtype=np.int32)
    max_err=0
    mse=0
    for i in range(len(new_indices)):
        if new_indices[i] is None:
            continue
        x_left=x[previous_indices[i]]
        x_mid=x[new_indices[i]]
        x_right=x[previous_indices[i+1]]
        y_left=element[np.ix_(np.array([previous_indices[i]]),np.array(next_indices),np.array(next_indices))]
        y_mid=element[np.ix_(np.array([new_indices[i]]),np.array(next_indices),np.array(next_indices))]
        y_right=element[np.ix_(np.array([previous_indices[i+1]]),np.array(next_indices),np.array(next_indices))]
        h=(y_left*(x_right-x_mid)+y_right*(x_mid-x_left))/(x_right-x_left)
        err=h-y_mid
        local_qb=np.round(err/(2*eb)).astype(int)
        qb=np.concatenate((qb,local_qb.flatten()))
        element[np.ix_(np.array([new_indices[i]]),np.array(next_indices),np.array(next_indices))]=h-2*eb*local_qb
        max_err=max(max_err,np.abs(err-2*eb*local_qb).max())
        mse+=((err-2*eb*local_qb)**2).sum()
    return qb,max_err,mse