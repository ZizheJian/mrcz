import numpy as np

def partition(x:np.ndarray,block_id:np.ndarray,indices:np.ndarray,block_num:int,block_element_num:int)->int:
    cur_x=x[indices]
    min_max=np.zeros((3,2))
    for i in range(3):
        min_max[i,0]=np.min(cur_x[:,i])
        min_max[i,1]=np.max(cur_x[:,i])
    # print(f"partitioning {min_max}")
    # print(f"element_num={cur_x.shape[0]}")
    # print(f"block_num={block_num}")
    # input()
    partition_dim=np.argmax(min_max[:,1]-min_max[:,0])
    mid=0.5*(min_max[partition_dim,0]+min_max[partition_dim,1])
    right_mask=cur_x[:,partition_dim]>mid
    right_indices=indices[right_mask]
    left_indices=indices[~right_mask]
    if len(left_indices)>block_element_num:
        block_num=partition(x,block_id,left_indices,block_num,block_element_num)
    block_id[right_indices]=block_num
    block_num+=1
    if len(right_indices)>block_element_num:
        block_num=partition(x,block_id,right_indices,block_num,block_element_num)
    # for i in range(block_num):
    #     print(f"block {i} has {(block_id==i).sum()} elements")
    # input()
    return block_num
    