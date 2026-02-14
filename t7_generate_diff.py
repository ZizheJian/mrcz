import numpy as np
from py_funcs.read_nekrs_header import read_nekrs_header

data_path1="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/to_transfer/pb0.f00001"
data_path2=data_path1+".sperr.bin"
diff_path=data_path2+".diff"

# f1=open(data_path1,"rb")
# f2=open(data_path2,"rb")
# f_diff=open(diff_path,"wb")
# header=read_nekrs_header(f1)
# f1.seek(0)
# header_binary=np.fromfile(f1,dtype=np.int8,count=136)
# np.fromfile(f2,dtype=np.int8,count=136)
# f_diff.write(header_binary.tobytes())

# element_map=np.fromfile(f1,dtype=np.int32,count=header.file_element_num)
# np.fromfile(f2,dtype=np.int32,count=header.file_element_num)
# f_diff.write(element_map.tobytes())

# X=np.fromfile(f1,dtype=np.float32,count=header.file_element_num*header.element_size*3)
# np.fromfile(f2,dtype=np.float32,count=header.file_element_num*header.element_size*3)
# f_diff.write(X.tobytes())

# fields1=np.fromfile(f1,dtype=np.float32,count=header.file_element_num*header.element_size*(3+1+1))
# fields2=np.fromfile(f2,dtype=np.float32,count=header.file_element_num*header.element_size*(3+1+1))
# diff=fields2-fields1
# diff=diff.reshape((5,header.file_element_num*header.element_size))
# print(np.max(np.abs(diff[4,:])))
# f_diff.write(diff.tobytes())

# f1.close()
# f2.close()
# f_diff.close()

f_diff=open(diff_path,"rb")
f_diff.seek(27922188496)
data=np.fromfile(f_diff,dtype=np.float32,count=996942848)
print(np.max(np.abs(data)))
f_diff.close()