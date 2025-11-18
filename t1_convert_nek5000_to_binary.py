import yaml,time
import numpy as np
from plotly import graph_objects as go
from concurrent.futures import ProcessPoolExecutor
from read_nek5000_header import read_nek5000_header
from partition import partition

multi_thread=False
file_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"
# If no coordinates are present in the main file, use a separate coordinate file
coordinate_file_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"
read_cache_size=128*1024*1024
block_size=128*1024*1024

f=open(file_path,"rb")
header=read_nek5000_header(f)
block_element_num=block_size//(header.element_size*header.word_size*4) # During each compression, using xyz and 1 variable
print(f"block_element_num={block_element_num}")

######## skip element map ########
element_map=f.seek(header.file_element_num*4,1)
######## partition ########
print("Generating partition...")
start_time=time.time()
x=np.zeros((header.file_element_num,3),dtype=header.np_dtype)
block_id=np.zeros((header.file_element_num,),dtype=np.int32)
if "X" in header.variables:
    if multi_thread:
        file_pos=f.tell()
        def read_chunk(args):
            start,end=args
            size=(end-start)*3*header.element_size*header.word_size
            f=open(file_path,"rb")
            f.seek(file_pos+start*3*header.element_size*header.word_size)
            cache=f.read(size)
            arr=np.frombuffer(cache,dtype=header.np_dtype).reshape((-1,3)+header.element_shape)
            return arr.mean(axis=(2,3,4))
        elements_per_read=read_cache_size//(3*header.element_size*header.word_size)
        chunks=[(i,min(i+elements_per_read,header.file_element_num))
                for i in range(0,header.file_element_num,elements_per_read)]
        with ProcessPoolExecutor() as exe:
            results=list(exe.map(read_chunk,chunks))
        x=np.concatenate(results,axis=0)
    else:
        elements_per_read=read_cache_size//(3*header.element_size*header.word_size)
        for i in range(0,header.file_element_num,elements_per_read):
            elements_in_this_read=min(elements_per_read,(header.file_element_num-i))
            cache=f.read(elements_in_this_read*3*header.element_size*header.word_size)
            temp_x=np.frombuffer(cache,dtype=header.np_dtype).reshape((-1,3)+header.element_shape)
            temp_x=np.mean(temp_x,axis=(2,3,4))
            x[i:i+temp_x.shape[0],:]=temp_x
    block_num=partition(x,block_id,np.arange(header.file_element_num),1,block_element_num)
else:
    f2=open(coordinate_file_path,"rb")
    coordinate_header=read_nek5000_header(f2)
    f2.read(coordinate_header.file_element_num*4)
    raise NotImplementedError("Coordinate file reading not implemented yet.")
print(f"time= {time.time()-start_time} s")
exit()
    
######## save partitioned element map ########
print(f"Saving element map...")
start_time=time.time()
f.seek(-header.file_element_num*(4+3*header.element_size*header.word_size),1)
element_map=f.read(header.file_element_num*4)
element_map=np.frombuffer(element_map,dtype=np.int32)
for i in range(block_num):
    indices=np.where(block_id==i)[0]
    f_element_map=open(f"{file_path}.block_{i}.element_map","wb")
    f_element_map.write(element_map[indices].tobytes())
    f_element_map.close()
print(f"time= {time.time()-start_time} s")
######## save partitioned coordinates ########
print(f"Saving coordinate...")
start_time=time.time()
if "X" in header.variables:
    header.data_range["X"]=[[float('inf'),float('-inf')] for _ in range(3)]
    for i in range(block_num):
        f_coordinate=open(f"{file_path}.block_{i}.coordinate","wb")
        f_coordinate.close()
    elements_per_read=read_cache_size//(3*header.element_size*header.word_size)
    for i in range(0,header.file_element_num,elements_per_read):
        elements_in_this_read=min(elements_per_read,(header.file_element_num-i))
        cache=f.read(elements_in_this_read*3*header.element_size*header.word_size)
        cache=np.frombuffer(cache,dtype=header.np_dtype).reshape((-1,3)+header.element_shape)
        for dim in range(3):
            header.data_range["X"][dim][0]=min(header.data_range["X"][dim][0],float(np.min(cache[:,dim])))
            header.data_range["X"][dim][1]=max(header.data_range["X"][dim][1],float(np.max(cache[:,dim])))
        for j in range(elements_in_this_read):
            f_coordinate=open(f"{file_path}.block_{block_id[i+j]}.coordinate","ab")
            f_coordinate.write(cache[j].tobytes())
            f_coordinate.close()
    for i in range(block_num):
        f_coordinate=open(f"{file_path}.block_{i}.coordinate","rb")
        cache=f_coordinate.read()
        f_coordinate.close()
        cache=np.frombuffer(cache,dtype=header.np_dtype).reshape((-1,3)+header.element_shape)
        cache=cache.transpose((1,0,2,3,4))
        f_coordinate=open(f"{file_path}.block_{i}.coordinate","wb")
        f_coordinate.write(cache.tobytes())
        f_coordinate.close()
else:
    raise NotImplementedError("Coordinate file saving not implemented yet.")
print(f"time= {time.time()-start_time} s")
######## save velocity ########
print(f"Saving velocity...")
start_time=time.time()
if "U" in header.variables:
    header.data_range["U"]=[[float('inf'),float('-inf')] for _ in range(3)]
    for i in range(block_num):
        f_velocityx=open(f"{file_path}.block_{i}.velocityx","wb")
        f_velocityy=open(f"{file_path}.block_{i}.velocityy","wb")
        f_velocityz=open(f"{file_path}.block_{i}.velocityz","wb")
        f_velocityx.close()
        f_velocityy.close()
        f_velocityz.close()
    elements_per_read=read_cache_size//(3*header.element_size*header.word_size)
    for i in range(0,header.file_element_num,elements_per_read):
        elements_in_this_read=min(elements_per_read,(header.file_element_num-i))
        cache=f.read(elements_in_this_read*3*header.element_size*header.word_size)
        cache=np.frombuffer(cache,dtype=header.np_dtype).reshape((-1,3)+header.element_shape)
        for dim in range(3):
            header.data_range["U"][dim][0]=min(header.data_range["U"][dim][0],float(np.min(cache[:,dim])))
            header.data_range["U"][dim][1]=max(header.data_range["U"][dim][1],float(np.max(cache[:,dim])))
        for j in range(elements_in_this_read):
            f_velocityx=open(f"{file_path}.block_{block_id[i+j]}.velocityx","ab")
            f_velocityx.write(cache[j,0].tobytes())
            f_velocityx.close()
            f_velocityy=open(f"{file_path}.block_{block_id[i+j]}.velocityy","ab")
            f_velocityy.write(cache[j,1].tobytes())
            f_velocityy.close()
            f_velocityz=open(f"{file_path}.block_{block_id[i+j]}.velocityz","ab")
            f_velocityz.write(cache[j,2].tobytes())
            f_velocityz.close()
print(f"time= {time.time()-start_time} s")
######## save P and T ########
print(f"Saving pressure and temperature...")
start_time=time.time()
if "P" in header.variables:
    header.data_range["P"]=[float('inf'),float('-inf')]
    for i in range(block_num):
        f_pressure=open(f"{file_path}.block_{i}.pressure","wb")
        f_pressure.close()
    elements_per_read=read_cache_size//(header.element_size*header.word_size)
    for i in range(0,header.file_element_num,elements_per_read):
        elements_in_this_read=min(elements_per_read,(header.file_element_num-i))
        cache=f.read(elements_in_this_read*header.element_size*header.word_size)
        cache=np.frombuffer(cache,dtype=header.np_dtype).reshape((-1,1)+header.element_shape)
        header.data_range["P"][0]=min(header.data_range["P"][0],float(np.min(cache)))
        header.data_range["P"][1]=max(header.data_range["P"][1],float(np.max(cache)))
        for j in range(elements_in_this_read):
            f_pressure=open(f"{file_path}.block_{block_id[i+j]}.pressure","ab")
            f_pressure.write(cache[j].tobytes())
            f_pressure.close()
if "T" in header.variables:
    header.data_range["T"]=[float('inf'),float('-inf')]
    for i in range(block_num):
        f_temperature=open(f"{file_path}.block_{i}.temperature","wb")
        f_temperature.close()
    elements_per_read=read_cache_size//(header.element_size*header.word_size)
    for i in range(0,header.file_element_num,elements_per_read):
        elements_in_this_read=min(elements_per_read,(header.file_element_num-i))
        cache=f.read(elements_in_this_read*header.element_size*header.word_size)
        cache=np.frombuffer(cache,dtype=header.np_dtype).reshape((-1,1)+header.element_shape)
        header.data_range["T"][0]=min(header.data_range["T"][0],float(np.min(cache)))
        header.data_range["T"][1]=max(header.data_range["T"][1],float(np.max(cache)))
        for j in range(elements_in_this_read):
            f_temperature=open(f"{file_path}.block_{block_id[i+j]}.temperature","ab")
            f_temperature.write(cache[j].tobytes())
            f_temperature.close()
print(f"time= {time.time()-start_time} s")
######## save passive scalars ########
print(f"Saving passive scalars...")
start_time=time.time()
header.data_range["S"]=[[float('inf'),float('-inf')] for _ in range(header.passive_scalar_num)]
for passive_scalar_id in range(header.passive_scalar_num):
    for i in range(block_num):
        f_passive_scalar=open(f"{file_path}.block_{i}.s{passive_scalar_id:02d}","wb")
        f_passive_scalar.close()
    elements_per_read=read_cache_size//(header.element_size*header.word_size)
    for i in range(0,header.file_element_num,elements_per_read):
        elements_in_this_read=min(elements_per_read,(header.file_element_num-i))
        cache=f.read(elements_in_this_read*header.element_size*header.word_size)
        cache=np.frombuffer(cache,dtype=header.np_dtype).reshape((-1,1)+header.element_shape)
        header.data_range["S"][passive_scalar_id][0]=min(header.data_range["S"][passive_scalar_id][0],float(np.min(cache)))
        header.data_range["S"][passive_scalar_id][1]=max(header.data_range["S"][passive_scalar_id][1],float(np.max(cache)))
        for j in range(elements_in_this_read):
            f_passive_scalar=open(f"{file_path}.block_{block_id[i+j]}.s{passive_scalar_id:02d}","ab")
            f_passive_scalar.write(cache[j].tobytes())
            f_passive_scalar.close()
print(f"time= {time.time()-start_time} s")
######## save header ########
print(f"Saving header...")
start_time=time.time()
header.data_range=dict(header.data_range)
header_and_block_info={}
header_and_block_info["header"]=header.__dict__
header_and_block_info["block_num"]=block_num
header_and_block_info["blocks"]=[]
for i in range(block_num):
    block_info={}
    block_info["block_id"]=i
    block_info["element_map_path"]=f"{file_path}.block_{i}.element_map"
    block_info["coordinate_path"]=f"{file_path}.block_{i}.coordinate"
    if "U" in header.variables:
        block_info["velocityx_path"]=f"{file_path}.block_{i}.velocityx"
        block_info["velocityy_path"]=f"{file_path}.block_{i}.velocityy"
        block_info["velocityz_path"]=f"{file_path}.block_{i}.velocityz"
    if "P" in header.variables:
        block_info["pressure_path"]=f"{file_path}.block_{i}.pressure"
    if "T" in header.variables:
        block_info["temperature_path"]=f"{file_path}.block_{i}.temperature"
    for ps_id in range(header.passive_scalar_num):
        block_info[f"s{ps_id:02d}_path"]=f"{file_path}.block_{i}.s{ps_id:02d}"
    header_and_block_info["blocks"].append(block_info)
with open(f"{file_path}.header.yaml","w") as f_header:
    yaml.dump(header_and_block_info,f_header,sort_keys=False)
print(f"time= {time.time()-start_time} s")

    