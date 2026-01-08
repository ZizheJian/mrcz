import yaml,time,os
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from py_funcs.read_nekrs_header import read_nekrs_header
from py_funcs.partition import partition

file_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"
# If no coordinates are present in the main file, use a separate coordinate file
coordinate_file_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"
read_cache_size=128*1024*1024
block_size=128*1024*1024

if os.environ.get("SLURM_NTASKS") is not None:
    thread_num=int(os.environ.get("SLURM_NTASKS"))
else:
    thread_num=1
f=open(file_path,"rb")
header=read_nekrs_header(f)
data_range={}
block_element_num=block_size//(header.element_size*header.word_size*4) # During each compression, using xyz and 1 variable
print(f"block_element_num={block_element_num}")

######## skip element map ########
element_map=f.seek(header.file_element_num*4,1)

######## partition ########
print("Generating partition...")
start_time=time.time()
x=np.zeros((header.file_element_num,3),dtype=header.np_dtype)
element2block=np.zeros((header.file_element_num,),dtype=np.int32)
if "X" in header.variables:
    file_pos=f.tell()
    elements_per_read=read_cache_size//(3*header.element_size*header.word_size)
    chunks=[(i,min(i+elements_per_read,header.file_element_num))
            for i in range(0,header.file_element_num,elements_per_read)]
    def read_coordinate_chunk_and_mean(args):
        start,end=args
        data_num=(end-start)*3*header.element_size
        with open(file_path,"rb") as f:
            f.seek(file_pos+start*3*header.element_size*header.word_size)
            cache=np.fromfile(f,dtype=header.np_dtype,count=data_num)
        return cache.reshape((-1,3)+header.element_shape).mean(axis=(2,3,4))
    with ProcessPoolExecutor(max_workers=thread_num) as exe:
        results=list(exe.map(read_coordinate_chunk_and_mean,chunks))
    x=np.concatenate(results,axis=0)
    block_num=partition(x,element2block,np.arange(header.file_element_num),1,block_element_num)
    block2element=[[] for _ in range(block_num)]
    for i in range(header.file_element_num):
        block2element[element2block[i]].append(i)
    f.seek(header.file_element_num*3*header.element_size*header.word_size,1)
else:
    raise NotImplementedError("Alternative coordinate file reading not implemented yet.")
print(f"time= {time.time()-start_time} s")
    
######## save partitioned element map ########
print(f"Saving element map...")
start_time=time.time()
f.seek(-header.file_element_num*(4+3*header.element_size*header.word_size),1)
element_map=np.fromfile(f,dtype=np.int32,count=header.file_element_num)
for i in range(block_num):
    indices=np.where(element2block==i)[0]
    with open(f"{file_path}.block_{i}.element_map","wb") as f_element_map:
        f_element_map.write(element_map[indices].tobytes())
print(f"time= {time.time()-start_time} s")

######## save partitioned coordinates ########
print(f"Saving coordinate...")
start_time=time.time()
if "X" in header.variables:
    data_range["X"]=[[float('inf'),float('-inf')] for _ in range(3)]
    file_pos=f.tell()
    block_per_round=thread_num*read_cache_size//(block_size*3//4)
    def read_and_save_coordinate_block(block_id):
        cache=np.zeros((len(block2element[block_id]),3)+header.element_shape,dtype=header.np_dtype)
        with open(file_path,"rb") as f:
            for cache_pos,element_id in enumerate(block2element[block_id]):
                f.seek(file_pos+element_id*3*header.element_size*header.word_size)
                cache2=np.fromfile(f,dtype=header.np_dtype,count=3*header.element_size)
                cache[cache_pos]=cache2.reshape((3,)+header.element_shape)
        cache=cache.transpose((1,0,2,3,4))
        with open(f"{file_path}.block_{block_id}.coordinate","wb") as f_coordinate:
            f_coordinate.write(cache.tobytes())
        data_range=[[float('inf'),float('-inf')] for _ in range(3)]
        for dim in range(3):
            data_range[dim][0]=float(np.min(cache[dim]))
            data_range[dim][1]=float(np.max(cache[dim]))
        return data_range
    for block_id_start in range(0,block_num,block_per_round):
        block_id_end=min(block_id_start+block_per_round,block_num)
        with ProcessPoolExecutor(max_workers=thread_num) as exe:
            results=list(exe.map(read_and_save_coordinate_block,range(block_id_start,block_id_end)))
        for dim in range(3):
            data_range["X"][dim][0]=min(data_range["X"][dim][0],min([r[dim][0] for r in results]))
            data_range["X"][dim][1]=max(data_range["X"][dim][1],max([r[dim][1] for r in results]))
    f.seek(header.file_element_num*3*header.element_size*header.word_size,1)
else:
    raise NotImplementedError("Alternative coordinate file saving not implemented yet.")
print(f"time= {time.time()-start_time} s")

######## save velocity ########
if "U" in header.variables:
    start_time=time.time()
    print(f"Saving velocity...")
    data_range["U"]=[[float('inf'),float('-inf')] for _ in range(3)]
    file_pos=f.tell()
    block_per_round=thread_num*read_cache_size//(block_size*3//4)
    def read_and_save_velocity_block(block_id):
        cache=np.zeros((len(block2element[block_id]),3)+header.element_shape,dtype=header.np_dtype)
        with open(file_path,"rb") as f:
            for cache_pos,element_id in enumerate(block2element[block_id]):
                f.seek(file_pos+element_id*3*header.element_size*header.word_size)
                cache2=np.fromfile(f,dtype=header.np_dtype,count=3*header.element_size)
                cache[cache_pos]=cache2.reshape((3,)+header.element_shape)
        cache=cache.transpose((1,0,2,3,4))
        with open(f"{file_path}.block_{block_id}.velocityx","wb") as f_coordinate:
            f_coordinate.write(cache[0].tobytes())
        with open(f"{file_path}.block_{block_id}.velocityy","wb") as f_coordinate:
            f_coordinate.write(cache[1].tobytes())
        with open(f"{file_path}.block_{block_id}.velocityz","wb") as  f_coordinate:
            f_coordinate.write(cache[2].tobytes())
        data_range=[[float('inf'),float('-inf')] for _ in range(3)]
        for dim in range(3):
            data_range[dim][0]=float(np.min(cache[dim]))
            data_range[dim][1]=float(np.max(cache[dim]))
        return data_range
    for block_id_start in range(0,block_num,block_per_round):
        block_id_end=min(block_id_start+block_per_round,block_num)
        with ProcessPoolExecutor(max_workers=thread_num) as exe:
            results=list(exe.map(read_and_save_velocity_block,range(block_id_start,block_id_end)))
        for dim in range(3):
            data_range["U"][dim][0]=min(data_range["U"][dim][0],min([r[dim][0] for r in results]))
            data_range["U"][dim][1]=max(data_range["U"][dim][1],max([r[dim][1] for r in results]))
    f.seek(header.file_element_num*3*header.element_size*header.word_size,1)
    print(f"time= {time.time()-start_time} s")

######## save P ########
if "P" in header.variables:
    print(f"Saving pressure...")
    start_time=time.time()
    data_range["P"]=[float('inf'),float('-inf')]
    file_pos=f.tell()
    block_per_round=thread_num*read_cache_size//(block_size//4)
    def read_and_save_pressure_block(block_id):
        cache=np.zeros((len(block2element[block_id]),1)+header.element_shape,dtype=header.np_dtype)
        with open(file_path,"rb") as f:
            for cache_pos,element_id in enumerate(block2element[block_id]):
                f.seek(file_pos+element_id*header.element_size*header.word_size)
                cache2=np.fromfile(f,dtype=header.np_dtype,count=header.element_size)
                cache[cache_pos]=cache2.reshape((1,)+header.element_shape)
        with open(f"{file_path}.block_{block_id}.pressure","wb") as f_pressure:
            f_pressure.write(cache.tobytes())
        return [float(np.min(cache)),float(np.max(cache))]
    for block_id_start in range(0,block_num,block_per_round):
        block_id_end=min(block_id_start+block_per_round,block_num)
        with ProcessPoolExecutor(max_workers=thread_num) as exe:
            results=list(exe.map(read_and_save_pressure_block,range(block_id_start,block_id_end)))
        data_range["P"][0]=min(data_range["P"][0],min([r[0] for r in results]))
        data_range["P"][1]=max(data_range["P"][1],max([r[1] for r in results]))
    f.seek(header.file_element_num*header.element_size*header.word_size,1)
    print(f"time= {time.time()-start_time} s")

######## save T ########
if "T" in header.variables:
    print(f"Saving temperature...")
    start_time=time.time()
    data_range["T"]=[float('inf'),float('-inf')]
    file_pos=f.tell()
    block_per_round=thread_num*read_cache_size//(block_size//4)
    def read_and_save_temperature_block(block_id):
        cache=np.zeros((len(block2element[block_id]),1)+header.element_shape,dtype=header.np_dtype)
        with open(file_path,"rb") as f:
            for cache_pos,element_id in enumerate(block2element[block_id]):
                f.seek(file_pos+element_id*header.element_size*header.word_size)
                cache2=np.fromfile(f,dtype=header.np_dtype,count=header.element_size)
                cache[cache_pos]=cache2.reshape((1,)+header.element_shape)
        with open(f"{file_path}.block_{block_id}.temperature","wb") as f_temperature:
            f_temperature.write(cache.tobytes())
        return [float(np.min(cache)),float(np.max(cache))]
    for block_id_start in range(0,block_num,block_per_round):
        block_id_end=min(block_id_start+block_per_round,block_num)
        with ProcessPoolExecutor(max_workers=thread_num) as exe:
            results=list(exe.map(read_and_save_temperature_block,range(block_id_start,block_id_end)))
        data_range["T"][0]=min(data_range["T"][0],min([r[0] for r in results]))
        data_range["T"][1]=max(data_range["T"][1],max([r[1] for r in results]))
    f.seek(header.file_element_num*header.element_size*header.word_size,1)
    print(f"time= {time.time()-start_time} s")

######## save passive scalars ########
print(f"Saving passive scalars...")
start_time=time.time()
data_range["S"]=[[float('inf'),float('-inf')] for _ in range(header.passive_scalar_num)]
for passive_scalar_id in range(header.passive_scalar_num):
    file_pos=f.tell()
    block_per_round=thread_num*read_cache_size//(block_size//4)
    def read_and_save_passive_scalar_block(block_id):
        cache=np.zeros((len(block2element[block_id]),1)+header.element_shape,dtype=header.np_dtype)
        with open(file_path,"rb") as f:
            for cache_pos,element_id in enumerate(block2element[block_id]):
                f.seek(file_pos+element_id*header.element_size*header.word_size)
                cache2=np.fromfile(f,dtype=header.np_dtype,count=header.element_size)
                cache[cache_pos]=cache2.reshape((1,)+header.element_shape)
        with open(f"{file_path}.block_{block_id}.s{passive_scalar_id:02d}","wb") as f_passive_scalar:
            f_passive_scalar.write(cache.tobytes())
        return [float(np.min(cache)),float(np.max(cache))]
    for block_id_start in range(0,block_num,block_per_round):
        block_id_end=min(block_id_start+block_per_round,block_num)
        with ProcessPoolExecutor(max_workers=thread_num) as exe:
            results=list(exe.map(read_and_save_passive_scalar_block,range(block_id_start,block_id_end)))
        data_range["S"][passive_scalar_id][0]=min(data_range["S"][passive_scalar_id][0],min([r[0] for r in results]))
        data_range["S"][passive_scalar_id][1]=max(data_range["S"][passive_scalar_id][1],max([r[1] for r in results]))
    f.seek(header.file_element_num*header.element_size*header.word_size,1)
print(f"time= {time.time()-start_time} s")

######## save header ########
print(f"Saving header...")
start_time=time.time()
data_range=dict(data_range)
header_and_block_info={}
header_and_block_info["header"]=header.__dict__
header_and_block_info["data_range"]=data_range
header_and_block_info["block_num"]=block_num
header_and_block_info["blocks"]=[]
for i in range(block_num):
    block_info={}
    block_info["block_id"]=i
    block_info["element_num"]=len(block2element[i])
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

    