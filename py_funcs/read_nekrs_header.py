import struct,re
import numpy as np
from typing import Tuple,Dict

class nekrs_header_c:
    word_size:int=4  # word size in bytes
    realtype:str="" # "f" for float32, "d" for float64
    np_dtype=np.float32
    element_shape:Tuple[int, ...]=(0,0,0)
    element_size:int=0
    element_num:int=0
    file_element_num:int=0 # number of elements in the file
    time:float=0.0 # current simulation time
    istep:int=0 # current time step
    file_id:int=0
    file_num:int=0 # number of files (only 1 file will be read at a time)
    variables:str=""
    passive_scalar_num:int=0
    endian_mode:str=""

def read_nek5000_header(f)->nekrs_header_c:
    header=f.read(132).split()
    if len(header)<12:
        raise IOError("Header of the file was too short.")
    etagb=f.read(4)
    etagL=int(struct.unpack("<f",etagb)[0]*1e5)/1e5
    etagB=int(struct.unpack(">f",etagb)[0]*1e5)/1e5
    if etagL==6.54321:
        endian_mode="<"
    elif etagB==6.54321:
        endian_mode=">"
    else:
        raise ValueError("Could not interpret endianness")
    ret=nekrs_header_c()
    ret.word_size=int(header[1])
    if ret.word_size==4:
        ret.realtype="f"
        ret.np_dtype=np.float32
    elif ret.word_size==8:
        ret.realtype="d"
        ret.np_dtype=np.float64
    else:
        raise ValueError("Unrecognized word size: {}".format(ret.word_size))
    ret.element_shape=(int(header[2]),int(header[3]),int(header[4]))
    ret.element_size=int(np.prod(ret.element_shape))
    ret.element_num=int(header[5])
    ret.file_element_num=int(header[6])
    ret.time=float(header[7])
    ret.istep=int(header[8])
    ret.file_id=int(header[9])
    ret.file_num=int(header[10])
    ret.variables=header[11].decode("utf-8")
    ret.passive_scalar_num=int(re.search(r"S(\d{2})",ret.variables).group(1))
    ret.endian_mode=endian_mode
    return ret