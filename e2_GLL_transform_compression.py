import yaml,struct,os,time
import numpy as np
import zstandard as zstd
from dahuffman import HuffmanCodec
from concurrent.futures import ProcessPoolExecutor
from matplotlib import pyplot as plt
from plotly import graph_objects as go
from read_nek5000_header import nek5000_header_c
from py_funcs.quantize_ndarray import quantize_ndarray1,quantize_ndarray2

file_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"
# rel_eb=1e-5
abs_eb=1e-6
param_eb_scale=0.1414
zstd_separate_compress=False

if os.environ.get("SLURM_NTASKS") is not None:
    thread_num=int(os.environ.get("SLURM_NTASKS"))
else:
    thread_num=1

with open(f"{file_path}.header.yaml","r") as f_header:
    header_and_block_info=yaml.load(f_header,Loader=yaml.UnsafeLoader)
    header=nek5000_header_c()
    header.__dict__.update(header_and_block_info["header"])
    block_num=header_and_block_info["block_num"]
    blocks=header_and_block_info["blocks"]

try:
    abs_eb
except:
    abs_eb=rel_eb*(header.data_range["P"][1]-header.data_range["P"][0]) # type: ignore
else:
    rel_eb=abs_eb/(header.data_range["P"][1]-header.data_range["P"][0])
param_eb=abs_eb*param_eb_scale
# param_eb=2*rel_eb*param_eb_scale
max_err=0
mse=0
print(f"presure_min= {header.data_range['P'][0]}, presure_max= {header.data_range['P'][1]}")
print(f"rel_eb= {rel_eb}")
print(f"abs_eb= {abs_eb}")
print(f"param_eb= {param_eb}")

p=np.zeros((header.element_shape[0],header.element_shape[0]),dtype=np.float64)
x=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(header.element_shape[0]-1).deriv().roots(),[1]))
p[0]=1
p[1]=x
for i in range(2,header.element_shape[0]):
    p[i]=((2*i-1)*x*p[i-1]-(i-1)*p[i-2])/i
w=2/(header.element_shape[0]*(header.element_shape[0]-1)*p[-1]**2)
d=np.einsum("ia,a,ia->i",p,w,p)
element_base_weight=np.ones((8,*header.element_shape),dtype=np.float64)
element_base_weight[0]=np.einsum("a,b,c->abc",1-x,1-x,1-x)/8
element_base_weight[1]=np.einsum("a,b,c->abc",1-x,1-x,x+1)/8
element_base_weight[2]=np.einsum("a,b,c->abc",1-x,x+1,1-x)/8
element_base_weight[3]=np.einsum("a,b,c->abc",1-x,x+1,x+1)/8
element_base_weight[4]=np.einsum("a,b,c->abc",x+1,1-x,1-x)/8
element_base_weight[5]=np.einsum("a,b,c->abc",x+1,1-x,x+1)/8
element_base_weight[6]=np.einsum("a,b,c->abc",x+1,x+1,1-x)/8
element_base_weight[7]=np.einsum("a,b,c->abc",x+1,x+1,x+1)/8

def compress_chunk(chunk):
    f_qb_list=[]
    qb_list=[]
    max_err=0
    mse=0
    for element_id in range(chunk.shape[1]):
        element=chunk[0,element_id]
        element_base=element[0,0,0]*element_base_weight[0]+element[0,0,-1]*element_base_weight[1]+\
                     element[0,-1,0]*element_base_weight[2]+element[0,-1,-1]*element_base_weight[3]+\
                     element[-1,0,0]*element_base_weight[4]+element[-1,0,-1]*element_base_weight[5]+\
                     element[-1,-1,0]*element_base_weight[6]+element[-1,-1,-1]*element_base_weight[7]
        element-=element_base
        f=np.einsum("abc,ia,a->ibc",element,p,w)/d.reshape(-1,1,1)
        f=np.einsum("ibc,jb,b->ijc",f,p,w)/d.reshape(1,-1,1)
        f=np.einsum("ijc,kc,c->ijk",f,p,w)/d.reshape(1,1,-1)
        f_qb,f_quantized=quantize_ndarray1(f,param_eb)
        f_qb[(-1<=f_qb) & (f_qb<=1)]=0
        f_quantized=f_qb*2*param_eb
        f_qb_list.append(f_qb)
        new_element=np.einsum("ijk,kc->ijc",f_quantized,p)
        new_element=np.einsum("ijc,jb->ibc",new_element,p)
        new_element=np.einsum("ibc,ia->abc",new_element,p)
        qb,new_element=quantize_ndarray2(new_element,element,abs_eb)
        qb_list.append(qb)
        max_err=max(max_err,np.abs(new_element-element).max())
        mse+=np.sum((new_element-element)**2)
    return f_qb_list,qb_list,max_err,mse

def lossless_encode_block(args):
    f_qb_list,qb_list=args
    f_qb=np.array(f_qb_list).transpose(1,2,3,0).flatten()
    qb=np.array(qb_list).flatten()
    ######## Huffman f ########
    f_hf_bs=b""
    f_values,f_counts=np.unique(f_qb,return_counts=True)
    f_freq=dict(zip(f_values,f_counts))
    f_hf_bs+=struct.pack("<Q",len(f_values))
    max_f=np.abs(f_values).max().item()
    if max_f<=2**15-1:
        f_hf_bs+=struct.pack("<b",15)
        f_hf_bs+=struct.pack(f"<{len(f_values)}h",*f_values.tolist())
    else:
        f_hf_bs+=struct.pack("<b",31)
        f_hf_bs+=struct.pack(f"<{len(f_values)}i",*f_values.tolist())
    f_hf_bs+=struct.pack(f"<{len(f_values)}Q",*f_counts.tolist())
    f_codec=HuffmanCodec.from_frequencies(f_freq)
    temp_bs=f_codec.encode(f_qb.flatten().tolist())
    f_hf_bs+=struct.pack("<Q",len(temp_bs))
    f_hf_bs+=temp_bs
    ######## Huffman qb ########
    hf_bs=b""
    qb_values,qb_counts=np.unique(qb,return_counts=True)
    qb_freq=dict(zip(qb_values,qb_counts))
    hf_bs+=struct.pack("<Q",len(qb_values))
    max_qb=np.abs(qb_values).max().item()
    if max_qb<=2**15-1:
        hf_bs+=struct.pack("<b",15)
        hf_bs+=struct.pack(f"<{len(qb_values)}h",*qb_values.tolist())
    else:
        hf_bs+=struct.pack("<b",31)
        hf_bs+=struct.pack(f"<{len(qb_values)}i",*qb_values.tolist())
    hf_bs+=struct.pack(f"<{len(qb_values)}Q",*qb_counts.tolist())
    qb_codec=HuffmanCodec.from_frequencies(qb_freq)
    temp_bs=qb_codec.encode(qb.flatten().tolist())
    hf_bs+=struct.pack("<Q",len(temp_bs))
    hf_bs+=temp_bs
    ######## Zstd ########
    cctx=zstd.ZstdCompressor()
    if zstd_separate_compress:
        zstd_f_qb_bs=cctx.compress(f_hf_bs)
        zstd_qb_bs=cctx.compress(hf_bs)
        return len(zstd_f_qb_bs),len(zstd_qb_bs)
    else:
        zstd_bs=cctx.compress(f_hf_bs+hf_bs)
        return len(zstd_bs)

qb_list=[]
with ProcessPoolExecutor(max_workers=thread_num) as exe:
    total_lossy_time=0
    for block in blocks:
        block_element_num=os.path.getsize(block["pressure_path"])//(header.word_size*header.element_size)
        with open(block["pressure_path"],"rb") as f_pressure:
            cache=np.fromfile(f_pressure,dtype=header.np_dtype,count=block_element_num*header.element_size).reshape(1,-1,*header.element_shape)
        start_element_id=[(block_element_num*i)//thread_num for i in range(thread_num)]
        end_element_id=[(block_element_num*(i+1))//thread_num for i in range(thread_num)]
        exe_input=[cache[:,start_element_id[i]:end_element_id[i]] for i in range(thread_num)]
        start_time=time.time()
        results=list(exe.map(compress_chunk,exe_input))
        total_lossy_time+=time.time()-start_time
        qb_list.append(([],[]))
        for local_f_qb_list,local_qb_list,local_max_err,local_mse in results:
            qb_list[-1][0].extend(local_f_qb_list)
            qb_list[-1][1].extend(local_qb_list)
            if local_max_err>max_err:
                max_err=local_max_err
            mse+=local_mse
    print(f"有损时间: {total_lossy_time} s")

    # for block_id,(f_qb_list,qb_list) in enumerate(qb_list):
    #     print(f"Block {block_id}")
    #     f_qb=np.array(f_qb_list)
    #     f_qb_0_freq=np.sum(f_qb==0,axis=0)/f_qb.shape[0]

    #     fig,axes=plt.subplots(2,5,figsize=(50,20))
    #     axes=axes.flatten()
    #     ims=[]
    #     for i in range(10):
    #         im=axes[i].imshow(f_qb_0_freq[i,:,:],origin="lower",cmap="gray",vmin=0,vmax=1)
    #         ims.append(im)
    #         axes[i].axis("off")
    #     cbar=fig.colorbar(ims[0],ax=axes,fraction=0.05,shrink=0.9,aspect=20,pad=0.05)
    #     cbar.set_ticks([0,1])
    #     cbar.set_ticklabels(["No Zero","All Zero"],fontsize=60)
    #     plt.savefig(f"fig/0_freq_{block_id}.png")
    #     plt.close()

    #     x=y=z=np.arange(header.element_shape[0])
    #     X,Y,Z=np.meshgrid(x,y,z,indexing="ij")
    #     fig=go.Figure()
    #     fig.add_trace(go.Volume(x=X.flatten(),y=Y.flatten(),z=Z.flatten(),value=f_qb_0_freq.flatten(),isomin=0,isomax=1,
    #                             opacity=0.1,surface_count=40,colorscale="Greys",showscale=True,
    #                             colorbar=dict(tickmode="array",tickvals=[0,1],ticktext=["No Zero","All Zero"])))
    #     fig.update_layout(scene=dict(aspectmode="cube"),width=800,height=800)
    #     fig.write_html(f"fig/0_freq_{block_id}.html")
    
    start_time=time.time()
    results=list(exe.map(lossless_encode_block,qb_list))
    print(f"无损时间: {time.time()-start_time} s")
    if zstd_separate_compress:
        zstd_f_qb_len=sum([res[0] for res in results])
        zstd_qb_len=sum([res[1] for res in results])
    else:
        zstd_len=sum(results)

if zstd_separate_compress:
    print(f"Zstd f_qb编码长度: {zstd_f_qb_len}")
    print(f"Zstd qb编码长度: {zstd_qb_len}")
    print(f"原始数据大小: {header.file_element_num*header.element_size*header.word_size} bytes")
    print(f"压缩比: {header.file_element_num*header.element_size*header.word_size/(zstd_f_qb_len+zstd_qb_len)}")
else:
    print(f"Zstd编码长度: {zstd_len}")
    print(f"原始数据大小: {header.file_element_num*header.element_size*header.word_size} bytes")
    print(f"压缩比: {header.file_element_num*header.element_size*header.word_size/zstd_len}")

print(f"max_err= {max_err}({max_err*100/abs_eb}%)")
mse/=(header.file_element_num*header.element_size)
psnr=10*np.log10((header.data_range["P"][1]-header.data_range["P"][0])**2/mse)
print(f"MSE= {mse} , PSNR= {psnr} dB")