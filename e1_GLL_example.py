import yaml
import numpy as np
from plotly import graph_objects as go
from matplotlib import pyplot as plt

file_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"

with open(f"{file_path}.header.yaml","r") as f_header:
    header_and_block_info=yaml.load(f_header,Loader=yaml.UnsafeLoader)
    header=header_and_block_info["header"]
    block_num=header_and_block_info["block_num"]
    blocks=header_and_block_info["blocks"]
block=blocks[0]

with open(block["coordinate_path"],"rb") as f:
    coordinate=np.fromfile(f,dtype=header["np_dtype"]).reshape((3,-1)+header["element_shape"])
    print(f"coordinate shape: {coordinate.shape}")
coordinate=coordinate[:,1:2]

with open(block["pressure_path"],"rb") as f:
    pressure=np.fromfile(f,dtype=header["np_dtype"]).reshape((1,-1)+header["element_shape"])
    print(f"pressure shape: {pressure.shape}")
pressure=pressure[:,1:2]
center=pressure[0,0,4:6,4:6,4:6].sum()/8
pressure-=center
print(f"pressure max={pressure.max()}, min={pressure.min()}")

fig=go.Figure()
fig.add_trace(go.Scatter3d(x=coordinate[0,0].flatten(),y=coordinate[1,0].flatten(),z=coordinate[2,0].flatten(),mode="markers",marker=dict(size=1)))
fig.update_layout(scene=dict(aspectmode="data"))
fig.write_html("temp.html")

p=np.zeros((header["element_shape"][0],header["element_shape"][0]),dtype=header["np_dtype"])
x=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(header["element_shape"][0]-1).deriv().roots(),[1]))
p[0]=1
p[1]=x
for i in range(2,header["element_shape"][0]):
    p[i]=((2*i-1)*x*p[i-1]-(i-1)*p[i-2])/i
w=2/(header["element_shape"][0]*(header["element_shape"][0]-1)*p[-1]**2)

f=np.einsum("ibc,ai,i->abc",pressure[0,0],p,w)
f=np.einsum("aic,bi,i->abc",f,p,w)
f=np.einsum("abi,ci,i->abc",f,p,w)
print(f.max(),f.min())

fig,axes=plt.subplots(2,5,figsize=(50,20))
axes=axes.flatten()
for i in range(10):
    im=axes[i].imshow(f[i,:,:],origin="lower",cmap="jet",vmin=-3,vmax=3)
    axes[i].axis("off")
    fig.colorbar(im,ax=axes[i])
plt.savefig("temp1.png")
