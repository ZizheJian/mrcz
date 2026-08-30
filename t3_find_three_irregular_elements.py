# 功能：寻找并绘制三个两两共享面的非平面 NekRS element。
#
# 使用方法：
# 1. 修改data_path，使其指向包含坐标场的NekRS field文件。
# 2. 如果已经知道三个element ID，将它们写入selected_element_ids，程序会跳过搜索。
# 3. 如果需要重新搜索，将selected_element_ids设为None。程序顺序读取element，用量化后
#    的面心代表面，并在每加入一个element后检查是否形成三角邻接关系，找到后立即停止。
#    面心相同不一定严格表示两个面完全相同，因此应通过输出图人工确认搜索结果。
# 4. 运行python3 t3_find_irregular_elements.py。图片写入fig目录：three_elements.html
#    显示三个element，three_elements2.html显示中间的element，three_elements_ref.html
#    显示黑色reference element。
#
# 资源估计：
# - 已指定element ID时使用np.memmap按需读取，通常1 CPU、2 GiB内存足够。
# - 自动搜索如果较早找到结果，通常8--16 GiB内存足够。
# - 如果必须检查接近全部约324万个element，Python面索引可能需要数十GiB内存，建议准备
#   约64 GiB。搜索是单进程的；Anvil shared分区可能为了满足内存请求自动连带分配更多
#   CPU，但这些额外CPU不会加速当前程序。

import os
import numpy as np
from plotly import graph_objects as go
from scipy.spatial import ConvexHull
from py_funcs.read_nekrs_header import read_nekrs_header

tol=0.1
face_match_tol=1e-8

edge_groups=(((0,1),(2,3),(4,5),(6,7)),
             ((0,2),(1,3),(4,6),(5,7)),
             ((0,4),(1,5),(2,6),(3,7)))

def get_edge_lengths(element):
    return np.array([np.linalg.norm(element[a]-element[b])
                     for edge_group in edge_groups for a,b in edge_group])

def element_is_planar(edge_lengths=None):
    element_scale=np.max(edge_lengths)
    if element_scale==0:
        return True
    zero_edge_threshold=tol*element_scale
    for group_id in range(3):
        group_edge_lengths=edge_lengths[group_id*4:(group_id+1)*4]
        if np.all(group_edge_lengths<=zero_edge_threshold):
            return True
    return False

face_corner_ids=((0,1,2,3),(0,1,4,5),(0,2,4,6),
                 (1,3,5,7),(2,3,6,7),(4,5,6,7))

def get_face_keys(element):
    face_keys=[]
    for face_ids in face_corner_ids:
        face_center=np.mean(element[list(face_ids)],axis=0,dtype=np.float64)
        quantized_center=np.rint(face_center/face_match_tol).astype(np.int64)
        face_keys.append(tuple(quantized_center))
    return tuple(face_keys)

def find_three_non_planar_elements(x,corner_pos):
    face_to_elements={}
    adjacency={}
    for element_id in range(len(x)):
        element_corner=np.asarray(x[element_id,corner_pos,:])
        if (element_id+1)%100000==0:
            print(f"已检查{element_id+1}个element")
        if element_is_planar(get_edge_lengths(element_corner)):
            continue
        adjacency[element_id]=set()
        for face_key in get_face_keys(element_corner):
            face_owners=face_to_elements.setdefault(face_key,[])
            for neighbor_id in face_owners:
                common_neighbors=adjacency[element_id] & adjacency[neighbor_id]
                if common_neighbors:
                    return tuple(sorted((element_id,neighbor_id,min(common_neighbors))))
                adjacency[element_id].add(neighbor_id)
                adjacency[neighbor_id].add(element_id)
            face_owners.append(element_id)
    return None

def add_element1(fig,x,element_id,corner_idx,color,with_hull=False):
    element=x[element_id]
    c=element[corner_idx]
    if with_hull:
        hull=ConvexHull(element)
        tri=hull.simplices
        fig.add_trace(go.Mesh3d(x=element[:,0],y=element[:,1],z=element[:,2],i=tri[:,0],j=tri[:,1],k=tri[:,2],opacity=0.1,color=color,showscale=False))
    edges=[(0,1),(1,3),(3,2),(2,0),
           (4,5),(5,7),(7,6),(6,4),
           (0,4),(1,5),(2,6),(3,7)]
    xs,ys,zs=[],[],[]
    for a,b in edges:
        xs+=[c[a,0],c[b,0],None]
        ys+=[c[a,1],c[b,1],None]
        zs+=[c[a,2],c[b,2],None]
    fig.add_trace(go.Scatter3d(x=xs,y=ys,z=zs,mode="lines",line=dict(color=color,width=6),showlegend=False))

data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/viral-flamespeed/jicf0.f00000"
selected_element_ids=(1120084,1120235,1120242)
os.makedirs("fig",exist_ok=True)

with open(data_path,"rb") as f:
    header=read_nekrs_header(f)
    print(f"总共有{header.file_element_num}个element")
    if len(header.c_element_shape)!=3 or len(set(header.c_element_shape))!=1 or min(header.c_element_shape)<4:
        raise ValueError(f"绘图要求每个element的三个维度相同且至少为4，实际为{header.c_element_shape}")
    f.seek(header.file_element_num*4,1)
    coordinate_offset=f.tell()
    coordinate_count=header.file_element_num*header.element_size*3
    file_dtype=np.dtype(header.np_dtype).newbyteorder(header.endian_mode)
coordinate_bytes=coordinate_count*file_dtype.itemsize
if os.path.getsize(data_path)<coordinate_offset+coordinate_bytes:
    raise ValueError("文件中的坐标数据不完整")
x=np.memmap(data_path,dtype=file_dtype,mode="r",offset=coordinate_offset,
            shape=(header.file_element_num,3,header.element_size))
x=x.transpose((0,2,1))  # (element_num,element_size,3)，按需读取所选element

n=header.c_element_shape[0]
n2=n*n
corner_pos=np.array([0,n-1,(n-1)*n,n2-1,(n-1)*n2,(n-1)*n2+n-1,(n-1)*n2+(n-1)*n,n**3-1],dtype=np.int64)
if selected_element_ids is None:
    selected_element_ids=find_three_non_planar_elements(x,corner_pos)
    if selected_element_ids is None:
        raise RuntimeError("没有找到三个两两共享面的非平面element")
    print(f"找到三个两两共享面的非平面element：{selected_element_ids}")
elif max(selected_element_ids)>=header.file_element_num:
    raise IndexError(f"指定的element id超出范围[0,{header.file_element_num})：{selected_element_ids}")

######## 第一张图，显示内部填充和边框 ########
fig=go.Figure()
add_element1(fig,x,selected_element_ids[0],corner_pos,"red",with_hull=True)
add_element1(fig,x,selected_element_ids[1],corner_pos,"green",with_hull=True)
add_element1(fig,x,selected_element_ids[2],corner_pos,"blue",with_hull=True)
fig.update_layout(scene_camera=dict(eye=dict(x=0.6,y=2,z=0.7)),
    scene=dict(aspectmode="auto",xaxis=dict(visible=False),yaxis=dict(visible=False),zaxis=dict(visible=False),bgcolor="rgba(0,0,0,0)"),
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    margin=dict(l=0,r=0,t=0,b=0)
)
fig.write_html("fig/three_elements.html")

######## 第二张图，主要显示一个element靠前的三个面的点 ########
fig=go.Figure()
temp_element=x[selected_element_ids[1]].reshape([n,n,n,3])
layer1=np.concatenate([temp_element[0,:,:].reshape(-1,3),temp_element[1:,0,:].reshape(-1,3),temp_element[1:,1:,0].reshape(-1,3)],axis=0)
layer2=np.concatenate([temp_element[1,1:,1:].reshape(-1,3),temp_element[2:,1,1:].reshape(-1,3),temp_element[2:,2:,1].reshape(-1,3)],axis=0)
layer3=np.concatenate([temp_element[2,2:,2:].reshape(-1,3),temp_element[3:,2,2:].reshape(-1,3),temp_element[3:,3:,2].reshape(-1,3)],axis=0)
fig.add_trace(go.Scatter3d(x=layer1[:,0],y=layer1[:,1],z=layer1[:,2],
                           mode="markers",name=f"element {selected_element_ids[1]}",marker=dict(size=2,color="green"),showlegend=False))
fig.add_trace(go.Scatter3d(x=layer2[:,0],y=layer2[:,1],z=layer2[:,2],opacity=0.2,
                           mode="markers",name=f"element {selected_element_ids[1]}",marker=dict(size=2,color="green"),showlegend=False))
fig.add_trace(go.Scatter3d(x=layer3[:,0],y=layer3[:,1],z=layer3[:,2],opacity=0.04,
                           mode="markers",name=f"element {selected_element_ids[1]}",marker=dict(size=2,color="green"),showlegend=False))
layer1_edge=[]
for dim in range(3):
    layer1_edge.append([])
    layer1_edge[dim]+=np.concatenate([[temp_element[0,0,i,dim],temp_element[0,-1,i,dim],None] for i in range(n)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[0,i,0,dim],temp_element[0,i,-1,dim],None] for i in range(n)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[0,i,0,dim],temp_element[-1,i,0,dim],None] for i in range(n)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[0,0,i,dim],temp_element[-1,0,i,dim],None] for i in range(n)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[i,0,0,dim],temp_element[i,-1,0,dim],None] for i in range(n)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[i,0,0,dim],temp_element[i,0,-1,dim],None] for i in range(n)]).tolist()
fig.add_trace(go.Scatter3d(x=layer1_edge[0],y=layer1_edge[1],z=layer1_edge[2],mode="lines",line=dict(color="green",width=1),showlegend=False))
fig.update_layout(scene_camera=dict(eye=dict(x=0.6,y=2,z=2)),
    scene=dict(aspectmode="data",xaxis=dict(visible=False),yaxis=dict(visible=False),zaxis=dict(visible=False),bgcolor="rgba(0,0,0,0)"),
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    margin=dict(l=0,r=0,t=0,b=0)
)
fig.write_html(f"fig/three_elements2.html")

######## 第三张图，同样的方法画reference element ########
gll=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(n-1).deriv().roots(),[1]))[::-1]
temp_element=np.stack(np.meshgrid(gll,gll,gll,indexing="ij"),axis=-1)
fig=go.Figure()
layer1=np.concatenate([temp_element[0,:,:].reshape(-1,3),temp_element[1:,0,:].reshape(-1,3),temp_element[1:,1:,0].reshape(-1,3)],axis=0)
layer2=np.concatenate([temp_element[1,1:,1:].reshape(-1,3),temp_element[2:,1,1:].reshape(-1,3),temp_element[2:,2:,1].reshape(-1,3)],axis=0)
layer3=np.concatenate([temp_element[2,2:,2:].reshape(-1,3),temp_element[3:,2,2:].reshape(-1,3),temp_element[3:,3:,2].reshape(-1,3)],axis=0)
fig.add_trace(go.Scatter3d(x=layer1[:,0],y=layer1[:,1],z=layer1[:,2],
                           mode="markers",name="reference element",marker=dict(size=2,color="black"),showlegend=False))
fig.add_trace(go.Scatter3d(x=layer2[:,0],y=layer2[:,1],z=layer2[:,2],opacity=0.2,
                           mode="markers",name="reference element",marker=dict(size=2,color="black"),showlegend=False))
fig.add_trace(go.Scatter3d(x=layer3[:,0],y=layer3[:,1],z=layer3[:,2],opacity=0.04,
                           mode="markers",name="reference element",marker=dict(size=2,color="black"),showlegend=False))
layer1_edge=[]
for dim in range(3):
    layer1_edge.append([])
    layer1_edge[dim]+=np.concatenate([[temp_element[0,0,i,dim],temp_element[0,-1,i,dim],None] for i in range(n)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[0,i,0,dim],temp_element[0,i,-1,dim],None] for i in range(n)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[0,i,0,dim],temp_element[-1,i,0,dim],None] for i in range(n)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[0,0,i,dim],temp_element[-1,0,i,dim],None] for i in range(n)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[i,0,0,dim],temp_element[i,-1,0,dim],None] for i in range(n)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[i,0,0,dim],temp_element[i,0,-1,dim],None] for i in range(n)]).tolist()
fig.add_trace(go.Scatter3d(x=layer1_edge[0],y=layer1_edge[1],z=layer1_edge[2],mode="lines",line=dict(color="black",width=1),showlegend=False))
fig.update_layout(scene_camera=dict(eye=dict(x=1.5,y=2,z=1.5)),
    scene=dict(aspectmode="data",xaxis=dict(visible=False),yaxis=dict(visible=False),zaxis=dict(visible=False),bgcolor="rgba(0,0,0,0)"),
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    margin=dict(l=0,r=0,t=0,b=0)
)
fig.write_html("fig/three_elements_ref.html")
