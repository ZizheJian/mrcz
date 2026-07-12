import numpy as np
from plotly import graph_objects as go
from scipy.spatial import ConvexHull
from py_funcs.read_nekrs_header import read_nekrs_header

tol=1e-8

def same_point(p1,p2):
    print(np.linalg.norm(p1-p2))
    return np.linalg.norm(p1-p2)<tol

def element_is_flat(element):
    if same_point(element[0],element[1]) and same_point(element[2],element[3]) and same_point(element[4],element[5]) and same_point(element[6],element[7]):
        return True
    if same_point(element[0],element[2]) and same_point(element[1],element[3]) and same_point(element[4],element[6]) and same_point(element[5],element[7]):
        return True
    if same_point(element[0],element[4]) and same_point(element[1],element[5]) and same_point(element[2],element[6]) and same_point(element[3],element[7]):
        return True
    return False

def shortest_and_longest_edge(element):
    max_ratio=0
    edges=[]
    for i in range(4):
        edges.append(np.linalg.norm(element[2*i]-element[2*i+1]))
    max_ratio=max(max_ratio,(max(edges)+tol)/(min(edges)+tol))
    edges=[]
    for i in range(2):
        edges.append(np.linalg.norm(element[4*i]-element[4*i+2]))
        edges.append(np.linalg.norm(element[4*i+1]-element[4*i+3]))
    max_ratio=max(max_ratio,(max(edges)+tol)/(min(edges)+tol))
    edges=[]
    for i in range(4):
        edges.append(np.linalg.norm(element[i]-element[i+4]))
    max_ratio=max(max_ratio,(max(edges)+tol)/(min(edges)+tol))
    return max_ratio

def find_and_write_three_elements(header,corner,x):
    surface_int=np.zeros((header.file_element_num,6,3),dtype=np.int64)
    surface_int[:,0,:]=np.round(0.25*(corner[:,0,:]+corner[:,1,:]+corner[:,2,:]+corner[:,3,:])/tol)
    surface_int[:,1,:]=np.round(0.25*(corner[:,0,:]+corner[:,1,:]+corner[:,4,:]+corner[:,5,:])/tol)
    surface_int[:,2,:]=np.round(0.25*(corner[:,0,:]+corner[:,2,:]+corner[:,4,:]+corner[:,6,:])/tol)
    surface_int[:,3,:]=np.round(0.25*(corner[:,1,:]+corner[:,3,:]+corner[:,5,:]+corner[:,7,:])/tol)
    surface_int[:,4,:]=np.round(0.25*(corner[:,2,:]+corner[:,3,:]+corner[:,6,:]+corner[:,7,:])/tol)
    surface_int[:,5,:]=np.round(0.25*(corner[:,4,:]+corner[:,5,:]+corner[:,6,:]+corner[:,7,:])/tol)

    surface2element={}
    for i in range(header.file_element_num):
        for j in range(6):
            key=tuple(surface_int[i,j,:])
            if key not in surface2element:
                surface2element[key]=[]
            surface2element[key].append(i)
    surface2element={k:v for k,v in surface2element.items() if len(v)==2}

    for surface1 in surface2element:
        elements1=surface2element[surface1][0:1] # 其中一个element
        surfaces2=set()
        for element1 in elements1:
            for surface_id in range(6):
                surface2=tuple(surface_int[element1,surface_id,:])
                if (surface1!=surface2) and (surface2 in surface2element):
                    surfaces2.add(surface2)
        elements2=set()
        for surface2 in surfaces2:
            for element2 in surface2element[surface2]:
                if element2 not in elements1:
                    elements2.add(element2)
        surfaces3=set()
        for element2 in elements2:
            for surface_id in range(6):
                surface3=tuple(surface_int[element2,surface_id,:])
                if surface3 not in surfaces2 and surface3 in surface2element:
                    surfaces3.add(surface3)
        elements3=set()
        for surface3 in surfaces3:
            for element3 in surface2element[surface3]:
                if element3 not in elements2:
                    elements3.add(element3)
        if (surface2element[surface1][1] in elements3):
            element_a=surface2element[surface1][0]
            element_b=surface2element[surface1][1]
            surfaces_a=set()
            for surface_id in range(6):
                surface_a=tuple(surface_int[element_a,surface_id,:])
                if surface_a!=surface1 and surface_a in surface2element:
                    surfaces_a.add(surface_a)
            surfaces_b=set()
            for surface_id in range(6):
                surface_b=tuple(surface_int[element_b,surface_id,:])
                if surface_b!=surface1 and surface_b in surface2element:
                    surfaces_b.add(surface_b)
            elements_a2=set()
            for surface_a in surfaces_a:
                for element_a2 in surface2element[surface_a]:
                    if element_a2!=element_a:
                        elements_a2.add(element_a2)
            elements_b2=set()
            for surface_b in surfaces_b:
                for element_b2 in surface2element[surface_b]:
                    if element_b2!=element_b:
                        elements_b2.add(element_b2)
            element_c=[element for element in elements_a2 if element in elements_b2][0]
            print(f"找到环形相邻结构: {element_a}, {element_b}, {element_c}")
            fig=go.Figure()
            fig.add_trace(go.Scatter3d(x=x[element_a,:,0],y=x[element_a,:,1],z=x[element_a,:,2],mode="markers",
                                        name=f"element {element_a}",marker=dict(size=2,color="red")))
            fig.add_trace(go.Scatter3d(x=x[element_b,:,0],y=x[element_b,:,1],z=x[element_b,:,2],mode="markers",
                                        name=f"element {element_b}",marker=dict(size=2,color="green")))
            fig.add_trace(go.Scatter3d(x=x[element_c,:,0],y=x[element_c,:,1],z=x[element_c,:,2],mode="markers",
                                        name=f"element {element_c}",marker=dict(size=2,color="blue")))
            fig.update_layout(scene=dict(aspectmode="data",xaxis_title="X",yaxis_title="Y",zaxis_title="Z"))
            fig.write_html(f"fig/irregular_elements_{element_a}_{element_b}_{element_c}.html")
            input()

def add_element1(fig,x,element_id,color,with_hull=False):
    element=x[element_id]
    corner_idx=[0,9,90,99,900,909,990,999]
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

f=open(data_path,"rb")
header=read_nekrs_header(f)
print(header.file_element_num)
f.seek(header.file_element_num*4,1)

x=np.fromfile(f,dtype=np.float32,count=header.file_element_num*header.element_size*3)
x=x.reshape((header.file_element_num,3,header.element_size))
x=x.transpose((0,2,1))  # (element_num,element_size,3)

corner_pos=np.array([0,9,90,99,900,909,990,999],dtype=np.int64)
corner=x[:,corner_pos,:]

# not_flat_elements=[i for i in range(header.file_element_num) if not element_is_flat(corner[i])]
# print(f"共有 {len(not_flat_elements)} 个非平面elements")

# best_element=0
# worst_element=0
# max_ratio=0
# min_ratio=1/tol
# for i in not_flat_elements:
#     ratio=shortest_and_longest_edge(corner[i])
#     if ratio<min_ratio:
#         min_ratio=ratio
#         best_element=i
#     if ratio>max_ratio:
#         max_ratio=ratio
#         worst_element=i
# print(f"最优元素编号: {best_element}, 最小边长比: {min_ratio}")
# print(f"最劣元素编号: {worst_element}, 最大边长比: {max_ratio}")

######## 第一张图，显示内部填充和边框 ########
fig=go.Figure()
add_element1(fig,x,1120084,"red",with_hull=True)
add_element1(fig,x,1120235,"green",with_hull=True)
add_element1(fig,x,1120242,"blue",with_hull=True)
fig.update_layout(scene_camera=dict(eye=dict(x=0.6,y=2,z=0.7)),
    scene=dict(xaxis=dict(visible=False),yaxis=dict(visible=False),zaxis=dict(visible=False),bgcolor="rgba(0,0,0,0)"),
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    margin=dict(l=0,r=0,t=0,b=0)
)
fig.write_html(f"fig/three_elements1.html")

######## 第二张图，主要显示一个element靠前的三个面的点 ########
fig=go.Figure()
temp_element=x[1120235].reshape([10,10,10,3])
layer1=np.concatenate([temp_element[0,:,:].reshape(-1,3),temp_element[1:,0,:].reshape(-1,3),temp_element[1:,1:,0].reshape(-1,3)],axis=0)
layer2=np.concatenate([temp_element[1,1:,1:].reshape(-1,3),temp_element[2:,1,1:].reshape(-1,3),temp_element[2:,2:,1].reshape(-1,3)],axis=0)
layer3=np.concatenate([temp_element[2,2:,2:].reshape(-1,3),temp_element[3:,2,2:].reshape(-1,3),temp_element[3:,3:,2].reshape(-1,3)],axis=0)
fig.add_trace(go.Scatter3d(x=layer1[:,0],y=layer1[:,1],z=layer1[:,2],
                           mode="markers",name=f"element {1120235}",marker=dict(size=2,color="green"),showlegend=False))
fig.add_trace(go.Scatter3d(x=layer2[:,0],y=layer2[:,1],z=layer2[:,2],opacity=0.2,
                           mode="markers",name=f"element {1120235}",marker=dict(size=2,color="green"),showlegend=False))
fig.add_trace(go.Scatter3d(x=layer3[:,0],y=layer3[:,1],z=layer3[:,2],opacity=0.04,
                           mode="markers",name=f"element {1120235}",marker=dict(size=2,color="green"),showlegend=False))
layer1_edge=[]
for dim in range(3):
    layer1_edge.append([])
    layer1_edge[dim]+=np.concatenate([[temp_element[0,0,i,dim],temp_element[0,-1,i,dim],None] for i in range(0,10)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[0,i,0,dim],temp_element[0,i,-1,dim],None] for i in range(0,10)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[0,i,0,dim],temp_element[-1,i,0,dim],None] for i in range(0,10)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[0,0,i,dim],temp_element[-1,0,i,dim],None] for i in range(0,10)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[i,0,0,dim],temp_element[i,-1,0,dim],None] for i in range(0,10)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[i,0,0,dim],temp_element[i,0,-1,dim],None] for i in range(0,10)]).tolist()
fig.add_trace(go.Scatter3d(x=layer1_edge[0],y=layer1_edge[1],z=layer1_edge[2],mode="lines",line=dict(color="green",width=1),showlegend=False))
fig.update_layout(scene_camera=dict(eye=dict(x=0.6,y=2,z=2)),
    scene=dict(xaxis=dict(visible=False),yaxis=dict(visible=False),zaxis=dict(visible=False),bgcolor="rgba(0,0,0,0)"),
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    margin=dict(l=0,r=0,t=0,b=0)
)
fig.write_html(f"fig/three_elements2.html")

######## 第三张图，同样的方法画reference element ########
gll=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(10-1).deriv().roots(),[1]))[::-1]
temp_element=np.stack(np.meshgrid(gll,gll,gll,indexing="ij"),axis=-1)
fig=go.Figure()
layer1=np.concatenate([temp_element[0,:,:].reshape(-1,3),temp_element[1:,0,:].reshape(-1,3),temp_element[1:,1:,0].reshape(-1,3)],axis=0)
layer2=np.concatenate([temp_element[1,1:,1:].reshape(-1,3),temp_element[2:,1,1:].reshape(-1,3),temp_element[2:,2:,1].reshape(-1,3)],axis=0)
layer3=np.concatenate([temp_element[2,2:,2:].reshape(-1,3),temp_element[3:,2,2:].reshape(-1,3),temp_element[3:,3:,2].reshape(-1,3)],axis=0)
fig.add_trace(go.Scatter3d(x=layer1[:,0],y=layer1[:,1],z=layer1[:,2],
                           mode="markers",name=f"element {1120235}",marker=dict(size=2,color="green"),showlegend=False))
fig.add_trace(go.Scatter3d(x=layer2[:,0],y=layer2[:,1],z=layer2[:,2],opacity=0.2,
                           mode="markers",name=f"element {1120235}",marker=dict(size=2,color="green"),showlegend=False))
fig.add_trace(go.Scatter3d(x=layer3[:,0],y=layer3[:,1],z=layer3[:,2],opacity=0.04,
                           mode="markers",name=f"element {1120235}",marker=dict(size=2,color="green"),showlegend=False))
layer1_edge=[]
for dim in range(3):
    layer1_edge.append([])
    layer1_edge[dim]+=np.concatenate([[temp_element[0,0,i,dim],temp_element[0,-1,i,dim],None] for i in range(0,10)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[0,i,0,dim],temp_element[0,i,-1,dim],None] for i in range(0,10)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[0,i,0,dim],temp_element[-1,i,0,dim],None] for i in range(0,10)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[0,0,i,dim],temp_element[-1,0,i,dim],None] for i in range(0,10)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[i,0,0,dim],temp_element[i,-1,0,dim],None] for i in range(0,10)]).tolist()
    layer1_edge[dim]+=np.concatenate([[temp_element[i,0,0,dim],temp_element[i,0,-1,dim],None] for i in range(0,10)]).tolist()
fig.add_trace(go.Scatter3d(x=layer1_edge[0],y=layer1_edge[1],z=layer1_edge[2],mode="lines",line=dict(color="green",width=1),showlegend=False))
fig.update_layout(scene_camera=dict(eye=dict(x=1.5,y=2,z=1.5)),
    scene=dict(xaxis=dict(visible=False),yaxis=dict(visible=False),zaxis=dict(visible=False),bgcolor="rgba(0,0,0,0)"),
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    margin=dict(l=0,r=0,t=0,b=0)
)
fig.write_html(f"fig/three_elements3.html")
