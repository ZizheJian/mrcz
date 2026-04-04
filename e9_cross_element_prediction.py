import numpy as np
import random
from py_funcs.read_nekrs_header import read_nekrs_header

class element_class:
    def __init__(self,X,U,P,T,corner_index):
        self.X=X
        self.U=U
        self.P=P
        self.T=T
        self.corner_index=np.array(corner_index).reshape(2,2,2)
    def rotate_a(self):
        self.X=np.rot90(self.X,k=1,axes=(2,1))
        self.U=np.rot90(self.U,k=1,axes=(2,1))
        self.P=np.rot90(self.P,k=1,axes=(2,1))
        self.T=np.rot90(self.T,k=1,axes=(2,1))
        self.corner_index=np.rot90(self.corner_index,k=1,axes=(1,0))
    def rotate_d(self):
        self.X=np.rot90(self.X,k=1,axes=(1,2))
        self.U=np.rot90(self.U,k=1,axes=(1,2))
        self.P=np.rot90(self.P,k=1,axes=(1,2))
        self.T=np.rot90(self.T,k=1,axes=(1,2))
        self.corner_index=np.rot90(self.corner_index,k=1,axes=(0,1))
    def rotate_w(self):
        self.X=np.rot90(self.X,k=1,axes=(3,2))
        self.U=np.rot90(self.U,k=1,axes=(3,2))
        self.P=np.rot90(self.P,k=1,axes=(3,2))
        self.T=np.rot90(self.T,k=1,axes=(3,2))
        self.corner_index=np.rot90(self.corner_index,k=1,axes=(2,1))
    def rotate_s(self):
        self.X=np.rot90(self.X,k=1,axes=(2,3))
        self.U=np.rot90(self.U,k=1,axes=(2,3))
        self.P=np.rot90(self.P,k=1,axes=(2,3))
        self.T=np.rot90(self.T,k=1,axes=(2,3))
        self.corner_index=np.rot90(self.corner_index,k=1,axes=(1,2))
    def rotate_q(self):
        self.X=np.rot90(self.X,k=1,axes=(1,3))
        self.U=np.rot90(self.U,k=1,axes=(1,3))
        self.P=np.rot90(self.P,k=1,axes=(1,3))
        self.T=np.rot90(self.T,k=1,axes=(1,3))
        self.corner_index=np.rot90(self.corner_index,k=1,axes=(0,2))
    def rotate_e(self):
        self.X=np.rot90(self.X,k=1,axes=(3,1))
        self.U=np.rot90(self.U,k=1,axes=(3,1))
        self.P=np.rot90(self.P,k=1,axes=(3,1))
        self.T=np.rot90(self.T,k=1,axes=(3,1))
        self.corner_index=np.rot90(self.corner_index,k=1,axes=(2,0))

# data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/fld_pyutils/demos/data/test0.f00001"
data_path="/anvil/projects/x-cis240192/x-zjian1-in-cis240192/unstructured-data/hyuan/pb0.f00001"
f=open(data_path,"rb")
header=read_nekrs_header(f)
element_map=np.fromfile(f,dtype=np.int32,count=header.file_element_num)
X=np.fromfile(f,dtype=np.float32,count=header.file_element_num*header.element_size*3)
X=X.reshape(header.file_element_num,3,header.element_shape[0],header.element_shape[1],header.element_shape[2]).transpose(1,0,2,3,4)
U=np.fromfile(f,dtype=np.float32,count=header.file_element_num*header.element_size*3)
U=U.reshape(header.file_element_num,3,header.element_shape[0],header.element_shape[1],header.element_shape[2]).transpose(1,0,2,3,4)
P=np.fromfile(f,dtype=np.float32,count=header.file_element_num*header.element_size)
P=P.reshape(1,header.file_element_num,header.element_shape[0],header.element_shape[1],header.element_shape[2])
T=np.fromfile(f,dtype=np.float32,count=header.file_element_num*header.element_size)
T=T.reshape(1,header.file_element_num,header.element_shape[0],header.element_shape[1],header.element_shape[2])
f.close()

assert(header.element_shape[0]==header.element_shape[1]==header.element_shape[2])
n=header.element_shape[0]
x_legendre=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(n-1).deriv().roots(),[1]))
gap_legendre=np.diff(x_legendre)
ei_list=random.sample(range(header.file_element_num),5)
for ei in ei_list:
    ei_corners=[X[:,ei,0,0,0],X[:,ei,n-1,0,0],X[:,ei,0,n-1,0],X[:,ei,n-1,n-1,0],X[:,ei,0,0,n-1],X[:,ei,n-1,0,n-1],X[:,ei,0,n-1,n-1],X[:,ei,n-1,n-1,n-1]]
    ei_center=np.mean(ei_corners,axis=0)
    dist=np.zeros(header.file_element_num)
    min_center_dist=np.inf
    for ej in range(header.file_element_num):
        if ej%10000==0:
            print("ej=",ej)
        if ej==ei:
            dist[ej]=np.inf
            continue
        ej_corners=[X[:,ej,0,0,0],X[:,ej,n-1,0,0],X[:,ej,0,n-1,0],X[:,ej,n-1,n-1,0],X[:,ej,0,0,n-1],X[:,ej,n-1,0,n-1],X[:,ej,0,n-1,n-1],X[:,ej,n-1,n-1,n-1]]
        ej_center=np.mean(ej_corners,axis=0)
        if np.linalg.norm(ei_center-ej_center)<min_center_dist:
            min_center_dist=np.linalg.norm(ei_center-ej_center)
        elif np.linalg.norm(ei_center-ej_center)>10*min_center_dist:
            dist[ej]=np.inf
            continue
        temp_dist=[np.linalg.norm(ei_corner-ej_corner) for ei_corner in ei_corners for ej_corner in ej_corners]
        dist[ej]=np.partition(temp_dist,3)[3]
    nearest_e=list(np.argpartition(dist,6)[:6])
    nearest_e=[ej for ej in nearest_e if dist[ej]<min_center_dist*1e-3]
    print("ei=",ei)
    print(ei_corners)
    for ej in nearest_e:
        print("ej=",ej)
        ej_coeners=[X[:,ej,0,0,0],X[:,ej,n-1,0,0],X[:,ej,0,n-1,0],X[:,ej,n-1,n-1,0],X[:,ej,0,0,n-1],X[:,ej,n-1,0,n-1],X[:,ej,0,n-1,n-1],X[:,ej,n-1,n-1,n-1]]
        print(ej_coeners)
        print("dist=",dist[ej])
    for ej in nearest_e:
        ej_coeners=[X[:,ej,0,0,0],X[:,ej,n-1,0,0],X[:,ej,0,n-1,0],X[:,ej,n-1,n-1,0],X[:,ej,0,0,n-1],X[:,ej,n-1,0,n-1],X[:,ej,0,n-1,n-1],X[:,ej,n-1,n-1,n-1]]
        pairs=[]
        for ii,cei in enumerate(ei_corners):
            for jj,cej in enumerate(ej_coeners):
                d=np.linalg.norm(cei-cej)
                pairs.append((d,ii,jj))
        print(pairs)
        pairs.sort(key=lambda x:x[0])
        nearest4pairs=pairs[:4]
        element_i=element_class(X[:,ei],U[:,ei],P[:,ei],T[:,ei],[0,4,2,6,1,5,3,7])
        while True:
            print("corner on shared face: ",[i for _,i,_ in nearest4pairs])
            print(f"  {element_i.corner_index[0,1,1]}-{element_i.corner_index[1,1,1]}")
            print(f" / /|")
            print(f"{element_i.corner_index[0,0,1]}-{element_i.corner_index[1,0,1]} {element_i.corner_index[1,1,0]}")
            print(f"| |/")
            print(f"{element_i.corner_index[0,0,0]}-{element_i.corner_index[1,0,0]}")
            ei_face=[element_i.X[:,n-1,0,0],element_i.X[:,n-1,n-1,0],element_i.X[:,n-1,0,n-1],element_i.X[:,n-1,n-1,n-1]]
            ej_corner=[X[:,ej,0,0,0],X[:,ej,n-1,0,0],X[:,ej,0,n-1,0],X[:,ej,n-1,n-1,0],X[:,ej,0,0,n-1],X[:,ej,n-1,0,n-1],X[:,ej,0,n-1,n-1],X[:,ej,n-1,n-1,n-1]]
            print("ei_face=",ei_face)
            print("ej_corner=",ej_corner)
            command=input()
            if command=="a":
                element_i.rotate_a()
            elif command=="d":
                element_i.rotate_d()
            elif command=="w":
                element_i.rotate_w()
            elif command=="s":
                element_i.rotate_s()
            elif command=="q":
                element_i.rotate_q()
            elif command=="e":
                element_i.rotate_e()
            elif command=="break":
                break
        element_j=element_class(X[:,ej],U[:,ej],P[:,ej],T[:,ej],[0,4,2,6,1,5,3,7])
        while True:
            print("corner on shared face: ",[j for _,_,j in nearest4pairs])
            print(f"  {element_j.corner_index[0,1,1]}-{element_j.corner_index[1,1,1]}")
            print(f" / /|")
            print(f"{element_j.corner_index[0,0,1]}-{element_j.corner_index[1,0,1]} {element_j.corner_index[1,1,0]}")
            print(f"| |/")
            print(f"{element_j.corner_index[0,0,0]}-{element_j.corner_index[1,0,0]}")
            ei_face=[element_i.X[:,n-1,0,0],element_i.X[:,n-1,n-1,0],element_i.X[:,n-1,0,n-1],element_i.X[:,n-1,n-1,n-1]]
            ej_face=[element_j.X[:,0,0,0],element_j.X[:,0,n-1,0],element_j.X[:,0,0,n-1],element_j.X[:,0,n-1,n-1]]
            sum_dist=np.sum([np.linalg.norm(ei_face[i]-ej_face[i]) for i in range(4)])
            print("ei_face=",ei_face)
            print("ej_face=",ej_face)
            print("dist=",sum_dist)
            command=input()
            if command=="a":
                element_j.rotate_a()
            elif command=="d":
                element_j.rotate_d()
            elif command=="w":
                element_j.rotate_w()
            elif command=="s":
                element_j.rotate_s()
            elif command=="q":
                element_j.rotate_q()
            elif command=="e":
                element_j.rotate_e()
            elif command=="break":
                break
        cross_element_err=((0.5*(element_i.U[0,-2]+element_j.U[0,1])-0.5*(element_i.U[0,-1]+element_j.U[0,0]))**2).mean()
        in_element_err=(((element_i.U[0,-3]*gap_legendre[0]+element_i.U[0,-1]*gap_legendre[1])/(gap_legendre[0]+gap_legendre[1])-element_i.U[0,-2])**2).mean()
        in_element_err+=(((element_j.U[0,2]*gap_legendre[0]+element_j.U[0,0]*gap_legendre[1])/(gap_legendre[0]+gap_legendre[1])-element_j.U[0,1])**2).mean()
        in_element_err/=2
        print("cross_element_err=",cross_element_err)
        print("in_element_err=",in_element_err)
                
        
