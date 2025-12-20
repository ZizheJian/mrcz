import numpy as np
import scipy

n=10

def max_non_diag(w,L):
    ret=np.einsum("mi,i,ni->mn",L,w,L)
    return np.max(np.abs(ret-np.diag(np.diag(ret))))

######## GLL + DLT ########
x=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(n-1).deriv().roots(),[1]))
L=np.zeros((n,n),dtype=np.float64)
L[0]=1
L[1]=x
for i in range(2,n):
    L[i]=((2*i-1)*x*L[i-1]-(i-1)*L[i-2])/i
w=2/(n*(n-1)*L[-1]**2)
print("max non-diag=",max_non_diag(w,L))
print("d=",np.einsum("ia,a,ia->i",L,w,L))

######## Uniform Grid + DLT ########
def max_non_diag(w,L):
    ret=np.einsum("mi,i,ni->mn",L,w,L)
    return np.max(np.abs(ret-np.diag(np.diag(ret))))

x=np.arange(0,n)/(n-1)*2-1
L=np.zeros((n,n),dtype=np.float64)
L[0]=1
L[1]=x
for i in range(2,n):
    L[i]=((2*i-1)*x*L[i-1]-(i-1)*L[i-2])/i
cons=({"type":"eq","fun":lambda w:np.sum(w)-2})
bounds=[(0,None)]*n
w0=np.full(n,2/n)
res=scipy.optimize.minimize(lambda w:max_non_diag(w,L),w0,method="SLSQP",bounds=bounds,constraints=cons)
print(res.x)
print("max non-diag=",max_non_diag(res.x,L))
print("d=",np.einsum("ia,a,ia->i",L,res.x,L))

######## GLL + DChT ########
x=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(n-1).deriv().roots(),[1]))
L=np.zeros((n,n),dtype=np.float64)
L[0]=1
L[1]=x
for i in range(2,n):
    L[i]=2*x*L[i-1]-L[i-2]
cons=({"type":"eq","fun":lambda w:np.sum(w)-2})
bounds=[(0,None)]*n
w0=np.full(n,2/n)
res=scipy.optimize.minimize(lambda w:max_non_diag(w,L),w0,method="SLSQP",bounds=bounds,constraints=cons)
print(res.x)
print("max non-diag=",max_non_diag(res.x,L))
print("d=",np.einsum("ia,a,ia->i",L,res.x,L))
