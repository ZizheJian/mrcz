import numpy as np

n=10
x_legendre=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(n-1).deriv().roots(),[1]))
x_chebyshev=np.sort(np.cos((np.arange(n))*np.pi/(n-1)))
x_uniform=np.linspace(-1,1,n)
print("Legendre nodes:",x_legendre)
print("Chebyshev nodes:",x_chebyshev)
print("Legendre gap:",np.diff(x_legendre))
print("Chebyshev gap:",np.diff(x_chebyshev))

dlt_gll=np.zeros((n,n),dtype=np.float64)
for i in range(n):
    if i==0:
        dlt_gll[:,i]=1.0
    elif i==1:
        dlt_gll[:,i]=x_legendre
    else:
        dlt_gll[:,i]=((2*i-1)*x_legendre*dlt_gll[:,i-1]-(i-1)*dlt_gll[:,i-2])/i
print(f"dlt_gll_vmax={np.max(np.abs(dlt_gll))}, dlt_gll_norm2={np.linalg.norm(dlt_gll,2)}")

dcht_gll=np.zeros((n,n),dtype=np.float64)
for i in range(n):
    if i==0:
        dcht_gll[:,i]=1.0
    elif i==1:
        dcht_gll[:,i]=x_legendre
    else:
        dcht_gll[:,i]=2*x_legendre*dcht_gll[:,i-1]-dcht_gll[:,i-2]
print(f"dcht_gll_vmax={np.max(np.abs(dcht_gll))}, dcht_gll_norm2={np.linalg.norm(dcht_gll,2)}")

dct_gll=np.zeros((n,n),dtype=np.float64)
for i in range(n):
    for j in range(n):
        dct_gll[j,i]=np.cos(i*np.pi*(x_legendre[j]+1)/2)
print(f"dct_gll_vmax={np.max(np.abs(dct_gll))}, dct_gll_norm2={np.linalg.norm(dct_gll,2)}")

dlt_uniform=np.zeros((n,n),dtype=np.float64)
for i in range(n):
    if i==0:
        dlt_uniform[:,i]=1.0
    elif i==1:
        dlt_uniform[:,i]=x_uniform
    else:
        dlt_uniform[:,i]=((2*i-1)*x_uniform*dlt_uniform[:,i-1]-(i-1)*dlt_uniform[:,i-2])/i
print(f"dlt_uniform_vmax={np.max(np.abs(dlt_uniform))}, dlt_uniform_norm2={np.linalg.norm(dlt_uniform,2)}")

dcht_uniform=np.zeros((n,n),dtype=np.float64)
for i in range(n):
    if i==0:
        dcht_uniform[:,i]=1.0
    elif i==1:
        dcht_uniform[:,i]=x_uniform
    else:
        dcht_uniform[:,i]=2*x_uniform*dcht_uniform[:,i-1]-dcht_uniform[:,i-2]
print(f"dcht_uniform_vmax={np.max(np.abs(dcht_uniform))}, dcht_uniform_norm2={np.linalg.norm(dcht_uniform,2)}")

dst_uniform=np.zeros((n,n),dtype=np.float64)
for i in range(n):
    for j in range(n):
        dst_uniform[j,i]=np.sin((i+1)*np.pi*(j+1)/(n+1))
print(f"dst_uniform_vmax={np.max(np.abs(dst_uniform))}, dst_uniform_norm2={np.linalg.norm(dst_uniform,2)}")

dct_uniform=np.zeros((n,n),dtype=np.float64)
for i in range(n):
    for j in range(n):
        dct_uniform[j,i]=np.cos(i*np.pi*(j+1)/(n+1))
print(f"dct_uniform_vmax={np.max(np.abs(dct_uniform))}, dct_uniform_norm2={np.linalg.norm(dct_uniform,2)}")

