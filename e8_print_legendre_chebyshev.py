import numpy as np

n=10
x_legendre=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(n-1).deriv().roots(),[1]))
x_chebyshev=np.sort(np.cos((np.arange(n))*np.pi/(n-1)))
print("Legendre nodes:",x_legendre)
print("Chebyshev nodes:",x_chebyshev)
print("Legendre gap:",np.diff(x_legendre))
print("Chebyshev gap:",np.diff(x_chebyshev))

# def normalize_columns(V):
#  norms=np.linalg.norm(V,axis=0)
#  return V/norms,norms

# def stability_report(name,V):
#     condV=np.linalg.cond(V)
#     G=V.T@V
#     condG=np.linalg.cond(G)
#     Vn,_=normalize_columns(V)
#     Gn=Vn.T@Vn
#     orth_err_2=np.linalg.norm(Gn-np.eye(Gn.shape[0]),2)
#     orth_err_f=np.linalg.norm(Gn-np.eye(Gn.shape[0]),'fro')
#     print(f"\n{name}")
#     print("shape:",V.shape)
#     print("cond(V):",condV)
#     print("cond(V^T V):",condG)
#     print("||G-I||_2 after column normalization:",orth_err_2)
#     print("||G-I||_F after column normalization:",orth_err_f)


# V_leg_on_gll=np.zeros((n,n),dtype=np.float64)
# for i in range(n):
#     if i==0:
#         V_leg_on_gll[:,i]=1.0
#     elif i==1:
#         V_leg_on_gll[:,i]=x_legendre
#     else:
#         V_leg_on_gll[:,i]=((2*i-1)*x_legendre*V_leg_on_gll[:,i-1]-(i-1)*V_leg_on_gll[:,i-2])/i

# V_cheb_on_gll=np.zeros((n,n),dtype=np.float64)
# for i in range(n):
#     if i==0:
#         V_cheb_on_gll[:,i]=1.0
#     elif i==1:
#         V_cheb_on_gll[:,i]=x_legendre
#     else:
#         V_cheb_on_gll[:,i]=2*x_legendre*V_cheb_on_gll[:,i-1]-V_cheb_on_gll[:,i-2]

# stability_report("Discrete Legendre matrix on GLL",V_leg_on_gll)
# stability_report("Discrete Chebyshev matrix on GLL",V_cheb_on_gll)