import numpy as np
from matplotlib import pyplot as plt

x=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(10-1).deriv().roots(),[1]))
X,Y=np.meshgrid(x, x)
plt.figure(figsize=(6, 6))
plt.scatter(X,Y,s=50,c="blue",marker="o")  # 或用 plt.plot(X, Y, 'o-') 连线
plt.axis("equal")
plt.savefig("fig/2D_GLL_points.png",dpi=300)