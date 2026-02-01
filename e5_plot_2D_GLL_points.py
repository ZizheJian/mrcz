import numpy as np
from matplotlib import pyplot as plt

x=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(10-1).deriv().roots(),[1]))
X,Y=np.meshgrid(x, x)

fig, ax = plt.subplots(figsize=(6, 6))

ax.scatter(X, Y, s=50, c="blue", marker="o")

# 1️⃣ 坐标比例一致
ax.set_aspect("equal")

# 2️⃣ 去掉外边框（top / right）
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# 4️⃣ 可选：轴线加粗一点，更清楚
# ax.spines["bottom"].set_linewidth(1.5)
# ax.spines["left"].set_linewidth(1.5)

# 5️⃣ 增大刻度数字字号
ax.tick_params(axis="both", which="major", labelsize=20)

ax.set_xticks(np.linspace(-1, 1, 5))
ax.set_yticks(np.linspace(-1, 1, 5))

# 6️⃣ 设置坐标轴名称和字号
ax.set_xlabel("x", fontsize=20)
ax.set_ylabel("y", fontsize=20)
plt.savefig("fig/2D_GLL_points.png", dpi=300, bbox_inches="tight")