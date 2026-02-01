# TODO: 3阶样条插值、改变预测顺序
# TODO: 对element按照连续性进行分类
# TODO: tune各种参数
# TODO: 所有huffman共用一个字典
# TODO: 根据error bound选择数据类型

import numpy as np

n=6
x_legendre=np.concatenate(([-1],np.polynomial.legendre.Legendre.basis(n-1).deriv().roots(),[1]))
x_chebyshev=np.sort(np.cos((np.arange(n))*np.pi/(n-1)))
print("Legendre nodes:",x_legendre)
print("Chebyshev nodes:",x_chebyshev)
print("Legendre gap:",np.diff(x_legendre))
print("Chebyshev gap:",np.diff(x_chebyshev))