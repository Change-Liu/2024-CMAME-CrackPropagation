import vtk
from vtk.util.numpy_support import vtk_to_numpy
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy as np

def read_vtu_file(file_name):
    # 创建一个读取器并加载VTU文件
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()

    # 获取数据和点坐标
    output = reader.GetOutput()
    points = vtk_to_numpy(output.GetPoints().GetData())
    sigma_yy = vtk_to_numpy(output.GetPointData().GetArray('sigma_yy'))

    # 定义新的规则网格，限定在x:[30,50], y:[30,50]
    grid_x, grid_y = np.mgrid[41:90:1000j, 45:55:1000j]

    # 插值到新的网格
    grid_sigma_yy = griddata(points[:,0:2], sigma_yy, (grid_x, grid_y), method='cubic')

    y_start = 45
    y_end = 55
    y_points = 1000

    # 计算间隔
    dy = (y_end - y_start) / (y_points - 1)

    # 定义 x_val 和计算对应的索引
    y_val = 50
    idy = int((y_val - y_start) / dy)

    # 提取对应于 x = 42 的 y 和 v 值
    x_slice = grid_x.T[idy, :]
    sigma_yy_slice = grid_sigma_yy.T[idy, :]

    return x_slice, sigma_yy_slice

# 启用LaTeX渲染和设置全局字体大小
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 30
plt.rcParams['font.family'] = 'Times New Roman'

# 绘制 v 随 y 变化的图
plt.figure(figsize=(10,8))
colors=['#008a9d',"#eeb0b0","#b1c44e","#243258","#8e5629","#631f66"]
markers=['o','s','^','D','v','<','>','p','h','H','d','P','X','*','+','|','_']
labels=[r'$P\rightarrow$', r'$P\leftarrow$']

files=["flex1_eyy.vtu", "flex2_eyy.vtu"]
cidx=0
for i in files:
    x_slice, sigma_yy_slice = read_vtu_file(i)

    plt.plot(x_slice-40, sigma_yy_slice, color=colors[cidx], marker=markers[cidx], markevery=100,  markersize=16, label=labels[cidx], lw=3)

    cidx = cidx + 1
    
plt.yticks([-0.02, -0.01, 0, 0.01, 0.02])
plt.xlabel(r'Distance from crack tip $x$')
plt.ylabel(r'Eigenstrain $\varepsilon_{22}^{\mathrm{eige}}$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('fig7b.svg')
plt.show()