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
    sigma_xx = vtk_to_numpy(output.GetPointData().GetArray('sigma_xx'))
    sigma_xy = vtk_to_numpy(output.GetPointData().GetArray('sigma_xy'))
    v = vtk_to_numpy(output.GetPointData().GetArray('v'))

    # 定义新的规则网格，限定在x:[30,50], y:[30,50]
    grid_x, grid_y = np.mgrid[35:45:1000j, 45:55:1000j]

    # 插值到新的网格
    grid_sigma_yy = griddata(points[:,0:2], sigma_yy, (grid_x, grid_y), method='cubic')
    grid_sigma_xx = griddata(points[:,0:2], sigma_xx, (grid_x, grid_y), method='cubic')
    grid_sigma_xy = griddata(points[:,0:2], sigma_xy, (grid_x, grid_y), method='cubic')

    x_start = 35
    x_end = 45
    y_start = 45
    y_end = 55
    x_points = y_points = 1000

    # 计算间隔
    dx = (x_end - x_start) / (x_points - 1)
    dy = (y_end - y_start) / (y_points - 1)

    # 定义 x_val 和计算对应的索引
    x_val = 42
    idx = int((x_val - x_start) / dx)

    # 提取对应于 x = 42 的 y 和 v 值
    y_slice = grid_y.T[:, idx]
    sigma_yy_slice = grid_sigma_yy.T[:, idx]
    sigma_xx_slice = grid_sigma_xx.T[:, idx]
    sigma_xy_slice = grid_sigma_xy.T[:, idx]

    return y_slice, sigma_yy_slice, sigma_xx_slice, sigma_xy_slice

# 启用LaTeX渲染和设置全局字体大小
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 24
plt.rcParams['font.family'] = 'Times New Roman'

# 绘制 v 随 y 变化的图
plt.figure(figsize=(10,8))
colors=['#008a9d',"#eeb0b0","#b1c44e","#243258","#8e5629","#631f66"]
markers=['o','s','^','D','v','<','>','p','h','H','d','P','X','*','+','|','_']
labels=[r'$P\downarrow$', r'$P\uparrow$']

files=["flex3-1500.vtu", "flex4-1500.vtu"]

y_slice, sigma_yy_slice, sigma_xx_slice, sigma_xy_slice = read_vtu_file(files[1])

plt.plot(y_slice, sigma_xx_slice, color=colors[0], marker=markers[0], markevery=75, markersize=12, label=r'$\varepsilon_{11}^{\mathrm{eige}}$', lw=3)
plt.plot(y_slice, sigma_yy_slice, color=colors[1], marker=markers[1], markevery=75, markersize=12, label=r'$\varepsilon_{22}^{\mathrm{eige}}$', lw=3)
# plt.plot(y_slice, sigma_xy_slice, color=colors[2], marker=markers[2], markevery=100, label=r'$\varepsilon_{xy}^{\mathrm{eige}}$', lw=3)
    
plt.xlabel(r'$y$ coordinate')
plt.ylabel(r'Eigenstrain $\varepsilon_{ij}^{\mathrm{eige}}$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Fig14d.svg')
plt.show()