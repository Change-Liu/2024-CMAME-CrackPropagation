import vtk
from vtk.util.numpy_support import vtk_to_numpy
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy as np

# 启用LaTeX渲染和设置全局字体大小
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 24
plt.rcParams['font.family'] = 'Times New Roman'

# 创建一个读取器并加载VTU文件
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("flex2_eyy.vtu")
reader.Update()

# 获取数据和点坐标
output = reader.GetOutput()
points = vtk_to_numpy(output.GetPoints().GetData())
sigma_yy = vtk_to_numpy(output.GetPointData().GetArray('sigma_yy'))
v = vtk_to_numpy(output.GetPointData().GetArray('v'))

# 定义新的规则网格，限定在x:[30,50], y:[30,50]
grid_x, grid_y = np.mgrid[35:45:1000j, 45:55:1000j]

# 插值到新的网格
grid_sigma_yy = griddata(points[:,0:2], sigma_yy, (grid_x, grid_y), method='cubic')
grid_v = griddata(points[:,0:2], v, (grid_x, grid_y), method='cubic')

# 绘制，设置显示范围为x:[30,50], y:[30,50]
plt.figure(figsize=(10, 8))
im = plt.imshow(grid_sigma_yy.T, extent=(35, 45, 45, 55), origin='lower', cmap='rainbow', vmin=-0.02, vmax=0.02)
cp = plt.contour(grid_x, grid_y, grid_v, levels=[0.9], colors='black',linestyles="--")  # 设置等高线级别为0.9

# 设置x和y的坐标轴标记
plt.xticks([ 35, 40, 45])
plt.yticks([ 45, 50, 55])

# 添加颜色条和设置颜色条的刻度及标题
cbar = plt.colorbar(im, ticks=[-0.02, -0.01, 0, 0.01, 0.02])
# cbar.set_label(r'$\varepsilon^{\mathrm{eige}}$', labelpad=20)

plt.tight_layout()
plt.savefig('fig8b.svg')
plt.show()