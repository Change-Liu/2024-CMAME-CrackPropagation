import vtk
from vtk.util.numpy_support import vtk_to_numpy
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import colorcet as cc
import numpy as np

def plotfig1(datafile):
    # 启用LaTeX渲染和设置全局字体大小
    # plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size'] = 24
    plt.rcParams['font.family'] = 'Times New Roman'

    # 创建一个读取器并加载VTU文件
    base_path = 'T{}.vtu'
    paths = base_path.format(datafile)
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(paths)
    reader.Update()

    # 获取数据和点坐标
    output = reader.GetOutput()
    points = vtk_to_numpy(output.GetPoints().GetData())
    v_array = vtk_to_numpy(output.GetPointData().GetArray('v'))
    P_array = vtk_to_numpy(output.GetPointData().GetArray('P'))

    # 定义新的规则网格，限定在x:[30,50], y:[30,50]
    grid_x, grid_y = np.mgrid[0:100:1000j, 0:100:1000j]

    # 插值到新的网格
    grid_v = griddata(points[:,0:2], v_array, (grid_x, grid_y), method='cubic')
    grid_Px = griddata(points[:,0:2], P_array[:,0], (grid_x, grid_y), method='cubic')  # X 分量
    grid_Py = griddata(points[:,0:2], P_array[:,1], (grid_x, grid_y), method='cubic')  # Y 分量

    # 计算极化方向角
    angle_matrix = np.arctan2(grid_Py, grid_Px) * (180 / np.pi)  # 转换为度数
    # 将负角度转换为正角度（0°到360°）
    angle_matrix[angle_matrix < 0] += 360


    # 绘制，设置显示范围为x:[30,50], y:[30,50]
    fig,ax=plt.subplots(figsize=(8, 8))
    im = ax.imshow(grid_v.T, extent=(0, 100, 0, 100), origin='lower', cmap='binary', vmin=0, vmax=1)

    # 筛选v<0.6的区域并叠加矢量场P
    mask = grid_v < 0.4
    sampling_stride = 50  # 采样步长
    quiver_mask = mask[::sampling_stride, ::sampling_stride]  # 应用同样的采样步长来创建掩码

    qi = ax.quiver(grid_x[::sampling_stride, ::sampling_stride][quiver_mask], 
            grid_y[::sampling_stride, ::sampling_stride][quiver_mask],
            grid_Px[::sampling_stride, ::sampling_stride][quiver_mask],
            grid_Py[::sampling_stride, ::sampling_stride][quiver_mask],
            angle_matrix[::sampling_stride, ::sampling_stride][quiver_mask],  # 使用角度的值设置颜色
            scale=20, width=0.006, cmap=cc.cm.CET_C6s,norm=plt.Normalize(vmin=0, vmax=360))

    ax.set_xticks([])
    ax.set_yticks([])

    # 添加颜色条和设置颜色条的刻度及标题
    # cbar1 = plt.colorbar(im, ticks=[0, 0.5, 1])
    # cbar1.set_label(r'$\nu$', labelpad=20)
    # cbar2 = plt.colorbar(qi)
    # cbar2.set_label(r'$P$', labelpad=20)

    plt.tight_layout()
    plt.savefig(f'frame_{datafile:04d}.svg')
    plt.close(fig)

    return

for i in [200,700,900,1600,2000,2200]:
    plotfig1(i)