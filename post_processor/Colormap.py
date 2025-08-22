import numpy as np
import matplotlib.pyplot as plt
import colorcet as cc
from matplotlib.colors import hsv_to_rgb

def create_annulus_with_angle(inner_radius, outer_radius, image_size):
    """创建一个包含圆环的二维数组，其中圆环的值为每个像素的角度"""
    x = np.linspace(-1, 1, image_size)
    y = np.linspace(-1, 1, image_size)
    X, Y = np.meshgrid(x, y)
    
    distances = np.sqrt(X**2 + Y**2)  # 计算每个点到中心的距离
    angles = np.arctan2(Y, X)  # 计算每个点的角度

    annulus = np.logical_and(distances > inner_radius, distances < outer_radius)
    angle_map = np.where(annulus, angles, np.nan)  # 圆环外设置为NaN
    return angle_map

# 设置圆环的内外半径和图像尺寸
inner_radius = 0.2
outer_radius = 0.5
image_size = 500  # 图像尺寸

# 生成带角度的圆环数据
annulus_image = create_annulus_with_angle(inner_radius, outer_radius, image_size)

# 角度转换为颜色
angle_colors = (annulus_image + np.pi) / (2 * np.pi)  # 归一化到0-1

import colorcet as cc
# 创建colormap，并设置NaN值显示为白色
cmap = cc.cm.CET_C6s
cmap.set_bad('white')

# 应用colormap
color_image = cmap(angle_colors)

# 使用imshow绘制图像
plt.axis('off')  # 隐藏坐标轴
plt.imshow(color_image, origin='lower', extent=[-1, 1, -1, 1])
plt.tight_layout()
plt.savefig('colormap.svg')
plt.show()