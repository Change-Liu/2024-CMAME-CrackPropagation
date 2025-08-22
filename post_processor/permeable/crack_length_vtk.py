from paraview.simple import *
import numpy as np
import matplotlib.pyplot as plt

# 获取活动源，即当前已经加载的文件
source = GetActiveSource()

# 获取时间步长信息
timeKeeper = GetTimeKeeper()
timeSteps = timeKeeper.TimestepValues

# 保存每个时间步的长度
lengths = []

# 遍历所有时间步
for time in timeSteps:
    # 设置当前时间步
    timeKeeper.Time = time
    Render()
    
    # 应用 Threshold 过滤器，筛选 v > 0.9 的区域
    threshold = Threshold(Input=source)
    threshold.Scalars = ['POINTS', 'v']  # 修改为你数据中的标量名称
    threshold.ThresholdMethod = "Above Upper Threshold"  # Uses only upper value
    threshold.UpperThreshold = 0.9
    
    # 获取 Threshold 后的输出数据的包围盒
    bounds = threshold.GetDataInformation().GetBounds()
    length_x = bounds[1] - bounds[0]
    length_y = bounds[3] - bounds[2]
    length_z = bounds[5] - bounds[4]

    # 计算对角线长度
    length = np.sqrt(length_x**2 + length_y**2 + length_z**2)
    lengths.append(length)

# 绘制长度随时间变化的图
plt.figure()
plt.plot(timeSteps, lengths, marker='o')
plt.xlabel('Time Step')
plt.ylabel('Length of v > 0.9 Region')
plt.title('Length Change Over Time')
plt.grid(True)
plt.show()
