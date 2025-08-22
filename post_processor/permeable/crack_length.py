import time
import pytesseract
from PIL import ImageGrab, Image, ImageOps, ImageDraw

# 设置截图区域的像素坐标，需要根据实际情况调整
# 这里假设你需要捕获屏幕上两个不同位置的数字
bbox1 = (800, 380, 1100, 490)  # (left, top, right, bottom)
bbox2 = (700, 500, 1300, 600)  # (left, top, right, bottom)

# 保存数据的文件路径
output_file = 'screen_data.txt'

# 主循环，每隔0.1秒读取屏幕数据并保存
while True:
    
    # 截取屏幕区域并转换为灰度图像
    screen = ImageGrab.grab()
    gray_image = ImageOps.grayscale(screen)  # 转换为灰度图像
    binary_image = gray_image.point(lambda p: p > 100 and 255)  # 二值化处理
    
    # 截取屏幕区域并转换为灰度图像
    screen = ImageGrab.grab()
    
    # 在截图上标记所选定的区域
    draw = ImageDraw.Draw(screen)
    draw.rectangle(bbox1, outline='red', width=2)
    draw.rectangle(bbox2, outline='blue', width=2)
    
    # 使用 pytesseract 进行文字识别
    number1 = pytesseract.image_to_string(screen.crop(bbox2),  config='--psm 7 --oem 2')
    number2 = pytesseract.image_to_string(screen.crop(bbox1), config='--psm 7 --oem 2')

    print(number1, number2)
    
    # 将识别结果写入文件
    with open(output_file, 'a') as f:
        f.write(f'{number1.strip()}, {number2.strip()}\n')
    
    # 显示带有标记的截图（可选）
    screen.show()
    
    # 每隔0.1秒执行一次
    time.sleep(30)
