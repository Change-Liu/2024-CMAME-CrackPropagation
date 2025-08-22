import cairosvg
import os

def svg_to_pdf(source_folder, target_folder):
    # 确保目标文件夹存在
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    # 定义需要转换的文件名列表
    filelist=[2,3,4,5,6,7,9,10,11,12,13,14,16]
    valid_filenames = [f"Fig{i}.svg" for i in filelist]  # 生成文件名列表 Fig1.svg 到 Fig10.svg

    # 遍历源文件夹中的所有文件
    for filename in os.listdir(source_folder):
        if filename in valid_filenames:  # 检查文件名是否在需要转换的列表中
            # 定义源文件和目标文件的完整路径
            source_file = os.path.join(source_folder, filename)
            target_file = os.path.join(target_folder, filename.replace('.svg', '.pdf'))
            
            # 转换SVG到PDF
            cairosvg.svg2pdf(url=source_file, write_to=target_file)
            print(f"Converted {source_file} to {target_file}")

# 使用示例
source_folder = './'
target_folder = './'
svg_to_pdf(source_folder, target_folder)
