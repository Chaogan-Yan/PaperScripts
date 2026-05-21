# 卡阈值 排序

import os
import pandas as pd

def process_csv_files(directory):
    # 遍历指定文件夹中的所有CSV文件
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            filepath = os.path.join(directory, filename)
            
            # 读取CSV文件内容
            df = pd.read_csv(filepath)
            
            # 删除'Vertex'列中值小于100的所有行
            df = df[df['Vertex'] >= 100]
            
            # 根据'Peak_Intensity'列的绝对值对表格排序
            df = df.reindex(df['Peak_Intensity'].abs().sort_values(ascending=False).index)
            
            # 为新文件创建一个名字
            new_filename = os.path.splitext(filename)[0] + '_processed.csv'
            new_filepath = os.path.join(directory, new_filename)
            
            # 保存处理后的文件为新文件
            df.to_csv(new_filepath, index=False)

def process_single_csv_files(directory, filename):
    filepath = os.path.join(directory, filename)
    # 读取CSV文件内容
    df = pd.read_csv(filepath)

    # 删除'Vertex'列中值小于100的所有行
    df = df[df['Vertex'] >= 100]

    # 根据'Peak_Intensity'列的绝对值对表格排序
    df = df.reindex(df['Peak_Intensity'].abs().sort_values(ascending=False).index)

    # 为新文件创建一个名字
    new_filename = os.path.splitext(filename)[0] + '_processed.csv'
    new_filepath = os.path.join(directory, new_filename)

    # 保存处理后的文件为新文件
    df.to_csv(new_filepath, index=False)


# 调用函数并提供您的文件夹路径  
process_single_csv_files('', '')

