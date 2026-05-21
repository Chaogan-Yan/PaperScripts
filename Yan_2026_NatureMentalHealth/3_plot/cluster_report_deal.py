# 提取峰值坐标
import re
import csv

def write_to_csv(filename, columns, headers):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(headers)
        rows = zip(*columns)
        writer.writerows(rows)
    print(f"Data has been written to {filename}.")

def split_text(text):
    clusters = text.split("Cluster")
    return clusters


def extract_numbers(text):
    numbers = re.findall(r'-?\d+\.\d+|-?\d+', text) #负数 小数 整数
    # numbers = re.findall(r'\d+\.\d+|\d+', text)
    # numbers = re.findall(r'\d+', text)
    return numbers

def extract_words(text):
    pattern = r'\b(\w+)\s*\(Currently Displayed\)'
    matches = re.findall(pattern, text)
    return matches

import re

def extract_region(text):
    pattern = r'\[\d+\]\s*(\w+)'
    matches = re.findall(pattern, text)
    return matches

import re

def extract_and_join_words(text):
    pattern = r'\[\d+\]\s*(\w+)'
    matches = re.findall(pattern, text)
    
    if len(matches) > 3:
        joined_words = '/'.join(matches[1:4])  # 只连接前三个单词
        return joined_words
    elif len(matches) == 3:
        joined_words = '/'.join(matches[1:3])  # 只连接前三个单词
        return joined_words
    elif len(matches) == 2:
        return matches[1]
    
    return None





Vertex_Num = []
Peak_Index = []
Peak_Coord_X = []
Peak_Coord_Y = []
Peak_Coord_Z = []
Peak_Intensity = []
Peak_Label = []
other_info = []
# 替换 'filename.txt' 为您实际的文件名
path = '/' 
name = '' 
filename = path + name+'.txt'
with open(filename, 'r') as file:
    text = file.read()
result = split_text(text)
for i in range(int((len(result)-1)/3)):
    j = int(i*3)+2
    Vertex_Num.append(int(extract_numbers(result[j])[0]))
    vertex_info = extract_numbers(result[j+1])
    Peak_Index.append(vertex_info[2])
    Peak_Coord_X.append(vertex_info[3])
    Peak_Coord_Y.append(vertex_info[4])
    Peak_Coord_Z.append(vertex_info[5])
    Peak_Intensity.append(vertex_info[6])
    # label_info = extract_words(result[j+1])
    other = extract_and_join_words(result[j+1])
    # Peak_Label.append(label_info[0])
    other_info.append(other)




# 替换 'filename.csv' 为您想要保存的实际文件名
filename = path + name+'_result.csv'
# 替换 'Header 1', 'Header 2', 'Header 3' 为每列的标题
headers = ['Region', 'Peak Coord X', 'Peak Coord Y', 'Peak Coord Z', 'Peak Index', 'Vertex', 'Peak_Intensity']
# 将所有列表存储为一个二维列表
columns = [other_info, Peak_Coord_X, Peak_Coord_Y, Peak_Coord_Z, Peak_Index, Vertex_Num, Peak_Intensity]

write_to_csv(filename, columns, headers)

