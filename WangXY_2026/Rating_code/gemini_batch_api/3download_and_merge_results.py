
import json
import pandas as pd
from google import genai

# ================= 配置区 =================
API_KEY = "YOUR_API_KEY"
BATCH_JOB_NAME = "batches/8wvkx54os94ss098uhccq88c3z2fbtdodgf5"
ORIGINAL_EXCEL = "D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/AllRowData_2_2.xlsx"
FINAL_OUTPUT_EXCEL = "D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/gemini_2_2.xlsx"

client = genai.Client(api_key=API_KEY)

# 1. 状态检查
print(f"正在查询任务状态: {BATCH_JOB_NAME}...")
job = client.batches.get(name=BATCH_JOB_NAME)
if "SUCCEEDED" not in job.state.name:
    print(f"❌ 任务未完成，当前状态为: {job.state.name}")
    exit()

# 2. 获取输出文件名称并下载 (核心修正点)
# 真正的结果是存储在 job.dest.file_name 中的一个云端文件
try:
    output_file_name = job.dest.file_name
    print(f"任务已成功！结果存储在云端文件: {output_file_name}")
    print("正在将结果文件下载到本地...")
    
    # 将文件以字节形式下载，并解码为字符串
    content_bytes = client.files.download(file=output_file_name)
    jsonl_text = content_bytes.decode('utf-8')
except AttributeError:
    print("❌ 无法获取到输出文件。任务可能已过期被清理，或者 API 发生变动。")
    exit()
except Exception as e:
    print(f"❌ 下载文件时发生错误: {e}")
    exit()

# 3. 按行解析下载下来的文件内容
results = []
error_count = 0

# 去除空行并将文本切分为按行的列表
lines = [line for line in jsonl_text.split('\n') if line.strip()]
print(f"下载完毕，共获取到 {len(lines)} 条请求记录，正在解析并对齐 Excel...")

for i, line in enumerate(lines):
    # 【防御机制】预设这一行的默认值为 None
    row_data = {f"q{j}": None for j in range(1, 27)} 
    
    try:
        # 将当前行解析为 Python 字典
        response_data = json.loads(line)
        
        # 检查是否成功生成了回复 (拦截 Safety Filter)
        if 'response' not in response_data or 'candidates' not in response_data['response']:
            print(f"⚠️ 第 {i+1} 行未生成有效内容 (大概率触发了安全拦截)。")
            error_count += 1
            results.append(row_data)
            continue
            
        # 提取模型真正回答的文本
        raw_text = response_data['response']['candidates'][0]['content']['parts'][0]['text']
        clean_json = raw_text.strip().replace("```json", "").replace("```", "")
        
        # 解析分数字典
        scores = json.loads(clean_json)
        
        # 将分数填入对应行
        for key, value in scores.items():
            if key in row_data:
                row_data[key] = value
                
    except json.JSONDecodeError:
        print(f"⚠️ 第 {i+1} 行数据或评分 JSON 解析失败。")
        error_count += 1
    except Exception as e:
        print(f"⚠️ 第 {i+1} 行发生未知解析错误: {e}")
        error_count += 1
        
    results.append(row_data)

print(f"✅ 数据解析完毕。其中 {error_count} 行异常（将显示为留空）。")

# 4. 对齐并合并 Excel
print("正在与原 Excel 文件合并...")
try:
    df_original = pd.read_excel(ORIGINAL_EXCEL)
    df_scores = pd.DataFrame(results)

    if len(df_scores) != len(df_original):
        print(f"🚨 致命错误拦截：抓取到的结果数 ({len(df_scores)}) 与 Excel 原行数 ({len(df_original)}) 不符！")
        print("请检查原始数据量。")
    else:
        df_final = pd.concat([df_original.reset_index(drop=True), df_scores], axis=1)
        df_final.to_excel(FINAL_OUTPUT_EXCEL, index=False)
        print(f"🎉 完美合并！带评分的结果已安全保存至: {FINAL_OUTPUT_EXCEL}")
except Exception as e:
    print(f"❌ 读写 Excel 发生错误: {e}")