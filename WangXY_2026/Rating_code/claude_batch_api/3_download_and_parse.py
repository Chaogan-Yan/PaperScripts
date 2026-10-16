# 3_download_and_parse.py
import os
import json
import pandas as pd
import anthropic

client = anthropic.Anthropic(api_key="YOUR_API_KEY")

batch_id = "msgbatch_01XeDua3YfB7kodruZhERUHm"   # 换成你的 batch id
input_excel  = "D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/claude_fail.xlsx"
final_excel  = "D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/claude_add.xlsx"

def clean_score(x):
    try:
        x = int(x)
        return x if 1 <= x <= 9 else None
    except Exception:
        return None

# 读入原始 Excel
df = pd.read_excel(input_excel, engine="openpyxl")
scores_map = {}

# 流式获取结果（Claude 不需要下载文件）
for entry in client.messages.batches.results(batch_id):
    custom_id = entry.custom_id

    if entry.result.type == "succeeded":
        # 提取文本内容
        content = entry.result.message.content[0].text
        try:
            data = json.loads(content)
            scores_map[custom_id] = {
                f"q{i}": clean_score(data.get(f"q{i}")) for i in range(1, 27)
            }
        except Exception:
            scores_map[custom_id] = {f"q{i}": None for i in range(1, 27)}

    elif entry.result.type == "errored":
        print(f"[错误] {custom_id}: {entry.result.error}")
        scores_map[custom_id] = {f"q{i}": None for i in range(1, 27)}

# 按原始 Excel 行号组装
score_rows = []
for idx in df.index:
    key = f"row_{idx}"
    score_rows.append(scores_map.get(key, {f"q{i}": None for i in range(1, 27)}))

score_df = pd.DataFrame(score_rows)
df_out = pd.concat([df, score_df], axis=1)
df_out.to_excel(final_excel, index=False)

print(f"已保存到: {final_excel}")
print(f"原始行数: {len(df)}")
print(f"成功解析: {sum(k in scores_map for k in [f'row_{i}' for i in df.index])}")