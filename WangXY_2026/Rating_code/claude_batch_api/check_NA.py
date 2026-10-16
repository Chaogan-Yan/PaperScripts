# 3_download_and_parse.py
import os
import json
import pandas as pd
import anthropic

client = anthropic.Anthropic(api_key="YOUR_API_KEY")

batch_id = "msgbatch_011UGea3yDtbEgQ8U6WAqhVv"
input_excel  = "D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/progress/claude_2.xlsx"
final_excel  = "D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/progress/claude_out2.xlsx"

def clean_score(x):
    try:
        x = int(x)
        return x if 1 <= x <= 9 else None
    except Exception:
        return None

df = pd.read_excel(input_excel, engine="openpyxl")
scores_map = {}

# 诊断计数器
json_success = 0
json_fail = 0
first_entry_printed = False  # 只打印第一条的原始内容

for entry in client.messages.batches.results(batch_id):
    custom_id = entry.custom_id

    if entry.result.type == "succeeded":
        content = entry.result.message.content[0].text

        # 打印第一条原始内容，看格式是否正常
        if not first_entry_printed:
            print(f"=== 第一条原始返回内容 ===")
            print(f"custom_id: {custom_id}")
            print(f"content: {content[:300]}")
            print(f"==========================")
            first_entry_printed = True

        try:
            data = json.loads(content)
            scores_map[custom_id] = {
                f"q{i}": clean_score(data.get(f"q{i}")) for i in range(1, 27)
            }
            json_success += 1
        except Exception:
            # 打印解析失败的原始内容
            print(f"[JSON解析失败] {custom_id}: {content[:200]}")
            scores_map[custom_id] = {f"q{i}": None for i in range(1, 27)}
            json_fail += 1

    elif entry.result.type == "errored":
        print(f"[API错误] {custom_id}: {entry.result.error}")
        scores_map[custom_id] = {f"q{i}": None for i in range(1, 27)}

print(f"\n=== 诊断结果 ===")
print(f"JSON解析成功: {json_success}")
print(f"JSON解析失败: {json_fail}")

score_rows = []
for idx in df.index:
    key = f"row_{idx}"
    score_rows.append(scores_map.get(key, {f"q{i}": None for i in range(1, 27)}))

score_df = pd.DataFrame(score_rows)
df_out = pd.concat([df, score_df], axis=1)
df_out.to_excel(final_excel, index=False)

print(f"\n=== 保存结果 ===")
print(f"已保存到: {final_excel}")
print(f"原始行数: {len(df)}")
print(f"成功写入评分的行数: {sum(k in scores_map for k in [f'row_{i}' for i in df.index])}")