import json
import pandas as pd

# ===== 1. 改成你的路径 =====
input_excel = "D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/AllRowData_last2.xlsx"
output_jsonl = "D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/batch_output_last2.jsonl"
final_excel = "D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/gpt5.4_last2.xlsx"

def clean_score(x):
    try:
        x = int(x)
        return x if 1 <= x <= 9 else None
    except Exception:
        return None

# 读入原始 Excel
df = pd.read_excel(input_excel, engine="openpyxl")

# 用来存每一行的评分结果
scores_map = {}

# 读取 batch 输出
with open(output_jsonl, "r", encoding="utf-8") as f:
    for line in f:
        item = json.loads(line)

        custom_id = item.get("custom_id")
        response_body = item.get("response", {}).get("body", {})
        choices = response_body.get("choices", [])

        if not custom_id or not choices:
            continue

        content = choices[0]["message"]["content"]

        try:
            data = json.loads(content)
            scores_map[custom_id] = {
                f"q{i}": clean_score(data.get(f"q{i}")) for i in range(1, 27)
            }
        except Exception:
            scores_map[custom_id] = {
                f"q{i}": None for i in range(1, 27)
            }

# 按原始 Excel 的行号组装结果
score_rows = []
for idx in df.index:
    key = f"row_{idx}"
    score_rows.append(scores_map.get(key, {f"q{i}": None for i in range(1, 27)}))

# 合并成你熟悉的 Excel 形式
score_df = pd.DataFrame(score_rows)
df_out = pd.concat([df, score_df], axis=1)

# 保存
df_out.to_excel(final_excel, index=False)

print(f"已保存到: {final_excel}")
print(f"原始行数: {len(df)}")
print(f"成功解析的评分行数: {sum(k in scores_map for k in [f'row_{i}' for i in df.index])}")