import re
import json
import pandas as pd
import anthropic

client = anthropic.Anthropic(api_key="YOUR_API_KEY")

# 每次batch改这三个变量
batch_id    = "msgbatch_01MCRM6CCPa9wdA2DTHKxFqR"
input_excel = "D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/progress/claude_2.xlsx"       # 对应这次batch的原始输入文件
final_excel = "D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/claude/claude_out2_new.xlsx" # 输出

def clean_score(x):
    try:
        x = int(x)
        return x if 1 <= x <= 9 else None
    except Exception:
        return None

def parse_json_content(content):
    content = re.sub(r'```json\s*', '', content)
    content = re.sub(r'```\s*', '', content)
    content = content.strip()
    content = re.sub(r',\s*}', '}', content)
    content = re.sub(r',\s*$', '', content)
    if not content.endswith('}'):
        content = content + '}'
    return json.loads(content)

df = pd.read_excel(input_excel, engine="openpyxl")
scores_map = {}

for entry in client.messages.batches.results(batch_id):
    custom_id = entry.custom_id
    if entry.result.type == "succeeded":
        content = entry.result.message.content[0].text
        try:
            data = parse_json_content(content)
            scores_map[custom_id] = {
                f"q{i}": clean_score(data.get(f"q{i}")) for i in range(1, 27)
            }
        except Exception:
            print(f"[仍然失败] {custom_id}: {content[:200]}")
            scores_map[custom_id] = {f"q{i}": None for i in range(1, 27)}

score_rows = []
for idx in df.index:
    key = f"row_{idx}"
    score_rows.append(scores_map.get(key, {f"q{i}": None for i in range(1, 27)}))

score_df = pd.DataFrame(score_rows)
df_out = pd.concat([df, score_df], axis=1)
df_out.to_excel(final_excel, index=False)

# 统计q26缺失情况
print(f"总行数: {len(df_out)}")
print(f"q26有值: {df_out['q26'].notna().sum()}")
print(f"q26缺失: {df_out['q26'].isna().sum()}")