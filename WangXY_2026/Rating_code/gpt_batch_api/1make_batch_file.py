import os
import json
import pandas as pd
from openai import OpenAI

client = OpenAI(api_key=os.environ["OPENAI_API_KEY"])

QUESTIONS = [
    "1.您认为该思维片段多大可能是【抑郁症患者】说的?",
    "2.您认为该思维片段多大程度上是反映说话者【内部感受】的（躯体，内心感觉）?",
    "3.您认为该思维片段多大程度上是【内部导向】的（独立于外部刺激，自发产生的想法/画面）?",
    "4.您认为该思维片段多大程度上是【外部导向】的（由外界刺激诱发的想法，例如噪声、设备、实验任务等）?",
    "5.您认为该思维片段多大程度上是关注【自我】的?",
    "6.您认为该思维片段多大程度上是关注【他人】的?",
    "7.您认为该思维片段多大程度上说的是【正性的、积极】的事情？",
    "8.您认为该思维片段多大程度上说的是【负性的、消极】的事情？",
    "9.您认为该思维片段涉及【问题解决】的程度。（试图解决某个问题，或思考如何达到某个目标）？",
    "10.您认为该思维片段反映说话者【自责、内疚】的程度？",
    "11.您认为该思维片段反映说话者【担忧】的程度？",
    "12.您认为该思维片段反映说话者【生气】的程度？",
    "13.您认为该思维片段反映说话者【厌恶】的程度？",
    "14.您认为该思维片段反映说话者【恐惧】的程度？",
    "15.您认为该思维片段反映说话者【悲伤】的程度？",
    "16.您认为该思维片段反映说话者【惊喜】的程度？",
    "17.您认为该思维片段反映说话者【幸福】的程度？",
    "18.您认为该思维片段反映说话者【喜欢】的程度？",
    "19.您认为该思维片段反映说话者【愉悦】的程度？",
    "20.您认为该思维片段反映说话者【自卑】的程度？",
    "21.您认为该思维片段反映说话者【失望】的程度？",
    "22.您认为该思维片段反映说话者【反省】的程度？",
    "23.您认为该思维片段多大程度上是关注【过去】的事情？",
    "24.您认为该思维片段多大程度上是关注【未来】的事情？",
    "25.您认为该思维片段多大程度上是关注【当下】的事情？",
    "26.您认为该思维片段反映说话者【反刍】的程度？（反刍：对负性情绪本身及其原因和可能后果的反复思考）",
]

QUESTION_IDS = list(range(1, 27))

SYSTEM_PROMPT = """你是一个严格的心理学专家，请根据给定思维片段完成Q1-Q26评分。
除Q1外，所有题目均按以下标准评分：1=完全没有，9=几乎都是。
Q1按以下标准评分：1=完全不是，9=绝对是。
必须只输出JSON对象，键为q1...q26，值为1-9整数。
禁止输出任何解释或额外文本。"""

Q_BLOCK = "\n".join([f"Q{i}: {q}" for i, q in zip(QUESTION_IDS, QUESTIONS)])

file_path = "D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/AllRowData_last2.xlsx"
out_jsonl = "D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/batch_input_last2.jsonl"

df = pd.read_excel(file_path, engine="openpyxl")

with open(out_jsonl, "w", encoding="utf-8") as f:
    for idx, text in enumerate(df["Response"]):
        if not isinstance(text, str) or not text.strip():
            continue

        user_prompt = f"""请按Q1-Q26评分：

{Q_BLOCK}

文本：{text}

只输出JSON对象。"""

        req = {
            "custom_id": f"row_{idx}",
            "method": "POST",
            "url": "/v1/chat/completions",
            "body": {
                "model": "gpt-5.4",
                "messages": [
                    {"role": "system", "content": SYSTEM_PROMPT},
                    {"role": "user", "content": user_prompt}
                ],
                "temperature": 0,
                "response_format": {"type": "json_object"}
            }
        }

        f.write(json.dumps(req, ensure_ascii=False) + "\n")

print(f"Done. JSONL saved to: {out_jsonl}")
print(f"Requests written: {sum(df['Response'].apply(lambda x: isinstance(x, str) and x.strip() != ''))}")