import time
import json
import pandas as pd
from google import genai
from google.genai import types

# ================= 配置区 =================
API_KEY = "YOUR_API_KEY"
INPUT_EXCEL = "D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/AllRowData.xlsx"
OUTPUT_JSONL = "D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/gemini_batch_input.jsonl"
MODEL_NAME = "models/gemini-3.1-pro-preview"  # 推荐先用 flash 测试，性价比最高

client = genai.Client(api_key=API_KEY)

# ================= 1. 生成符合规范的 JSONL =================
QUESTIONS = [
    f"Q{i+1}: {q}" for i, q in enumerate([
        "您认为该思维片段多大可能是【抑郁症患者】说的?",
        "您认为该思维片段多大程度上是反映说话者【内部感受】的（躯体，内心感觉）?",
        "您认为该思维片段多大程度上是【内部导向】的（独立于外部刺激，自发产生的想法/画面）?",
        "您认为该思维片段多大程度上是【外部导向】的（由外界刺激诱发的想法，例如噪声、设备、实验任务等）?",
        "您认为该思维片段多大程度上是关注【自我】的?",
        "您认为该思维片段多大程度上是关注【他人】的?",
        "您认为该思维片段多大程度上说的是【正性的、积极】的事情？",
        "您认为该思维片段多大程度上说的是【负性的、消极】的事情？",
        "您认为该思维片段涉及【问题解决】的程度。（试图解决某个问题，或思考如何达到某个目标）？",
        "您认为该思维片段反映说话者【自责、内疚】的程度？",
        "您认为该思维片段反映说话者【担忧】的程度？",
        "您认为该思维片段反映说话者【生气】的程度？",
        "您认为该思维片段反映说话者【厌恶】的程度？",
        "您认为该思维片段反映说话者【恐惧】的程度？",
        "您认为该思维片段反映说话者【悲伤】的程度？",
        "您认为该思维片段反映说话者【惊喜】的程度？",
        "您认为该思维片段反映说话者【幸福】的程度？",
        "您认为该思维片段反映说话者【喜欢】的程度？",
        "您认为该思维片段反映说话者【愉悦】的程度？",
        "您认为该思维片段反映说话者【自卑】的程度？",
        "您认为该思维片段反映说话者【失望】的程度？",
        "您认为该思维片段反映说话者【反省】的程度？",
        "您认为该思维片段多大程度上是关注【过去】的事情？",
        "您认为该思维片段多大程度上是关注【未来】的事情？",
        "您认为该思维片段多大程度上是关注【当下】的事情？",
        "您认为该思维片段反映说话者【反刍】的程度？",
    ])
]

Q_BLOCK = "\n".join(QUESTIONS)

SYSTEM_PROMPT = """你是一个严格的心理学专家，请根据给定思维片段完成Q1-Q26评分。
除Q1外，所有题目均按以下标准评分：1=完全没有，9=几乎都是。
Q1按以下标准评分：1=完全不是，9=绝对是。
必须只输出JSON对象，键为q1...q26，值为1-9整数。
禁止输出任何解释或额外文本。"""

print("Step 1: Generating JSONL file...")
df = pd.read_excel(INPUT_EXCEL, engine="openpyxl")
written_count = 0

with open(OUTPUT_JSONL, "w", encoding="utf-8") as f:
    for idx, text in enumerate(df["Response"]):
        if not isinstance(text, str) or not text.strip():
            continue

        user_content = f"请按Q1-Q26评分：\n\n{Q_BLOCK}\n\n文本：{text}"
        
        # 核心修正：使用 generation_config 字段
        record = {
            "request": {
                "contents": [
                    {
                        "role": "user",
                        "parts": [{"text": f"{SYSTEM_PROMPT}\n\n{user_content}"}]
                    }
                ],
                "generation_config": {
                    "temperature": 0.0,
                    "response_mime_type": "application/json"
                }
            }
        }
        f.write(json.dumps(record, ensure_ascii=False) + "\n")
        written_count += 1

print(f"Success: {written_count} requests written to {OUTPUT_JSONL}")

# ================= 2. 上传文件 =================
print("\nStep 2: Uploading file to Gemini File API...")
uploaded_file = client.files.upload(
    file=OUTPUT_JSONL,
    config=types.UploadFileConfig(
        display_name="llm-rating-task-data",
        mime_type="application/jsonl"
    )
)
print(f"File uploaded: {uploaded_file.name}")

# ================= 3. 等待文件生效 (重要) =================
print("Waiting for file to become ACTIVE...")
while True:
    file_info = client.files.get(name=uploaded_file.name)
    if file_info.state.name == "ACTIVE":
        print("File is now ACTIVE.")
        break
    elif file_info.state.name == "FAILED":
        raise Exception("File processing failed.")
    time.sleep(2)

# ================= 4. 创建 Batch Job =================
print("\nStep 3: Creating Batch Job...")
try:
    batch_job = client.batches.create(
        model=MODEL_NAME,
        src=uploaded_file.name,
        config={
            "display_name": "psychology-rating-batch-001"
        }
    )
    print("-" * 30)
    print(f"Batch Job Created!")
    print(f"Job Name: {batch_job.name}")
    print(f"Initial State: {batch_job.state.name}")
    print("-" * 30)
    print("你可以通过 client.batches.get(name='{}') 来查看进度。".format(batch_job.name))

except Exception as e:
    print(f"Error creating job: {e}")