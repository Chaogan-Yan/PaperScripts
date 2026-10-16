import os
from openai import OpenAI

client = OpenAI(api_key=os.environ["OPENAI_API_KEY"])

batch_job = client.batches.create(
    input_file_id="file-PbhE8hsisjTheiaFrX3dvU",   # 这里放你刚才的 file id
    endpoint="/v1/chat/completions",
    completion_window="24h"
)

print("Batch id:", batch_job.id)
print("Status:", batch_job.status)