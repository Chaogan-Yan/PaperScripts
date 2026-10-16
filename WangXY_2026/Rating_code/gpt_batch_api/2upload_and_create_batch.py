import os
from openai import OpenAI

client = OpenAI(api_key=os.environ["OPENAI_API_KEY"])

batch_input_file = client.files.create(
    file=open("D:/2025THU/SELF_OTHER/LLMRating_Task/Raw/batch_input_last2.jsonl", "rb"),
    purpose="batch"
)

print("Uploaded file id:", batch_input_file.id)

