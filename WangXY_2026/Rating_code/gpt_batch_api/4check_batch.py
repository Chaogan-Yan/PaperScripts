import os
from openai import OpenAI

client = OpenAI(api_key=os.environ["OPENAI_API_KEY"])

batch = client.batches.retrieve("batch_69c78a24502c81909ab24832913e4d39")  # 换成你的 batch id

print("Status:", batch.status)
print("Output file id:", batch.output_file_id)
print("Error file id:", batch.error_file_id)
print("Errors:", batch.errors)