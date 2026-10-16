import os
from openai import OpenAI

client = OpenAI(api_key=os.environ["OPENAI_API_KEY"])

file_id = "file-CYNH54faeTyE5ynucE4QeJ"  # output file id

result = client.files.content(file_id)

with open("batch_output_rest2.jsonl", "wb") as f:
    f.write(result.read())

print("Downloaded batch_output_rest2.jsonl")