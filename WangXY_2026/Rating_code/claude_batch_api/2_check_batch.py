# 2_check_batch.py
import os
import anthropic

client = anthropic.Anthropic(api_key="YOUR_API_KEY")

batch_id = "msgbatch_01XeDua3YfB7kodruZhERUHm"  # 换成你的 batch id

batch = client.messages.batches.retrieve(batch_id)

print(f"Status: {batch.processing_status}")   # in_progress / ended / canceling
print(f"成功: {batch.request_counts.succeeded}")
print(f"失败: {batch.request_counts.errored}")
print(f"处理中: {batch.request_counts.processing}")
print(f"结果可下载: {batch.results_url}")