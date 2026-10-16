from google import genai

client = genai.Client(api_key="YOUR_API_KEY")

# 填入你刚才创建任务后返回的 Job Name
#job_name = "batches/nldwg50490ga11s0qmz64p07hoem973t2tsc" 
job_name = "batches/8wvkx54os94ss098uhccq88c3z2fbtdodgf5"

try:
    batch_job = client.batches.get(name=job_name)
    
    print("-" * 30)
    print(f"任务名称: {batch_job.name}")
    print(f"当前状态: {batch_job.state.name}") # 常见状态: IN_PROGRESS, SUCCEEDED, FAILED
    print(f"创建时间: {batch_job.create_time}")
    
    if batch_job.state.name == "SUCCEEDED":
        print(f"恭喜！任务已完成。")
        print(f"结果文件位置: {batch_job.output_config.gcs_dest.output_uri_prefix if hasattr(batch_job.output_config, 'gcs_dest') else '检查 API 响应获取下载链接'}")
    elif batch_job.state.name == "FAILED":
        print(f"任务失败，原因: {batch_job.error}")
    else:
        print("任务仍在处理中，请稍后再试...")
    print("-" * 30)

except Exception as e:
    print(f"查询出错: {e}")