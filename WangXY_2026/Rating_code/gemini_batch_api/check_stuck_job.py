from google import genai

# 请填入你之前使用的真实 API KEY
client = genai.Client(api_key="YOUR_API_KEY")
stuck_job_name = "batches/nldwg50490ga11s0qmz64p07hoem973t2tsc"

try:
    print(f"正在向服务器请求任务 {stuck_job_name} 的详细诊断信息...\n")
    job = client.batches.get(name=stuck_job_name)
    
    print("-" * 40)
    print(f"【任务名称】: {job.name}")
    print(f"【当前状态】: {job.state.name}")
    print(f"【创建时间】: {job.create_time}")
    print(f"【更新时间】: {job.update_time}")
    print("-" * 40)
    
    if job.state.name == "FAILED":
        print("【诊断结果】: 任务已失败！")
        print(f"【具体报错代码/信息】: {job.error}")
    elif job.state.name == "IN_PROGRESS":
        print("【诊断结果】: 任务仍在运行中（或者已卡死）。")
    elif job.state.name == "SUCCEEDED":
        print("【诊断结果】: 任务其实已经成功了！可能是之前的脚本没有正确捕捉到状态。")
    else:
        print(f"【其他状态】: {job.state.name}")
        
except Exception as e:
    print(f"查询失败，可能是 API 密钥错误或该任务已被服务器清理。报错详情: {e}")