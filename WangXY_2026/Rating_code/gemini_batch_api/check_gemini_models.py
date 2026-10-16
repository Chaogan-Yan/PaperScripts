from google import genai

client = genai.Client(api_key="YOUR_API_KEY")

for m in client.models.list():
    try:
        name = getattr(m, "name", None)
        actions = getattr(m, "supported_actions", None) or getattr(m, "supported_generation_methods", None)
        print(name, actions)
    except Exception as e:
        print("error:", e)