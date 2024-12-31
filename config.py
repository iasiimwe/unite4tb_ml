import os

OPENAI_API_KEY = "xxx" # Add own API keys - only the OpenAI key required for this analysis
REPLICATE_API_TOKEN = "xxx"
WOLFRAM_ALPHA_APPID = "xxx"
GOOGLE_API_KEY = "xxx"
GOOGLE_CSE_ID = "xxx"
HuggingFace = "xxx"

def set_environment():
    variable_dict = globals().items()
    for key, value in variable_dict:
        if ("API" in key or "ID" in key) and isinstance(value, str):
            # Set the environment variable
            os.environ[key] = value
