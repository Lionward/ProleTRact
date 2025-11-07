import subprocess
from pathlib import Path


def main():
    app_path = Path(__file__).resolve().parent / "app.py"
    subprocess.run(["streamlit", "run", str(app_path)], check=True)
