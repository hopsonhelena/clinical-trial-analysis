make setup:
pip install -r requirements.txt

make pipeline:
python load_data.py
python analysis.py

make dashboard:
streamlit run dashboard.py