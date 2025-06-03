import pandas as pd
import streamlit as st
from io import BytesIO
from utils.utils import *
from utils.metabolomics import *
from utils.cardiotoxicity import *
from utils.neurotoxicity import *
from utils.hepatotoxicity import *

# Инициализация session state для fc_mode, если его еще нет
if "fc_mode" not in st.session_state:
    st.session_state["fc_mode"] = 'ratio (B/A)'  # Значение по умолчанию

st.title("Органотоксичность")

panel = main_radio_button_study()



if panel == "Метаболомика":
   metabolomika_app()
elif panel == "Кардиотоксичность":
   st.subheader("Кардиотоксичность")
elif panel == "Нейротоксичность":
   st.subheader("Нейротоксичность")
else:
   st.subheader("Гепатотоксичность")
