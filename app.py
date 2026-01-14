import streamlit as st
from utils.utils import *
from utils.metabolomics import *
from utils.cardiotoxicity import *
from utils.neurotoxicity import *
from utils.hepatotoxicity import *

st.set_page_config(page_title="Органотоксичность", layout="wide")
st.title("Органотоксичность")

panel = main_radio_button_study()



if panel == "Метаболомика":
   metabolomika_app()
elif panel == "Кардиотоксичность":
   cardiotoxicity_app()
elif panel == "Нейротоксичность":
   neurotoxicity_app()
else:
   st.subheader("Гепатотоксичность")