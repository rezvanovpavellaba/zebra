import streamlit as st

# Опции для радиокнопок
time_units = ['с', 'мин', 'ч', 'д', 'неделя', 'другое (ввести единицы измерения)']

concentration_units = [
    'пг/мл', 'нг/мл', 'мкг/мл', 'мг/мл', 'пмоль/мл', 'нмоль/мл', 'мкмоль/мл', 'ммоль/мл',
    'пг/л', 'нг/л', 'мкг/л', 'мг/л', 'пмоль/л', 'нмоль/л', 'мкмоль/л', 'ммоль/л',
    'другое (ввести единицы измерения)'
]
dose_units = [
    'пг', 'нг', 'мкг', 'мг', 'г', 'мкг/кг', 'мг/кг', 'г/кг', 'мкг/м²', 'мг/м²', 'г/м²',
    'пмоль', 'нмоль', 'мкмоль', 'ммоль', 'моль', 'мкмоль/кг', 'ммоль/кг', 'моль/кг',
    'мкмоль/м²', 'ммоль/м²', 'моль/м²', 'другое (ввести единицы измерения)'
]

organ_concentration_units = [
    'мг/г', 'мкг/г', 'нг/г', 'мг/мл', 'мкг/мл', 'нг/мл', 'ммоль/л', 'моль/л',
    'другое (ввести единицы измерения)'
]

# Функция для отображения радиокнопок с пользовательским вводом
def radio_with_custom_input(label, options, session_key, selector_research, key):
    unique_key = f"{session_key}_{selector_research}"  # Объединяем сессии с уникальным идентификатором
    custom_selected_key = f"custom_{unique_key}_selected"

    if f"unit_choice{unique_key}" not in st.session_state:
        st.session_state[f"unit_choice{unique_key}"] = 0  # Значение по умолчанию
    
    if custom_selected_key not in st.session_state:
        st.session_state[custom_selected_key] = False

    if st.session_state[custom_selected_key]:
        selected_option = 'другое (ввести единицы измерения)'
    else:
        selected_option = st.session_state.get(unique_key, options[0])
    
    if key is not None:
       selected_option = st.radio(label, options, index=options.index(selected_option),key = selector_research)
    else:
       selected_option = st.radio(label, options, index=options.index(selected_option))
    
    if selected_option == 'другое (ввести единицы измерения)':
        st.session_state[custom_selected_key] = True

        custom_value = st.text_input(f"Введите единицы измерения", st.session_state.get(f"custom_{unique_key}", ""),key = key)

        st.session_state[f"custom_{unique_key}"] = custom_value
        return custom_value
    else:
        st.session_state[unique_key] = selected_option
        st.session_state[custom_selected_key] = False
        return selected_option
    
# Отдельные функции для каждой радиокнопки

def select_time_unit(selector_research,key = None):
    with st.expander("Выбрать единицы измерения времени"):
        return radio_with_custom_input("Выберите единицу времени", time_units, 'selected_time',selector_research,key)

def select_concentration_unit(selector_research,key,name_drug):
    with st.expander(f"Выбрать единицы измерения концентрации для ({name_drug})"):
        return radio_with_custom_input("Выберите единицу концентрации", concentration_units, 'selected_concentration',selector_research,key)

def select_dose_unit(selector_research,key = None):
    with st.expander("Выбрать единицы измерения дозы"):
        return radio_with_custom_input("Выберите единицу дозы", dose_units, 'selected_dose',selector_research,key)

def select_organ_concentration_unit(selector_research,key = None):
    with st.expander("Выбрать единицы измерения концентрации в органах"):
        return radio_with_custom_input("Выберите единицу концентрации в органах", organ_concentration_units, 'selected_organ_concentration',selector_research,key)

