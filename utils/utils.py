from streamlit_option_menu import option_menu

#основная радиокнопка
def main_radio_button_study():
    panel = option_menu(
        "",
        ["Метаболомика", "Кардиотоксичность","Нейротоксичность","Гепатотоксичность"],
        icons=["droplet-half", "heart-pulse","lightning-charge","recycle"],
        menu_icon="cast",
        default_index=0,
        orientation="horizontal",
        key=f"menu",
        styles = {
    "container": {"padding": "0px", "background-color": "#73b5f2"},  # Светло-голубой фон
    "icon": {"color": "#ffff", "font-size": "18px"},  # Голубые иконки
    "nav-link": {
        "font-size": "16px",
        "text-align": "center",
        "margin": "0px",
        "--hover-color": "#138abd",
        "color": "#ffff",
    },
    "nav-link-selected": {"background-color": "#4985c1", 'color': '#ffff',"font-weight": "normal","font-size": "18px"},  # Голубой активный пункт
}
    )
    return panel