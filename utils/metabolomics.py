import streamlit as st
import pandas as pd
from datetime import datetime
import io
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px


def plot_fold_change(fc_df, selected_drugs, selected_metabolites):
    if not selected_drugs or not selected_metabolites:
        st.info("Выберите хотя бы один препарат и метаболит для отображения графика.")
        return
    
    data = fc_df[
        (fc_df['Drug'].isin(selected_drugs)) & (fc_df['Group'] == 'test')
    ].copy()

    long_df = data.melt(
        id_vars=['Drug', 'Concentration'],
        value_vars=selected_metabolites,
        var_name='Metabolite',
        value_name='FoldChange'
    )

    long_df['Drug_Conc'] = long_df['Drug'].astype(str) + ' (' + long_df['Concentration'].astype(str) + ')'

    st.write(long_df)

    fig = px.bar(
        long_df,
        x='Metabolite',
        y='FoldChange',
        color='Drug_Conc',
        barmode='group',
        title="Fold Change по метаболитам",
        labels={
            "FoldChange": "Fold Change (Test / Control)",
            "Metabolite": "Метаболиты"
        },
        height=600
    )

    # Добавим горизонтальную линию (y=1.0)
    fig.add_shape(
        type="line",
        x0=-0.5,
        x1=len(selected_metabolites)-0.5,
        y0=1.0,
        y1=1.0,
        line=dict(color="red", dash="dash"),
    )

    st.plotly_chart(fig, use_container_width=True)

def calculate_fold_change(df, mode='ratio'):
    """
    Рассчитывает Fold Change для метаболитов между тестовыми и контрольными группами.
    Добавляет строку для контроля (Control vs Control), где FC = 1 (ratio) или 0 (difference),
    со средней концентрацией и диапазоном дат.
    :param df: DataFrame с исходными данными
    :param mode: 'ratio' для B/A или 'difference' для (B-A)/A
    :return: DataFrame с рассчитанными Fold Change
    """
    metabolite_cols = [col for col in df.columns 
                      if col not in ['Drug', 'Experiment date', 'Group', 'Concentration']]
    
    results = []
    
    for drug in df['Drug'].unique():
        drug_data = df[df['Drug'] == drug]
        control_data = drug_data[drug_data['Group'] == 'control_neg']
        control_means = control_data[metabolite_cols].mean()

        # Средняя концентрация контролей
        avg_concentration = control_data['Concentration'].mean()

        # Диапазон дат эксперимента
        date_min = control_data['Experiment date'].min()
        date_max = control_data['Experiment date'].max()
        if pd.notnull(date_min) and pd.notnull(date_max):
            date_range_str = f"{date_min} — {date_max}"
        else:
            date_range_str = None

        # Добавляем строку "контроль против контроля"
        control_row = {
            'Drug': drug,
            'Experiment date': date_range_str,
            'Group': 'control_vs_control',
            'Concentration': avg_concentration
        }
        for metabolite in metabolite_cols:
            control_row[metabolite] = 1.0 if mode == 'ratio' else 0.0
        results.append(control_row)
        
        # Обрабатываем тестовые строки
        test_data = drug_data[drug_data['Group'] == 'test']
        
        for _, row in test_data.iterrows():
            result_row = {
                'Drug': drug,
                'Experiment date': row['Experiment date'],
                'Group': 'test',
                'Concentration': row['Concentration']
            }
            for metabolite in metabolite_cols:
                control_val = control_means[metabolite]
                test_val = row[metabolite]
                
                if mode == 'ratio' and control_val != 0:
                    fc = test_val / control_val
                elif mode == 'difference' and control_val != 0:
                    fc = (test_val - control_val) / control_val
                else:
                    fc = None
                result_row[metabolite] = fc
            results.append(result_row)
    
    return pd.DataFrame(results)

def load_and_preprocess_data(uploaded_file):
    """Загружает и предобрабатывает данные из Excel файла."""
    try:
        df = pd.read_excel(uploaded_file)
        
        # Проверка обязательных колонок
        required_columns = ['Drug', 'Experiment date', 'Group', 'Concentration']
        for col in required_columns:
            if col not in df.columns:
                st.error(f"Отсутствует обязательная колонка: {col}")
                return None
        
        # Приведение названий колонок к нижнему регистру для удобства
        #df.columns = df.columns.str.strip().str.lower()
        
        # Предобработка даты
        if 'Experiment date' in df.columns:
            df['Experiment date'] = pd.to_datetime(df['Experiment date'], errors='coerce').dt.date
        
        # Приведение группы к нижнему регистру и проверка допустимых значений
        if 'Group' in df.columns:
            #df['Group'] = df['Group'].astype(str).str.lower().str.strip()
            valid_groups = ['control_neg', 'test']
            if not df['Group'].isin(valid_groups).all():
                st.warning("Обнаружены нестандартные значения в колонке Group")
        
        # Преобразование концентрации в числовой формат
        if 'Concentration' in df.columns:
            df['Concentration'] = pd.to_numeric(df['Concentration'], errors='coerce')
        
        return df
    
    except Exception as e:
        st.error(f"Ошибка при загрузке файла: {str(e)}")
        return None

def display_metabolite_data(df):
    """Отображает данные метаболитов в интерактивной таблице."""
    if df is None:
        return
    
    st.subheader("Данные метаболитов")
    
    # Определяем колонки с метаболитами (все, кроме обязательных)
    metabolite_cols = [col for col in df.columns 
                      if col.lower() not in ['Drug', 'Experiment date', 'Group', 'Concentration']]
    
    if not metabolite_cols:
        st.warning("Не обнаружены колонки с метаболитами")
        return
    
    # Показываем весь датафрейм
    st.dataframe(df)

def metabolomika_app():
    """Основная функция приложения для анализа данных метаболомики."""
    st.write("Загрузите Excel-файл с данными метаболитов для анализа")
    
    # Загрузка файла
    uploaded_file = st.file_uploader("Выберите Excel файл", type=["xlsx", "xls"])
    
    if uploaded_file is not None:
        # Загрузка и предобработка данных
        df = load_and_preprocess_data(uploaded_file)
        
        if df is not None:
            # Отображение основной информации о данных
            with st.sidebar:
                st.subheader("Общая информация о загруженных данных")
                
                # Уникальные препараты
                with st.expander(f"Уникальные препараты ({len(df['Drug'].unique())})"):
                    st.write(", ".join(sorted(df['Drug'].unique())))
                
                # Диапазоны дат и концентраций по препаратам
                with st.expander("Детали по препаратам"):
                    for drug in sorted(df['Drug'].unique()):
                        drug_data = df[df['Drug'] == drug]
                        drug_dates = drug_data['Experiment date']
                        
                        control_counts = drug_data['Group'].value_counts().to_dict()
                        control_info = ", ".join([f"{k}: {v}" for k, v in control_counts.items()])
                        
                        test_concentrations = drug_data[drug_data['Group'] == 'test']['Concentration']
                        conc_range = f"{test_concentrations.min()} - {test_concentrations.max()}" if not test_concentrations.empty else "N/A"
                        
                        st.write(f"""
                        **{drug}**  
                        • Даты: {drug_dates.min()} - {drug_dates.max()}  
                        • Контрольные группы: {control_info}  
                        • Диапазон концентраций (Test): {conc_range}
                        """)
                
                # Выбор метаболита через selectbox
                metabolite_cols = [col for col in df.columns 
                                if col.lower() not in ['Drug', 'Experiment date', 'Group', 'Concentration']]
                
                st.selectbox(
                    f"Выберите метаболит ({len(metabolite_cols)} доступно)",
                    sorted(metabolite_cols),
                    index=0,
                    key="metabolite_selector"
                )
                
                # Добавляем интерфейс для расчета Fold Change
                st.subheader("Анализ Fold Change")
                fc_mode = st.radio(
                    "Режим расчета Fold Change",
                    ('ratio (B/A)', 'difference ((B-A)/A)'),
                    index=0
                )
                calculate_fc = st.button("Рассчитать Fold Change")
            
            # Отображение данных метаболитов
            display_metabolite_data(df)
            
            # Расчет и отображение Fold Change по запросу
            if calculate_fc and df is not None:
                st.subheader("Результаты расчета Fold Change")
                
                # Определяем режим расчета
                mode = 'ratio' if fc_mode == 'ratio (B/A)' else 'difference'
                
                # Рассчитываем Fold Change
                fc_df = calculate_fold_change(df, mode=mode)

                # ✅ Сохраняем в сессию
                st.session_state['fc_df'] = fc_df
                
                # Отображаем результаты
                st.dataframe(fc_df)
                
                # Добавляем возможность скачать результаты
                output_fc = io.BytesIO()
                with pd.ExcelWriter(output_fc, engine='openpyxl') as writer:
                    fc_df.to_excel(writer, index=False, sheet_name='Fold_Change')
                output_fc.seek(0)
                
                st.download_button(
                    label="Скачать результаты Fold Change",
                    data=output_fc,
                    file_name="fold_change_results.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )

            if 'fc_df' in st.session_state:
                fc_df = st.session_state['fc_df']
                
                available_drugs = sorted(fc_df['Drug'].unique())
                available_metabolites = [col for col in fc_df.columns if col not in ['Drug', 'Experiment date', 'Group', 'Concentration']]

                if 'show_plot' not in st.session_state:
                    st.session_state['show_plot'] = False

                with st.form("plot_form"):
                    selected_drugs = st.multiselect("Выберите препарат(ы)", available_drugs, default=available_drugs, key="drug_select")
                    selected_metabolites = st.multiselect("Выберите метаболиты", available_metabolites, default=available_metabolites, key="met_select")
                    submitted = st.form_submit_button("Показать график")

                    if submitted:
                        st.session_state['show_plot'] = True

                if st.session_state['show_plot']:
                    plot_fold_change(fc_df, selected_drugs, selected_metabolites)    

# Запуск приложения
if __name__ == "__main__":
    metabolomika_app()