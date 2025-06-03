import streamlit as st
import pandas as pd
import numpy as np
import io
import plotly.express as px
from scipy.stats import ttest_ind

# Инициализация session state для fc_mode, если его еще нет
if 'fc_mode' not in st.session_state:
    st.session_state["fc_mode"] = 'ratio (B/A)'  # Значение по умолчанию


def plot_volcano(fc_df, selected_drugs, p_value_threshold=0.05, log2fc_threshold=1.0):
    """
    Строит Volcano Plot с реальными p-values.
    Данные берутся из fc_df (должны быть колонки '*_pvalue').
    """
    if not selected_drugs:
        st.warning("Выберите препараты для Volcano Plot.")
        return
    
    # Фильтруем данные: только тестовая группа и максимальная концентрация
    volcano_data = []
    for drug in selected_drugs:
        drug_data = fc_df[(fc_df['Drug'] == drug) & (fc_df['Group'] == 'test')]
        if not drug_data.empty:
            max_conc = drug_data['Concentration'].max()
            max_conc_data = drug_data[drug_data['Concentration'] == max_conc].copy()
            volcano_data.append(max_conc_data)
    
    if not volcano_data:
        st.error("Нет данных для Volcano Plot.")
        return
    
    volcano_df = pd.concat(volcano_data)
    
    # Собираем метаболиты и их p-values
    metabolite_cols = [col for col in volcano_df.columns 
                      if col not in ['Drug', 'Group', 'Concentration'] 
                      and not col.endswith('(pvalue)')]
    
    # Преобразуем в "длинный" формат
    long_data = []
    for _, row in volcano_df.iterrows():
        for metabolite in metabolite_cols:
            log2fc = row[metabolite]
            p_value = row[f'{metabolite} (pvalue)']
            
            long_data.append({
                'Drug': row['Drug'],
                'Concentration': row['Concentration'],
                'Metabolite': metabolite,
                'log2FC': log2fc,
                'p_value': p_value,
                '-log10(p_value)': -np.log10(p_value) if p_value > 0 else 10  # избегаем деления на 0
            })
    
    long_df = pd.DataFrame(long_data)
    
    # Строим график
    fig = px.scatter(
        long_df,
        x='log2FC',
        y='-log10(p_value)',
        color='Drug',
        hover_data=['Metabolite', 'Concentration', 'p_value'],
        title=f"Volcano Plot (макс. концентрация, p < {p_value_threshold}, |log2FC| > {log2fc_threshold})",
        labels={
            'log2FC': 'log₂(Fold Change)',
            '-log10(p_value)': '-log₁₀(p-value)'
        },
        height=600
    )
    
    # Добавляем пороговые линии
    fig.add_shape(
        type='line',
        x0=-log2fc_threshold,
        x1=log2fc_threshold,
        y0=-np.log10(p_value_threshold),
        y1=-np.log10(p_value_threshold),
        line=dict(color='red', dash='dash'),
    )
    
    fig.add_shape(
        type='line',
        x0=-log2fc_threshold,
        x1=-log2fc_threshold,
        y0=0,
        y1=long_df['-log10(p_value)'].max() + 1,
        line=dict(color='gray', dash='dot'),
    )
    
    fig.add_shape(
        type='line',
        x0=log2fc_threshold,
        x1=log2fc_threshold,
        y0=0,
        y1=long_df['-log10(p_value)'].max() + 1,
        line=dict(color='gray', dash='dot'),
    )
    
    st.plotly_chart(fig, use_container_width=True)


def calculate_descriptive_stats_new(df, group_cols, value_cols):
    """
    Функция расчета описательной статистики в требуемом формате:
    - Drug, Group, Concentration - заполняются только при смене группы
    - Parameter - содержит названия статистик (count, mean, std и т.д.)
    - Остальные колонки - метаболиты
    """
    # Определяем агрегационные функции
    def q1(x):
        return np.percentile(x, 25)
    
    def q3(x):
        return np.percentile(x, 75)
    
    # Основные статистики
    stats = df.groupby(group_cols)[value_cols].agg(
        ['count', 'mean', 'std', 'min', 'median', 'max', q1, q3]
    )
    
    # Переименовываем колонки (убираем MultiIndex)
    stats.columns = ['_'.join(col).strip() for col in stats.columns.values]
    
    # Сбрасываем индекс и преобразуем в нужный формат
    stats = stats.reset_index()
    
    # Создаем длинный формат таблицы
    stats = stats.melt(
        id_vars=group_cols,
        var_name='temp',
        value_name='value'
    )
    
    # Разделяем метаболит и параметр
    stats[['Metabolite', 'Parameter']] = stats['temp'].str.split('_', n=1, expand=True)
    stats.drop(columns=['temp'], inplace=True)
    
    # Преобразуем обратно в широкий формат
    stats = stats.pivot_table(
        index=group_cols + ['Parameter'],
        columns='Metabolite',
        values='value'
    ).reset_index()
    
    # Возвращаем к нормальному порядку колонок
    stats.columns.name = None
    column_order = group_cols + ['Parameter'] + value_cols
    stats = stats[column_order]
    
    return stats


def plot_fold_change(fc_df, selected_drugs, selected_metabolites):
    if not selected_drugs or not selected_metabolites:
        st.info("Выберите хотя бы один препарат и метаболит для отображения графика.")
        return
    
    # Определяем базовое значение линии в зависимости от режима
    if st.session_state["fc_mode"] == 'ratio (B/A)':
        baseline = 1.0
    elif st.session_state["fc_mode"] in ['difference ((B-A)/A)', 'log₂(B/A)']:
        baseline = 0.0
    else:
        baseline = 1.0  # fallback
    
    # Фильтруем данные и исключаем колонки с p-value
    data = fc_df[
        (fc_df['Drug'].isin(selected_drugs)) & 
        (fc_df['Group'] == 'test')
    ].copy()
    
    # Исключаем колонки с p-value из выбранных метаболитов
    selected_metabolites = [m for m in selected_metabolites if not m.endswith('(pvalue)')]
    
    # Проверяем, что остались метаболиты для отображения
    if not selected_metabolites:
        st.warning("Нет метаболитов для отображения после фильтрации p-value.")
        return
    
    long_df = data.melt(
        id_vars=['Drug', 'Concentration'],
        value_vars=selected_metabolites,
        var_name='Metabolite',
        value_name='FoldChange'
    )

    long_df['Drug_Conc'] = long_df['Drug'].astype(str) + ' (' + long_df['Concentration'].astype(str) + ')'

    fig = px.bar(
        long_df,
        x='Metabolite',
        y='FoldChange',
        color='Drug_Conc',
        barmode='group',
        title=f"Fold Change по метаболитам ({st.session_state["fc_mode"]})",
        labels={
            "FoldChange": f"Fold Change ({st.session_state["fc_mode"]})",
            "Metabolite": "Метаболиты"
        },
        height=600
    )

    # Добавляем горизонтальную линию с учетом режима
    fig.add_shape(
        type="line",
        x0=-0.5,
        x1=len(selected_metabolites)-0.5,
        y0=baseline,
        y1=baseline,
        line=dict(color="red", dash="dash"),
    )

    st.plotly_chart(fig, use_container_width=True)

def calculate_fold_change_with_pvalues(df, mode='ratio'):
    """
    Расчёт Fold Change и p-values с использованием t-теста Стьюдента.
    Предполагает равенство дисперсий в группах.
    """
    metabolite_cols = [col for col in df.columns 
                      if col not in ['Drug', 'Group', 'Concentration']]
    
    results = []
    warnings = []
    
    for drug in df['Drug'].unique():
        drug_data = df[df['Drug'] == drug]
        control_data = drug_data[drug_data['Group'] == 'control_neg']
        test_data = drug_data[drug_data['Group'] == 'test']
        
        # Добавляем строку "контроль против контроля"
        control_row = {
            'Drug': drug,
            'Group': 'control_vs_control',
            'Concentration': control_data['Concentration'].mean()
        }
        for metabolite in metabolite_cols:
            if mode == 'ratio':
                control_row[metabolite] = 1.0 #(проверено в Excel) т.к число деленное на само себя
            elif mode == 'difference':
                control_row[metabolite] = 0.0 #(проверено в Excel) т.к 0 в числителе
            elif mode == 'log2_ratio':
                control_row[metabolite] = 0.0 #(проверено в Excel) т.к правила вычисления логарифма
            control_row[f'{metabolite} (pvalue)'] = 1.0  # p-value для контроля = 1 (проверено в Excel)
        results.append(control_row)
        
        # Обрабатываем тестовые группы
        for conc in test_data['Concentration'].unique():
            test_subset = test_data[test_data['Concentration'] == conc]
            
            result_row = {
                'Drug': drug,
                'Group': 'test',
                'Concentration': conc
            }
            
            for metabolite in metabolite_cols:
                control_vals = control_data[metabolite].values
                test_vals = test_subset[metabolite].values
                
                # Проверка количества данных
                if len(control_vals) < 2 or len(test_vals) < 2:  # Хотя бы 2 точки в тесте для FC
                    p_value = None
                    if len(test_vals) < 2:
                        warnings.append(f"{metabolite} (Drug: {drug}, Conc: {conc})")
                    continue
                
                # Рассчитываем Fold Change
                control_mean = np.mean(control_vals)
                if control_mean == 0:
                    fc = None
                else:
                    ratio = np.mean(test_vals) / control_mean
                    if mode == 'ratio':
                        fc = ratio
                    elif mode == 'difference':
                        fc = (np.mean(test_vals) - control_mean) / control_mean
                    elif mode == 'log2_ratio':
                        fc = np.log2(ratio) if ratio > 0 else None #(проверено в Excel) т.к нельзя вычислить логарифм нуля
                
                # Классический t-тест Стьюдента (равные дисперсии)
                try:
                    _, p_value = ttest_ind(control_vals, test_vals)
                except:
                    p_value = None  # в случае ошибки
                
                result_row[metabolite] = fc
                result_row[f'{metabolite} (pvalue)'] = p_value
            
            results.append(result_row)
    
    # Создаём DataFrame
    result_df = pd.DataFrame(results)
    
    # Переупорядочиваем столбцы: сначала все обычные, затем все p-value
    non_pvalue_cols = [col for col in result_df.columns if '(pvalue)' not in col]
    pvalue_cols = [col for col in result_df.columns if '(pvalue)' in col]
    
    # Упорядочиваем столбцы
    ordered_cols = non_pvalue_cols + pvalue_cols
    result_df = result_df[ordered_cols]
    
    return result_df, warnings


def load_and_preprocess_data(uploaded_file):
    """Загружает и предобрабатывает данные из Excel файла."""
    try:
        df = pd.read_excel(uploaded_file)
        
        # Проверка обязательных колонок
        required_columns = ['Drug', 'Group', 'Concentration']
        for col in required_columns:
            if col not in df.columns:
                st.error(f"Отсутствует обязательная колонка: {col}")
                return None
        
        # Приведение группы к нижнему регистру и проверка допустимых значений
        if 'Group' in df.columns:
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
    st.dataframe(df)


def metabolomika_app():
    """Основная функция приложения для анализа данных метаболомики."""
    st.title("Анализ данных метаболомики")
    st.write("Загрузите Excel-файл с данными метаболитов для анализа")
    
    uploaded_file = st.file_uploader("Выберите Excel файл", type=["xlsx", "xls"])
    
    if uploaded_file is not None:
        df = load_and_preprocess_data(uploaded_file)
        
        if df is not None:
            with st.sidebar:
                st.subheader("Общая информация о загруженных данных")
                
                with st.expander(f"Уникальные препараты ({len(df['Drug'].unique())})"):
                    st.write(", ".join(df['Drug'].unique()))
                
                with st.expander("Детали по препаратам"):
                    for drug in df['Drug'].unique():
                        drug_data = df[df['Drug'] == drug]
                        
                        control_counts = drug_data['Group'].value_counts().to_dict()
                        control_info = ", ".join([f"{k}: {v}" for k, v in control_counts.items()])
                        
                        test_concentrations = drug_data[drug_data['Group'] == 'test']['Concentration']
                        conc_range = f"{test_concentrations.min()} - {test_concentrations.max()}" if not test_concentrations.empty else "N/A"
                        
                        st.write(f"""
                        **{drug}**  
                        • Контрольные группы: {control_info}  
                        • Диапазон концентраций (Test): {conc_range}
                        """)
                
                metabolite_cols = [col for col in df.columns 
                                if col not in ['Drug', 'Group', 'Concentration']]
                
                st.selectbox(
                    f"Выберите метаболит ({len(metabolite_cols)} доступно)",
                    metabolite_cols,
                    index=0,
                    key="metabolite_selector"
                )
                
                st.subheader("Анализ Fold Change")
                fc_mode = st.radio(
                    "Режим расчета Fold Change",
                    ('ratio (B/A)', 'difference ((B-A)/A)', 'log₂(B/A)'),
                    index=['ratio (B/A)', 'difference ((B-A)/A)', 'log₂(B/A)'].index(st.session_state["fc_mode"])
                )
                calculate_fc = st.button("Рассчитать Fold Change")
            
            display_metabolite_data(df)
            
            # Отображаем описательную статистику исходных данных
            st.subheader("Описательная статистика исходных концентраций")
            
            # Используйте один вызов для всех данных:
            all_stats = calculate_descriptive_stats_new(
                df,
                group_cols=['Drug', 'Group', 'Concentration'],
                value_cols=metabolite_cols
            )

            st.dataframe(all_stats)
            
            # Кнопка скачивания статистики исходных данных
            output_original_stats = io.BytesIO()
            with pd.ExcelWriter(output_original_stats, engine='openpyxl') as writer:
                all_stats.to_excel(writer, index=False, sheet_name='Stats')
            output_original_stats.seek(0)
            
            st.download_button(
                label="Скачать статистику исходных данных",
                data=output_original_stats,
                file_name="original_metabolites_stats.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
            
            if calculate_fc and df is not None:
                if fc_mode != st.session_state["fc_mode"]:
                    st.session_state["fc_mode"] = fc_mode
                
                mode = 'ratio' if fc_mode == 'ratio (B/A)' else 'difference' if fc_mode == 'difference ((B-A)/A)' else 'log2_ratio'

                fc_df,warnings_fc = calculate_fold_change_with_pvalues(df, mode=mode)

                st.session_state['fc_df'] = fc_df
                st.session_state['warnings_fc'] = warnings_fc
                
            if 'fc_df' in st.session_state:
                fc_df = st.session_state['fc_df']
                warnings_fc = st.session_state['warnings_fc']
                
                st.subheader(f"Результаты расчета Fold Change ({st.session_state["fc_mode"]})")
                st.dataframe(fc_df)

                with st.sidebar:

                    if len(warnings_fc) !=0:

                        selected_warnings_fc = st.selectbox(
                            f"Отчет об ошибках в расчётах ({len(warnings_fc)} всего)",
                            warnings_fc,
                            index=0,
                            key="warnings_fc_select"
                        )

                        st.error(f"Недостаточно данных для {selected_warnings_fc}, проверьте исходные данные на количество повторностей")
                        

                # Кнопка скачивания полных результатов Fold Change
                output_fc = io.BytesIO()
                with pd.ExcelWriter(output_fc, engine='openpyxl') as writer:
                    fc_df.to_excel(writer, index=False, sheet_name='Fold_Change')
                output_fc.seek(0)
                
                st.download_button(
                    label="Скачать полные результаты Fold Change",
                    data=output_fc,
                    file_name=f"fold_change_results_{st.session_state["fc_mode"]}.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )
                
                # Графики Fold Change
                available_drugs = fc_df['Drug'].unique()

                # Получаем список метаболитов, исключая p-value и служебные колонки
                available_metabolites = [
                    col for col in fc_df.columns 
                    if col not in ['Drug', 'Group', 'Concentration'] 
                    and not col.endswith('(pvalue)')
                ]
                
                if 'show_plot' not in st.session_state:
                    st.session_state['show_plot'] = False
                
                with st.form("plot_form"):
                    selected_drugs = st.multiselect(
                        "Выберите препарат(ы)", 
                        available_drugs, 
                        default=available_drugs, 
                        key="drug_select"
                    )
                    selected_metabolites = st.multiselect(
                        "Выберите метаболиты", 
                        available_metabolites, 
                        default=available_metabolites, 
                        key="met_select"
                    )
                    submitted = st.form_submit_button("Перерисовать график")
                    
                    if submitted:
                        st.session_state['show_plot'] = True
                

                if st.session_state['show_plot']:
                    plot_fold_change(fc_df, selected_drugs, selected_metabolites)

                    # Добавляем Volcano Plot только для режима log₂(B/A)
                    if st.session_state["fc_mode"] == 'log₂(B/A)':
                        st.subheader("Volcano Plot")
                        st.write("""
                        **Интерпретация Volcano Plot:**
                        - Точки в верхних правом/левом углах — значимые изменения (большой |log2FC| и низкий p-value).
                        - Горизонтальная линия — порог значимости (p < 0.05).
                        - Вертикальные линии — порог изменения (|log2FC| > 1).
                        """)
                        
                        # Выбираем препараты для Volcano Plot (можно ограничить выбор)
                        volcano_drugs = st.multiselect(
                            "Выберите препараты для Volcano Plot",
                            selected_drugs,
                            default=selected_drugs,
                            key="volcano_drugs"
                        )
                        
                        plot_volcano(fc_df, volcano_drugs)
