import streamlit as st
import pandas as pd
import numpy as np
import io
import plotly.express as px
from scipy.stats import ttest_ind
from utils.radio_unit import *
import re
from collections import defaultdict


def plot_fold_change_horizontal(fc_df, selected_drugs, selected_metabolites, list_measure_unit_concentration,apply_colors,
                              legend_position="Снизу", show_legend=True, opacity=0.8,
                              baseline_color="black", custom_colors=None):  # Добавлен новый параметр custom_colors
    
    # Предварительный расчет количества концентраций у каждого препарата в виде фрейма
    fc_df_for_count_conc = fc_df[fc_df['Concentration'] != 0]
    result_fc_df_for_count_conc = fc_df_for_count_conc.groupby('Drug')['Concentration'].nunique().reset_index()
    result_fc_df_for_count_conc.columns = ['Drug', 'Number_of_Concentrations']

    if not selected_drugs or not selected_metabolites:
        st.info("Выберите хотя бы один препарат и метаболит для отображения графика.")
        return
    
    # Определяем базовое значение линии
    baseline = 1.0 if st.session_state["fc_mode"] == 'ratio (B/A)' else 0.0
    
    # Фильтруем данные
    data = fc_df[
        (fc_df['Drug'].isin(selected_drugs)) & 
        (fc_df['Group'] == 'test')
    ].copy()
    
    # Исключаем p-value
    selected_metabolites = [m for m in selected_metabolites if not m.endswith('(pvalue)')]
    
    if not selected_metabolites:
        st.warning("Нет метаболитов для отображения после фильтрации p-value.")
        return
    
    # Подготовка данных
    long_df = data.melt(
        id_vars=['Drug', 'Concentration'],
        value_vars=selected_metabolites,
        var_name='Metabolite',
        value_name='FoldChange'
    )
    
    # Расчет для формирования списка с единицами измерения
    number_of_concentrations = result_fc_df_for_count_conc['Number_of_Concentrations'].tolist()

    result_measure_unit = []
    for measure_unit, number in zip(list_measure_unit_concentration, number_of_concentrations):
        for _ in range(number):
            result_measure_unit.append(measure_unit)

    result_measure_unit_final = result_measure_unit * len(selected_metabolites)

    long_df['Препарат (Концентрация)'] = [
        f"{drug} ({conc} {unit})" 
        for drug, conc, unit in zip(
            long_df['Drug'], 
            long_df['Concentration'], 
            result_measure_unit_final
        )
    ]

    # Рассчитываем оптимальные параметры
    num_metabolites = len(selected_metabolites)
    base_width = 600
    width_per_metabolite = max(100, base_width / max(1, num_metabolites*0.7))
    total_width = width_per_metabolite * num_metabolites
    
    # Определяем цветовую схему
    if custom_colors and len(custom_colors) == len(long_df['Препарат (Концентрация)'].unique()) and apply_colors:
        # Используем пользовательские цвета, если они предоставлены и их количество совпадает
        colors = custom_colors
    else:  # Custom
        colors = ["#0037ff", "#ff0000", "#3cff00", "#dbf75f", "#046f22",
                 "#1de0fa", "#e928af", "#11e490", "#f7ff01", "#e16d07"]

    combo_order = (
        fc_df[fc_df['Concentration'] != 0]
          .assign(
              combo=lambda d: d['Drug'] + ' (' 
                             + d['Concentration'].astype(str) 
                             + ' ' 
                             + result_measure_unit_final[:len(d)]  # или ваш unit-список
                             + ')'
          )['combo']
          .drop_duplicates()
          .tolist()
    )
        
    long_df = long_df.dropna(subset=['FoldChange']).reset_index(drop=True)

    long_df['Препарат (Концентрация)'] = pd.Categorical(
        long_df['Препарат (Концентрация)'],
        categories=combo_order,
        ordered=True
    )

    # Создаем график
    fig = px.bar(
        long_df,
        x='Metabolite',
        y='FoldChange',
        color='Препарат (Концентрация)',
        barmode='group',
        category_orders={'Препарат (Концентрация)': combo_order},
        color_discrete_sequence=colors[:len(combo_order)],
        labels={
            "FoldChange": f'Fold Change ({st.session_state["fc_mode"]})',
            "Metabolite": "Метаболиты",
            "color": f"Препарат (концентрация, единицы измерения)"
        },
        height=600,
        width=total_width
    )
    
    # Настройки положения легенды
    legend_orientation = 'h' if legend_position in ["Снизу", "Сверху"] else 'v'
    legend_y = -0.5 if legend_position == "Снизу" else 1.0 if legend_position == "Сверху" else 0.9
    legend_x = 0.5 if legend_position in ["Снизу", "Сверху"] else 1.1 if legend_position == "Справа" else -0.2
    
    fig.update_layout(legend_traceorder='normal')

    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(color='black'),
        margin=dict(l=50, r=50, t=50, b=50),
        title={
            'text': f'Изменение уровня метаболитов',
            'y':0.95,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {'size': 16}
        },
        xaxis=dict(
            title_font=dict(color='black'),
            tickfont=dict(color='black'),
            showline=False,
            mirror=False,
            tickangle=270
        ),
        yaxis=dict(
            title_font=dict(color='black'),
            tickfont=dict(color='black'),
            showline=False,
            mirror=False
        ),
        legend=dict(
            title_text=f'Препарат (концентрация, единицы измерения)',
            title_font=dict(color='black'),
            font=dict(color='black'),
            itemsizing='constant',
            bgcolor='rgba(255,255,255,0.5)',
            orientation=legend_orientation,
            yanchor="bottom" if legend_position == "Снизу" else "top" if legend_position == "Сверху" else "middle",
            xanchor="center" if legend_position in ["Снизу", "Сверху"] else "left" if legend_position == "Справа" else "right",
            y=legend_y,
            x=legend_x
        ) if show_legend else dict(visible=False),
        autosize=True
    )
    
    # Настраиваем прозрачность
    fig.update_traces(opacity=opacity)
    
    # Горизонтальная линия
    fig.add_shape(
        type="line",
        x0=-0.5,
        x1=num_metabolites-0.5,
        y0=baseline,
        y1=baseline,
        line=dict(color=baseline_color, dash="dash"),  # Используем новый параметр
    )

    fig.update_xaxes(
        categoryorder='array',
        categoryarray=selected_metabolites
    )

    # Конфигурация для скачивания
    config = {
        'toImageButtonOptions': {
            'format': 'jpeg',
            'filename': 'fold_change_plot',
            'scale': 4,
            'dpi': 600
        }
    }

    st.plotly_chart(fig, use_container_width=True, config=config)

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


def plot_volcano(fc_df, selected_drugs, p_value_threshold=0.05, log2fc_threshold=1.0, 
                custom_colors=None, show_legend=True, 
                hline_color='red', vline_color='gray'):
    """
    Улучшенный Volcano Plot с фильтрацией None/NaN значений и раздельной настройкой цветов линий
    Возвращает DataFrame с значимыми метаболитами
    """
    if not selected_drugs:
        st.warning("Выберите препараты для Volcano Plot.")
        return None
    
    # Фильтрация данных
    volcano_data = []
    for drug in selected_drugs:
        drug_data = fc_df[(fc_df['Drug'] == drug) & (fc_df['Group'] == 'test')]
        if not drug_data.empty:
            max_conc = drug_data['Concentration'].max()
            max_conc_data = drug_data[drug_data['Concentration'] == max_conc].copy()
            volcano_data.append(max_conc_data)
    
    if not volcano_data:
        st.error("Нет данных для Volcano Plot.")
        return None
    
    volcano_df = pd.concat(volcano_data)
    
    # Подготовка данных с фильтрацией None/NaN
    metabolite_cols = [col for col in volcano_df.columns 
                     if col not in ['Drug', 'Group', 'Concentration'] 
                     and not col.endswith('(pvalue)')]
    
    long_data = []
    significant_metabolites = []  # Для хранения значимых метаболитов
    
    for _, row in volcano_df.iterrows():
        for metabolite in metabolite_cols:
            log2fc = row[metabolite]
            p_value = row.get(f'{metabolite} (pvalue)', None)
            
            # Пропускаем точки с отсутствующими значениями
            if pd.isna(log2fc) or pd.isna(p_value) or p_value is None or log2fc is None:
                continue
                
            # Пропускаем нулевые p-value (чтобы избежать деления на ноль в логарифме)
            if p_value <= 0:
                continue
                
            # Проверяем, является ли метаболит значимым
            is_significant = (abs(log2fc) >= log2fc_threshold) and (p_value <= p_value_threshold)
            
            long_data.append({
                'Drug': row['Drug'],
                'Concentration': row['Concentration'],
                'Metabolite': metabolite,
                'log2FC': log2fc,
                'p_value': p_value,
                '-log10(p_value)': -np.log10(p_value),
                'Significant': is_significant
            })
            
            # Добавляем в список значимых метаболитов
            if is_significant:
                significant_metabolites.append({
                    'Drug': row['Drug'],
                    'Concentration': row['Concentration'],
                    'Metabolite': metabolite,
                    'log2FC': log2fc,
                    'p_value': p_value,
                    'Fold Change (ratio)': 2**log2fc if log2fc is not None else None,
                    'Change Direction': 'Up' if log2fc > 0 else 'Down'
                })
    
    if not long_data:
        st.error("Нет данных после фильтрации None/NaN значений.")
        return None
    
    long_df = pd.DataFrame(long_data)
    significant_df = pd.DataFrame(significant_metabolites)
    
    # Настройка цветов
    color_discrete_map = None
    if custom_colors:
        color_discrete_map = {drug: color for drug, color in custom_colors.items()}
    
    # Рассчитываем позиции для делений на оси Y
    p_value_ticks = {
        0.001: -np.log10(0.001),  # 3.0
        0.01: -np.log10(0.01),    # 2.0
        0.05: -np.log10(0.05),     # ~1.3
        0.1: -np.log10(0.1)        # 1.0
    }
    
    # Создание графика
    fig = px.scatter(
        long_df,
        x='log2FC',
        y='-log10(p_value)',
        color='Drug',
        color_discrete_map=color_discrete_map,
        hover_data=['Metabolite', 'Concentration', 'p_value', 'Significant'],
        labels={
            'log2FC': 'log₂(Fold Change)',
            '-log10(p_value)': '-log₁₀(p-value)'
        },
        height=600
    )
    
    # Настройка осей с указанием порогов
    fig.update_xaxes(
        tickvals=[-3, -2, -log2fc_threshold, -1, 0, 1, log2fc_threshold, 2, 3],
        ticktext=[
            '-3', '-2', f'-{log2fc_threshold}', '-1', '0', 
            '1', f'{log2fc_threshold}', '2', '3'
        ],
        tickfont=dict(color='black'),
        title_font=dict(color='black')
    )
    
    # Интеллектуальное размещение делений для оси Y
    y_tickvals = [0, 1, 2, 3, 4, 5]
    y_ticktext = ['0', '1', '2', '3', '4', '5']
    
    # Добавляем точное значение для выбранного p-value threshold
    if p_value_threshold in p_value_ticks:
        threshold_pos = p_value_ticks[p_value_threshold]
        if threshold_pos not in y_tickvals:
            y_tickvals.append(threshold_pos)
            y_ticktext.append(f"{round(threshold_pos, 2)}")
    
    fig.update_yaxes(
        tickvals=sorted(y_tickvals),
        ticktext=[y_ticktext[y_tickvals.index(x)] for x in sorted(y_tickvals)],
        tickfont=dict(color='black'),
        title_font=dict(color='black')
    )
    
    # Горизонтальная линия (p-value threshold)
    fig.add_shape(
        type='line',
        x0=long_df['log2FC'].min() - 0.5,  # Динамический минимум по X
        x1=long_df['log2FC'].max() + 0.5,  # Динамический максимум по X
        y0=-np.log10(p_value_threshold),
        y1=-np.log10(p_value_threshold),
        line=dict(color=hline_color, dash='dash', width=2),
        name=f'p-value threshold ({p_value_threshold})'
    )
    
    # Вертикальные линии (log2FC thresholds)
    fig.add_shape(
        type='line',
        x0=-log2fc_threshold,
        x1=-log2fc_threshold,
        y0=0,
        y1=long_df['-log10(p_value)'].max() + 1,
        line=dict(color=vline_color, dash='dot', width=2),
        name=f'log2FC threshold (-{log2fc_threshold})'
    )
    
    fig.add_shape(
        type='line',
        x0=log2fc_threshold,
        x1=log2fc_threshold,
        y0=0,
        y1=long_df['-log10(p_value)'].max() + 1,
        line=dict(color=vline_color, dash='dot', width=2),
        name=f'log2FC threshold ({log2fc_threshold})'
    )
    
    # Настройка легенды (без рамки)
    fig.update_layout(
        font=dict(color='black'),
        showlegend=show_legend,
        legend=dict(
            font=dict(color='black'),
            title_font=dict(color='black'),
            title_text='Препараты',
            bgcolor='rgba(255,255,255,0.7)',
            bordercolor='rgba(0,0,0,0)',
            borderwidth=0
        ),
        margin=dict(
        l=100,   # отступ слева
        r=100,   # отступ справа
        t=100,  # отступ сверху (под заголовок)
        b=80,   # отступ снизу (для подписи оси X)
        pad=20  # дополнительный «буфер»
        ),
        title={
            'text': f"Volcano Plot (p < {p_value_threshold}, |log2FC| > {log2fc_threshold})",
            'y':0.95,
            'x':0.475,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {'size': 16}
        },

    )
    
    # Конфигурация для скачивания в высоком качестве
    config = {
        'toImageButtonOptions': {
            'format': 'jpeg',
            'filename': 'volcano_plot',
            'scale': 4,
            'dpi': 600
        }
    }
    
    st.plotly_chart(fig, use_container_width=True, config=config)
    
    return significant_df


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

    # Инициализация session state для fc_mode, если его еще нет
    if "fc_mode" not in st.session_state:
        st.session_state["fc_mode"] = 'log₂(B/A)'  # Значение по умолчанию
     
    
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
                
                list_measure_unit_concentration = []
                for drug in df['Drug'].unique().tolist():
                     measure_unit_concentration = select_concentration_unit(f"select_concentration_unit{drug}",f"key_select_concentration_unit{drug}",f"{drug}")
                     list_measure_unit_concentration.append(measure_unit_concentration)
                
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
                
                st.subheader(f"Результаты расчета Fold Change ({st.session_state['fc_mode']})")
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
                    file_name=f"fold_change_results_{st.session_state['fc_mode']}.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )

                # Добавляем фильтрацию по log2FC только для режима log2

                fc_df_for_diagram = fc_df.copy()

                if st.session_state["fc_mode"] == 'log₂(B/A)':
                    log2fc_threshold = st.radio(
                        "Фильтрация по log2FC для графика столбчатых диаграмм:",
                        options=[None, 1.0, 0.58],
                        index=0,  # None по умолчанию (без фильтрации)
                        format_func=lambda x: "Без фильтрации" if x is None else f"{x} ({'2-кратное' if x == 1.0 else '1.5-кратное'} изменение)",
                        horizontal=True,
                        key="horizontal_log2fc_threshold"
                    )
                else:
                    log2fc_threshold = None


                if log2fc_threshold is not None:
                        # Определяем колонки, которые НЕ нужно фильтровать
                        exclude_cols = ['Drug', 'Group', 'Concentration']
                        
                        # Получаем список колонок для фильтрации
                        cols_to_filter = [col for col in fc_df_for_diagram.columns if col not in exclude_cols]
                        
                        # Применяем фильтрацию: заменяем значения, которые по модулю меньше порога, на None
                        fc_df_for_diagram[cols_to_filter] = fc_df_for_diagram[cols_to_filter].apply(
                            lambda x: x.where(abs(x) >= log2fc_threshold)
                        )

                        # Удаляем строки, где все значения (кроме exclude_cols) равны None
                        # Создаем маску: True для строк, которые нужно сохранить (хотя бы одно не-None значение в cols_to_filter)
                        mask = fc_df_for_diagram[cols_to_filter].notna().any(axis=1)
                        fc_df_for_diagram = fc_df_for_diagram[mask]
                        
                        # Удаляем столбцы, которые полностью состоят из None (после фильтрации)
                        fc_df_for_diagram.dropna(axis=1, how='all', inplace=True)
                
                if st.session_state["fc_mode"] == 'log₂(B/A)':
                    st.subheader("Данные после фильтрации по log2FC")

                    # Шаблон для поиска "(p-value)" в имени колонки
                    pval_pattern = re.compile(r'\(pvalue\)\s*$', flags=re.IGNORECASE)

                    # список колонок, которые надо удалить
                    cols_to_drop = [col for col in fc_df_for_diagram.columns if pval_pattern.search(col)]

                    # удаляем их
                    fc_df_for_diagram = fc_df_for_diagram.drop(columns=cols_to_drop)

                    # Отбираем только группу test
                    df_test = fc_df_for_diagram[fc_df_for_diagram['Group'] == 'test']

                    df_test.drop('Group', axis=1, inplace=True)
                    
                    coords = (
                        df_test
                        .set_index(['Drug', 'Concentration'])    # делаем составной индекс
                        .stack()                                 # «разворачиваем» столбцы метаболитов
                        .reset_index(name='Value')               # приводим обратно в датафрейм
                        .rename(columns={'level_2': 'Metabolite'})
                    )

                    # оставляем только нужные колонки
                    df_result_for_download = coords[['Drug', 'Concentration', 'Metabolite']]

                    # Для ускорения: сделаем индекс по Drug и Concentration
                    fc_indexed = fc_df.set_index(['Drug', 'Concentration'])

                    def fetch_values(row):
                        key = (row['Drug'], row['Concentration'])
                        met = row['Metabolite']
                        try:
                            log2fc = fc_indexed.loc[key, met]
                        except KeyError:
                            log2fc = np.nan
                        try:
                            pval = fc_indexed.loc[key, f"{met} (pvalue)"]
                        except KeyError:
                            pval = np.nan
                        return pd.Series({'log2FC': log2fc, 'pvalue': pval})

                    # Применяем к каждой строке df_result
                    df_result_for_download[['log2FC', 'pvalue']] = df_result_for_download.apply(fetch_values, axis=1)

                    # 1) Fold Change (ratio) из log2FC
                    #    просто обратно: ratio = 2**log2FC
                    df_result_for_download['Fold Change (ratio)'] = 2 ** df_result_for_download['log2FC']

                    # 2) Change Direction: «Up» если log2FC > 0, «Down» если < 0, иначе «No change»
                    df_result_for_download['Change Direction'] = np.where(
                        df_result_for_download['log2FC'] > 0, 'Up',
                        np.where(df_result_for_download['log2FC'] < 0, 'Down', 'No change')
                    )

                    num_unique = df_result_for_download['Metabolite'].nunique()


                    st.dataframe(df_result_for_download)
                    # Кнопка скачивания значимых метаболитов
                    output_df_result_for_download = io.BytesIO()
                    with pd.ExcelWriter(output_df_result_for_download, engine='openpyxl') as writer:
                        df_result_for_download.to_excel(writer, index=False, sheet_name=f'Данные после фильтрации по log2FC {log2fc_threshold}')
                    output_df_result_for_download.seek(0)
                    
                    st.download_button(
                        label="Скачать таблицу данных после фильтрации по log2FC",
                        data=output_df_result_for_download,
                        file_name=f"Данные после фильтрации по log2FC {log2fc_threshold}.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )           
                
                # Графики Fold Change
                available_drugs = fc_df_for_diagram['Drug'].unique()

                # Получаем список метаболитов, исключая p-value и служебные колонки
                available_metabolites = [
                    col for col in fc_df_for_diagram.columns 
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
                    
                    # Создаем словарь для хранения выбранных концентраций для каждого препарата
                    selected_concentrations = {}
                    
                    # Для каждого выбранного препарата создаем мультиселектор концентраций
                    for drug in selected_drugs:
                        # Получаем уникальные концентрации для данного препарата, исключая нулевые
                        drug_concentrations = [
                            conc for conc in fc_df_for_diagram[fc_df_for_diagram['Drug'] == drug]['Concentration'].unique() 
                            if conc != 0  # Исключаем нулевые концентрации
                        ]

                        # Проверяем, есть ли ненулевые концентрации
                        if drug_concentrations:
                            selected_concentrations[drug] = st.multiselect(
                                f"Концентрации для {drug}",
                                options=drug_concentrations,
                                default=drug_concentrations,
                                key=f"conc_{drug}"
                            )
                        else:
                            st.warning(f"Для препарата {drug} нет доступных ненулевых концентраций")
                            selected_concentrations[drug] = []
                    
                    selected_metabolites = st.multiselect(
                        "Выберите метаболиты", 
                        available_metabolites, 
                        default=available_metabolites, 
                        key="met_select"
                    )
                    
                    # Добавляем настройки отображения
                    with st.expander("Настройки отображения"):
                        col1, col2 = st.columns(2)
                        with col1:
                            legend_position = st.radio(
                                "Положение легенды",
                                ["Снизу", "Справа"],
                                index=0,
                                horizontal=True
                            )
                            show_legend = st.checkbox("Показывать легенду", value=True)
                        with col2:
  
                            opacity = st.slider("Прозрачность колонок", 0.1, 1.0, 0.8)

                            
                        with st.container(border=True):
                            # Добавляем возможность задать цвета для каждого столбца
                            if selected_drugs and selected_metabolites:
                                unique_combinations = fc_df_for_diagram[
                                    (fc_df_for_diagram['Drug'].isin(selected_drugs)) & 
                                    (fc_df_for_diagram['Concentration'].isin([conc for drug in selected_drugs for conc in selected_concentrations[drug]]))
                                ].groupby(['Drug', 'Concentration']).size().reset_index()
                                
                                custom_colors = []
                                st.write("Настройте цвета для каждого столбца:")

                                num_cols = 4  # Количество колонок
                                cols = st.columns(num_cols)

                                custom_colors = []
                                for idx, (_, row) in enumerate(unique_combinations.iterrows()):
                                    with cols[idx % num_cols]:
                                        color = st.color_picker(
                                            f"Цвет для {row['Drug']} ({row['Concentration']})",
                                            value=px.colors.qualitative.Plotly[idx % len(px.colors.qualitative.Plotly)],
                                            key=f"color_{row['Drug']}_{row['Concentration']}"
                                        )
                                        custom_colors.append(color)
                            else:
                                custom_colors = None

                            apply_colors = st.checkbox("Применить цвета")

                        with col2:    
                            baseline_color = st.color_picker("Цвет контрольной линии", "#000000")
                    
                    submitted = st.form_submit_button("Перерисовать график")
                    
                    if submitted:
                        st.session_state['show_plot'] = True

                if st.session_state['show_plot']:
                    # Фильтруем данные по выбранным концентрациям
                    filtered_fc_df = fc_df_for_diagram[
                        (fc_df_for_diagram['Drug'].isin(selected_drugs)) & 
                        (fc_df_for_diagram['Concentration'].isin([conc for drug in selected_drugs for conc in selected_concentrations[drug]]))
                    ]
                
                    
                    # В вызове функции добавьте параметр custom_colors:
                    plot_fold_change_horizontal(
                        filtered_fc_df, 
                        selected_drugs, 
                        selected_metabolites, 
                        list_measure_unit_concentration,
                        apply_colors,
                        legend_position=legend_position,
                        show_legend=show_legend,
                        opacity=opacity,
                        baseline_color=baseline_color,
                        custom_colors=custom_colors  # Добавляем пользовательские цвета
                    )

                    # Добавляем Volcano Plot только для режима log₂(B/A)
                    if st.session_state["fc_mode"] == 'log₂(B/A)':
                        st.subheader("Volcano Plot (для максимальной концентрации)")
                        st.write("""
                        **Интерпретация Volcano Plot:**
                        - Точки в верхних правом/левом углах — значимые изменения (большой |log2FC| и низкий p-value)
                        - Горизонтальная линия — порог значимости (p < выбранное значение)
                        - Вертикальные линии — порог изменения (|log2FC| > выбранное значение)
                        - Доступные пороги:
                        - log2FC: 0.58 (1.5x), 1.0 (2x)
                        - p-value: 0.001 (0.1%), 0.01 (1%), 0.05 (5%), 0.1 (10%)
                        """)

                        # Графики Vulcano
                        available_drugs = fc_df['Drug'].unique()


                        if 'show_vulcano_plot' not in st.session_state:
                            st.session_state['show_vulcano_plot'] = False
                        
                        # Создаем отдельную форму для Volcano Plot
                        with st.form("volcano_form"):
                            # Независимый выбор препаратов
                            volcano_drugs = st.multiselect(
                                "Выберите препараты для Volcano Plot",
                                available_drugs,
                                default=available_drugs[:min(5, len(available_drugs))],  # Первые 5 по умолчанию
                                key="volcano_drugs"
                            )
                            
                            # Настройки отображения
                            with st.expander("Настройки Volcano Plot"):
                                col1, col2 = st.columns(2)
                                with col1:
                                    p_value_threshold = st.selectbox(
                                        "Порог p-value (статистическая значимость)",
                                        options=[0.001, 0.01, 0.05, 0.1],
                                        index=2,  # 0.05 по умолчанию
                                        format_func=lambda x: f"{x} ({'***' if x == 0.001 else '**' if x == 0.01 else '*' if x == 0.05 else 'ns'})",
                                        help="Уровни значимости: *** - 0.1%, ** - 1%, * - 5%, ns - не значимо"
                                    )
                                    
                                    log2fc_threshold = st.radio(
                                        "Порог log2FC (кратность изменения)",
                                        options=[1.0, 0.58],
                                        index=0,  # 1.0 по умолчанию
                                        format_func=lambda x: f"{x} ({'2-кратное' if x == 1.0 else '1.5-кратное'} изменение)",
                                        horizontal=True
                                    )
                                
                                with col2:
                                    # Настройки отображения
                                    show_legend = st.checkbox(
                                        "Показывать легенду", 
                                        value=True,
                                        key='volcano_show_legend'
                                    )
                                    
                                    # Раздельные настройки цветов линий
                                    hline_color = st.color_picker(
                                        "Цвет горизонтальной линии (p-value)",
                                        value='#FF0000',
                                        key='volcano_hline_color'
                                    )
                                    
                                    vline_color = st.color_picker(
                                        "Цвет вертикальных линий (log2FC)",
                                        value='#808080',
                                        key='volcano_vline_color'
                                    )
                                    
                                    # Настройки цветов препаратов
                                    volcano_colors = {}
                                    if volcano_drugs:
                                        st.write("Настройте цвета для препаратов:")
                                        cols = st.columns(4)
                                        for idx, drug in enumerate(volcano_drugs):
                                            with cols[idx % 4]:
                                                volcano_colors[drug] = st.color_picker(
                                                    f"Цвет для {drug}",
                                                    value=px.colors.qualitative.Plotly[idx % len(px.colors.qualitative.Plotly)],
                                                    key=f"volcano_color_{drug}"
                                                )
                            
                            submitted_volcano = st.form_submit_button("Построить Volcano Plot")

                            if submitted_volcano:
                               st.session_state['show_vulcano_plot'] = True

                        if st.session_state['show_vulcano_plot'] and volcano_drugs:
                            significant_df = plot_volcano(
                                fc_df, 
                                volcano_drugs, 
                                p_value_threshold=p_value_threshold,
                                log2fc_threshold=log2fc_threshold,
                                custom_colors=volcano_colors,
                                show_legend=show_legend,
                                hline_color=hline_color,
                                vline_color=vline_color
                            )
                            
                            if significant_df is not None and not significant_df.empty:
                                st.subheader("Значимые метаболиты")
                                st.write(f"Найдено {len(significant_df)} значимых метаболитов (p < {p_value_threshold}, |log2FC| > {log2fc_threshold})")
                                st.dataframe(significant_df)
                                
                                # Кнопка скачивания значимых метаболитов
                                output_significant = io.BytesIO()
                                with pd.ExcelWriter(output_significant, engine='openpyxl') as writer:
                                    significant_df.to_excel(writer, index=False, sheet_name='Significant_Metabolites')
                                output_significant.seek(0)
                                
                                st.download_button(
                                    label="Скачать таблицу значимых метаболитов",
                                    data=output_significant,
                                    file_name=f"significant_metabolites_p{p_value_threshold}_fc{log2fc_threshold}.xlsx",
                                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                                )
                            elif significant_df is not None and significant_df.empty:
                                st.warning("Нет значимых метаболитов при заданных параметрах.")