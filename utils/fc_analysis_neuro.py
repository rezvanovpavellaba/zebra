import streamlit as st
from scipy.stats import ttest_ind, mannwhitneyu, shapiro
import plotly.express as px
import pandas as pd
import numpy as np

def plot_volcano(fc_df, selected_drugs,selected_conc, p_value_threshold=0.05, log2fc_threshold=1.0, 
                custom_colors=None, show_legend=True, 
                hline_color='red', vline_color='gray'):
    """
    Улучшенный Volcano Plot с фильтрацией None/NaN значений и раздельной настройкой цветов линий
    Возвращает DataFrame с значимыми метаболитами
    """
    if not selected_drugs:
        st.warning("Выберите препараты для Volcano Plot.")
        return None
    
    if not selected_conc:
        st.warning("Выберите концентрацию для Volcano Plot.")
        return None
    
    # Фильтрация данных
    volcano_data = []
    for drug in selected_drugs:
        drug_data = fc_df[(fc_df['Compound'] == drug) & (fc_df['Test/control'] == 'Test')]
        if not drug_data.empty:
            current_conc = selected_conc
            current_conc_data = drug_data[drug_data['Concentration'] == current_conc].copy()
            volcano_data.append(current_conc_data)
    
    if not volcano_data:
        st.error("Нет данных для Volcano Plot.")
        return None
    
    volcano_df = pd.concat(volcano_data)
    
    # Подготовка данных с фильтрацией None/NaN
    metabolite_cols = [col for col in volcano_df.columns 
                     if col not in ['Compound', 'Test/control', 'Concentration'] 
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
                'Compound': row['Compound'],
                'Concentration': row['Concentration'],
                'Metric': metabolite,
                'log2FC': log2fc,
                'p_value': p_value,
                '-log10(p_value)': -np.log10(p_value),
                'Significant': is_significant
            })
            
            # Добавляем в список значимых метаболитов
            if is_significant:
                significant_metabolites.append({
                    'Compound': row['Compound'],
                    'Concentration': row['Concentration'],
                    'Metric': metabolite,
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
        color='Compound',
        color_discrete_map=color_discrete_map,
        hover_data=['Metric', 'Concentration', 'p_value', 'Significant'],
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
            title_text='Drugs',
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

def calculate_fold_change_with_pvalues(df, mode='ratio'):
    """
    Расчёт Fold Change и p-values с использованием t-теста Стьюдента.
    Предполагает равенство дисперсий в группах.
    """
    metabolite_cols = [col for col in df.columns 
                      if col not in ['Compound', 'Test/control', 'Concentration',"well_id"]]
    
    results = []
    warnings = []
    
    for drug in df['Compound'].unique():
        drug_data = df[df['Compound'] == drug]
        control_data = drug_data[drug_data['Test/control'] == 'Control']
        test_data = drug_data[drug_data['Test/control'] == 'Test']
        
        # Добавляем строку "контроль против контроля"
        control_row = {
            'Compound': drug,
            'Test/control': 'control_vs_control',
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
                'Compound': drug,
                'Test/control': 'Test',
                'Concentration': conc
            }
            
            for metabolite in metabolite_cols:
                control_vals = control_data[metabolite].values
                test_vals = test_subset[metabolite].values
                
                # Проверка количества данных
                if len(control_vals) < 2 or len(test_vals) < 2:  # Хотя бы 2 точки в тесте для FC
                    p_value = None
                    if len(test_vals) < 2:
                        warnings.append(f"{metabolite} (Compound: {drug}, Conc: {conc})")
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
                
                # Тест на нормальность
                try:
                    p_norm_control = shapiro(control_vals).pvalue if len(control_vals) >= 3 else 1
                    p_norm_test = shapiro(test_vals).pvalue if len(test_vals) >= 3 else 1
                except:
                    p_norm_control = p_norm_test = 1  # безопасное значение

                # Выбор теста
                try:
                    if p_norm_control < 0.05 or p_norm_test < 0.05:
                        _, p_value = mannwhitneyu(control_vals, test_vals, alternative='two-sided')
                        print("mann")
                    else:
                        _, p_value = ttest_ind(control_vals, test_vals, equal_var=False)
                        print("student")
                except:
                    p_value = None
                
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