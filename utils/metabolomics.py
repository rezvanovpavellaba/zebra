import streamlit as st
import pandas as pd
import numpy as np
import io
import plotly.express as px
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, shapiro, mannwhitneyu
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


def plot_volcano(
    fc_df,
    selected_drugs,
    p_value_threshold=0.05,
    log2fc_threshold=1.0,
    custom_colors=None,
    show_legend=True,
    hline_color='red',
    vline_color='gray',
    legend_title='Drugs',
    metabolite_label_map=None,
    drug_label_map=None
):
    """
    Улучшенный Volcano Plot с фильтрацией None/NaN значений и раздельной настройкой цветов линий
    Возвращает DataFrame с значимыми метаболитами
    """
    if not selected_drugs:
        st.warning("Выберите препараты для Volcano Plot.")
        return None
    
    # Фильтрация данных
    metabolite_label_map = metabolite_label_map or {}
    drug_label_map = drug_label_map or {}
    legend_title = legend_title.strip() if isinstance(legend_title, str) and legend_title.strip() else 'Drugs'

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
    long_df['Display Drug'] = long_df['Drug'].apply(
        lambda drug: drug_label_map.get(drug, drug)
    )
    long_df['Display Metabolite'] = long_df['Metabolite'].apply(
        lambda metabolite: metabolite_label_map.get(metabolite, metabolite)
    )
    if not significant_df.empty and metabolite_label_map:
        significant_df['Display Metabolite'] = significant_df['Metabolite'].apply(
            lambda metabolite: metabolite_label_map.get(metabolite, metabolite)
        )
    if not significant_df.empty and drug_label_map:
        significant_df['Display Drug'] = significant_df['Drug'].apply(
            lambda drug: drug_label_map.get(drug, drug)
        )
    
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
        custom_data=['Display Drug', 'Display Metabolite', 'Concentration', 'p_value', 'Significant'],
        labels={
            'Display Metabolite': 'Metabolite',
            'log2FC': 'log₂(Fold Change)',
            '-log10(p_value)': '-log₁₀(p-value)'
        },
        height=600
    )

    fig.update_traces(
        hovertemplate=(
            f"{legend_title}=%{{customdata[0]}}<br>"
            "Metabolite=%{customdata[1]}<br>"
            "Concentration=%{customdata[2]}<br>"
            "p_value=%{customdata[3]:.4g}<br>"
            "Significant=%{customdata[4]}<br>"
            "log2FC=%{x}<br>"
            "-log10(p-value)=%{y}<extra></extra>"
        ),
        selector=dict(mode='markers')
    )
    for trace in fig.data:
        if getattr(trace, 'mode', None) == 'markers':
            display_name = drug_label_map.get(trace.name, trace.name)
            trace.name = display_name
            trace.legendgroup = display_name

    labeled_points = long_df[long_df['Significant']].copy()
    if not labeled_points.empty:
        labeled_points['LabelPosition'] = np.where(
            labeled_points['log2FC'] >= 0,
            'top right',
            'top left'
        )
        fig.add_scatter(
            x=labeled_points['log2FC'],
            y=labeled_points['-log10(p_value)'],
            mode='text',
            text=labeled_points['Display Metabolite'],
            textposition=labeled_points['LabelPosition'],
            textfont=dict(color='black', size=11),
            showlegend=False,
            hoverinfo='skip',
            cliponaxis=False
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
            title_text=legend_title,
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


def _format_numeric_label(value):
    """Компактное текстовое представление чисел для подписей."""
    if pd.isna(value):
        return ""
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        if float(value).is_integer():
            return str(int(value))
        return f"{float(value):g}"
    return str(value)


def _format_template_label(template, context, fallback=""):
    """Безопасно форматирует шаблоны подписей вида '{drug} | {metabolite}'."""
    if not isinstance(template, str):
        return fallback

    template = template.strip()
    if not template:
        return fallback

    try:
        return template.format(**context)
    except (KeyError, IndexError, ValueError):
        return template


def _cardio_like_color_cycle():
    """Возвращает цветовой цикл, визуально совпадающий с кардио-боксплотами."""
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]
    return cycle


def _format_metabolomics_significance_label(p_value, selected_threshold=0.05):
    """
    Подпись значимости для boxplot.
    Для p-value между 0.05 и выбранным порогом показываем точное p,
    чтобы сохранить логику текущего метаболомического пайплайна.
    """
    if p_value is None or pd.isna(p_value):
        return None

    p_value = float(p_value)
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    if selected_threshold is not None and p_value <= selected_threshold:
        return f"p={p_value:.2g}"
    return "ns"


def _metabolomics_test_group_key(concentration):
    """Стабильный ключ тестовой группы для конкретной концентрации."""
    return f"test::{repr(concentration)}"


def prepare_metabolomics_boxplot_data(
    df,
    drug,
    concentration,
    metabolite,
    control_tick="Control",
    test_tick=None
):
    """
    Готовит raw-данные для boxplot по одной паре Drug/Metabolite на максимальной концентрации.
    Контроль берётся только у того же препарата, что и в основном статистическом пайплайне.
    """
    control_values = pd.to_numeric(
        df[
            (df["Drug"] == drug) &
            (df["Group"] == "control_neg")
        ][metabolite],
        errors="coerce"
    ).dropna()

    concentrations = concentration if isinstance(concentration, (list, tuple, set, np.ndarray, pd.Series)) else [concentration]
    concentrations = [conc for conc in concentrations if pd.notna(conc)]

    rows = []
    rows.extend(
        {
            "GroupKey": "control",
            "TickLabel": control_tick,
            "GroupOrder": 0,
            "Value": value
        }
        for value in control_values
    )

    for idx, conc in enumerate(concentrations, start=1):
        current_test_tick = test_tick(conc) if callable(test_tick) else test_tick
        test_values = pd.to_numeric(
            df[
                (df["Drug"] == drug) &
                (df["Group"] == "test") &
                (df["Concentration"] == conc)
            ][metabolite],
            errors="coerce"
        ).dropna()

        rows.extend(
            {
                "GroupKey": _metabolomics_test_group_key(conc),
                "TickLabel": current_test_tick,
                "GroupOrder": idx,
                "Value": value
            }
            for value in test_values
        )

    return pd.DataFrame(rows)


def build_metabolomics_boxplot(
    plot_df,
    title=None,
    ylabel=None,
    xlabel=None,
    show_points=True,
    points_jitter=0.08,
    points_alpha=0.6,
    figsize=(7.0, 4.5),
    signif_by_group=None,
    signif_fontsize=12,
    signif_y_pad_frac=0.06,
):
    """
    Boxplot в визуальном стиле кардиотоксичности:
    чёрные границы/медианы, разные заливки боксов из цикла matplotlib.
    """
    if plot_df is None or plot_df.empty:
        return None

    group_columns = ["GroupKey", "TickLabel"]
    if "GroupOrder" in plot_df.columns:
        group_columns.append("GroupOrder")
        ordered_groups = (
            plot_df[group_columns]
            .drop_duplicates()
            .sort_values(["GroupOrder", "GroupKey"])
        )
    else:
        ordered_groups = (
            plot_df[group_columns]
            .drop_duplicates()
            .sort_values("GroupKey", key=lambda x: x.map({"control": 0, "test": 1}).fillna(99))
        )
    groups = ordered_groups["GroupKey"].tolist()
    tick_map = dict(zip(ordered_groups["GroupKey"], ordered_groups["TickLabel"]))

    data = [
        plot_df.loc[plot_df["GroupKey"] == group_key, "Value"].to_numpy()
        for group_key in groups
    ]
    tick_labels = [tick_map.get(group_key, group_key) for group_key in groups]

    fig, ax = plt.subplots(figsize=figsize)

    bp = ax.boxplot(
        data,
        labels=tick_labels,
        patch_artist=True,
        widths=0.55,
        showfliers=False,
        medianprops=dict(color="black", linewidth=2.0),
        boxprops=dict(edgecolor="black", linewidth=1.8),
        whiskerprops=dict(color="black", linewidth=1.6),
        capprops=dict(color="black", linewidth=1.6),
    )

    cycle = _cardio_like_color_cycle()
    for idx, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(cycle[idx % len(cycle)])
        patch.set_alpha(0.55)
        patch.set_edgecolor("black")
        patch.set_linewidth(1.8)

    if show_points:
        rng = np.random.default_rng(0)
        for idx, values in enumerate(data, start=1):
            if len(values) == 0:
                continue
            x = rng.normal(loc=idx, scale=points_jitter, size=len(values))
            ax.scatter(x, values, s=28, alpha=points_alpha, edgecolors="none")

    if signif_by_group:
        all_values = np.concatenate([values for values in data if len(values) > 0]) if any(len(values) > 0 for values in data) else np.array([0.0])
        y_max = float(np.nanmax(all_values))
        y_min = float(np.nanmin(all_values))
        y_range = max(1e-9, y_max - y_min)
        y_text = y_max + signif_y_pad_frac * y_range

        ax.set_ylim(top=y_text + 0.10 * y_range)

        for idx, group_key in enumerate(groups, start=1):
            signif_label = signif_by_group.get(group_key)
            if not signif_label:
                continue
            ax.text(
                idx,
                y_text,
                signif_label,
                ha="center",
                va="bottom",
                fontsize=signif_fontsize,
                fontweight="bold",
            )

    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.tick_params(axis="x", rotation=0)

    fig.tight_layout()
    return fig


def fig_to_png_bytes(fig):
    """Подготавливает matplotlib figure к скачиванию в PNG."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=600, bbox_inches="tight")
    buf.seek(0)
    return buf.getvalue()


def render_metabolomics_fc_barplot_section(fc_df, list_measure_unit_concentration):
    """Рендер вкладки с диаграммой Fold Change и связанными только с ней настройками."""
    fc_df_for_diagram = fc_df.copy()

    if st.session_state["fc_mode"] == 'log₂(B/A)':
        log2fc_threshold = st.radio(
            "Фильтрация по log2FC для графика столбчатых диаграмм:",
            options=[None, 1.0, 0.58],
            index=0,
            format_func=lambda x: "Без фильтрации" if x is None else f"{x} ({'2-кратное' if x == 1.0 else '1.5-кратное'} изменение)",
            horizontal=True,
            key="horizontal_log2fc_threshold"
        )
    else:
        log2fc_threshold = None

    if log2fc_threshold is not None:
        exclude_cols = ['Drug', 'Group', 'Concentration']
        cols_to_filter = [col for col in fc_df_for_diagram.columns if col not in exclude_cols]
        fc_df_for_diagram[cols_to_filter] = fc_df_for_diagram[cols_to_filter].apply(
            lambda x: x.where(abs(x) >= log2fc_threshold)
        )
        mask = fc_df_for_diagram[cols_to_filter].notna().any(axis=1)
        fc_df_for_diagram = fc_df_for_diagram[mask]
        fc_df_for_diagram.dropna(axis=1, how='all', inplace=True)

    if st.session_state["fc_mode"] == 'log₂(B/A)':
        st.subheader("Данные после фильтрации по log2FC")
        preview_df = fc_df_for_diagram.copy()

        pval_pattern = re.compile(r'\(pvalue\)\s*$', flags=re.IGNORECASE)
        cols_to_drop = [col for col in preview_df.columns if pval_pattern.search(col)]
        preview_df = preview_df.drop(columns=cols_to_drop)

        df_test = preview_df[preview_df['Group'] == 'test'].copy()
        if 'Group' in df_test.columns:
            df_test.drop('Group', axis=1, inplace=True)

        coords = (
            df_test
            .set_index(['Drug', 'Concentration'])
            .stack()
            .reset_index(name='Value')
            .rename(columns={'level_2': 'Metabolite'})
        ) if not df_test.empty else pd.DataFrame(columns=['Drug', 'Concentration', 'Metabolite', 'Value'])

        df_result_for_download = coords[['Drug', 'Concentration', 'Metabolite']].copy() if not coords.empty else pd.DataFrame(columns=['Drug', 'Concentration', 'Metabolite'])
        fc_indexed = fc_df.set_index(['Drug', 'Concentration'])

        if not df_result_for_download.empty:
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

            df_result_for_download[['log2FC', 'pvalue']] = df_result_for_download.apply(fetch_values, axis=1)
            df_result_for_download['Fold Change (ratio)'] = 2 ** df_result_for_download['log2FC']
            df_result_for_download['Change Direction'] = np.where(
                df_result_for_download['log2FC'] > 0, 'Up',
                np.where(df_result_for_download['log2FC'] < 0, 'Down', 'No change')
            )
        else:
            df_result_for_download = pd.DataFrame(
                columns=['Drug', 'Concentration', 'Metabolite', 'log2FC', 'pvalue', 'Fold Change (ratio)', 'Change Direction']
            )

        st.dataframe(df_result_for_download)

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

    available_drugs = fc_df_for_diagram['Drug'].unique()
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

        selected_concentrations = {}
        for drug in selected_drugs:
            drug_concentrations = [
                conc for conc in fc_df_for_diagram[fc_df_for_diagram['Drug'] == drug]['Concentration'].unique()
                if conc != 0
            ]

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
                if selected_drugs and selected_metabolites:
                    unique_combinations = fc_df_for_diagram[
                        (fc_df_for_diagram['Drug'].isin(selected_drugs)) &
                        (fc_df_for_diagram['Concentration'].isin([conc for drug in selected_drugs for conc in selected_concentrations[drug]]))
                    ].groupby(['Drug', 'Concentration']).size().reset_index()

                    st.write("Настройте цвета для каждого столбца:")
                    num_cols = 4
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
        filtered_fc_df = fc_df_for_diagram[
            (fc_df_for_diagram['Drug'].isin(selected_drugs)) &
            (fc_df_for_diagram['Concentration'].isin([conc for drug in selected_drugs for conc in selected_concentrations[drug]]))
        ]

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
            custom_colors=custom_colors
        )
    else:
        st.info("Настройте параметры и нажмите «Перерисовать график», чтобы построить диаграмму.")


def render_metabolomics_volcano_section(fc_df, available_metabolites):
    """Отдельный рендер Volcano Plot, независимый от barplot-секции."""
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

    available_drugs = fc_df['Drug'].unique()

    if 'show_vulcano_plot' not in st.session_state:
        st.session_state['show_vulcano_plot'] = False

    with st.form("volcano_form"):
        volcano_drugs = st.multiselect(
            "Выберите препараты для Volcano Plot",
            available_drugs,
            default=available_drugs[:min(5, len(available_drugs))],
            key="volcano_drugs"
        )

        with st.expander("Настройки Volcano Plot"):
            col1, col2 = st.columns(2)
            with col1:
                p_value_threshold = st.selectbox(
                    "Порог p-value (статистическая значимость)",
                    options=[0.001, 0.01, 0.05, 0.1],
                    index=2,
                    format_func=lambda x: f"{x} ({'***' if x == 0.001 else '**' if x == 0.01 else '*' if x == 0.05 else 'ns'})",
                    help="Уровни значимости: *** - 0.1%, ** - 1%, * - 5%, ns - не значимо"
                )

                log2fc_threshold = st.radio(
                    "Порог log2FC (кратность изменения)",
                    options=[1.0, 0.58],
                    index=0,
                    format_func=lambda x: f"{x} ({'2-кратное' if x == 1.0 else '1.5-кратное'} изменение)",
                    horizontal=True
                )

            with col2:
                show_legend = st.checkbox(
                    "Показывать легенду",
                    value=True,
                    key='volcano_show_legend'
                )
                legend_title = st.text_input(
                    "Заголовок легенды",
                    value="Drugs",
                    key='volcano_legend_title'
                )

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

            st.write("Подписи препаратов на графике:")
            st.caption("Изменённые подписи будут использоваться в легенде Volcano Plot и в таблице значимых метаболитов.")
            drug_labels_df = pd.DataFrame({
                'Drug': volcano_drugs,
                'Display Name': volcano_drugs
            })
            edited_drug_labels = st.data_editor(
                drug_labels_df,
                hide_index=True,
                use_container_width=True,
                disabled=['Drug'],
                column_config={
                    'Drug': st.column_config.TextColumn('Препарат'),
                    'Display Name': st.column_config.TextColumn('Подпись на графике')
                },
                key='volcano_drug_labels_editor'
            )
            drug_label_map = {
                row['Drug']: row['Display Name'].strip()
                for _, row in edited_drug_labels.iterrows()
                if isinstance(row['Display Name'], str)
                and row['Display Name'].strip()
                and row['Display Name'].strip() != row['Drug']
            }

            st.write("Подписи метаболитов на графике:")
            st.caption("Изменённые подписи будут использоваться для значимых точек Volcano Plot и в таблице значимых метаболитов.")
            metabolite_labels_df = pd.DataFrame({
                'Metabolite': available_metabolites,
                'Display Name': available_metabolites
            })
            edited_metabolite_labels = st.data_editor(
                metabolite_labels_df,
                hide_index=True,
                use_container_width=True,
                disabled=['Metabolite'],
                column_config={
                    'Metabolite': st.column_config.TextColumn('Метаболит'),
                    'Display Name': st.column_config.TextColumn('Подпись на графике')
                },
                key='volcano_metabolite_labels_editor'
            )
            metabolite_label_map = {
                row['Metabolite']: row['Display Name'].strip()
                for _, row in edited_metabolite_labels.iterrows()
                if isinstance(row['Display Name'], str)
                and row['Display Name'].strip()
                and row['Display Name'].strip() != row['Metabolite']
            }

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
            vline_color=vline_color,
            legend_title=legend_title,
            metabolite_label_map=metabolite_label_map,
            drug_label_map=drug_label_map
        )

        if significant_df is not None and not significant_df.empty:
            st.subheader("Значимые метаболиты")
            st.write(f"Найдено {len(significant_df)} значимых метаболитов (p < {p_value_threshold}, |log2FC| > {log2fc_threshold})")
            st.dataframe(significant_df)

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


def render_metabolomics_boxplot_section(df, fc_df, measure_unit_by_drug, available_metabolites):
    """Отдельный boxplot-раздел по всем метаболитам с использованием полного fc_df."""
    st.subheader("📊 Boxplots")
    st.caption("Секция работает независимо от Volcano Plot и использует полный набор метаболитов из расчёта p-value.")

    test_fc_df = fc_df[fc_df['Group'] == 'test'].copy()
    if test_fc_df.empty:
        st.info("Нет test-данных для построения boxplot.")
        return

    if 'metabolomics_boxplot_request' not in st.session_state:
        st.session_state['metabolomics_boxplot_request'] = None

    drug_options = sorted(test_fc_df['Drug'].dropna().unique().tolist())

    selected_boxplot_drug = st.selectbox(
        "Выберите препарат",
        options=[""] + drug_options,
        index=0,
        format_func=lambda value: "—" if value == "" else value,
        key="metabolomics_boxplot_drug"
    )

    concentration_options = []
    if selected_boxplot_drug:
        concentration_options = sorted(
            test_fc_df.loc[test_fc_df['Drug'] == selected_boxplot_drug, 'Concentration']
            .dropna()
            .unique()
            .tolist()
        )

    selected_boxplot_concentrations = st.multiselect(
        "Выберите концентрацию",
        options=concentration_options,
        default=[],
        format_func=lambda value: _format_numeric_label(value),
        key="metabolomics_boxplot_concentration"
    )

    selected_boxplot_metabolite = st.selectbox(
        "Выберите метаболит",
        options=[""] + list(available_metabolites),
        index=0,
        format_func=lambda value: "—" if value == "" else value,
        key="metabolomics_boxplot_metabolite"
    )

    with st.expander("Настройки Boxplot", expanded=True):
        st.caption("Шаблоны поддерживают переменные: {drug}, {metabolite}, {concentration}, {unit}")

        col1, col2, col3 = st.columns(3)
        with col1:
            boxplot_title_template = st.text_input(
                "Заголовок",
                value="{drug}",
                key="metabolomics_boxplot_title"
            )
        with col2:
            boxplot_xlabel_template = st.text_input(
                "X-label",
                value="Concentration",
                key="metabolomics_boxplot_xlabel",
                help="Можно оставить literal-текст или использовать шаблон с переменными."
            )
        with col3:
            boxplot_pvalue_threshold = st.selectbox(
                "Порог p-value для подписи",
                options=[0.001, 0.01, 0.05, 0.1],
                index=2,
                format_func=lambda x: f"{x} ({'***' if x == 0.001 else '**' if x == 0.01 else '*' if x == 0.05 else 'ns'})",
                key="metabolomics_boxplot_pvalue_threshold"
            )

        col4, col5, col6 = st.columns(3)
        with col4:
            boxplot_control_tick = st.text_input(
                "Подпись тика Control",
                value="Control",
                key="metabolomics_boxplot_control_tick"
            )
        with col5:
            boxplot_test_tick_template = st.text_input(
                "Подпись тика Test",
                value="{concentration}",
                key="metabolomics_boxplot_test_tick"
            )
        with col6:
            show_boxplot_points = st.checkbox(
                "Показывать точки",
                value=True,
                key="metabolomics_boxplot_show_points"
            )

    build_boxplot = st.button("Построить Boxplot", key="metabolomics_boxplot_build")
    if build_boxplot:
        if selected_boxplot_drug and selected_boxplot_concentrations and selected_boxplot_metabolite:
            st.session_state['metabolomics_boxplot_request'] = {
                'drug': selected_boxplot_drug,
                'concentrations': [conc for conc in concentration_options if conc in selected_boxplot_concentrations],
                'metabolite': selected_boxplot_metabolite,
            }
        else:
            st.session_state['metabolomics_boxplot_request'] = None
            st.warning("Выберите препарат, хотя бы одну концентрацию и метаболит перед построением boxplot.")

    request = st.session_state.get('metabolomics_boxplot_request')
    if not request:
        st.info("Выберите параметры выше и нажмите «Построить Boxplot», чтобы отрисовать график.")
        return

    concentration_labels = [_format_numeric_label(conc) for conc in request['concentrations']]
    concentration_label = ", ".join(concentration_labels)
    unit_label = measure_unit_by_drug.get(request['drug'], '')
    context = {
        "drug": request['drug'],
        "metabolite": request['metabolite'],
        "concentration": concentration_label,
        "unit": unit_label,
    }

    title_label = _format_template_label(
        boxplot_title_template,
        context,
        fallback=request['drug']
    )
    xlabel_label = _format_template_label(
        boxplot_xlabel_template,
        context,
        fallback="Concentration"
    )
    def format_test_tick(concentration):
        conc_context = {
            **context,
            "concentration": _format_numeric_label(concentration),
        }
        fallback_tick = _format_numeric_label(concentration)
        return _format_template_label(
            boxplot_test_tick_template,
            conc_context,
            fallback=fallback_tick
        )

    safe_drug_key = re.sub(r'[^\w.-]+', '_', str(request['drug']), flags=re.UNICODE).strip('_') or "drug"
    safe_metabolite_key = re.sub(r'[^\w.-]+', '_', str(request['metabolite']), flags=re.UNICODE).strip('_') or "metabolite"
    safe_conc_key = re.sub(r'[^\w.-]+', '_', concentration_label, flags=re.UNICODE).strip('_') or "conc"

    ylabel_label = st.text_input(
        "Y-label",
        value=request['metabolite'],
        key=f"metabolomics_boxplot_ylabel_{safe_drug_key}_{safe_metabolite_key}_{safe_conc_key}"
    )

    boxplot_source_df = prepare_metabolomics_boxplot_data(
        df,
        drug=request['drug'],
        concentration=request['concentrations'],
        metabolite=request['metabolite'],
        control_tick=boxplot_control_tick,
        test_tick=format_test_tick
    )

    if boxplot_source_df.empty or boxplot_source_df['GroupKey'].nunique() < 2:
        st.warning("Не удалось построить boxplot: недостаточно raw-данных в control/test.")
        return

    summary_rows = []
    signif_by_group = {}
    for conc in request['concentrations']:
        stats_row = test_fc_df[
            (test_fc_df['Drug'] == request['drug']) &
            (test_fc_df['Concentration'] == conc)
        ]
        if stats_row.empty:
            effect_value = np.nan
            p_value = np.nan
        else:
            stats_row = stats_row.iloc[0]
            effect_value = stats_row.get(request['metabolite'], np.nan)
            p_value = stats_row.get(f"{request['metabolite']} (pvalue)", np.nan)

        group_key = _metabolomics_test_group_key(conc)
        signif_by_group[group_key] = _format_metabolomics_significance_label(
            p_value,
            selected_threshold=boxplot_pvalue_threshold
        )
        summary_rows.append({
            "Concentration": _format_numeric_label(conc),
            "p-value": p_value,
            st.session_state['fc_mode']: effect_value
        })

    figure_width = max(6.5, 2.8 * boxplot_source_df['GroupKey'].nunique())
    fig_boxplot = build_metabolomics_boxplot(
        boxplot_source_df,
        title=title_label,
        ylabel=ylabel_label,
        xlabel=xlabel_label,
        show_points=show_boxplot_points,
        figsize=(figure_width, 4.5),
        signif_by_group=signif_by_group
    )

    if fig_boxplot is None:
        st.warning("Не удалось построить boxplot для выбранной комбинации.")
        return

    st.pyplot(fig_boxplot, use_container_width=True)
    summary_df = pd.DataFrame(summary_rows)
    if not summary_df.empty:
        st.dataframe(summary_df, hide_index=True, use_container_width=True)

    st.download_button(
        "Скачать boxplot (PNG)",
        data=fig_to_png_bytes(fig_boxplot),
        file_name=f"boxplot_{safe_drug_key}_{safe_metabolite_key}_{safe_conc_key}.png",
        mime="image/png",
        key=f"download_metabolomics_boxplot_{safe_drug_key}_{safe_metabolite_key}_{safe_conc_key}"
    )
    plt.close(fig_boxplot)


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
                measure_unit_by_drug = dict(zip(df['Drug'].unique().tolist(), list_measure_unit_concentration))
                
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

                all_available_metabolites = [
                    col for col in fc_df.columns
                    if col not in ['Drug', 'Group', 'Concentration']
                    and not col.endswith('(pvalue)')
                ]

                fc_barplot_tab, volcano_tab, boxplots_tab = st.tabs([
                    "Диаграмма по Fold Change",
                    "Volcano Plot",
                    "Boxplots"
                ])

                with fc_barplot_tab:
                    render_metabolomics_fc_barplot_section(
                        fc_df,
                        list_measure_unit_concentration
                    )

                with volcano_tab:
                    if st.session_state["fc_mode"] == 'log₂(B/A)':
                        render_metabolomics_volcano_section(fc_df, all_available_metabolites)
                    else:
                        st.info("Volcano Plot доступен только для режима log₂(B/A).")

                with boxplots_tab:
                    render_metabolomics_boxplot_section(
                        df,
                        fc_df,
                        measure_unit_by_drug,
                        all_available_metabolites
                    )

                return

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
                all_available_metabolites = [
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

                barplot_tab, volcano_tab, boxplots_tab = st.tabs([
                    "Fold Change",
                    "Volcano Plot",
                    "Boxplots"
                ])

                with barplot_tab:
                    if st.session_state['show_plot']:
                        filtered_fc_df = fc_df_for_diagram[
                            (fc_df_for_diagram['Drug'].isin(selected_drugs)) &
                            (fc_df_for_diagram['Concentration'].isin([conc for drug in selected_drugs for conc in selected_concentrations[drug]]))
                        ]

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
                            custom_colors=custom_colors
                        )
                    else:
                        st.info("Настройте параметры и нажмите «Перерисовать график», чтобы построить Fold Change chart.")

                with volcano_tab:
                    if st.session_state["fc_mode"] == 'log₂(B/A)':
                        render_metabolomics_volcano_section(fc_df, all_available_metabolites)
                    else:
                        st.info("Volcano Plot доступен только для режима log₂(B/A).")

                with boxplots_tab:
                    render_metabolomics_boxplot_section(
                        df,
                        fc_df,
                        measure_unit_by_drug,
                        all_available_metabolites
                    )

                return

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
                                    legend_title = st.text_input(
                                        "Заголовок легенды",
                                        value="Drugs",
                                        key='volcano_legend_title'
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

                                st.write("Подписи препаратов на графике:")
                                st.caption("Изменённые подписи будут использоваться в легенде Volcano Plot и в таблице значимых метаболитов.")
                                drug_labels_df = pd.DataFrame({
                                    'Drug': volcano_drugs,
                                    'Display Name': volcano_drugs
                                })
                                edited_drug_labels = st.data_editor(
                                    drug_labels_df,
                                    hide_index=True,
                                    use_container_width=True,
                                    disabled=['Drug'],
                                    column_config={
                                        'Drug': st.column_config.TextColumn('Препарат'),
                                        'Display Name': st.column_config.TextColumn('Подпись на графике')
                                    },
                                    key='volcano_drug_labels_editor'
                                )
                                drug_label_map = {
                                    row['Drug']: row['Display Name'].strip()
                                    for _, row in edited_drug_labels.iterrows()
                                    if isinstance(row['Display Name'], str)
                                    and row['Display Name'].strip()
                                    and row['Display Name'].strip() != row['Drug']
                                }

                                st.write("Подписи метаболитов на графике:")
                                st.caption("Изменённые подписи будут использоваться для значимых точек Volcano Plot и в таблице значимых метаболитов.")
                                metabolite_labels_df = pd.DataFrame({
                                    'Metabolite': available_metabolites,
                                    'Display Name': available_metabolites
                                })
                                edited_metabolite_labels = st.data_editor(
                                    metabolite_labels_df,
                                    hide_index=True,
                                    use_container_width=True,
                                    disabled=['Metabolite'],
                                    column_config={
                                        'Metabolite': st.column_config.TextColumn('Метаболит'),
                                        'Display Name': st.column_config.TextColumn('Подпись на графике')
                                    },
                                    key='volcano_metabolite_labels_editor'
                                )
                                metabolite_label_map = {
                                    row['Metabolite']: row['Display Name'].strip()
                                    for _, row in edited_metabolite_labels.iterrows()
                                    if isinstance(row['Display Name'], str)
                                    and row['Display Name'].strip()
                                    and row['Display Name'].strip() != row['Metabolite']
                                }
                            
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
                                vline_color=vline_color,
                                legend_title=legend_title,
                                metabolite_label_map=metabolite_label_map,
                                drug_label_map=drug_label_map
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
                                st.subheader("📊 Boxplots по значимым метаболитам")

                                significant_boxplot_df = (
                                    significant_df.copy()
                                    .drop_duplicates(subset=['Drug', 'Concentration', 'Metabolite'])
                                    .sort_values(['Drug', 'Concentration', 'Metabolite'])
                                    .reset_index(drop=True)
                                )

                                if 'Display Drug' not in significant_boxplot_df.columns:
                                    significant_boxplot_df['Display Drug'] = significant_boxplot_df['Drug']
                                if 'Display Metabolite' not in significant_boxplot_df.columns:
                                    significant_boxplot_df['Display Metabolite'] = significant_boxplot_df['Metabolite']

                                significant_boxplot_df['Concentration Label'] = significant_boxplot_df['Concentration'].apply(_format_numeric_label)
                                significant_boxplot_df['Unit'] = significant_boxplot_df['Drug'].map(measure_unit_by_drug).fillna('')
                                significant_boxplot_df['Pair Label'] = significant_boxplot_df.apply(
                                    lambda row: (
                                        f"{row['Display Drug']} | "
                                        f"{row['Display Metabolite']} | "
                                        f"{row['Concentration Label']}"
                                        f"{(' ' + row['Unit']) if row['Unit'] else ''}"
                                    ),
                                    axis=1
                                )

                                drug_options = significant_boxplot_df['Drug'].drop_duplicates().tolist()
                                drug_display_map = (
                                    significant_boxplot_df[['Drug', 'Display Drug']]
                                    .drop_duplicates()
                                    .set_index('Drug')['Display Drug']
                                    .to_dict()
                                )

                                selected_boxplot_drug = st.selectbox(
                                    "Выберите препарат для boxplot",
                                    options=[""] + drug_options,
                                    index=0,
                                    format_func=lambda value: "—" if value == "" else drug_display_map.get(value, value),
                                    key="metabolomics_boxplot_drug"
                                )

                                metabolite_options = []
                                metabolite_display_map = {}
                                if selected_boxplot_drug:
                                    metabolite_subset = significant_boxplot_df[
                                        significant_boxplot_df['Drug'] == selected_boxplot_drug
                                    ][['Metabolite', 'Display Metabolite']].drop_duplicates()
                                    metabolite_options = metabolite_subset['Metabolite'].tolist()
                                    metabolite_display_map = (
                                        metabolite_subset
                                        .set_index('Metabolite')['Display Metabolite']
                                        .to_dict()
                                    )

                                selected_boxplot_metabolite = st.selectbox(
                                    "Выберите метаболит для boxplot",
                                    options=[""] + metabolite_options,
                                    index=0,
                                    format_func=lambda value: "—" if value == "" else metabolite_display_map.get(value, value),
                                    key="metabolomics_boxplot_metabolite"
                                )

                                with st.expander("Настройки Boxplots", expanded=True):
                                    st.caption("Шаблоны поддерживают переменные: {drug}, {metabolite}, {concentration}, {unit}")

                                    col1, col2 = st.columns(2)
                                    with col1:
                                        boxplot_title_template = st.text_input(
                                            "Заголовок",
                                            value="{drug}",
                                            key="metabolomics_boxplot_title"
                                        )
                                    with col2:
                                        boxplot_xlabel_template = st.text_input(
                                            "X-label",
                                            value="Concentration",
                                            key="metabolomics_boxplot_xlabel",
                                            help="Можно оставить literal-текст или использовать шаблон с переменными."
                                        )

                                    col4, col5, col6 = st.columns(3)
                                    with col4:
                                        boxplot_control_tick = st.text_input(
                                            "Подпись тика Control",
                                            value="Control",
                                            key="metabolomics_boxplot_control_tick"
                                        )
                                    with col5:
                                        boxplot_test_tick_template = st.text_input(
                                            "Подпись тика Test",
                                            value="{concentration}",
                                            key="metabolomics_boxplot_test_tick"
                                        )
                                    with col6:
                                        show_boxplot_points = st.checkbox(
                                            "Показывать точки",
                                            value=True,
                                            key="metabolomics_boxplot_show_points"
                                        )

                                boxplot_targets = significant_boxplot_df[
                                    (significant_boxplot_df['Drug'] == selected_boxplot_drug) &
                                    (significant_boxplot_df['Metabolite'] == selected_boxplot_metabolite)
                                ].reset_index(drop=True)

                                if not selected_boxplot_drug or not selected_boxplot_metabolite:
                                    st.info("Выберите препарат и метаболит, чтобы построить boxplot.")
                                elif boxplot_targets.empty:
                                    st.info("Для выбранной комбинации нет значимого метаболита для построения boxplot.")
                                else:
                                    boxplot_columns = st.columns(2)

                                    for idx, row in boxplot_targets.iterrows():
                                        concentration_label = row['Concentration Label']
                                        unit_label = row['Unit']
                                        default_test_tick = (
                                            f"{concentration_label}{(' ' + unit_label) if unit_label else ''}"
                                        ).strip()
                                        context = {
                                            "drug": row['Display Drug'],
                                            "metabolite": row['Display Metabolite'],
                                            "concentration": concentration_label,
                                            "unit": unit_label,
                                        }

                                        title_label = _format_template_label(
                                            boxplot_title_template,
                                            context,
                                            fallback=row['Display Drug']
                                        )
                                        xlabel_label = _format_template_label(
                                            boxplot_xlabel_template,
                                            context,
                                            fallback="Concentration"
                                        )
                                        test_tick_label = _format_template_label(
                                            boxplot_test_tick_template,
                                            context,
                                            fallback=default_test_tick
                                        )

                                        boxplot_source_df = prepare_metabolomics_boxplot_data(
                                            df,
                                            drug=row['Drug'],
                                            concentration=row['Concentration'],
                                            metabolite=row['Metabolite'],
                                            control_tick=boxplot_control_tick,
                                            test_tick=test_tick_label
                                        )

                                        with boxplot_columns[idx % 2]:
                                            if boxplot_source_df.empty or boxplot_source_df['GroupKey'].nunique() < 2:
                                                st.warning(
                                                    f"Не удалось построить boxplot для {row['Pair Label']}: недостаточно raw-данных в control/test."
                                                )
                                                continue

                                            safe_drug_key = re.sub(r'[^\w.-]+', '_', str(row['Drug']), flags=re.UNICODE).strip('_') or "drug"
                                            safe_metabolite_key = re.sub(r'[^\w.-]+', '_', str(row['Metabolite']), flags=re.UNICODE).strip('_') or "metabolite"

                                            ylabel_label = st.text_input(
                                                "Y-label",
                                                value=row['Display Metabolite'],
                                                key=f"metabolomics_boxplot_ylabel_{safe_drug_key}_{safe_metabolite_key}"
                                            )

                                            figure_width = max(6.5, 2.8 * boxplot_source_df['GroupKey'].nunique())
                                            fig_boxplot = build_metabolomics_boxplot(
                                                boxplot_source_df,
                                                title=title_label,
                                                ylabel=ylabel_label,
                                                xlabel=xlabel_label,
                                                show_points=show_boxplot_points,
                                                figsize=(figure_width, 4.5),
                                                signif_by_group={
                                                    "test": _format_metabolomics_significance_label(
                                                        row['p_value'],
                                                        selected_threshold=p_value_threshold
                                                    )
                                                }
                                            )

                                            if fig_boxplot is None:
                                                st.warning(f"Не удалось построить boxplot для {row['Pair Label']}.")
                                                continue

                                            st.pyplot(fig_boxplot, use_container_width=True)
                                            st.caption(
                                                f"p-value = {row['p_value']:.4g}; log2FC = {row['log2FC']:.4g}"
                                            )

                                            safe_drug = safe_drug_key
                                            safe_metabolite = safe_metabolite_key
                                            safe_conc = re.sub(r'[^\w.-]+', '_', concentration_label, flags=re.UNICODE).strip('_') or "conc"

                                            st.download_button(
                                                "Скачать boxplot (PNG)",
                                                data=fig_to_png_bytes(fig_boxplot),
                                                file_name=f"boxplot_{safe_drug}_{safe_metabolite}_{safe_conc}.png",
                                                mime="image/png",
                                                key=f"download_metabolomics_boxplot_{idx}"
                                            )
                                            plt.close(fig_boxplot)
                            elif significant_df is not None and significant_df.empty:
                                st.warning("Нет значимых метаболитов при заданных параметрах.")
