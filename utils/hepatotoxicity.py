import re
from io import BytesIO

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from utils.hr_cardio import (
    build_decision_and_global_table,
    build_hr_boxplot,
    build_levene_report_table,
    build_normality_report_table,
    build_posthoc_reports,
    dataframe_to_excel_bytes,
)


REQUIRED_HEPATO_COLUMNS = ["Sample", "Compound", "Concentration"]
ALPHA_OPTIONS = [0.05, 0.01, 0.001]


def _normalize_column_name(value: str) -> str:
    return re.sub(r"\s+", " ", str(value).strip()).lower()


def _find_required_columns(df: pd.DataFrame) -> dict[str, str]:
    normalized_map = {_normalize_column_name(col): col for col in df.columns}
    matched: dict[str, str] = {}

    for required in REQUIRED_HEPATO_COLUMNS:
        source_col = normalized_map.get(_normalize_column_name(required))
        if source_col is None:
            raise ValueError(
                "В файле отсутствуют обязательные столбцы: "
                + ", ".join(REQUIRED_HEPATO_COLUMNS)
            )
        matched[required] = source_col

    return matched


def _format_concentration(value) -> str:
    conc = pd.to_numeric(pd.Series([value]), errors="coerce").iat[0]
    if pd.isna(conc):
        return ""
    if float(conc).is_integer():
        return str(int(conc))
    return f"{float(conc):g}"


def _build_group_label(compound, concentration) -> str:
    compound_text = "" if pd.isna(compound) else str(compound).strip()
    if not compound_text:
        return ""
    if compound_text.lower() == "control":
        return "Control"

    concentration_text = _format_concentration(concentration)
    return f"{compound_text} {concentration_text}".strip() if concentration_text else compound_text


def _dose_from_group(group_name: str) -> float:
    group_text = str(group_name).strip()
    if group_text.lower() == "control":
        return 0.0

    match = re.search(r"(\d+(\.\d+)?)", group_text)
    if match:
        return float(match.group(1))

    return np.nan


def _metric_slug(metric: str) -> str:
    slug = re.sub(r"[^0-9A-Za-z]+", "_", str(metric)).strip("_").lower()
    return slug or "metric"


def _fig_to_png_bytes(fig) -> bytes:
    buffer = BytesIO()
    fig.savefig(buffer, format="png", dpi=600, bbox_inches="tight")
    buffer.seek(0)
    return buffer.getvalue()


def _build_signif_map_for_metric(
    analysis_df: pd.DataFrame,
    metric: str,
    posthoc_details: dict[str, pd.DataFrame],
    control_label: str = "Control",
) -> dict[str, str]:
    groups = sorted(analysis_df["Compound"].astype(str).unique().tolist(), key=_dose_from_group)
    result = {
        group: "ns"
        for group in groups
        if group.strip().lower() != control_label.strip().lower()
    }

    posthoc_df = posthoc_details.get(metric)
    if posthoc_df is None or posthoc_df.empty:
        return result

    if "Group 2" not in posthoc_df.columns or "p_value_adj" not in posthoc_df.columns:
        return result

    for _, row in posthoc_df.iterrows():
        group_name = str(row["Group 2"])
        p_value = pd.to_numeric(pd.Series([row["p_value_adj"]]), errors="coerce").iat[0]
        if group_name.strip().lower() == control_label.strip().lower():
            continue
        if pd.isna(p_value):
            result[group_name] = "ns"
        elif p_value < 0.001:
            result[group_name] = "***"
        elif p_value < 0.01:
            result[group_name] = "**"
        elif p_value < 0.05:
            result[group_name] = "*"
        else:
            result[group_name] = "ns"

    return result


def load_hepatotoxicity_tables(file_like) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    raw_df = pd.read_excel(file_like)
    required_columns = _find_required_columns(raw_df)
    raw_df = raw_df.rename(columns={source: target for target, source in required_columns.items()})

    metrics = [col for col in raw_df.columns if col not in REQUIRED_HEPATO_COLUMNS]
    if not metrics:
        raise ValueError("В файле не найдены числовые метрики для анализа.")

    prepared_df = raw_df.copy()
    prepared_df["Sample"] = prepared_df["Sample"].astype(str).str.strip()
    prepared_df["Compound"] = prepared_df["Compound"].astype(str).str.strip()
    prepared_df["Concentration"] = pd.to_numeric(prepared_df["Concentration"], errors="coerce")

    valid_metrics: list[str] = []
    for metric in metrics:
        prepared_df[metric] = pd.to_numeric(prepared_df[metric], errors="coerce")
        if prepared_df[metric].notna().any():
            valid_metrics.append(metric)

    if not valid_metrics:
        raise ValueError("В файле не найдены числовые столбцы с измеряемыми показателями.")

    prepared_df = prepared_df.dropna(subset=["Compound"]).copy()
    prepared_df["Analysis group"] = prepared_df.apply(
        lambda row: _build_group_label(row["Compound"], row["Concentration"]),
        axis=1,
    )
    prepared_df = prepared_df[prepared_df["Analysis group"] != ""].copy()
    prepared_df = prepared_df.dropna(subset=valid_metrics, how="all").reset_index(drop=True)

    analysis_df = prepared_df[["Analysis group", *valid_metrics]].rename(
        columns={"Analysis group": "Compound"}
    )

    raw_display_df = prepared_df[["Sample", "Compound", "Concentration", *valid_metrics, "Analysis group"]]
    return raw_display_df, analysis_df, valid_metrics


def build_hepatotoxicity_summary_table(
    analysis_df: pd.DataFrame,
    metrics: list[str],
) -> pd.DataFrame:
    grouped = analysis_df.groupby("Compound", dropna=False)
    summary_df = pd.DataFrame({"Compound": list(grouped.groups.keys())})

    for metric in metrics:
        aggregated = grouped[metric].agg(["count", "mean", "std"]).reset_index()
        aggregated[f"Standard Error | {metric}"] = aggregated["std"] / np.sqrt(aggregated["count"])
        aggregated = aggregated.rename(
            columns={
                "count": f"N | {metric}",
                "mean": f"Mean | {metric}",
            }
        )
        aggregated = aggregated.drop(columns=["std"])
        summary_df = summary_df.merge(aggregated, on="Compound", how="left")

    summary_df["_dose"] = summary_df["Compound"].apply(_dose_from_group)
    summary_df = (
        summary_df
        .sort_values(["_dose", "Compound"], na_position="last")
        .drop(columns="_dose")
        .reset_index(drop=True)
    )
    return summary_df


def hepatotoxicity_app():
    st.title("🧪 Гепатотоксичность")

    uploaded_file = st.file_uploader(
        label="Загрузите Excel файл (morphology)",
        type=["xlsx", "xls"],
        key="hepato_uploader",
    )

    if uploaded_file is None:
        return

    try:
        raw_df, analysis_df, metrics = load_hepatotoxicity_tables(uploaded_file)
    except Exception as exc:
        st.error(f"Ошибка при обработке файла гепатотоксичности: {exc}")
        st.stop()

    st.subheader("📋 Исходная таблица морфометрии")
    st.dataframe(raw_df, use_container_width=True)
    st.download_button(
        label="⬇️ Скачать исходную таблицу (Excel)",
        data=dataframe_to_excel_bytes(raw_df, sheet_name="Hepato_raw"),
        file_name="hepatotoxicity_raw.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

    summary_df = build_hepatotoxicity_summary_table(analysis_df, metrics)
    st.subheader("📊 Агрегация по Compound (N / Mean / SE)")
    st.dataframe(summary_df, use_container_width=True)
    st.download_button(
        label="⬇️ Скачать агрегированную таблицу (Excel)",
        data=dataframe_to_excel_bytes(summary_df, sheet_name="Hepato_summary"),
        file_name="hepatotoxicity_summary.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

    st.subheader("Статистический анализ")

    st.subheader("🧪 Проверка нормальности по группам (Shapiro–Wilk)")
    alpha_norm = st.selectbox(
        "α (нормальность, Shapiro)",
        ALPHA_OPTIONS,
        index=0,
        key="hepato_alpha_norm",
    )
    normality_df = build_normality_report_table(analysis_df, metrics=metrics, alpha=alpha_norm)
    st.dataframe(normality_df, use_container_width=True)
    st.download_button(
        label="⬇️ Скачать таблицу нормальности (Excel)",
        data=dataframe_to_excel_bytes(normality_df, sheet_name="Normality"),
        file_name="hepato_normality_shapiro.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

    st.subheader("🧪 Проверка гомогенности дисперсий (Levene)")
    alpha_levene = st.selectbox(
        "α (дисперсии, Levene)",
        ALPHA_OPTIONS,
        index=0,
        key="hepato_alpha_levene",
    )
    levene_df = build_levene_report_table(analysis_df, metrics=metrics, alpha=alpha_levene)
    st.dataframe(levene_df, use_container_width=True)
    st.download_button(
        label="⬇️ Скачать таблицу Levene (Excel)",
        data=dataframe_to_excel_bytes(levene_df, sheet_name="Levene"),
        file_name="hepato_levene.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

    st.subheader("📌 Выбор теста + глобальный тест (ANOVA / Kruskal–Wallis)")
    alpha_global = st.selectbox(
        "α (глобальный тест)",
        ALPHA_OPTIONS,
        index=0,
        key="hepato_alpha_global",
    )
    decision_global_df = build_decision_and_global_table(
        hr_df=analysis_df,
        normality_df=normality_df,
        levene_df=levene_df,
        metrics=metrics,
        alpha=alpha_global,
    )
    st.dataframe(decision_global_df, use_container_width=True)
    st.download_button(
        label="⬇️ Скачать выбор теста и глобальный тест (Excel)",
        data=dataframe_to_excel_bytes(decision_global_df, sheet_name="Decision_Global"),
        file_name="hepato_decision_global.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

    st.subheader("🔍 Post-hoc анализ (vs Control)")
    alpha_posthoc = st.selectbox(
        "α (post-hoc)",
        ALPHA_OPTIONS,
        index=0,
        key="hepato_alpha_posthoc",
    )
    posthoc_summary_df, posthoc_details = build_posthoc_reports(
        hr_df=analysis_df,
        decision_global_df=decision_global_df,
        metrics=metrics,
        alpha_posthoc=alpha_posthoc,
        alpha_global=alpha_global,
        control_label="Control",
        only_vs_control=True,
        dunn_p_adjust="holm",
    )
    st.dataframe(posthoc_summary_df, use_container_width=True)
    st.download_button(
        label="⬇️ Скачать сводку post-hoc (Excel)",
        data=dataframe_to_excel_bytes(posthoc_summary_df, sheet_name="Posthoc_summary"),
        file_name="hepato_posthoc_summary.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

    if posthoc_details:
        for metric, metric_df in posthoc_details.items():
            st.markdown(f"**{metric}**")
            st.dataframe(metric_df, use_container_width=True)
            st.download_button(
                label=f"⬇️ Скачать post-hoc для {metric} (Excel)",
                data=dataframe_to_excel_bytes(metric_df, sheet_name="Posthoc"),
                file_name=f"hepato_posthoc_{_metric_slug(metric)}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                key=f"download_hepato_posthoc_{_metric_slug(metric)}",
            )

    st.subheader("📊 Boxplots")
    default_title = uploaded_file.name.rsplit(".", 1)[0] if getattr(uploaded_file, "name", None) else "Hepatotoxicity"
    plot_title = st.text_input("Заголовок", value=default_title, key="hepato_plot_title")
    xlabel_common = st.text_input(
        "X-label",
        value="",
        key="hepato_xlabel",
        help="Оставь пустым для автоподстановки названия соединения из групп Compound.",
    )
    control_tick = st.text_input(
        "Подпись тика Control",
        value="Control",
        key="hepato_control_tick",
    )
    show_points = st.checkbox(
        "Показывать точки (individual values)",
        value=True,
        key="hepato_show_points",
    )

    ylabels: dict[str, str] = {}
    for metric in metrics:
        metric_slug = _metric_slug(metric)
        ylabels[metric] = st.text_input(
            f"Y-label для {metric}",
            value=metric,
            key=f"hepato_ylabel_{metric_slug}",
        )

    for metric in metrics:
        metric_slug = _metric_slug(metric)
        significance_map = _build_signif_map_for_metric(
            analysis_df=analysis_df,
            metric=metric,
            posthoc_details=posthoc_details,
            control_label="Control",
        )
        fig = build_hr_boxplot(
            hr_df=analysis_df,
            metric=metric,
            title=plot_title,
            ylabel=ylabels[metric],
            xlabel=xlabel_common,
            show_points=show_points,
            control_tick=control_tick,
            signif_by_group=significance_map,
        )
        st.markdown(f"**{metric}**")
        st.pyplot(fig, use_container_width=True)
        st.download_button(
            label=f"⬇️ Скачать boxplot для {metric} (PNG)",
            data=_fig_to_png_bytes(fig),
            file_name=f"hepato_boxplot_{metric_slug}.png",
            mime="image/png",
            key=f"download_hepato_boxplot_{metric_slug}",
        )
        plt.close(fig)
