import re
from io import BytesIO

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from utils.hr_cardio import (
    build_decision_and_global_table,
    build_levene_report_table,
    build_normality_report_table,
    build_posthoc_reports,
)


REQUIRED_HEPATO_COLUMNS = ["compound", "concentration", "liver_area", "yolk_sac_area"]
HEPATO_METRICS = ["liver_size_percent", "yolk_sac_retention_percent"]
ALPHA_OPTIONS = [0.05, 0.01, 0.001]

COLUMN_ALIASES = {
    "compound": {
        "compound", "drug", "preparation", "препарат", "вещество", "соединение",
        "лекарство", "drug name",
    },
    "sample": {
        "sample", "subject", "id", "sample id", "образец", "номер", "рыба", "fish",
    },
    "group": {
        "group", "группа",
    },
    "concentration": {
        "concentration", "conc", "dose", "доза", "концентрация",
    },
    "liver_area": {
        "liver area", "liver_area", "площадь печени", "печень", "liver",
    },
    "yolk_sac_area": {
        "yolk sac area", "yolk_sac_area", "площадь желточного мешка",
        "желточный мешок", "yolk sac", "yolk",
    },
}


def _clean_text(value) -> str:
    return re.sub(r"\s+", " ", str(value).strip())


def _normalize_key(value) -> str:
    value = _clean_text(value).lower().replace("ё", "е")
    value = value.replace("%", " percent ")
    value = re.sub(r"[_\-./(),;:]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def _format_concentration(value) -> str:
    if pd.isna(value):
        return ""
    value = float(value)
    if value.is_integer():
        return str(int(value))
    return f"{value:g}"


def _format_mean_se(mean_value, se_value) -> str:
    if pd.isna(mean_value):
        return ""
    if pd.isna(se_value):
        return f"{mean_value:.1f} ± —"
    return f"{mean_value:.1f} ± {se_value:.1f}"


def _format_p_value(value) -> str:
    if pd.isna(value):
        return "ns"
    value = float(value)
    if value < 0.001:
        return "<0.001"
    return f"{value:.4f}".rstrip("0").rstrip(".")


def _p_to_significance(value) -> str:
    if pd.isna(value):
        return "ns"
    value = float(value)
    if value < 0.001:
        return "***"
    if value < 0.01:
        return "**"
    if value < 0.05:
        return "*"
    return "ns"


def _metric_title(metric: str) -> str:
    if metric == "liver_size_percent":
        return "Liver size (% of control)"
    if metric == "yolk_sac_retention_percent":
        return "Yolk sac size (% of control)"
    return metric


def _metric_slug(metric: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", str(metric)).strip("_").lower()


def _soft_numeric(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.strip()
    text = text.str.replace("\u00a0", "", regex=False)
    text = text.str.replace(" ", "", regex=False)
    text = text.str.replace(",", ".", regex=False)
    extracted = text.str.extract(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", expand=False)
    return pd.to_numeric(extracted, errors="coerce")


def _sort_summary(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["_compound_sort"] = out["compound"].astype(str).str.lower()
    out = out.sort_values(["_compound_sort", "concentration"]).drop(columns="_compound_sort")
    return out.reset_index(drop=True)


def _fig_to_png_bytes(fig, dpi: int = 600) -> bytes:
    buffer = BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi, bbox_inches="tight")
    buffer.seek(0)
    return buffer.getvalue()


def _safe_filename(value: str) -> str:
    filename = re.sub(r"[^0-9A-Za-zА-Яа-я_.-]+", "_", str(value)).strip("_")
    return filename or "plot"


def normalize_hepatotoxicity_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Нормализует названия колонок."""
    out = df.copy()
    rename_map = {}
    normalized_columns = {_normalize_key(col): col for col in out.columns}

    for target, aliases in COLUMN_ALIASES.items():
        alias_keys = {_normalize_key(alias) for alias in aliases}
        for alias in alias_keys:
            source_col = normalized_columns.get(alias)
            if source_col is not None:
                rename_map[source_col] = target
                break

    out = out.rename(columns=rename_map)
    return out


def validate_hepatotoxicity_dataframe(df: pd.DataFrame) -> list[str]:
    """Проверяет обязательные колонки."""
    missing = [col for col in REQUIRED_HEPATO_COLUMNS if col not in df.columns]
    if not missing:
        return []

    labels = {
        "compound": "Drug/Compound",
        "concentration": "Concentration",
        "liver_area": "Площадь печени / Liver area",
        "yolk_sac_area": "Площадь желточного мешка / Yolk sac area",
    }
    return [f"Отсутствуют обязательные колонки: {', '.join(labels[col] for col in missing)}"]


def prepare_hepatotoxicity_dataframe(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Готовит данные к расчетам."""
    warnings: list[str] = []
    out = normalize_hepatotoxicity_columns(df)

    errors = validate_hepatotoxicity_dataframe(out)
    if errors:
        raise ValueError("; ".join(errors))

    keep_cols = [col for col in ["sample", "compound", "group", "concentration", "liver_area", "yolk_sac_area"] if col in out.columns]
    out = out[keep_cols].copy()

    if "sample" not in out.columns:
        out.insert(0, "sample", np.arange(1, len(out) + 1))
    if "group" not in out.columns:
        out.insert(2, "group", "")

    out["sample"] = out["sample"].astype(str).str.strip()
    out["compound"] = out["compound"].astype(str).map(_clean_text)
    out["group"] = out["group"].astype(str).map(_clean_text)
    out["concentration"] = _soft_numeric(out["concentration"])
    out["liver_area"] = _soft_numeric(out["liver_area"])
    out["yolk_sac_area"] = _soft_numeric(out["yolk_sac_area"])

    invalid_concentration = out["concentration"].isna().sum()
    if invalid_concentration:
        warnings.append(
            f"Удалено строк с нечисловой концентрацией: {invalid_concentration}."
        )

    invalid_metrics = out[["liver_area", "yolk_sac_area"]].isna().any(axis=1).sum()
    if invalid_metrics:
        warnings.append(
            f"Удалено строк с пустыми или нечисловыми площадями печени/желточного мешка: {invalid_metrics}."
        )

    empty_compound = out["compound"].eq("").sum()
    if empty_compound:
        warnings.append(f"Удалено строк без названия препарата: {empty_compound}.")

    out = out.dropna(subset=["concentration", "liver_area", "yolk_sac_area"])
    out = out[out["compound"] != ""].copy()
    out = out.reset_index(drop=True)

    if out.empty:
        raise ValueError("После очистки не осталось строк для анализа.")

    return out, warnings


def calculate_hepatotoxicity_metrics(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Считает проценты от контроля."""
    warnings: list[str] = []
    parts: list[pd.DataFrame] = []

    for compound, sub_df in df.groupby("compound", sort=False):
        control_df = sub_df[np.isclose(sub_df["concentration"], 0)]
        if control_df.empty:
            warnings.append(
                f"Для препарата {compound} нет контроля с concentration = 0. Препарат исключен из расчета."
            )
            continue

        mean_liver_control = control_df["liver_area"].mean()
        mean_yolk_control = control_df["yolk_sac_area"].mean()
        if not np.isfinite(mean_liver_control) or mean_liver_control == 0:
            warnings.append(
                f"Для препарата {compound} средняя площадь печени в контроле равна 0 или нечисловая. Препарат исключен."
            )
            continue
        if not np.isfinite(mean_yolk_control) or mean_yolk_control == 0:
            warnings.append(
                f"Для препарата {compound} средняя площадь желточного мешка в контроле равна 0 или нечисловая. Препарат исключен."
            )
            continue

        metric_df = sub_df.copy()
        metric_df["control_liver_area_mean"] = mean_liver_control
        metric_df["control_yolk_sac_area_mean"] = mean_yolk_control
        metric_df["liver_size_percent"] = metric_df["liver_area"] / mean_liver_control * 100
        metric_df["yolk_sac_retention_percent"] = metric_df["yolk_sac_area"] / mean_yolk_control * 100
        parts.append(metric_df)

    if not parts:
        raise ValueError("Нет препаратов, пригодных для расчета: у всех отсутствует контроль или контроль некорректен.")

    return pd.concat(parts, ignore_index=True), warnings


def _aggregate_metric(df: pd.DataFrame, metric: str, prefix: str) -> pd.DataFrame:
    agg = (
        df.groupby(["compound", "concentration"], dropna=False)[metric]
        .agg(["count", "mean", "std"])
        .reset_index()
    )
    agg[f"{prefix}_se"] = agg["std"] / np.sqrt(agg["count"])
    agg = agg.rename(
        columns={
            "count": "n",
            "mean": f"{prefix}_mean",
            "std": f"{prefix}_sd",
        }
    )
    return agg


def _analysis_df_for_compound(metrics_df: pd.DataFrame, compound: str, metric: str) -> pd.DataFrame:
    sub_df = metrics_df.loc[metrics_df["compound"] == compound, ["concentration", metric]].copy()
    sub_df["Compound"] = sub_df["concentration"].apply(
        lambda value: "Control" if np.isclose(value, 0) else _format_concentration(value)
    )
    return sub_df.rename(columns={metric: _metric_title(metric)})[["Compound", _metric_title(metric)]]


def _posthoc_p_values_by_compound(
    metrics_df: pd.DataFrame,
    metric: str,
    alpha_norm: float,
    alpha_levene: float,
    alpha_global: float,
    alpha_posthoc: float,
) -> tuple[dict[tuple[str, float], float], pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    p_values: dict[tuple[str, float], float] = {}
    normality_tables = []
    levene_tables = []
    decision_tables = []
    posthoc_tables = []
    metric_title = _metric_title(metric)

    for compound in metrics_df["compound"].drop_duplicates():
        analysis_df = _analysis_df_for_compound(metrics_df, compound, metric)
        concentrations = sorted(metrics_df.loc[metrics_df["compound"] == compound, "concentration"].dropna().unique())
        if len(concentrations) < 2:
            continue

        normality_df = build_normality_report_table(analysis_df, metrics=[metric_title], alpha=alpha_norm)
        levene_df = build_levene_report_table(analysis_df, metrics=[metric_title], alpha=alpha_levene)
        decision_df = build_decision_and_global_table(
            hr_df=analysis_df,
            normality_df=normality_df,
            levene_df=levene_df,
            metrics=[metric_title],
            alpha=alpha_global,
        )
        posthoc_summary_df, posthoc_details = build_posthoc_reports(
            hr_df=analysis_df,
            decision_global_df=decision_df,
            metrics=[metric_title],
            alpha_posthoc=alpha_posthoc,
            alpha_global=alpha_global,
            control_label="Control",
            only_vs_control=True,
            dunn_p_adjust="holm",
        )

        normality_tables.append(normality_df.assign(Drug=compound))
        levene_tables.append(levene_df.assign(Drug=compound))
        decision_tables.append(decision_df.assign(Drug=compound))
        posthoc_summary_df = posthoc_summary_df.assign(Drug=compound)
        posthoc_tables.append(posthoc_summary_df)

        detail_df = posthoc_details.get(metric_title)
        if detail_df is None or detail_df.empty:
            continue

        for _, row in detail_df.iterrows():
            concentration = _soft_numeric(pd.Series([row.get("Group 2")])).iat[0]
            p_value = pd.to_numeric(pd.Series([row.get("p_value_adj")]), errors="coerce").iat[0]
            if pd.notna(concentration) and pd.notna(p_value):
                p_values[(compound, float(concentration))] = float(p_value)

    empty_normality = pd.DataFrame(columns=["Drug", "Compound", "Metric", "N", "Shapiro_W", "p_value", "Normal (p>α)", "Note"])
    empty_levene = pd.DataFrame(columns=["Drug", "Metric", "Levene_stat", "p_value", "Equal_variance", "Note"])
    empty_decision = pd.DataFrame(columns=["Drug", "Metric", "Normality_OK", "Equal_variance", "Selected_test", "Statistic", "p_value", "Reject_H0", "Note"])
    empty_posthoc = pd.DataFrame(columns=["Drug", "Metric", "Selected_test", "Global_p_value", "Posthoc_test", "Posthoc_applied", "Note"])

    return (
        p_values,
        pd.concat(normality_tables, ignore_index=True) if normality_tables else empty_normality,
        pd.concat(levene_tables, ignore_index=True) if levene_tables else empty_levene,
        pd.concat(decision_tables, ignore_index=True) if decision_tables else empty_decision,
        pd.concat(posthoc_tables, ignore_index=True) if posthoc_tables else empty_posthoc,
    )


def build_hepatotoxicity_summary_tables(
    df_metrics: pd.DataFrame,
    alpha_norm: float = 0.05,
    alpha_levene: float = 0.05,
    alpha_global: float = 0.05,
    alpha_posthoc: float = 0.05,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame]]:
    """Строит итоговые таблицы."""
    liver_agg = _aggregate_metric(df_metrics, "liver_size_percent", "liver_size")
    yolk_agg = _aggregate_metric(df_metrics, "yolk_sac_retention_percent", "yolk_sac_retention")

    summary = liver_agg.merge(
        yolk_agg.drop(columns=["n"]),
        on=["compound", "concentration"],
        how="outer",
    )

    liver_p, liver_norm, liver_lev, liver_decision, liver_posthoc = _posthoc_p_values_by_compound(
        df_metrics, "liver_size_percent", alpha_norm, alpha_levene, alpha_global, alpha_posthoc
    )
    yolk_p, yolk_norm, yolk_lev, yolk_decision, yolk_posthoc = _posthoc_p_values_by_compound(
        df_metrics, "yolk_sac_retention_percent", alpha_norm, alpha_levene, alpha_global, alpha_posthoc
    )

    summary["liver_size_p_value"] = summary.apply(
        lambda row: np.nan if np.isclose(row["concentration"], 0) else liver_p.get((row["compound"], float(row["concentration"])), np.nan),
        axis=1,
    )
    summary["yolk_sac_retention_p_value"] = summary.apply(
        lambda row: np.nan if np.isclose(row["concentration"], 0) else yolk_p.get((row["compound"], float(row["concentration"])), np.nan),
        axis=1,
    )
    summary["liver_size_significance"] = summary["liver_size_p_value"].apply(_p_to_significance)
    summary["yolk_sac_retention_significance"] = summary["yolk_sac_retention_p_value"].apply(_p_to_significance)

    detailed_columns = [
        "compound",
        "concentration",
        "n",
        "liver_size_mean",
        "liver_size_sd",
        "liver_size_se",
        "liver_size_p_value",
        "liver_size_significance",
        "yolk_sac_retention_mean",
        "yolk_sac_retention_sd",
        "yolk_sac_retention_se",
        "yolk_sac_retention_p_value",
        "yolk_sac_retention_significance",
    ]
    summary = _sort_summary(summary[detailed_columns])

    journal = pd.DataFrame({
        "Drug": summary["compound"],
        "Concentration": summary["concentration"].apply(_format_concentration),
        "Liver size, % of control (mean ± SE)": summary.apply(
            lambda row: _format_mean_se(row["liver_size_mean"], row["liver_size_se"]),
            axis=1,
        ),
        "p-value": summary["liver_size_p_value"].apply(_format_p_value),
        "Yolk sac size, % of control (mean ± SE)": summary.apply(
            lambda row: _format_mean_se(row["yolk_sac_retention_mean"], row["yolk_sac_retention_se"]),
            axis=1,
        ),
        "p-value ": summary["yolk_sac_retention_p_value"].apply(_format_p_value),
    })

    stats_tables = {
        "normality": pd.concat([liver_norm, yolk_norm], ignore_index=True),
        "levene": pd.concat([liver_lev, yolk_lev], ignore_index=True),
        "decision_global": pd.concat([liver_decision, yolk_decision], ignore_index=True),
        "posthoc_summary": pd.concat([liver_posthoc, yolk_posthoc], ignore_index=True),
    }

    return summary, journal, stats_tables


def _plot_compound_bar(
    summary_df: pd.DataFrame,
    compound: str,
    metric_prefix: str,
    title: str,
    xlabel: str,
    ylabel: str,
):
    plot_df = (
        summary_df[summary_df["compound"] == compound]
        .sort_values("concentration")
        .reset_index(drop=True)
    )
    x_positions = np.arange(len(plot_df))
    means = pd.to_numeric(plot_df[f"{metric_prefix}_mean"], errors="coerce").to_numpy(dtype=float)
    se_values = pd.to_numeric(plot_df[f"{metric_prefix}_se"], errors="coerce").fillna(0).to_numpy(dtype=float)
    labels = plot_df["concentration"].apply(_format_concentration).tolist()

    fig, ax = plt.subplots(figsize=(5.8, 4.2))
    ax.bar(
        x_positions,
        means,
        yerr=se_values,
        width=0.68,
        color="#d9e6e2",
        edgecolor="black",
        linewidth=1,
        error_kw={
            "ecolor": "black",
            "elinewidth": 1.2,
            "capsize": 4,
            "capthick": 1.2,
        },
    )

    ax.axhline(100, linestyle="--", color="0.35", linewidth=1)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, rotation=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.18)

    y_top = np.nanmax(means + se_values) if len(means) else 100
    y_bottom = min(0, np.nanmin(means - se_values)) if len(means) else 0
    y_span = max(y_top - y_bottom, 1)
    label_offset = y_span * 0.045

    for idx, row in plot_df.iterrows():
        if np.isclose(row["concentration"], 0):
            continue
        mean_value = pd.to_numeric(pd.Series([row.get(f"{metric_prefix}_mean")]), errors="coerce").iat[0]
        se_value = pd.to_numeric(pd.Series([row.get(f"{metric_prefix}_se")]), errors="coerce").fillna(0).iat[0]
        significance = str(row.get(f"{metric_prefix}_significance", "ns"))
        if pd.isna(mean_value):
            continue
        ax.text(
            idx,
            float(mean_value) + float(se_value) + label_offset,
            significance,
            ha="center",
            va="bottom",
            fontsize=11,
            color="black",
        )

    ax.set_ylim(y_bottom, y_top + label_offset * 4)
    fig.tight_layout()
    return fig


def plot_liver_size(
    summary_df: pd.DataFrame,
    compound: str,
    title: str,
    xlabel: str = "Concentration",
    ylabel: str = "Liver size (% of control)",
):
    """График печени."""
    return _plot_compound_bar(summary_df, compound, "liver_size", title, xlabel, ylabel)


def plot_yolk_sac_retention(
    summary_df: pd.DataFrame,
    compound: str,
    title: str,
    xlabel: str = "Concentration",
    ylabel: str = "Yolk sac size (% of control)",
):
    """График желточного мешка."""
    return _plot_compound_bar(summary_df, compound, "yolk_sac_retention", title, xlabel, ylabel)


def export_hepatotoxicity_results(
    raw_normalized: pd.DataFrame,
    individual_metrics: pd.DataFrame,
    summary_detailed: pd.DataFrame,
    summary_journal: pd.DataFrame,
    stats_tables: dict[str, pd.DataFrame] | None = None,
) -> bytes:
    """Экспортирует Excel."""
    buffer = BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        raw_normalized.to_excel(writer, index=False, sheet_name="raw_normalized")
        individual_metrics.to_excel(writer, index=False, sheet_name="individual_metrics")
        summary_detailed.to_excel(writer, index=False, sheet_name="summary_detailed")
        summary_journal.to_excel(writer, index=False, sheet_name="summary_journal")
        for name, table in (stats_tables or {}).items():
            table.to_excel(writer, index=False, sheet_name=name[:31])
    buffer.seek(0)
    return buffer.getvalue()


def _dataframe_to_excel_bytes(df: pd.DataFrame, sheet_name: str) -> bytes:
    """Экспортирует один лист."""
    buffer = BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name[:31])
    buffer.seek(0)
    return buffer.getvalue()


def build_individual_hepatotoxicity_table(df_metrics: pd.DataFrame) -> pd.DataFrame:
    """Готовит индивидуальные проценты."""
    columns = [
        "sample",
        "compound",
        "group",
        "concentration",
        "liver_area",
        "control_liver_area_mean",
        "liver_size_percent",
        "yolk_sac_area",
        "control_yolk_sac_area_mean",
        "yolk_sac_retention_percent",
    ]
    available_columns = [col for col in columns if col in df_metrics.columns]
    out = df_metrics[available_columns].copy()
    out = out.rename(columns={"yolk_sac_retention_percent": "yolk_sac_size_percent"})
    return _sort_summary(out)


def display_hepatotoxicity_summary_table(summary_df: pd.DataFrame) -> pd.DataFrame:
    """Переименовывает желточный мешок для вывода."""
    return summary_df.rename(
        columns={
            "yolk_sac_retention_mean": "yolk_sac_size_mean",
            "yolk_sac_retention_sd": "yolk_sac_size_sd",
            "yolk_sac_retention_se": "yolk_sac_size_se",
            "yolk_sac_retention_p_value": "yolk_sac_size_p_value",
            "yolk_sac_retention_significance": "yolk_sac_size_significance",
        }
    )


def _read_uploaded_excel(uploaded_file) -> pd.DataFrame:
    return pd.read_excel(uploaded_file)


def hepatotoxicity_app():
    st.title("Гепатотоксичность")
    st.caption("Расчет Liver size и Yolk sac size как % от контроля внутри каждого препарата.")

    uploaded_file = st.file_uploader(
        label="Загрузите Excel-файл с данными гепатотоксичности (.xlsx)",
        type=["xlsx"],
        key="hepato_uploader",
    )

    if uploaded_file is None:
        st.info("Загрузите файл с колонками Drug, Concentration, Площадь печени и Площадь желточного мешка.")
        return

    try:
        raw_df = _read_uploaded_excel(uploaded_file)
        prepared_df, prepare_warnings = prepare_hepatotoxicity_dataframe(raw_df)
        metrics_df, metric_warnings = calculate_hepatotoxicity_metrics(prepared_df)
    except Exception as exc:
        st.error(f"Ошибка при обработке файла гепатотоксичности: {exc}")
        return

    all_warnings = prepare_warnings + metric_warnings
    for warning in all_warnings:
        st.warning(warning)

    st.subheader("Предпросмотр данных")
    preview_df = prepared_df.drop(columns=["sample"], errors="ignore")
    st.dataframe(preview_df, use_container_width=True)
    st.download_button(
        label="Скачать нормализованные исходные данные .xlsx",
        data=_dataframe_to_excel_bytes(prepared_df, "raw_normalized"),
        file_name="hepatotoxicity_raw_normalized.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key="download_hepato_raw_normalized",
    )
    col1, col2, col3 = st.columns(3)
    col1.metric("Строк после очистки", len(prepared_df))
    col2.metric("Препаратов", prepared_df["compound"].nunique())
    col3.metric("Концентраций", prepared_df["concentration"].nunique())

    st.subheader("Расчеты")
    stat_col1, stat_col2, stat_col3, stat_col4 = st.columns(4)
    alpha_norm = stat_col1.selectbox("α Shapiro", ALPHA_OPTIONS, index=0, key="hepato_alpha_norm")
    alpha_levene = stat_col2.selectbox("α Levene", ALPHA_OPTIONS, index=0, key="hepato_alpha_levene")
    alpha_global = stat_col3.selectbox("α ANOVA/Kruskal", ALPHA_OPTIONS, index=0, key="hepato_alpha_global")
    alpha_posthoc = stat_col4.selectbox("α post-hoc", ALPHA_OPTIONS, index=0, key="hepato_alpha_posthoc")

    summary_detailed, summary_journal, stats_tables = build_hepatotoxicity_summary_tables(
        metrics_df,
        alpha_norm=alpha_norm,
        alpha_levene=alpha_levene,
        alpha_global=alpha_global,
        alpha_posthoc=alpha_posthoc,
    )

    individual_metrics_table = build_individual_hepatotoxicity_table(metrics_df)
    summary_detailed_display = display_hepatotoxicity_summary_table(summary_detailed)
    individual_excel = _dataframe_to_excel_bytes(individual_metrics_table, "individual_metrics")
    detailed_excel = _dataframe_to_excel_bytes(summary_detailed_display, "summary_detailed")
    journal_excel = _dataframe_to_excel_bytes(summary_journal, "summary_journal")
    all_excel = export_hepatotoxicity_results(
        raw_normalized=prepared_df,
        individual_metrics=individual_metrics_table,
        summary_detailed=summary_detailed_display,
        summary_journal=summary_journal,
        stats_tables=stats_tables,
    )

    st.markdown("**Индивидуальные расчетные значения (% от контроля)**")
    st.dataframe(individual_metrics_table, use_container_width=True)
    st.download_button(
        label="Скачать индивидуальные расчетные значения .xlsx",
        data=individual_excel,
        file_name="hepatotoxicity_individual_metrics.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key="download_hepato_individual_metrics",
    )

    st.markdown("**Рабочая подробная таблица**")
    st.dataframe(summary_detailed_display, use_container_width=True)
    st.download_button(
        label="Скачать подробную таблицу .xlsx",
        data=detailed_excel,
        file_name="hepatotoxicity_summary_detailed.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key="download_hepato_summary_detailed",
    )

    st.markdown("**Таблица в журнальном стиле**")
    st.dataframe(summary_journal, use_container_width=True)
    st.download_button(
        label="Скачать журнальную таблицу .xlsx",
        data=journal_excel,
        file_name="hepatotoxicity_summary_journal.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key="download_hepato_summary_journal",
    )

    with st.expander("Статистический пайплайн ANOVA/Kruskal и post-hoc"):
        st.markdown("**Выбор теста и глобальный тест**")
        st.dataframe(stats_tables["decision_global"], use_container_width=True)
        st.download_button(
            label="Скачать выбор теста и глобальный тест .xlsx",
            data=_dataframe_to_excel_bytes(stats_tables["decision_global"], "decision_global"),
            file_name="hepatotoxicity_decision_global.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            key="download_hepato_decision_global",
        )
        st.markdown("**Post-hoc сводка**")
        st.dataframe(stats_tables["posthoc_summary"], use_container_width=True)
        st.download_button(
            label="Скачать post-hoc сводку .xlsx",
            data=_dataframe_to_excel_bytes(stats_tables["posthoc_summary"], "posthoc_summary"),
            file_name="hepatotoxicity_posthoc_summary.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            key="download_hepato_posthoc_summary",
        )

    st.download_button(
        label="Скачать все результаты одним Excel .xlsx",
        data=all_excel,
        file_name="hepatotoxicity_results.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key="download_hepato_all_results",
    )

    st.subheader("Графики")
    compounds = summary_detailed["compound"].drop_duplicates().tolist()

    with st.expander("Настройки подписей графиков", expanded=True):
        xlabel_common = st.text_input(
            "Общая подпись оси X",
            value="Concentration",
            key="hepato_plot_xlabel",
        )
        liver_ylabel = st.text_input(
            "Подпись оси Y для Liver size",
            value="Liver size (% of control)",
            key="hepato_liver_ylabel",
        )
        yolk_ylabel = st.text_input(
            "Подпись оси Y для Yolk sac size",
            value="Yolk sac size (% of control)",
            key="hepato_yolk_ylabel",
        )

        st.markdown("**Индивидуальные названия препаратов**")
        title_by_compound = {}
        for compound in compounds:
            slug = _safe_filename(compound)
            title_by_compound[compound] = st.text_input(
                f"Название для {compound}",
                value=str(compound).title(),
                key=f"hepato_plot_title_{slug}",
            )

    st.markdown("**Liver size (% of control)**")
    for row_start in range(0, len(compounds), 2):
        chart_cols = st.columns(2)
        for chart_col, compound in zip(chart_cols, compounds[row_start:row_start + 2]):
            with chart_col:
                fig = plot_liver_size(
                    summary_detailed,
                    compound=compound,
                    title=title_by_compound[compound],
                    xlabel=xlabel_common,
                    ylabel=liver_ylabel,
                )
                st.pyplot(fig, use_container_width=True)
                st.download_button(
                    label="Скачать PNG 600 dpi",
                    data=_fig_to_png_bytes(fig, dpi=600),
                    file_name=f"hepatotoxicity_liver_size_{_safe_filename(compound)}.png",
                    mime="image/png",
                    key=f"download_hepato_liver_{_safe_filename(compound)}",
                )
                plt.close(fig)

    st.markdown("**Yolk sac size (% of control)**")
    for row_start in range(0, len(compounds), 2):
        chart_cols = st.columns(2)
        for chart_col, compound in zip(chart_cols, compounds[row_start:row_start + 2]):
            with chart_col:
                fig = plot_yolk_sac_retention(
                    summary_detailed,
                    compound=compound,
                    title=title_by_compound[compound],
                    xlabel=xlabel_common,
                    ylabel=yolk_ylabel,
                )
                st.pyplot(fig, use_container_width=True)
                st.download_button(
                    label="Скачать PNG 600 dpi",
                    data=_fig_to_png_bytes(fig, dpi=600),
                    file_name=f"hepatotoxicity_yolk_sac_size_{_safe_filename(compound)}.png",
                    mime="image/png",
                    key=f"download_hepato_yolk_{_safe_filename(compound)}",
                )
                plt.close(fig)
