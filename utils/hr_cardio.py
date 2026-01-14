import pandas as pd
from io import BytesIO
import numpy as np
import re
from scipy.stats import shapiro, levene, f_oneway, kruskal
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import scikit_posthocs as sp
from itertools import combinations
from statsmodels.stats.multitest import multipletests

#Загрузка и чтение Excel
def load_hr_analysis_sheet(
    file_like,
    sheet_name: str = "Analysis"
) -> pd.DataFrame:
    """
    Загружает лист Analysis с двухуровневым заголовком.
    """
    df = pd.read_excel(
        file_like,
        sheet_name=sheet_name,
        header=[0, 1]
    )
    return df

#Функция поиска колонки по префиксу
def find_column_by_prefix(
    columns: list[str],
    prefix: str
) -> str:
    """
    Ищет колонку, начинающуюся с prefix.
    """
    for col in columns:
        if col == prefix or col.startswith(prefix + "|"):
            return col
    raise ValueError(f"Column with prefix '{prefix}' not found")

#Приведение мульти-заголовков к плоскому виду
def flatten_multilevel_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Превращает multi-index колонки в плоские строки вида 'BPS|Mean'.
    """
    df = df.copy()

    df.columns = [
        f"{a}|{b}" if b and not str(b).startswith("Unnamed") else a
        for a, b in df.columns
    ]

    return df

#Извлечение нужных колонок (Compound, BPS Mean, BPM Mean)
def extract_hr_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Извлекает Compound, BPS Mean и BPM Mean,
    устойчиво к многоуровневым заголовкам.
    """
    cols = df.columns.tolist()

    compound_col = find_column_by_prefix(cols, "Compound")
    bps_col = find_column_by_prefix(cols, "BPS|Mean")
    bpm_col = find_column_by_prefix(cols, "BPM|Mean")

    out = df[[compound_col, bps_col, bpm_col]].copy()
    out.columns = ["Compound", "BPS Mean", "BPM Mean"]

    return out

#Очистка и приведение типов
def clean_hr_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Убирает пустые строки и приводит BPS/BPM к числам.
    """
    df = df.copy()

    df = df.dropna(subset=["Compound", "BPS Mean", "BPM Mean"])

    df["BPS Mean"] = pd.to_numeric(df["BPS Mean"], errors="coerce")
    df["BPM Mean"] = pd.to_numeric(df["BPM Mean"], errors="coerce")

    df = df.dropna(subset=["BPS Mean", "BPM Mean"])

    return df.reset_index(drop=True)

#Формирование Subject внутри каждой группы Compound
def assign_subjects_within_compound(
    df: pd.DataFrame,
    subject_col: str = "Subject"
) -> pd.DataFrame:
    """
    Формирует Subject (1, 2, 3, ...) внутри каждой группы Compound.
    """
    df = df.copy()

    df[subject_col] = (
        df
        .groupby("Compound")
        .cumcount()
        .add(1)
        .astype(str)
    )

    return df

#Финальная сборка таблицы ЧСС
def build_hr_table(
    file_like,
    sheet_name: str = "Analysis"
) -> pd.DataFrame:
    df = load_hr_analysis_sheet(file_like, sheet_name)
    df = flatten_multilevel_columns(df)
    df = extract_hr_columns(df)
    df = clean_hr_table(df)
    df = assign_subjects_within_compound(df)

    return df[["Subject", "Compound", "BPS Mean", "BPM Mean"]]

#Универсальная выгрузка в Excel
def dataframe_to_excel_bytes(
    df: pd.DataFrame,
    sheet_name: str = "Data"
) -> bytes:
    """
    Преобразует DataFrame в Excel (bytes).
    """
    buffer = BytesIO()
    with pd.ExcelWriter(buffer) as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)

    return buffer.getvalue()

#Чистая функция агрегации (основная)
def aggregate_hr_by_compound(
    df: pd.DataFrame
) -> pd.DataFrame:
    """
    Агрегация ЧСС по Compound:
    N, Mean, Standard Error отдельно для BPS и BPM.
    """
    grouped = df.groupby("Compound")

    agg = grouped.agg(
        N_BPS=("BPS Mean", "count"),
        Mean_BPS=("BPS Mean", "mean"),
        Std_BPS=("BPS Mean", "std"),

        N_BPM=("BPM Mean", "count"),
        Mean_BPM=("BPM Mean", "mean"),
        Std_BPM=("BPM Mean", "std"),
    ).reset_index()

    # Standard Error
    agg["SE_BPS"] = agg["Std_BPS"] / np.sqrt(agg["N_BPS"])
    agg["SE_BPM"] = agg["Std_BPM"] / np.sqrt(agg["N_BPM"])

    # Убираем std (не нужен в финальной таблице)
    agg = agg.drop(columns=["Std_BPS", "Std_BPM"])

    return agg

#Функция извлечения дозы из Compound
def extract_dose_from_compound(compound: str) -> float:
    """
    Control -> 0
    Propranolol 10 -> 10
    Propranolol 100 -> 100
    """
    if compound.lower() == "control":
        return 0.0

    match = re.search(r"(\d+(\.\d+)?)", compound)
    if match:
        return float(match.group(1))

    return np.nan

############## для сортировки от меньшей к большой дозе

def _dose_from_compound(compound: str) -> float:
    """
    Control -> 0
    Propranolol 10 -> 10
    Propranolol 100 -> 100
    """
    s = str(compound).strip()
    if s.lower() == "control":
        return 0.0

    m = re.search(r"(\d+(\.\d+)?)", s)
    return float(m.group(1)) if m else np.nan

def format_hr_summary_table(
    agg_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Формирует итоговую таблицу ЧСС
    с явными суффиксами BPS / BPM
    и корректной сортировкой по дозе.
    """
    out = agg_df.copy()

    # --- извлекаем дозу для сортировки ---
    out["_dose"] = out["Compound"].apply(extract_dose_from_compound)
    out = out.sort_values("_dose").reset_index(drop=True)

    # --- переименование колонок ---
    out = out.rename(columns={
        "N_BPS": "N_BPS",
        "Mean_BPS": "Mean_BPS",
        "SE_BPS": "Standard Error_BPS",

        "N_BPM": "N_BPM",
        "Mean_BPM": "Mean_BPM",
        "SE_BPM": "Standard Error_BPM",
    })

    # --- финальный порядок колонок ---
    out = out[
        [
            "Compound",
            "N_BPS", "Mean_BPS", "Standard Error_BPS",
            "N_BPM", "Mean_BPM", "Standard Error_BPM",
        ]
    ]

    return out


#####Статистика

def split_groups_for_metric(
    df: pd.DataFrame,
    metric: str,
    group_col: str = "Compound",
    min_n: int = 1,
) -> dict[str, np.ndarray]:
    """
    Возвращает {group: np.array(values)} для выбранной метрики.
    min_n — минимальный размер группы, чтобы она попала в результат.
    """
    groups: dict[str, np.ndarray] = {}
    for group, sub_df in df.groupby(group_col):
        vals = pd.to_numeric(sub_df[metric], errors="coerce").dropna().to_numpy()
        if vals.size >= min_n:
            groups[str(group)] = vals
    return groups

#Проверка на нормальность тестом Shapiro–Wilk

SHAPIRO_MIN_N = 3

def shapiro_by_group(df: pd.DataFrame, metric: str, alpha: float = 0.05) -> pd.DataFrame:
    """
    Shapiro–Wilk по каждой группе Compound для одной метрики.
    Возвращает таблицу для UI.
    """
    groups = split_groups_for_metric(df, metric)

    rows = []
    for compound, values in groups.items():
        n = int(values.size)

        if n < SHAPIRO_MIN_N:
            rows.append({
                "Compound": compound,
                "Metric": metric,
                "N": n,
                "Shapiro_W": np.nan,
                "p_value": np.nan,
                "Normal (p>α)": np.nan,
                "Note": f"n < {SHAPIRO_MIN_N}, Shapiro not applied",
            })
            continue

        w, p = shapiro(values)
        rows.append({
            "Compound": compound,
            "Metric": metric,
            "N": n,
            "Shapiro_W": float(w),
            "p_value": float(p),
            "Normal (p>α)": bool(p > alpha),
            "Note": "applied",
        })

    out = pd.DataFrame(rows)
    if not out.empty:
        out["_dose"] = out["Compound"].apply(_dose_from_compound)
        out = out.sort_values(["Metric", "_dose", "Compound"], na_position="last").drop(columns="_dose").reset_index(drop=True)

    return out


def build_normality_report_table(
    df: pd.DataFrame,
    metrics: list[str] | None = None,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """
    Общая таблица нормальности по всем метрикам (по умолчанию BPS/BPM).
    """
    if metrics is None:
        metrics = ["BPS Mean", "BPM Mean"]

    tables = [shapiro_by_group(df, m, alpha=alpha) for m in metrics]
    tables = [t for t in tables if not t.empty]

    if not tables:
        return pd.DataFrame(columns=[
            "Compound", "Metric", "N", "Shapiro_W", "p_value", "Normal (p>α)", "Note"
        ])

    return pd.concat(tables, ignore_index=True)

#Проверка на гомогенность дисперсий тестом Levene

LEVENE_MIN_N = 2          # минимально допустимое n в группе для дисперсии
LEVENE_MIN_GROUPS = 2     # нужно минимум 2 группы


def levene_by_metric(df: pd.DataFrame, metric: str, alpha: float = 0.05) -> pd.DataFrame:
    """
    Levene по всем группам Compound для одной метрики.
    Возвращает 1 строку:
    Metric | Levene_stat | p_value | Equal_variance | Note
    """
    groups = split_groups_for_metric(df, metric)

    # оставляем только группы, где хватает точек для дисперсии
    valid = {g: v for g, v in groups.items() if v.size >= LEVENE_MIN_N}

    # если групп мало — тест не применяем
    if len(valid) < LEVENE_MIN_GROUPS:
        return pd.DataFrame([{
            "Metric": metric,
            "Levene_stat": np.nan,
            "p_value": np.nan,
            "Equal_variance": np.nan,
            "Note": f"not applied: need >= {LEVENE_MIN_GROUPS} groups with n>={LEVENE_MIN_N}",
        }])

    stat, p = levene(*valid.values(), center="median")  # робастнее к выбросам
    return pd.DataFrame([{
        "Metric": metric,
        "Levene_stat": float(stat),
        "p_value": float(p),
        "Equal_variance": bool(p > alpha),
        "Note": "applied",
    }])


def build_levene_report_table(df: pd.DataFrame, metrics: list[str] | None = None, alpha: float = 0.05) -> pd.DataFrame:
    """
    Таблица Levene для UI:
    Metric | Levene_stat | p_value | Equal_variance | Note
    """
    if metrics is None:
        metrics = ["BPS Mean", "BPM Mean"]

    rows = [levene_by_metric(df, m, alpha=alpha) for m in metrics]
    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["Metric", "Levene_stat", "p_value", "Equal_variance", "Note"]
    )
    
    # 🔹 ЯВНЫЙ порядок метрик (а не алфавит)
    order_map = {m: i for i, m in enumerate(metrics)}
    out["_order"] = out["Metric"].map(order_map)

    out = (
        out
        .sort_values("_order")
        .drop(columns="_order")
        .reset_index(drop=True)
    )

    return out

#Выбор глобального теста (ANOVA/Kruskal-Wallis)

GLOBAL_MIN_GROUPS = 2
GLOBAL_MIN_N_PER_GROUP = 2

def _decide_selected_test_for_metric(
    normality_df: pd.DataFrame,
    levene_df: pd.DataFrame,
    metric: str,
) -> tuple[bool, bool, str, str]:
    """
    Возвращает:
      normality_ok, equal_var_ok, selected_test, note
    """
    # --- нормальность: все группы должны быть Normal (p>a) ---
    norm_subset = normality_df[normality_df["Metric"] == metric]
    if norm_subset.empty:
        normality_ok = False
        norm_note = "normality not evaluated"
    else:
        normality_ok = bool(norm_subset["Normal (p>α)"].all())
        norm_note = None

    # --- Levene ---
    lev_row = levene_df[levene_df["Metric"] == metric]
    if lev_row.empty:
        equal_var_ok = False
        lev_note = "levene not evaluated"
    else:
        equal_var_ok = bool(lev_row.iloc[0]["Equal_variance"])
        lev_note = None

    # --- выбор теста ---
    selected = "ANOVA" if (normality_ok and equal_var_ok) else "Kruskal-Wallis"

    notes = []
    if norm_note:
        notes.append(norm_note)
    if lev_note:
        notes.append(lev_note)

    note = "applied" if not notes else "; ".join(notes)

    return normality_ok, equal_var_ok, selected, note


def _run_global_test(metric: str, selected_test: str, groups: dict[str, np.ndarray]) -> tuple[float, float, str]:
    """
    Возвращает: statistic, p_value, note_suffix
    """
    st = (selected_test or "").strip().lower()

    if st == "anova":
        stat, p = f_oneway(*groups.values())
        return float(stat), float(p), "applied (F-statistic)"

    if st in ("kruskal-wallis", "kruskal wallis", "kruskal–wallis", "kruskal"):
        stat, p = kruskal(*groups.values())
        return float(stat), float(p), "applied (H-statistic)"

    return np.nan, np.nan, f"not applied: unknown Selected_test='{selected_test}'"


def build_decision_and_global_table(
    hr_df: pd.DataFrame,
    normality_df: pd.DataFrame,
    levene_df: pd.DataFrame,
    metrics: list[str] | None = None,
    alpha: float = 0.05,   # ← ВАЖНО
) -> pd.DataFrame:
    """
    ЕДИНАЯ таблица для UI/Excel:
    Metric | Normality_OK | Equal_variance | Selected_test | Statistic | p_value | Reject_H0 | Note
    """
    if metrics is None:
        metrics = ["BPS Mean", "BPM Mean"]

    rows = []

    for metric in metrics:
        normality_ok, equal_var_ok, selected_test, decision_note = _decide_selected_test_for_metric(
            normality_df=normality_df,
            levene_df=levene_df,
            metric=metric,
        )

        groups_all = split_groups_for_metric(hr_df, metric)
        groups = {g: v for g, v in groups_all.items() if v.size >= GLOBAL_MIN_N_PER_GROUP}

        if len(groups) < GLOBAL_MIN_GROUPS:
            stat, p = np.nan, np.nan
            global_note = f"not applied: need >= {GLOBAL_MIN_GROUPS} groups with n>={GLOBAL_MIN_N_PER_GROUP}"
        else:
            stat, p, global_note = _run_global_test(metric, selected_test, groups)

        # NOTE: если decision_note уже "not evaluated", а глобальный тест всё равно пошёл,
        # то decision_note сохраняем (чтобы было видно, что проверка предпосылок не проводилась)
        # А глобальный note показывает, что именно посчитали.
        note = global_note if decision_note == "applied" else f"{decision_note}; {global_note}"

        reject_h0 = bool(np.isfinite(p) and p < alpha)

        rows.append({
            "Metric": metric,
            "Normality_OK": normality_ok,
            "Equal_variance": equal_var_ok,
            "Selected_test": selected_test,
            "Statistic": stat,
            "p_value": p,
            "Reject_H0": reject_h0,
            "Note": note,
        })

    out = pd.DataFrame(rows)

    # порядок строк ровно как metrics (чтобы BPS/BPM не прыгали)
    order_map = {m: i for i, m in enumerate(metrics)}
    out["_order"] = out["Metric"].map(order_map)
    out = out.sort_values("_order").drop(columns="_order").reset_index(drop=True)

    return out

#Post-hoc тесты

def _ensure_control_in_group1_vs_control(
    df: pd.DataFrame,
    control_label: str = "Control",
) -> pd.DataFrame:
    """
    Делает так, чтобы во всех строках:
    Group 1 == Control, Group 2 == НЕ Control.
    Если пара была (Drug, Control) — разворачиваем.
    Для Tukey также корректно переворачиваем meandiff и CI.
    """
    out = df.copy()
    c = _norm_compound_name(control_label)

    g1 = out["Group 1"].astype(str).str.strip().str.lower()
    g2 = out["Group 2"].astype(str).str.strip().str.lower()

    swap_mask = (g1 != c) & (g2 == c)
    if not swap_mask.any():
        return out

    # swap Group 1/2
    out.loc[swap_mask, ["Group 1", "Group 2"]] = out.loc[swap_mask, ["Group 2", "Group 1"]].to_numpy()

    # если это Tukey-таблица — переворачиваем знак эффекта и CI
    if "meandiff" in out.columns:
        out.loc[swap_mask, "meandiff"] = -pd.to_numeric(out.loc[swap_mask, "meandiff"], errors="coerce")

    if "lower" in out.columns and "upper" in out.columns:
        lower = pd.to_numeric(out.loc[swap_mask, "lower"], errors="coerce").to_numpy()
        upper = pd.to_numeric(out.loc[swap_mask, "upper"], errors="coerce").to_numpy()
        out.loc[swap_mask, "lower"] = -upper
        out.loc[swap_mask, "upper"] = -lower

    return out


def _sort_vs_control_by_dose(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Сортирует строки post-hoc (vs control) по дозе Group 2.
    """
    out = df.copy()
    out["_dose2"] = out["Group 2"].apply(_dose_from_compound)
    out = out.sort_values(["_dose2", "Group 2"], na_position="last").drop(columns=["_dose2"]).reset_index(drop=True)
    return out

########################

def _norm_compound_name(x: str) -> str:
    return str(x).strip().lower()


def _filter_pairs_vs_control(df: pd.DataFrame, col1: str, col2: str, control_label: str = "Control") -> pd.DataFrame:
    """
    Оставляет только строки, где одна из групп = Control, а другая != Control.
    """
    c = _norm_compound_name(control_label)

    g1 = df[col1].astype(str).str.strip().str.lower()
    g2 = df[col2].astype(str).str.strip().str.lower()

    mask = ((g1 == c) & (g2 != c)) | ((g2 == c) & (g1 != c))
    return df.loc[mask].reset_index(drop=True)


def run_tukey_posthoc_table(
    df: pd.DataFrame,
    metric: str,
    alpha: float=0.05,
    control_label: str = "Control",
    only_vs_control: bool = True,
) -> pd.DataFrame:
    """
    Tukey HSD (all-pairs), затем опционально фильтруем сравнения vs Control.
    Возвращает "tidy" таблицу.
    """
    sub = df[["Compound", metric]].copy()
    sub[metric] = pd.to_numeric(sub[metric], errors="coerce")
    sub = sub.dropna(subset=["Compound", metric])

    tukey = pairwise_tukeyhsd(
        endog=sub[metric].to_numpy(),
        groups=sub["Compound"].astype(str).to_numpy(),
        alpha=alpha,
    )
    
    #берем сырые значения теста без округлений
    groups = list(tukey.groupsunique)
    pairs = list(combinations(range(len(groups)), 2))  # i<j

    rows = []
    for k, (i, j) in enumerate(pairs):
        g1 = groups[i]
        g2 = groups[j]

        rows.append({
            "Group 1": g1,
            "Group 2": g2,
            "meandiff": float(tukey.meandiffs[k]),
            "p_value_adj": float(tukey.pvalues[k]),
            "lower": float(tukey.confint[k][0]),
            "upper": float(tukey.confint[k][1]),
            "Reject_H0": bool(tukey.reject[k]),
        })

    out = pd.DataFrame(rows)

    # унифицируем имена колонок
    # statsmodels обычно: group1 group2 meandiff p-adj lower upper reject
    rename_map = {
        "group1": "Group 1",
        "group2": "Group 2",
        "p-adj": "p_value_adj",
        "reject": "Reject_H0",
    }
    out = out.rename(columns=rename_map)

    # типы
    for c in ["meandiff", "lower", "upper", "p_value_adj"]:
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce")

    if "Reject_H0" in out.columns:
        out["Reject_H0"] = out["Reject_H0"].astype(str).str.lower().isin(["true", "1", "yes"])

    if only_vs_control:
        out = _filter_pairs_vs_control(out, "Group 1", "Group 2", control_label=control_label)
        out = _ensure_control_in_group1_vs_control(out, control_label=control_label)
        out = _sort_vs_control_by_dose(out)

    # удобный порядок
    cols = [c for c in ["Group 1", "Group 2", "meandiff", "p_value_adj", "lower", "upper", "Reject_H0"] if c in out.columns]
    out = out[cols].reset_index(drop=True)

    return out


def run_dunn_posthoc_table(
    df: pd.DataFrame,
    metric: str,
    p_adjust: str = "holm",
    control_label: str = "Control",
    only_vs_control: bool = True,
    alpha: float=0.05,
) -> pd.DataFrame:
    """
    Dunn post-hoc (после Kruskal), p_adjust по умолчанию holm.
    Возвращает tidy таблицу: Group 1 | Group 2 | p_value_adj | Reject_H0
    """
    sub = df[["Compound", metric]].copy()
    sub[metric] = pd.to_numeric(sub[metric], errors="coerce")
    sub = sub.dropna(subset=["Compound", metric])

    # 1) сначала считаем НЕскорректированные p-values по Dunn
    mat = sp.posthoc_dunn(
        sub,
        val_col=metric,
        group_col="Compound",
        p_adjust=None,   # <-- ВАЖНО: без коррекции здесь
    )

    out = (
        mat.reset_index()
        .melt(id_vars="index", var_name="Group 2", value_name="p_value")  # <-- сырое p
        .rename(columns={"index": "Group 1"})
    )

    out = out[out["Group 1"] != out["Group 2"]].reset_index(drop=True)

    # 2) фильтруем пары vs Control (как у тебя)
    if only_vs_control:
        out = _filter_pairs_vs_control(out, "Group 1", "Group 2", control_label=control_label)
        out = _ensure_control_in_group1_vs_control(out, control_label=control_label)
        out = out.drop_duplicates(subset=["Group 1", "Group 2"], keep="first").reset_index(drop=True)
        out = _sort_vs_control_by_dose(out)

    out["p_value"] = pd.to_numeric(out["p_value"], errors="coerce")

    # 3) Holm ТОЛЬКО на comparisons vs control
    mask = out["p_value"].notna()
    pvals = out.loc[mask, "p_value"].to_numpy()

    if len(pvals) > 0:
        p_adj = multipletests(pvals, alpha=alpha, method="holm")[1]
        out.loc[mask, "p_value_adj"] = p_adj
    else:
        out["p_value_adj"] = pd.NA

    out["p_value_adj"] = pd.to_numeric(out["p_value_adj"], errors="coerce")

    # ✅ ключевая строка остаётся, но уже по p_value_adj после "vs control" Holm
    out["Reject_H0"] = out["p_value_adj"].lt(alpha)

    cols = [c for c in ["Group 1", "Group 2", "p_value_adj", "Reject_H0"] if c in out.columns]
    out = out[cols].reset_index(drop=True)
    return out


def build_posthoc_reports(
    hr_df: pd.DataFrame,
    decision_global_df: pd.DataFrame,
    metrics: list[str] | None = None,
    alpha_posthoc: float=0.05, 
    alpha_global: float=0.05,
    control_label: str = "Control",
    only_vs_control: bool = True,
    dunn_p_adjust: str = "holm",
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """
    Возвращает:
    1) summary-таблицу для UI: Metric | Selected_test | Global_p_value | Posthoc_test | Posthoc_applied | Note
    2) словарь {metric: posthoc_df} для детального вывода/экспорта
    """
    if metrics is None:
        metrics = ["BPS Mean", "BPM Mean"]

    rows = []
    details: dict[str, pd.DataFrame] = {}

    for metric in metrics:
        row = decision_global_df[decision_global_df["Metric"] == metric]
        if row.empty:
            rows.append({
                "Metric": metric,
                "Selected_test": None,
                "Global_p_value": np.nan,
                "Posthoc_test": None,
                "Posthoc_applied": False,
                "Note": "not applied: metric missing in decision/global table",
            })
            continue

        selected_test = str(row.iloc[0]["Selected_test"])
        gp = row.iloc[0].get("p_value", np.nan)
        gp = float(gp) if pd.notna(gp) else np.nan

        if not np.isfinite(gp):
            rows.append({
                "Metric": metric,
                "Selected_test": selected_test,
                "Global_p_value": gp,
                "Posthoc_test": None,
                "Posthoc_applied": False,
                "Note": "not applied: global p_value is NaN",
            })
            continue

        if gp >= alpha_global:
            rows.append({
                "Metric": metric,
                "Selected_test": selected_test,
                "Global_p_value": gp,
                "Posthoc_test": None,
                "Posthoc_applied": False,
                "Note": f"not applied: global p_value >= {alpha_global}",
            })
            continue

        st = selected_test.strip().lower()
        if st == "anova":
            posthoc_df = run_tukey_posthoc_table(
                hr_df, metric=metric, alpha=alpha_posthoc,
                control_label=control_label, only_vs_control=only_vs_control
            )
            posthoc_name = "Tukey HSD"
        else:
            posthoc_df = run_dunn_posthoc_table(
                hr_df,
                metric=metric,
                p_adjust=dunn_p_adjust,
                control_label=control_label,
                only_vs_control=only_vs_control,
                alpha=alpha_posthoc,
            )
            posthoc_name = f"Dunn ({dunn_p_adjust})"

        details[metric] = posthoc_df

        rows.append({
            "Metric": metric,
            "Selected_test": selected_test,
            "Global_p_value": gp,
            "Posthoc_test": posthoc_name,
            "Posthoc_applied": True,
            "Note": "applied",
        })

    summary = pd.DataFrame(rows)

    # порядок метрик как задано
    order_map = {m: i for i, m in enumerate(metrics)}
    summary["_order"] = summary["Metric"].map(order_map)
    summary = summary.sort_values("_order").drop(columns="_order").reset_index(drop=True)

    return summary, details


##############Графики
import matplotlib.pyplot as plt

def p_to_sig(p: float | None) -> str:
    """p -> stars/ns по твоим правилам."""
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"

def _drug_from_compound(compound: str) -> str | None:
    """
    Control -> None
    Propranolol 10 -> Propranolol
    """
    s = str(compound).strip()
    if s.lower() == "control":
        return None

    # берём всё до первого числа как "название вещества"
    m = re.match(r"^\s*([^\d]+?)\s*(\d+(\.\d+)?)\s*$", s)
    if m:
        return m.group(1).strip()

    # fallback: если формат неожиданный
    parts = s.split()
    return parts[0].strip() if parts else None


def _dose_label_from_compound(compound: str, control_as: str = "Control") -> str:
    """
    Control -> Control (или "0", если захочешь)
    Propranolol 10 -> 10
    """
    s = str(compound).strip()
    if s.lower() == "control":
        return control_as

    m = re.search(r"(\d+(\.\d+)?)", s)
    return m.group(1) if m else s

def _compound_order_by_dose(compounds: list[str]) -> list[str]:
    # Control -> 0, дальше по дозе
    tmp = pd.DataFrame({"Compound": list(compounds)})
    tmp["_dose"] = tmp["Compound"].apply(_dose_from_compound)  # у тебя уже есть
    tmp = tmp.sort_values(["_dose", "Compound"], na_position="last")
    return tmp["Compound"].tolist()

def build_hr_boxplot(
    hr_df: pd.DataFrame,
    metric: str,
    title: str | None = None,
    ylabel: str | None = None,
    xlabel: str | None = None,
    show_points: bool = True,
    points_jitter: float = 0.08,
    points_alpha: float = 0.6,
    figsize: tuple[float, float] = (10, 4.5),
    control_tick: str = "Control",  # <-- можно поменять на "0"
    # --- NEW ---
    signif_by_group: dict[str, str] | None = None,  # ключ = full group name (как в df["Compound"])
    show_signif: bool = True,
    signif_fontsize: int = 12,
    signif_y_pad_frac: float = 0.06,
) -> plt.Figure:
    """
    Строит boxplot по группам Compound для выбранной метрики.
    Границы боксов и медиана — чёрные, заливка боксов — разная (по дефолтному циклу matplotlib).
    """

    df = hr_df[["Compound", metric]].copy()
    df[metric] = pd.to_numeric(df[metric], errors="coerce")
    df = df.dropna(subset=["Compound", metric])

    # порядок групп
    groups = _compound_order_by_dose(sorted(df["Compound"].astype(str).unique().tolist()))
    data = [df.loc[df["Compound"].astype(str) == g, metric].to_numpy() for g in groups]

    # ✅ подписи тиков (только дозы)
    tick_labels = [_dose_label_from_compound(g, control_as=control_tick) for g in groups]

    # ✅ авто xlabel из названия вещества (если пусто)
    xlabel_clean = (xlabel or "").strip()
    if not xlabel_clean:
        drugs = [_drug_from_compound(g) for g in groups]
        drugs = [d for d in drugs if d]  # выкинули None/пустые
        xlabel_clean = drugs[0] if drugs else "Compound"

    fig, ax = plt.subplots(figsize=figsize)

    bp = ax.boxplot(
        data,
        labels=tick_labels,    # <-- ВАЖНО: показываем дозы, а не full Compound
        patch_artist=True,   # чтобы можно было красить заливку
        widths=0.55,
        showfliers=False,
        medianprops=dict(color="black", linewidth=2.0),
        boxprops=dict(edgecolor="black", linewidth=1.8),
        whiskerprops=dict(color="black", linewidth=1.6),
        capprops=dict(color="black", linewidth=1.6),
    )

    # заливка: разные цвета, но без “хардкода” — берём дефолтный цикл matplotlib
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]

    for i, patch in enumerate(bp["boxes"]):
        color = cycle[i % len(cycle)]
        patch.set_facecolor(color)
        patch.set_alpha(0.55)                # прозрачность ЗАЛИВКИ
        patch.set_edgecolor("black")          # граница
        patch.set_linewidth(1.8)

    # точки (по желанию) — дают ощущение n
    if show_points:
        rng = np.random.default_rng(0)
        for i, y in enumerate(data, start=1):
            if len(y) == 0:
                continue
            x = rng.normal(loc=i, scale=points_jitter, size=len(y))
            ax.scatter(x, y, s=28, alpha=points_alpha, edgecolors="none")

    # --- NEW: significance annotations ---
    if show_signif and signif_by_group:
        # общий y-уровень: чуть выше максимума по данным
        all_y = np.concatenate([y for y in data if len(y) > 0]) if any(len(y) > 0 for y in data) else np.array([0.0])
        y_max = float(np.nanmax(all_y))
        y_min = float(np.nanmin(all_y))
        yr = max(1e-9, (y_max - y_min))
        y_text = y_max + signif_y_pad_frac * yr

        # чтобы текст не обрезало
        ax.set_ylim(top=y_text + 0.10 * yr)

        for i, g in enumerate(groups, start=1):
            sig = signif_by_group.get(str(g))
            if not sig:
                continue
            ax.text(
                i,
                y_text,
                sig,
                ha="center",
                va="bottom",
                fontsize=signif_fontsize,
                fontweight="bold",
            )

    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel_clean)

    ax.tick_params(axis="x", rotation=0)

    fig.tight_layout()
    return fig