import matplotlib.pyplot as plt
import plotly.graph_objects as go
import pandas as pd
import re

##############################
#Индивидуальные графики
def get_available_subjects(df: pd.DataFrame):
    subjects = []

    for col in df.columns:
        if col.startswith("Time_"):
            subjects.append(col.replace("Time_", ""))

    subjects = list(set(subjects))

    def sort_key(s: str):
        s_lower = s.lower()

        # --- replicate (_1, _2, ...)
        rep_match = re.search(r"_(\d+)$", s_lower)
        replicate = int(rep_match.group(1)) if rep_match else 0

        # --- control
        if "control" in s_lower:
            return (0, 0, replicate)

        # --- dose
        dose_match = re.search(r"(\d+(\.\d+)?)", s_lower)
        if dose_match:
            dose = float(dose_match.group(1))
            return (1, dose, replicate)

        # --- fallback
        return (2, float("inf"), replicate)

    return sorted(subjects, key=sort_key)

def pretty_subject_label(
    suffix: str,
    control_name: str = "Control",
    subject_word: str = "Subject",
) -> str:
    """
    control_2 -> "Control (Subject №2)" (control_name может меняться)
    100_1     -> "100 (Subject №1)"
    """
    s = str(suffix).strip()

    m = re.match(r"^(.*)_(\d+)$", s)
    if m:
        base = m.group(1).strip()
        idx = int(m.group(2))
    else:
        base, idx = s, None

    is_control = "control" in base.lower()

    if idx is None:
        # fallback, если нет _n
        return control_name if is_control else base

    if is_control:
        return f"{control_name} ({subject_word} №{idx})"
    return f"{base} ({subject_word} №{idx})"

#статичный график
def plot_activity_curves(
    df: pd.DataFrame,
    selected_subjects: list[str],
    title: str = "Compound",
    xlabel: str = "Time (s)",
    ylabel: str = "Activity (%)",
    control_name: str = "Control",
    subject_word: str = "Subject",
):
    """
    Строит график Activity (%) vs Time (s) для всех валидных пар
    """
    fig, ax = plt.subplots(figsize=(10, 4))

    for suffix in selected_subjects:
        time_col = f"Time_{suffix}"
        percent_col = f"Percent_{suffix}"

        if time_col not in df.columns or percent_col not in df.columns:
            continue

        sub_df = df[[time_col, percent_col]].copy()
        sub_df = sub_df.dropna(how="all")
        if sub_df.empty:
            continue

        sub_df = sub_df[sub_df[percent_col] != "-"]
        if sub_df.empty:
            continue

        x = pd.to_numeric(sub_df[time_col], errors="coerce")
        y = pd.to_numeric(sub_df[percent_col], errors="coerce")

        mask = (~x.isna()) & (~y.isna())
        if mask.sum() == 0:
            continue

        ax.plot(
            x[mask],
            y[mask],
            linewidth=1,
            label = pretty_subject_label(
            suffix,
            control_name=control_name,
            subject_word=subject_word,
        ),
            alpha=0.6,
        )

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=8, ncol=1,loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    )

    fig.tight_layout()

    return fig

#динамический график
def plot_activity_curves_plotly(
    df: pd.DataFrame,
    selected_subjects: list[str],
    title: str = "Compound",
    xlabel: str = "Time (s)",
    ylabel: str = "Activity (%)",
    control_name: str = "Control",
    subject_word: str = "Subject",
):
    """
    Plotly-версия: Activity (%) vs Time (s) для выбранных пар Time_/Percent_.
    """
    fig = go.Figure()

    for suffix in selected_subjects:
        time_col = f"Time_{suffix}"
        percent_col = f"Percent_{suffix}"

        if time_col not in df.columns or percent_col not in df.columns:
            continue

        sub_df = df[[time_col, percent_col]].copy()
        sub_df = sub_df.dropna(how="all")
        if sub_df.empty:
            continue

        sub_df = sub_df[sub_df[percent_col] != "-"]
        if sub_df.empty:
            continue

        x = pd.to_numeric(sub_df[time_col], errors="coerce")
        y = pd.to_numeric(sub_df[percent_col], errors="coerce")

        mask = (~x.isna()) & (~y.isna())
        if mask.sum() == 0:
            continue

        label = pretty_subject_label(
            suffix,
            control_name=control_name,
            subject_word=subject_word,
        )

        fig.add_trace(
            go.Scatter(
                x=x[mask],
                y=y[mask],
                mode="lines",
                name=label,
            )
        )

    fig.update_layout(
        title=title,
        title_x=0.45,          # ← центрирование
        xaxis_title=xlabel,
        yaxis_title=ylabel,
        height=420,
        margin=dict(l=40, r=20, t=50, b=40),
    )

    return fig

#######################################

#######################################
#общие функции
def extract_base_name(sheet_name: str) -> str:
    """
    Берём текст после "_" и обрезаем по первому пробелу
    """
    if "_" not in sheet_name:
        return sheet_name.strip()

    after_underscore = sheet_name.split("_", 1)[1]
    return after_underscore.split(" ", 1)[0].strip()

def extract_time_percent(df: pd.DataFrame):
    """
    Ищет строку, где есть 'Time (s)' и '%', и парсит данные ниже
    """
    df = df.copy()

    header_row_idx = None

    for i in range(len(df)):
        row = df.iloc[i].astype(str)
        if "Time (s)" in row.values and "%" in row.values:
            header_row_idx = i
            break

    if header_row_idx is None:
        return None

    # Получаем индексы нужных колонок
    header_row = df.iloc[header_row_idx]
    time_col = header_row[header_row == "Time (s)"].index[0]
    percent_col = header_row[header_row == "%"].index[0]

    # Берём данные ниже заголовков
    out = df.loc[header_row_idx + 1:, [time_col, percent_col]].copy()
    out.columns = ["Time (s)", "%"]

    # Убираем пустые строки
    out = out.dropna(how="all")

    return out.reset_index(drop=True)