import io
import hashlib
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap


# === УТИЛИТЫ ===

def file_md5(file_bytes: bytes) -> str:
    return hashlib.md5(file_bytes).hexdigest()


@st.cache_data(show_spinner=False)
def parse_excel_to_df(file_bytes: bytes) -> pd.DataFrame:
    """Чтение Excel и сбор плоских заголовков как в исходном коде."""
    df_full = pd.read_excel(io.BytesIO(file_bytes), sheet_name=0, header=None)

    basic_columns = df_full.iloc[0, :8].tolist()
    metric = df_full.iloc[0, 8:]
    subtype = df_full.iloc[1, 8:]
    unit = df_full.iloc[2, 8:]
    extra = df_full.iloc[3, 8:]

    multi_columns = [
        f"{str(m).strip() if pd.notna(m) else ''} "
        f"({str(s).strip() if pd.notna(s) else ''}, "
        f"{str(u).strip() if pd.notna(u) else ''}, "
        f"{str(e).strip() if pd.notna(e) else ''})"
        for m, s, u, e in zip(metric, subtype, unit, extra)
    ]

    columns = basic_columns + multi_columns
    df_data = df_full.iloc[4:].copy().reset_index(drop=True)
    df_data.columns = columns
    return df_data


def show_table(df: pd.DataFrame, exposure_time: str, caption: str):
    st.markdown(f"**📄 {caption} — {exposure_time}**")
    st.dataframe(df, use_container_width=True)


def init_state_for_key(key: str):
    """Единоразовая инициализация служебных полей для конкретного времени экспозиции."""
    st.session_state.setdefault(f"ver_{key}", 0)              # версия мультиселекта
    st.session_state.setdefault(f"{key}_raw", None)           # оригинальный df
    st.session_state.setdefault(f"{key}_data", None)          # рабочая копия
    st.session_state.setdefault(f"{key}_file_hash", None)     # контроль смены файла


def get_well_col(df: pd.DataFrame) -> str | None:
    """Возвращает имя колонки с well_id (по договорённости — 3-я колонка)."""
    if df is None or df.empty or len(df.columns) < 3:
        return None
    return df.columns[2]


def compute_deletion_report(file_keys: list[str]) -> dict:
    """
    Считает для каждого файла:
      - исходное число строк и уникальных лунок (из RAW),
      - текущее число строк и уникальных лунок (из DATA),
      - удалённые строки и лунки как разницу,
      - список удалённых лунок с количеством строк.
    """
    report = {}
    for key in file_keys:
        raw_df = st.session_state.get(f"{key}_raw")
        cur_df = st.session_state.get(f"{key}_data")

        if raw_df is None or cur_df is None:
            report[key] = {
                "rows_raw": 0, "rows_cur": 0, "rows_removed": 0,
                "wells_raw": 0, "wells_cur": 0, "wells_removed": 0,
                "removed_wells_detail": {}
            }
            continue

        well_col = get_well_col(raw_df)
        rows_raw = len(raw_df)
        rows_cur = len(cur_df)
        rows_removed = max(rows_raw - rows_cur, 0)

        if well_col is None:
            wells_raw = wells_cur = wells_removed = 0
            removed_wells_detail = {}
        else:
            wells_raw = raw_df[well_col].astype(str).dropna().nunique()
            wells_cur = cur_df[well_col].astype(str).dropna().nunique()
            wells_removed = max(wells_raw - wells_cur, 0)

            # Определяем какие лунки были удалены
            wells_raw_set = set(raw_df[well_col].astype(str).dropna().unique())
            wells_cur_set = set(cur_df[well_col].astype(str).dropna().unique())
            removed_wells = wells_raw_set - wells_cur_set

            removed_wells_detail = {}
            for w in removed_wells:
                cnt = (raw_df[well_col].astype(str) == w).sum()
                removed_wells_detail[w] = int(cnt)

        report[key] = {
            "rows_raw": rows_raw,
            "rows_cur": rows_cur,
            "rows_removed": rows_removed,
            "wells_raw": wells_raw,
            "wells_cur": wells_cur,
            "wells_removed": wells_removed,
            "removed_wells_detail": removed_wells_detail,
        }
    return report


def render_sidebar_report(report: dict):
    """Рисует динамический отчёт в сайдбаре."""
    st.sidebar.header("🧾 Отчёт по удалению")

    total_rows_removed = sum(r["rows_removed"] for r in report.values())
    total_wells_removed = sum(r["wells_removed"] for r in report.values())

    st.sidebar.metric("Удалено строк (суммарно)", total_rows_removed)
    st.sidebar.metric("Удалено лунок (суммарно)", total_wells_removed)

    for key in ["1h", "4h", "12h", "24h"]:
        r = report.get(key, {})
        with st.sidebar.expander(f"{key}: детали", expanded=False):
            st.write(f"**Строки**: исходно {r.get('rows_raw',0)}, сейчас {r.get('rows_cur',0)}, удалено {r.get('rows_removed',0)}")
            st.write(f"**Лунки**: исходно {r.get('wells_raw',0)}, сейчас {r.get('wells_cur',0)}, удалено {r.get('wells_removed',0)}")

            if r.get("removed_wells_detail"):
                st.markdown("**Удалённые лунки:**")
                for w, cnt in r["removed_wells_detail"].items():
                    st.write(f"- {w}: {cnt} строк")


# === ПРИЛОЖЕНИЕ ===

def neurotoxicity_app():
    st.title("Анализ поведенческих данных DanioVision")

    tab_upload, tab_analysis, tab_agg = st.tabs(["📁 Загрузка и просмотр", "⚙️ Анализ данных", "📊 Графики"])
    file_keys = ["1h", "4h", "12h", "24h"]

    for key in file_keys:
        init_state_for_key(key)

    # =============================
    # 📁 Вкладка загрузки и просмотра
    # =============================
    with tab_upload:
        st.markdown("#### Загрузите 4 Excel-файла (1h, 4h, 12h, 24h)")

        uploaded_files = {}
        cols = st.columns(4)
        for i, key in enumerate(file_keys):
            with cols[i]:
                uploaded = st.file_uploader(f"Файл {key}", type="xlsx", key=f"uploader_{key}")
            if uploaded is not None:
                st.session_state[f"{key}_file"] = uploaded
            uploaded_files[key] = st.session_state.get(f"{key}_file")

        if all(uploaded_files.values()):
            st.success("Все 4 файла успешно загружены")

            for key, file in uploaded_files.items():
                file_bytes = file.getvalue()
                cur_hash = file_md5(file_bytes)

                if st.session_state[f"{key}_file_hash"] != cur_hash:
                    raw_df = parse_excel_to_df(file_bytes)
                    st.session_state[f"{key}_raw"] = raw_df
                    st.session_state[f"{key}_data"] = raw_df.copy()
                    st.session_state[f"{key}_file_hash"] = cur_hash
                    st.session_state[f"ver_{key}"] = 0

                col1, col2 = st.columns([3, 1])
                with col1:
                    show_table(st.session_state[f"{key}_data"], key, "Текущая рабочая копия")

                with col2:
                    st.markdown("##### ⚙️ Действия")
                    if st.button(f"↩️ Сбросить {key} к оригиналу", key=f"reset_{key}"):
                        st.session_state[f"{key}_data"] = st.session_state[f"{key}_raw"].copy()
                        st.session_state[f"ver_{key}"] += 1
                        st.success("Сброшено к оригиналу")
                        st.rerun()

                    towrite = io.BytesIO()
                    with pd.ExcelWriter(towrite, engine="openpyxl") as writer:
                        st.session_state[f"{key}_data"].to_excel(writer, index=False, sheet_name="Current")
                    towrite.seek(0)
                    st.download_button(
                        label=f"⬇️ Скачать {key} (текущая копия)",
                        data=towrite,
                        file_name=f"{key}_current.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        key=f"dl_current_{key}"
                    )
        else:
            st.warning("Пожалуйста, загрузите все 4 файла для продолжения.")

    # =================
    # ⚙️ Вкладка анализа
    # =================
    with tab_analysis:
        st.markdown("### ⚙️ Анализ по лункам с аномальным движением (Distance moved > 1000)")

        all_dfs_loaded = all(st.session_state.get(f"{key}_data") is not None for key in file_keys)
        if not all_dfs_loaded:
            st.warning("Сначала загрузите и обработайте таблицы на вкладке «Загрузка и просмотр».")
        else:
            for key in file_keys:
                df = st.session_state[f"{key}_data"].copy()

                if len(df.columns) < 3:
                    st.warning(f"Файл {key}: меньше 3 колонок — не могу определить well_id.")
                    continue
                well_col = df.columns[2]

                dist_cols = [col for col in df.columns if isinstance(col, str) and "Distance moved" in col]
                if not dist_cols:
                    st.warning(f"Файл {key}: не найдена колонка 'Distance moved'.")
                    continue
                dist_col = dist_cols[0]

                df[dist_col] = pd.to_numeric(df[dist_col], errors="coerce")

                high_wells = (
                    df[df[dist_col] > 1000][well_col]
                    .dropna()
                    .astype(str)
                    .unique()
                    .tolist()
                )

                st.markdown(f"### 📄 {key}")
                if not high_wells:
                    st.success("✅ Нет лунок с Distance moved > 1000")
                    st.dataframe(df, use_container_width=True)
                else:
                    st.error(f"⚠️ Обнаружены лунки с Distance moved > 1000: {high_wells}")

                    ms_key = f"delete_{key}_{st.session_state[f'ver_{key}']}"
                    selected_wells = st.multiselect(
                        f"Выберите лунки для удаления из {key}",
                        options=high_wells,
                        key=ms_key
                    )

                    if st.button(f"Удалить выбранные лунки из {key}", key=f"remove_btn_{key}") and selected_wells:
                        before = len(df)
                        new_df = df[~df[well_col].astype(str).isin([str(x) for x in selected_wells])].reset_index(drop=True)
                        st.session_state[f"{key}_data"] = new_df
                        st.session_state[f"ver_{key}"] += 1
                        removed = before - len(new_df)
                        st.success(f"✅ Удалено {removed} строк, связанных с {len(selected_wells)} лунками: {selected_wells}")
                        st.rerun()

                    st.dataframe(st.session_state[f"{key}_data"], use_container_width=True)

                towrite = io.BytesIO()
                with pd.ExcelWriter(towrite, engine="openpyxl") as writer:
                    st.session_state[f"{key}_data"].to_excel(writer, index=False, sheet_name="Filtered")
                towrite.seek(0)
                st.download_button(
                    label=f"⬇️ Скачать {key} как Excel",
                    data=towrite,
                    file_name=f"{key}_filtered.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    key=f"dl_filtered_{key}"
                )
    
    # ================
    # 📊 Вкладка агрегации
    # ================
    with tab_agg:
        st.markdown("### 📊 Агрегация Velocity по временным отрезкам и группам")

        key_selected = st.selectbox("Выберите файл:", file_keys, format_func=lambda x: f"{x}")

        df = st.session_state.get(f"{key_selected}_data")
        if df is None:
            st.warning("Сначала загрузите таблицу на вкладке «Загрузка и просмотр».")
        else:
            # === Определение ключевых колонок ===
            col_trial = df.columns[0]       # Trial
            col_exposure = df.columns[1]    # Exposure
            col_well = df.columns[2]        # well_id
            col_group = "Test/control"      # Подгруппа
            col_conc = "Concentration"      # Доза
            col_time = "Time"               # Время
            col_velocity = None             # Velocity колонка

            # Найдём колонку Velocity
            velocity_candidates = [col for col in df.columns if isinstance(col, str) and "Velocity" in col]
            if not velocity_candidates:
                st.error("Колонка с Velocity не найдена.")
                st.stop()
            col_velocity = velocity_candidates[0]

            # Проверим наличие всех нужных колонок
            missing_cols = [c for c in [col_trial, col_exposure, col_well, col_group, col_conc, col_time, col_velocity] if c not in df.columns]
            if missing_cols:
                st.error(f"Не найдены необходимые колонки: {missing_cols}")
                st.stop()

            # Приводим к нужным типам
            df[col_velocity] = pd.to_numeric(df[col_velocity], errors="coerce")
            df[col_time] = df[col_time].astype(str)
            df[col_conc] = df[col_conc].astype(str)

            # Заменяем NaN в Concentration для контроля на "0"
            df[col_conc] = df.apply(lambda row: "0" if row[col_group].strip().lower() == "control" else row[col_conc], axis=1)

            # Время → число минуты: Time = 00:00:00 → 2, 00:02:00 → 4, ...
            # Порядок должен сохраняться как в исходном файле
            time_order = df[col_time].dropna().unique().tolist()
            time_mapping = {t: str((i + 1) * 2) for i, t in enumerate(time_order)}
            df["TimeNumeric"] = df[col_time].map(time_mapping)

            # Группировка
            group_cols = ["TimeNumeric", col_group, col_conc]

            agg_df = (
                df.groupby(group_cols, sort=False)[col_velocity]
                .agg(['mean', 'std'])
                .reset_index()
                .rename(columns={
                    'mean': 'Velocity Mean',
                    'std': 'Velocity SD'
                })
            )

            st.dataframe(agg_df, use_container_width=True)

            # Скачивание
            towrite = io.BytesIO()
            with pd.ExcelWriter(towrite, engine="openpyxl") as writer:
                agg_df.to_excel(writer, index=False, sheet_name="Aggregated")
            towrite.seek(0)
            st.download_button(
                label="⬇️ Скачать агрегированные данные",
                data=towrite,
                file_name=f"{key_selected}_velocity_aggregated.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

            # === График ===
            st.markdown("### 📈 График: Velocity vs Time")

            agg_df["Time"] = pd.to_numeric(agg_df["TimeNumeric"], errors="coerce")
            agg_df["Concentration"] = pd.to_numeric(agg_df["Concentration"], errors="coerce")

            all_concentrations = agg_df["Concentration"].dropna().unique().tolist()
            all_concentrations.sort()

            selected_concs = st.multiselect(
                "Выберите концентрации для отображения (0 = контроль):",
                options=all_concentrations,
                default=all_concentrations
            )

            # === Настройки ===
            with st.expander("⚙️ Настройки графика", expanded=False):
                title = st.text_input("Заголовок графика", value=st.session_state.get("neurotoxicity_title", "Velocity vs Time"), key="key_neurotoxicity_title")
                st.session_state["neurotoxicity_title"] = title
                xlabel = st.text_input("Подпись оси X", value=st.session_state.get("neurotoxicity_xlabel", "Time (min)"), key="key_neurotoxicity_xlabel")
                st.session_state["neurotoxicity_xlabel"] = xlabel
                ylabel = st.text_input("Подпись оси Y", value=st.session_state.get("neurotoxicity_ylabel", "Velocity (mm/s)"), key="key_neurotoxicity_ylabel")
                st.session_state["neurotoxicity_ylabel"] = ylabel
                legend_title = st.text_input("Заголовок легенды", value=st.session_state.get("neurotoxicity_legend_title", "Concentration"), key="key_neurotoxicity_legend_title")
                st.session_state["neurotoxicity_legend_title"] = legend_title

                dark_label = st.text_input("Текст для тёмных сегментов", value=st.session_state.get("neurotoxicity_dark_label", "Dark"), key="key_neurotoxicity_dark_label")
                st.session_state["neurotoxicity_dark_label"] = dark_label

                light_label = st.text_input("Текст для светлых сегментов", value=st.session_state.get("neurotoxicity_light_label", "Light"), key="key_neurotoxicity_light_label")
                st.session_state["neurotoxicity_light_label"] = light_label

                legend_labels = {}
                st.markdown("Подписи для легенды:")
                for conc in selected_concs:
                    default_label = "Control" if float(conc) == 0 else f"{conc} µM"
                    legend_labels[conc] = st.text_input(
                        f"Концентрация {conc}",
                        value=st.session_state.get(f"neurotoxicity_legend_{conc}", default_label),
                        key=f"key_neurotoxicity_legend_{conc}"
                    )
                    st.session_state[f"neurotoxicity_legend_{conc}"] = legend_labels[conc]

                show_plot = st.checkbox("Показывать график", value=True)

            if show_plot:

                filtered_df = agg_df[agg_df["Concentration"].isin(selected_concs)]

                fig, ax = plt.subplots(figsize=(10, 6))

                # Палитра цветов
                cmap = get_cmap("Set2")
                color_list = cmap.colors
                color_cycle = color_list * (len(selected_concs) // len(color_list) + 1)

                for i, conc in enumerate(selected_concs):
                    df_line = filtered_df[filtered_df["Concentration"] == conc]
                    label = legend_labels.get(conc, str(conc))
                    ax.errorbar(
                        df_line["Time"],
                        df_line["Velocity Mean"],
                        yerr=df_line["Velocity SD"],
                        label=label,
                        marker='o',
                        capsize=3,
                        linewidth=1.5,
                        alpha=0.6,
                        color=color_cycle[i]
                    )

                ax.set_xlabel(xlabel, fontsize=12, labelpad=18)  # увеличенный отступ
                ax.set_ylabel(ylabel, fontsize=12)
                ax.set_title(title, fontsize=14)
                ax.legend(title=legend_title, fontsize=10, title_fontsize=11, loc='best')
                ax.get_legend().get_title().set_ha('left')  
                ax.grid(False)
                ax.set_xlim([0, 52])
                ax.set_xticks([0, 10, 20, 30, 40, 50])

                # Расчёт нижнего Y для подписей
                y_min, y_max = ax.get_ylim()
                text_y = y_min - 0.04 * (y_max - y_min)

                for i, label in enumerate([dark_label, light_label] * 3):
                    x_pos = i * 10 + 5
                    if x_pos > 50:
                        break
                    ax.text(
                        x_pos, text_y,
                        label,
                        ha='center',
                        va='top',
                        fontsize=10
                    )

                st.pyplot(fig)

                # Скачивание
                buf = io.BytesIO()
                fig.savefig(buf, format="png", dpi=600, bbox_inches="tight")
                st.download_button(
                    label="⬇️ Скачать график PNG (600 dpi)",
                    data=buf.getvalue(),
                    file_name=f"{key_selected}_velocity_plot.png",
                    mime="image/png"
                )


    # === САЙДБАР: динамический отчёт ===
    report = compute_deletion_report(file_keys)
    render_sidebar_report(report)


if __name__ == "__main__":
    neurotoxicity_app()
