import io
import hashlib
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import numpy as np
import plotly.graph_objects as go


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
                def well_sort_key(well: str):
                    row = well[0]
                    col = int(well[1:])
                    return (row, col)

                for w, cnt in sorted(r["removed_wells_detail"].items(), key=lambda x: well_sort_key(x[0])):
                    st.write(f"- {w}: {cnt} строк")


# === ПРИЛОЖЕНИЕ ===

def neurotoxicity_app():
    st.title("Анализ поведенческих данных DanioVision")

    tab_upload, tab_analysis, tab_agg, tab_calc = st.tabs(["📁 Загрузка и просмотр", "⚙️ Анализ данных", "📊 Графики","📐 Расчёты"])
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
        st.markdown("### ⚙️ Анализ по лункам с Distance moved > 1000")

        all_dfs_loaded = all(st.session_state.get(f"{key}_data") is not None for key in file_keys)
        if not all_dfs_loaded:
            st.warning("Сначала загрузите и обработайте таблицы на вкладке «Загрузка и просмотр».")
        else:
            selected_key = st.selectbox("Выберите временной интервал:", file_keys)

            df = st.session_state.get(f"{selected_key}_data")
            raw_df = st.session_state.get(f"{selected_key}_raw")
            if df is None or raw_df is None:
                st.warning(f"Файл {selected_key} не загружен.")
                st.stop()

            well_col = df.columns[2]
            dist_cols = [col for col in df.columns if isinstance(col, str) and "Distance moved" in col]
            if not dist_cols:
                st.warning(f"Файл {selected_key}: не найдена колонка 'Distance moved'.")
                st.stop()
            dist_col = dist_cols[0]

            dist_numeric = pd.to_numeric(df[dist_col], errors="coerce")

            wells_all = [f"{r}{c}" for r in "ABCDEFGH" for c in range(1, 13)]
            wells_present = raw_df[well_col].astype(str).dropna().unique().tolist()
            wells_current = df[well_col].astype(str).dropna().unique().tolist()
            wells_high = df[dist_numeric > 1000][well_col].astype(str).dropna().unique().tolist()
            removed_wells = list(set(wells_present) - set(wells_current))

            def well_sort_key(well: str):
                row = well[0]
                col = int(well[1:])
                return (row, col)

            distance_start_idx = df.columns.get_loc(next(col for col in df.columns if isinstance(col, str) and "Distance moved" in col))
            data_cols = df.columns[distance_start_idx:]
            dash_mask = df[data_cols].astype(str).applymap(lambda x: x.strip() == "-")
            wells_with_dash = df.loc[dash_mask.any(axis=1), well_col].astype(str).unique().tolist()

            well_colors = []
            for w in wells_all:
                if w in removed_wells:
                    color = "red"
                elif w in wells_high and w in wells_with_dash:
                    color = "purple"
                elif w in wells_with_dash:
                    color = "orange"
                elif w in wells_high:
                    color = "yellow"
                elif w in wells_current:
                    color = "green"
                else:
                    color = "lightgray"
                well_colors.append(color)

            x_vals = [int(w[1:]) for w in wells_all]
            y_vals = [8 - "ABCDEFGH".index(w[0]) for w in wells_all]

            fig = go.Figure(data=go.Scatter(
                x=x_vals,
                y=y_vals,
                mode='markers+text',
                text=wells_all,
                textposition='middle center',
                marker=dict(size=35, color=well_colors, line=dict(color='black', width=1)),
                textfont=dict(color='black'),
                hoverinfo='text'
            ))

            fig.update_layout(
                title=f"{selected_key}: интерактивный планшет",
                width=650,
                height=500,
                xaxis=dict(title=None, range=[0.5, 12.5], tickvals=list(range(1, 13))),
                yaxis=dict(title=None, range=[0.5, 8.5], tickvals=list(range(1, 9)), ticktext=list("HGFEDCBA")),
                margin=dict(l=20, r=20, t=40, b=20)
            )

            with st.expander("🧪 Планшет и удаление лунок", expanded=True):
                st.caption("🟢 валидная, 🟡 > 1000, 🟠 содержит '-', 🟣 и > 1000, и '-', 🔴 удалена, ⚪ нет данных")

                st.plotly_chart(fig, use_container_width=False, config={"displayModeBar": False, "scrollZoom": False, "staticPlot": True})

                col_yellow, col_orange, col_red = st.columns(3)

                # === ЖЁЛТЫЕ ===
                with col_yellow:
                    st.markdown("#### 🟡 Удаление лунок > 1000")
                    if wells_high:
                        well_to_remove = st.selectbox("Одиночный выбор:", sorted(wells_high, key=well_sort_key), key=f"select_remove_{selected_key}")
                        if st.button("Удалить выбранную", key=f"remove_btn_{selected_key}"):
                            df_new = df[df[well_col].astype(str) != well_to_remove].reset_index(drop=True)
                            st.session_state[f"{selected_key}_data"] = df_new
                            st.session_state[f"ver_{selected_key}"] += 1
                            st.success(f"Удалена лунка {well_to_remove}")
                            st.rerun()

                        rows_to_remove = df[df[well_col].astype(str) == well_to_remove]
                        st.dataframe(rows_to_remove.style.applymap(
                            lambda val: "background-color: #f0ff47" if isinstance(val, (int, float)) and val > 1000 else "",
                            subset=[dist_col]
                        ), use_container_width=True)

                        multi_remove = st.multiselect("Множественный выбор:", sorted(wells_high, key=well_sort_key), key=f"multi_remove_{selected_key}")
                        if multi_remove and st.button("✅ Удалить выбранные", key=f"multi_remove_btn_{selected_key}"):
                            df_new = df[~df[well_col].astype(str).isin(multi_remove)].reset_index(drop=True)
                            st.session_state[f"{selected_key}_data"] = df_new
                            st.session_state[f"ver_{selected_key}"] += 1
                            st.success(f"Удалено {len(multi_remove)} лунок")
                            st.rerun()

                        if st.button("❌ Удалить все жёлтые лунки", key=f"remove_all_high_{selected_key}"):
                            df_new = df[~df[well_col].astype(str).isin(wells_high)].reset_index(drop=True)
                            st.session_state[f"{selected_key}_data"] = df_new
                            st.session_state[f"ver_{selected_key}"] += 1
                            st.success(f"Удалены все жёлтые лунки: {len(wells_high)}")
                            st.rerun()
                    else:
                        st.info("Нет жёлтых лунок")

                # === ОРАНЖЕВЫЕ ===
                with col_orange:
                    st.markdown("#### 🟠 Обработка лунок с '-'")
                    if wells_with_dash:
                        well_dash = st.selectbox("Одиночный выбор:", sorted(wells_with_dash, key=well_sort_key), key=f"dash_select_{selected_key}")
                        action = st.radio("Действие:", ["Удалить лунку", "Заменить '-' на NaN"], key=f"dash_action_{selected_key}")
                        if st.button("Применить", key=f"dash_apply_{selected_key}"):
                            if action == "Удалить лунку":
                                df = df[df[well_col].astype(str) != well_dash].reset_index(drop=True)
                            else:
                                mask = df[well_col].astype(str) == well_dash
                                df.loc[mask] = df.loc[mask].replace("-", np.nan)
                            st.session_state[f"{selected_key}_data"] = df
                            st.session_state[f"ver_{selected_key}"] += 1
                            st.rerun()

                        st.dataframe(df[df[well_col].astype(str) == well_dash].style.applymap(
                            lambda val: "background-color: #ffd699" if isinstance(val, str) and val.strip() == "-" else ""
                        ), use_container_width=True)

                        multi_dash = st.multiselect("Множественный выбор:", sorted(wells_with_dash, key=well_sort_key), key=f"dash_multi_select_{selected_key}")
                        if multi_dash:
                            multi_action = st.radio("Действие для выбранных:", ["Удалить", "Заменить '-' на NaN"], key=f"dash_multi_action_{selected_key}")
                            if st.button("⚙️ Применить к выбранным", key=f"dash_multi_apply_{selected_key}"):
                                if multi_action == "Удалить":
                                    df = df[~df[well_col].astype(str).isin(multi_dash)].reset_index(drop=True)
                                else:
                                    mask = df[well_col].astype(str).isin(multi_dash)
                                    df.loc[mask] = df.loc[mask].replace("-", np.nan)
                                st.session_state[f"{selected_key}_data"] = df
                                st.session_state[f"ver_{selected_key}"] += 1
                                st.rerun()

                        action_all = st.radio("Действие для всех:", ["Удалить все", "Заменить все '-' на NaN"], key=f"dash_all_action_{selected_key}")
                        if st.button("⚙️ Применить ко всем", key=f"dash_all_apply_{selected_key}"):
                            if action_all == "Удалить все":
                                df = df[~df[well_col].astype(str).isin(wells_with_dash)].reset_index(drop=True)
                            else:
                                mask = df[well_col].astype(str).isin(wells_with_dash)
                                df.loc[mask] = df.loc[mask].replace("-", np.nan)
                            st.session_state[f"{selected_key}_data"] = df
                            st.session_state[f"ver_{selected_key}"] += 1
                            st.success("Обработка завершена")
                            st.rerun()
                    else:
                        st.info("Нет лунок с '-'")

                # === КРАСНЫЕ ===
                with col_red:
                    st.markdown("#### 🔴 Восстановление лунок")
                    if removed_wells:
                        # Одиночный выбор
                        well_to_restore = st.selectbox("Одиночный выбор:", sorted(removed_wells, key=well_sort_key), key=f"restore_select_{selected_key}")
                        if st.button("Восстановить", key=f"restore_btn_{selected_key}"):
                            df_clean = df[df[well_col].astype(str) != well_to_restore].copy()
                            df_restore = raw_df[raw_df[well_col].astype(str) == well_to_restore].copy()
                            df_new = pd.concat([df_clean, df_restore], ignore_index=True)
                            st.session_state[f"{selected_key}_data"] = df_new.reset_index(drop=True)
                            st.session_state[f"ver_{selected_key}"] += 1
                            st.success(f"Лунка {well_to_restore} восстановлена")
                            st.rerun()

                        # Множественный выбор
                        multi_restore = st.multiselect("Множественный выбор:", sorted(removed_wells, key=well_sort_key), key=f"multi_restore_select_{selected_key}")
                        if multi_restore and st.button("♻️ Восстановить выбранные", key=f"multi_restore_btn_{selected_key}"):
                            df_clean = df[~df[well_col].astype(str).isin(multi_restore)].copy()
                            df_restore = raw_df[raw_df[well_col].astype(str).isin(multi_restore)].copy()
                            df_new = pd.concat([df_clean, df_restore], ignore_index=True)
                            st.session_state[f"{selected_key}_data"] = df_new.reset_index(drop=True)
                            st.session_state[f"ver_{selected_key}"] += 1
                            st.success(f"Восстановлены: {', '.join(multi_restore)}")
                            st.rerun()

                        # Восстановить все
                        if st.button("♻️ Восстановить все лунки", key=f"restore_all_{selected_key}"):
                            df_clean = df.copy()
                            df_restore = raw_df[raw_df[well_col].astype(str).isin(removed_wells)].copy()
                            df_new = pd.concat([df_clean, df_restore], ignore_index=True)
                            st.session_state[f"{selected_key}_data"] = df_new.reset_index(drop=True)
                            st.session_state[f"ver_{selected_key}"] += 1
                            st.success("Восстановлены все удалённые лунки")
                            st.rerun()
                    else:
                        st.info("Нет удалённых лунок")

            with st.expander("📋 Показать таблицу данных по всем лункам"):
                st.dataframe(df, use_container_width=True)

        # === САЙДБАР: динамический отчёт ===
        report = compute_deletion_report(file_keys)
        render_sidebar_report(report)

    # ================
    # 📊 Вкладка агрегации
    # ================
    with tab_agg:
        st.markdown("### 📊 Агрегация Velocity по временным отрезкам и группам")

        key_selected = st.selectbox("Выберите файл:", file_keys, format_func=lambda x: f"{x}")

        df = st.session_state.get(f"{key_selected}_data")
        if df is None:
            st.warning("Сначала загрузите таблицы на вкладке «Загрузка и просмотр».")
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
            try:
               df[col_velocity] = pd.to_numeric(df[col_velocity])#, errors="coerce")
 
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
               
               agg_df["Time"] = pd.to_numeric(agg_df["TimeNumeric"])#, errors="coerce")
               agg_df["Concentration"] = pd.to_numeric(agg_df["Concentration"])#, errors="coerce")

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

            except ValueError as e:
               st.error(f"❌ В столбце `{col_velocity}` есть необработанные значения ('-'). Перейдите во вкладку «⚙️ Анализ данных» и обработайте их.")

    # ================
    # 📐 Вкладка расчётов
    # ================
    with tab_calc:
        st.markdown("### 📐 Расчёт Distance и Velocity по каждой лунке и сегментам (10-мин)")

        key_selected_calc = st.selectbox("Выберите файл для расчётов:", file_keys, key="select_calc")

        df = st.session_state.get(f"{key_selected_calc}_data")
        if df is None:
            st.warning("Сначала загрузите таблицы на вкладке «Загрузка и просмотр».")
        else:
            col_well = df.columns[2]
            col_time = "Time"

            dist_cols = [col for col in df.columns if isinstance(col, str) and "Distance moved" in col]
            velocity_cols = [col for col in df.columns if isinstance(col, str) and "Velocity" in col]
            moving_cols = [col for col in df.columns if isinstance(col, str) and "Moving" in col]
            not_moving_cols = [col for col in df.columns if isinstance(col, str) and "Not Moving" in col]
            turn_angle_cols = [col for col in df.columns if isinstance(col, str) and "Turn angle" in col]

            if not dist_cols or not velocity_cols or not turn_angle_cols or not moving_cols or not not_moving_cols:
                st.error("Отсутствуют необходимые колонки.")
                st.stop()

            col_dist = dist_cols[0]
            col_velocity = velocity_cols[0]
            col_moving = moving_cols[0]
            col_not_moving = not_moving_cols[0]
            col_turn_angle = turn_angle_cols[0]
            
            try:
              df[col_dist] = pd.to_numeric(df[col_dist])
              df[col_velocity] = pd.to_numeric(df[col_velocity])
              df[col_moving] = pd.to_numeric(df[col_moving])
              df[col_not_moving] = pd.to_numeric(df[col_not_moving])
              df[col_turn_angle] = pd.to_numeric(df[col_turn_angle])

              time_order = df[col_time].dropna().unique().tolist()
              time_mapping = {t: str((i + 1) * 2) for i, t in enumerate(time_order)}
              df["TimeNumeric"] = df[col_time].map(time_mapping).astype(float)

              original_well_order = df[col_well].dropna().astype(str).drop_duplicates().tolist()

              def get_segment(tn: float) -> str:
                  try:
                      tn = int(tn)
                      if 2 <= tn <= 10:
                          return "1"
                      elif 12 <= tn <= 20:
                          return "2"
                      elif 22 <= tn <= 30:
                          return "3"
                      elif 32 <= tn <= 40:
                          return "4"
                      elif 42 <= tn <= 50:
                          return "5"
                      else:
                          return "Other"
                  except:
                      return "Invalid"

              df["Segment"] = df["TimeNumeric"].apply(get_segment)

              # === MOVING / NOT MOVING RATIO ===
              moving_col = "Movement (Moving / Center-point, Cumulative Duration, s)"
              notmoving_col = "Movement (Not Moving / Center-point, Cumulative Duration, s)"

              if moving_col not in df.columns or notmoving_col not in df.columns:
                  st.error("❌ В таблице не найдены нужные колонки для расчёта отношения Moving / Not Moving.")
                  st.stop()

              mov_avg = (
                  df.groupby([col_well, "Segment"], sort=False)[moving_col]
                  .mean()
                  .unstack()
                  .add_prefix("mov_avg_")
              )

              notmov_avg = (
                  df.groupby([col_well, "Segment"], sort=False)[notmoving_col]
                  .mean()
                  .unstack()
                  .add_prefix("notmov_avg_")
              )

              moving_ratio_df = pd.DataFrame(index=mov_avg.index)

              for seg in ["1", "2", "3", "4", "5"]:
                  mov = mov_avg.get(f"mov_avg_{seg}", np.nan)
                  notmov = notmov_avg.get(f"notmov_avg_{seg}", np.nan).replace(0, np.nan)
                  moving_ratio_df[f"M_moving_not_moving_ratio_{seg}"] = mov / notmov

              moving_ratio_df = moving_ratio_df.reset_index()

              # === TURN ANGLE MEAN ±abs по сегментам ===
              angle_col = "Turn angle (Center-point / relative, Mean, deg)"
              if angle_col not in df.columns:
                  st.error("❌ Колонка с Turn angle не найдена.")
                  st.stop()

              angle_mean_df = (
                  df.groupby([col_well, "Segment"], sort=False)[angle_col]
                  .agg(['mean', lambda x: x.abs().mean()])
                  .rename(columns={"mean": "signed", "<lambda_0>": "abs_mean"})
                  .reset_index()
              )

              angle_pivot_signed = angle_mean_df.pivot(index=col_well, columns="Segment", values="signed").add_prefix("T_angle_mean_")
              angle_pivot_abs = angle_mean_df.pivot(index=col_well, columns="Segment", values="abs_mean").add_prefix("T_angle_mean_abs_")

              # Собрать в нужном порядке: abs_1, mean_1, abs_2, mean_2, ...
              angle_combined = angle_pivot_abs.copy()
              for s in ["1", "2", "3", "4", "5"]:
                  angle_combined[f"T_angle_mean_{s}"] = angle_pivot_signed[f"T_angle_mean_{s}"]

              combined_cols = []
              for s in ["1", "2", "3", "4", "5"]:
                  combined_cols += [f"T_angle_mean_abs_{s}", f"T_angle_mean_{s}"]

              angle_combined = angle_combined[combined_cols].reset_index()

              # === DISTANCE ===
              total_df = df.groupby(col_well, sort=False)[col_dist].sum().reset_index()
              total_df = total_df.rename(columns={col_dist: "Distance Total"})

              seg_df = df.groupby([col_well, "Segment"], sort=False)[col_dist].sum().reset_index()
              seg_pivot = seg_df.pivot(index=col_well, columns="Segment", values=col_dist).reset_index()
              seg_pivot.columns.name = None
              seg_pivot = seg_pivot.rename(columns={seg: f"D_{seg}" for seg in seg_pivot.columns if seg not in [col_well]})
              result_dist = pd.merge(total_df, seg_pivot, on=col_well, how="left").fillna(0)

              # === VELOCITY MEAN ===
              vel_mean_df = (
                  df.groupby([col_well, "Segment"], sort=False)[col_velocity]
                  .mean().reset_index()
                  .pivot(index=col_well, columns="Segment", values=col_velocity)
                  .reset_index()
                  .fillna(0)
              )
              vel_mean_df.columns.name = None
              vel_mean_df = vel_mean_df.rename(columns=lambda x: f"V_mean_{x}" if x != col_well else x)

              # === VELOCITY VARIANCE ===
              vel_var_df = (
                  df.groupby([col_well, "Segment"], sort=False)[col_velocity]
                  .var(ddof=1).reset_index()
                  .pivot(index=col_well, columns="Segment", values=col_velocity)
                  .reset_index()
                  .fillna(0)
              )
              vel_var_df.columns.name = None
              vel_var_df = vel_var_df.rename(columns=lambda x: f"V_var_{x}" if x != col_well else x)

              # === ОБЪЕДИНЕНИЕ ВСЕХ ===
              merged_df = result_dist.merge(vel_mean_df, on=col_well, how="left")
              merged_df = merged_df.merge(vel_var_df, on=col_well, how="left")
              merged_df = merged_df.merge(moving_ratio_df, on=col_well, how="left")
              merged_df = merged_df.merge(angle_combined, on=col_well, how="left")

              # === ДОБАВИМ ОТНОШЕНИЯ V_mean ===
              def compute_ratio(df, num, denom):
                  col_n = f"V_mean_{num}"
                  col_d = f"V_mean_{denom}"
                  new_col = f"V_ratio_{num}_{denom}"
                  if col_n in df.columns and col_d in df.columns:
                      df[new_col] = df[col_n] / df[col_d].replace(0, np.nan)
                  return df

              for a, b in [("1", "2"), ("3", "2"), ("3", "4"), ("5", "4")]:
                  merged_df = compute_ratio(merged_df, a, b)

              # === ПОРЯДОК КОЛОНОК ===
              base_cols = [col_well, "Distance Total"]
              dist_seg_cols = [f"D_{i}" for i in ["1", "2", "3", "4", "5"]]
              vel_cols = []
              for seg in ["1", "2", "3", "4", "5"]:
                  m = f"V_mean_{seg}"
                  v = f"V_var_{seg}"
                  if m in merged_df.columns and v in merged_df.columns:
                      vel_cols += [m, v]
              ratio_cols = [col for col in merged_df.columns if col.startswith("V_ratio_")]
              moving_ratio_cols = [col for col in merged_df.columns if col.startswith("M_moving_not_moving_ratio_")]
              angle_cols = []
              for s in ["1", "2", "3", "4", "5"]:
                  angle_cols += [f"T_angle_mean_abs_{s}", f"T_angle_mean_{s}"]

              final_cols = base_cols + dist_seg_cols + vel_cols + ratio_cols + moving_ratio_cols + angle_cols
              merged_df = merged_df[[c for c in final_cols if c in merged_df.columns]]

              # Порядок лунок
              merged_df[col_well] = pd.Categorical(merged_df[col_well], categories=original_well_order, ordered=True)
              merged_df = merged_df.sort_values(col_well).reset_index(drop=True)

              st.dataframe(merged_df, use_container_width=True)

              # Скачивание
              towrite = io.BytesIO()
              with pd.ExcelWriter(towrite, engine="openpyxl") as writer:
                  merged_df.to_excel(writer, index=False, sheet_name="Distance_Velocity_Summary")
              towrite.seek(0)
              st.download_button(
                  label="⬇️ Скачать таблицу расчётов",
                  data=towrite,
                  file_name=f"{key_selected_calc}_summary.xlsx",
                  mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                  key="dl_summary_calc"
              )
            except ValueError as e:
               st.error(f"❌ В данных есть необработанные значения ('-'). Перейдите во вкладку «⚙️ Анализ данных» и обработайте их.")
