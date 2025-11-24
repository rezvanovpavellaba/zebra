import io
import hashlib
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import numpy as np
import plotly.graph_objects as go
import pathlib
from utils.fc_analysis_neuro import *


# === УТИЛИТЫ ===

# Если у Control пустой Compound -> взять значение Compound из Test для того же эксперимента
def fill_control_compound(df):
    # Словарь: концентрация -> compound (по Test)
    mapping = (
        df.loc[df["Test/control"] == "Test", ["Concentration", "Compound"]]
        .drop_duplicates(subset=["Concentration"])
        .set_index("Concentration")["Compound"]
        .to_dict()
    )
    # Заполняем Control тем же Compound, что и у Test (берём первую попавшуюся запись)
    mask = (df["Test/control"] == "Control") & (df["Compound"].isna() | (df["Compound"] == ""))
    if mask.any() and mapping:
        common_compound = list(mapping.values())[0]  # берём первое (Compound один для всех Test)
        df.loc[mask, "Compound"] = common_compound
    return df

def sort_wells(wells):
    def key_fn(w):
        letter = w[0]                       # буква (A–H)
        number = int(w[1:]) if w[1:].isdigit() else 0  # число после буквы
        return (letter, number)
    return ", ".join(sorted(wells, key=key_fn))

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

    selected_key = st.sidebar.selectbox("Выберите временной интервал:", file_keys)

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

            dist_numeric = pd.to_numeric(df[dist_col])#, errors="coerce")

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
                    st.markdown("#### 🟡 Обработка лунок с Distance > 1000")
                    if wells_high:

                        data_cols_dist = df.columns[distance_start_idx:]

                        well_high = st.selectbox(
                            "Одиночный выбор:",
                            sorted(wells_high, key=well_sort_key),
                            key=f"high_select_{selected_key}"
                        )

                        action = st.radio(
                            "Действие:",
                            ["Удалить лунку", "Заменить >1000 на NaN", "Заменить >1000 на 0"],
                            key=f"high_action_{selected_key}"
                        )

                        if st.button("Применить", key=f"high_apply_{selected_key}"):
                            mask = df[well_col].astype(str) == well_high
                            

                            if action == "Удалить лунку":
                                df = df[~mask].reset_index(drop=True)
                            else:
                                # Определяем строки, где Distance > 1000
                                distance_numeric = pd.to_numeric(df.loc[mask, dist_col])#, errors="coerce")
                                target_rows = mask & (distance_numeric > 1000)

                                # Применяем замену ко всем колонкам после Distance
                                if action == "Заменить >1000 на NaN":
                                    df.loc[target_rows, data_cols_dist] = np.nan
                                else:
                                    # Заполняем значения с логикой "Not Moving → 120, остальные → 0"
                                    for col in data_cols_dist:
                                        if "Not Moving" in col:
                                            df.loc[target_rows, col] = 120
                                        else:
                                            df.loc[target_rows, col] = 0

                            st.session_state[f"{selected_key}_data"] = df
                            st.session_state[f"ver_{selected_key}"] += 1
                            st.rerun()

                        # === Редактирование вручную
                        st.markdown("##### 📝 Редактируемая таблица")
                        highlight_df = df[df[well_col].astype(str) == well_high].copy()
                        editable_df = highlight_df.copy()
                        editable_df[data_cols] = editable_df[data_cols].astype(str)
                        editable_df.reset_index(drop=True, inplace=True)

                        editable_df["⚠️"] = editable_df[dist_col].apply(
                            lambda x: "⚠️" if pd.to_numeric(x) > 1000 else "" #, errors="coerce")
                        )

                        edited_df = st.data_editor(
                            editable_df,
                            use_container_width=True,
                            key=f"editor_high_{selected_key}_{well_high}"
                        )

                        if st.button("💾 Сохранить изменения", key=f"save_edits_high_{selected_key}_{well_high}"):
                            mask = df[well_col].astype(str) == well_high
                            well_indices = df[mask].index
                            edited_clean = edited_df.drop(columns=["⚠️"])
                            if len(well_indices) == len(edited_clean):
                                df_updated = df.copy()
                                for orig_idx, (_, row) in zip(well_indices, edited_clean.iterrows()):
                                    df_updated.loc[orig_idx] = row
                                st.session_state[f"{selected_key}_data"] = df_updated
                                st.session_state[f"ver_{selected_key}"] += 1
                                st.success(f"Изменения для лунки {well_high} сохранены.")
                                st.rerun()
                            else:
                                st.error("❌ Количество строк изменилось. Возможно, вы удалили строку — сохранение невозможно.")

                        # === Множественный выбор
                        multi_high = st.multiselect(
                            "Множественный выбор:",
                            sorted(wells_high, key=well_sort_key),
                            key=f"high_multi_select_{selected_key}"
                        )
                        if multi_high:
                            multi_action = st.radio(
                                "Действие для выбранных:",
                                ["Удалить", "Заменить >1000 на NaN", "Заменить >1000 на 0"],
                                key=f"high_multi_action_{selected_key}"
                            )
                            if st.button("⚙️ Применить к выбранным", key=f"high_multi_apply_{selected_key}"):
                                mask = df[well_col].astype(str).isin(multi_high)
                                if multi_action == "Удалить":
                                    df = df[~mask].reset_index(drop=True)
                                else:
                                    # Определяем строки, где Distance > 1000
                                    distance_numeric = pd.to_numeric(df.loc[mask, dist_col])#, errors="coerce")
                                    target_rows = mask & (distance_numeric > 1000)

                                    # Применяем замену ко всем колонкам после Distance
                                    if action == "Заменить >1000 на NaN":
                                        df.loc[target_rows, data_cols_dist] = np.nan
                                    else:
                                        # Заполняем значения с логикой "Not Moving → 120, остальные → 0"
                                        for col in data_cols_dist:
                                            if "Not Moving" in col:
                                                df.loc[target_rows, col] = 120
                                            else:
                                                df.loc[target_rows, col] = 0

                                st.session_state[f"{selected_key}_data"] = df
                                st.session_state[f"ver_{selected_key}"] += 1
                                st.rerun()

                        # === Массовая обработка
                        action_all = st.radio(
                            "Действие для всех:",
                            ["Удалить все", "Заменить все >1000 на NaN", "Заменить все >1000 на 0"],
                            key=f"high_all_action_{selected_key}"
                        )
                        if st.button("⚙️ Применить ко всем", key=f"high_all_apply_{selected_key}"):
                            mask = df[well_col].astype(str).isin(wells_high)
                            if action_all == "Удалить все":
                                df = df[~mask].reset_index(drop=True)
                            else:
                                # Определяем строки, где Distance > 1000
                                distance_numeric = pd.to_numeric(df.loc[mask, dist_col])#, errors="coerce")
                                target_rows = mask & (distance_numeric > 1000)

                                # Применяем замену ко всем колонкам после Distance
                                if action == "Заменить >1000 на NaN":
                                    df.loc[target_rows, data_cols_dist] = np.nan
                                else:
                                    # Заполняем значения с логикой "Not Moving → 120, остальные → 0"
                                    for col in data_cols_dist:
                                        if "Not Moving" in col:
                                            df.loc[target_rows, col] = 120
                                        else:
                                            df.loc[target_rows, col] = 0

                            st.session_state[f"{selected_key}_data"] = df
                            st.session_state[f"ver_{selected_key}"] += 1
                            st.success("Обработка завершена")
                            st.rerun()
                    else:
                        st.info("Нет лунок с Distance > 1000")


                # === ОРАНЖЕВЫЕ ===
                with col_orange:
                    st.markdown("#### 🟠 Обработка лунок с '-'")
                    if wells_with_dash:
                        well_dash = st.selectbox(
                            "Одиночный выбор:",
                            sorted(wells_with_dash, key=well_sort_key),
                            key=f"dash_select_{selected_key}"
                        )

                        # ➕ ДОБАВИЛИ третий режим: "Заменить '-' на 0"
                        action = st.radio(
                            "Действие:",
                            ["Удалить лунку", "Заменить '-' на NaN", "Заменить '-' на 0"],
                            key=f"dash_action_{selected_key}"
                        )

                        if st.button("Применить", key=f"dash_apply_{selected_key}"):
                            if action == "Удалить лунку":
                                df = df[df[well_col].astype(str) != well_dash].reset_index(drop=True)
                            else:
                                # заменяем только в метриках, а не в мета-колонках
                                mask = df[well_col].astype(str) == well_dash
                                if action == "Заменить '-' на NaN":
                                    df.loc[mask, data_cols] = df.loc[mask, data_cols].replace("-", np.nan)
                                else:  # "Заменить '-' на 0"
                                    # сразу ставим числовой 0, чтобы потом .to_numeric проходил без ошибок
                                    df.loc[mask, data_cols] = df.loc[mask, data_cols].replace("-", 0)

                            st.session_state[f"{selected_key}_data"] = df
                            st.session_state[f"ver_{selected_key}"] += 1
                            st.rerun()

                        # # Подсветка ячеек с '-'
                        # # 📋 Просмотр с подсветкой
                        # st.markdown("##### 📋 Таблица с подсветкой '-' (только просмотр)")
                        highlight_df = df[df[well_col].astype(str) == well_dash].copy()
                        # st.dataframe(
                        #     highlight_df.style.applymap(
                        #         lambda val: "background-color: #ffd699" if isinstance(val, str) and val.strip() == "-" else ""
                        #     ),
                        #     use_container_width=True
                        # )

                        # 🧠 Создаём копию для редактирования и добавим колонку ⚠️
                        st.markdown("##### 📝 Редактируемая таблица")
                        editable_df = highlight_df.copy()

                        # Преобразуем все метрики в строки, чтобы точно редактировалось
                        editable_df[data_cols] = editable_df[data_cols].astype(str)
                        editable_df.reset_index(drop=True, inplace=True)

                        # Добавим колонку "⚠️", если есть хотя бы одно "-"
                        editable_df["⚠️"] = editable_df[data_cols].astype(str).apply(
                            lambda row: "⚠️" if any(cell.strip() == "-" for cell in row) else "",
                            axis=1
                        )

                        # Показываем редактор
                        edited_df = st.data_editor(
                            editable_df,
                            use_container_width=True,
                            key=f"editor_{selected_key}_{well_dash}"
                        )

                        # Кнопка сохранения изменений
                        if st.button("💾 Сохранить изменения", key=f"save_edits_{selected_key}_{well_dash}"):
                           # Определим строки, соответствующие этой лунке
                           mask = df[well_col].astype(str) == well_dash
                           well_indices = df[mask].index

                           # Удалим колонку "⚠️" перед вставкой
                           edited_clean = edited_df.drop(columns=["⚠️"])

                           if len(well_indices) == len(edited_clean):
                               df_updated = df.copy()
                               for orig_idx, (_, row) in zip(well_indices, edited_clean.iterrows()):
                                   df_updated.loc[orig_idx] = row

                               st.session_state[f"{selected_key}_data"] = df_updated
                               st.session_state[f"ver_{selected_key}"] += 1
                               st.success(f"Изменения для лунки {well_dash} сохранены.")
                               st.rerun()
                           else:
                               st.error("❌ Количество строк изменилось. Возможно, вы удалили строку — сохранение невозможно.")

                        # === Множественный выбор ===
                        multi_dash = st.multiselect(
                            "Множественный выбор:",
                            sorted(wells_with_dash, key=well_sort_key),
                            key=f"dash_multi_select_{selected_key}"
                        )
                        if multi_dash:
                            multi_action = st.radio(
                                "Действие для выбранных:",
                                ["Удалить", "Заменить '-' на NaN", "Заменить '-' на 0"],
                                key=f"dash_multi_action_{selected_key}"
                            )
                            if st.button("⚙️ Применить к выбранным", key=f"dash_multi_apply_{selected_key}"):
                                if multi_action == "Удалить":
                                    df = df[~df[well_col].astype(str).isin(multi_dash)].reset_index(drop=True)
                                elif multi_action == "Заменить '-' на NaN":
                                    mask = df[well_col].astype(str).isin(multi_dash)
                                    df.loc[mask, data_cols] = df.loc[mask, data_cols].replace("-", np.nan)
                                else:  # "Заменить '-' на 0"
                                    mask = df[well_col].astype(str).isin(multi_dash)
                                    df.loc[mask, data_cols] = df.loc[mask, data_cols].replace("-", 0)

                                st.session_state[f"{selected_key}_data"] = df
                                st.session_state[f"ver_{selected_key}"] += 1
                                st.rerun()

                        # === Массовая обработка всех оранжевых ===
                        action_all = st.radio(
                            "Действие для всех:",
                            ["Удалить все", "Заменить все '-' на NaN", "Заменить все '-' на 0"],
                            key=f"dash_all_action_{selected_key}"
                        )
                        if st.button("⚙️ Применить ко всем", key=f"dash_all_apply_{selected_key}"):
                            if action_all == "Удалить все":
                                df = df[~df[well_col].astype(str).isin(wells_with_dash)].reset_index(drop=True)
                            elif action_all == "Заменить все '-' на NaN":
                                mask = df[well_col].astype(str).isin(wells_with_dash)
                                df.loc[mask, data_cols] = df.loc[mask, data_cols].replace("-", np.nan)
                            else:  # "Заменить все '-' на 0"
                                mask = df[well_col].astype(str).isin(wells_with_dash)
                                df.loc[mask, data_cols] = df.loc[mask, data_cols].replace("-", 0)

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
                            
                            #при восстановлении сохраняем прежний порядок строк
                            # Соединяем с сохранением всех строк
                            df_new = pd.concat([df_clean, df_restore], ignore_index=False)

                            # Получаем порядок индексов исходного файла
                            raw_index_order = list(raw_df.index)

                            # Убираем дубликаты индексов в порядке, где они впервые встретились
                            raw_index_order = pd.Index(raw_index_order).drop_duplicates().tolist()

                            # Чтобы восстановленные строки попали на свои места:
                            # 1. Берём все строки из raw_df в исходном порядке
                            # 2. Отфильтровываем только те, которые есть в df_new по well_id
                            well_col = get_well_col(raw_df)
                            valid_wells = set(df_new[well_col].astype(str))

                            # Формируем новый DataFrame, следуя точному порядку строк из raw_df
                            df_new_ordered = raw_df[raw_df[well_col].astype(str).isin(valid_wells)].copy().reset_index(drop=True)
                            #######################################################################################################

                            st.session_state[f"{selected_key}_data"] = df_new_ordered
                            st.session_state[f"ver_{selected_key}"] += 1
                            st.success(f"Лунка {well_to_restore} восстановлена")
                            st.rerun()

                        # Множественный выбор
                        multi_restore = st.multiselect("Множественный выбор:", sorted(removed_wells, key=well_sort_key), key=f"multi_restore_select_{selected_key}")
                        if multi_restore and st.button("♻️ Восстановить выбранные", key=f"multi_restore_btn_{selected_key}"):
                            df_clean = df[~df[well_col].astype(str).isin(multi_restore)].copy()
                            df_restore = raw_df[raw_df[well_col].astype(str).isin(multi_restore)].copy()
                            
                            #при восстановлении сохраняем прежний порядок строк
                            # Соединяем с сохранением всех строк
                            df_new = pd.concat([df_clean, df_restore], ignore_index=False)

                            # Получаем порядок индексов исходного файла
                            raw_index_order = list(raw_df.index)

                            # Убираем дубликаты индексов в порядке, где они впервые встретились
                            raw_index_order = pd.Index(raw_index_order).drop_duplicates().tolist()

                            # Чтобы восстановленные строки попали на свои места:
                            # 1. Берём все строки из raw_df в исходном порядке
                            # 2. Отфильтровываем только те, которые есть в df_new по well_id
                            well_col = get_well_col(raw_df)
                            valid_wells = set(df_new[well_col].astype(str))

                            # Формируем новый DataFrame, следуя точному порядку строк из raw_df
                            df_new_ordered = raw_df[raw_df[well_col].astype(str).isin(valid_wells)].copy().reset_index(drop=True)
                            #######################################################################################################

                            st.session_state[f"{selected_key}_data"] = df_new_ordered
                            st.session_state[f"ver_{selected_key}"] += 1
                            st.success(f"Восстановлены: {', '.join(multi_restore)}")
                            st.rerun()

                        # Восстановить все
                        if st.button("♻️ Восстановить все лунки", key=f"restore_all_{selected_key}"):
                            df_clean = df.copy()
                            df_restore = raw_df[raw_df[well_col].astype(str).isin(removed_wells)].copy()
                            
                            #при восстановлении сохраняем прежний порядок строк
                            # Соединяем с сохранением всех строк
                            df_new = pd.concat([df_clean, df_restore], ignore_index=False)

                            # Получаем порядок индексов исходного файла
                            raw_index_order = list(raw_df.index)

                            # Убираем дубликаты индексов в порядке, где они впервые встретились
                            raw_index_order = pd.Index(raw_index_order).drop_duplicates().tolist()

                            # Чтобы восстановленные строки попали на свои места:
                            # 1. Берём все строки из raw_df в исходном порядке
                            # 2. Отфильтровываем только те, которые есть в df_new по well_id
                            well_col = get_well_col(raw_df)
                            valid_wells = set(df_new[well_col].astype(str))

                            # Формируем новый DataFrame, следуя точному порядку строк из raw_df
                            df_new_ordered = raw_df[raw_df[well_col].astype(str).isin(valid_wells)].copy().reset_index(drop=True)
                            #######################################################################################################

                            st.session_state[f"{selected_key}_data"] = df_new_ordered
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

        key_selected = selected_key

        df = st.session_state.get(f"{key_selected}_data")
        if df is None:
            st.warning("Сначала загрузите таблицы на вкладке «Загрузка и просмотр».")
        else:
            df = df.copy()  # 🔧 вот эта строка — ключ к изоляции изменений
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
               st.write(agg_df["Concentration"])
               st.error(f"❌ В столбце `{col_velocity}` есть необработанные значения ('-'). Перейдите во вкладку «⚙️ Анализ данных» и обработайте их.")

    # ================
    # 📐 Вкладка расчётов
    # ================
    with tab_calc:
        st.markdown("### 📐 Расчёт поведенческих метрик по каждой лунке и сегментам (10-мин)")

        key_selected_calc = selected_key

        df = st.session_state.get(f"{key_selected_calc}_data")
        if df is None:
            st.warning("Сначала загрузите таблицы на вкладке «Загрузка и просмотр».")
        else:
            df = df.copy()  # 🔧 вот эта строка — ключ к изоляции изменений
            col_well = df.columns[2]
            col_time = "Time"

            dist_cols = [col for col in df.columns if isinstance(col, str) and "Distance moved" in col]
            velocity_cols = [col for col in df.columns if isinstance(col, str) and "Velocity" in col]
            moving_cols = [col for col in df.columns if isinstance(col, str) and "Moving" in col]
            not_moving_cols = [col for col in df.columns if isinstance(col, str) and "Not Moving" in col]
            turn_angle_cols = [col for col in df.columns if isinstance(col, str) and "Turn angle" in col]
            angular_velocity_cols = [col for col in df.columns if isinstance(col, str) and "Angular velocity" in col]
            meander_cols = [col for col in df.columns if isinstance(col, str) and "Meander" in col and "Mean" in col and "Total" not in col]
            meander_total_cols = [col for col in df.columns if isinstance(col, str) and "Meander" in col and "Total" in col]
            cw_cols = [col for col in df.columns if isinstance(col, str) and col.startswith("CW Rotation")]
            ccw_cols = [col for col in df.columns if isinstance(col, str) and col.startswith("CCW Rotation")]
       

            if not dist_cols or not velocity_cols or not turn_angle_cols or not moving_cols or not not_moving_cols or not angular_velocity_cols or not meander_cols or not meander_total_cols or not cw_cols or not ccw_cols:
                st.error("Отсутствуют необходимые колонки.")
                st.stop()

            col_dist = dist_cols[0]
            col_velocity = velocity_cols[0]
            col_moving = moving_cols[0]
            col_not_moving = not_moving_cols[0]
            col_turn_angle = turn_angle_cols[0]
            col_ang_velocity = angular_velocity_cols[0]
            col_meander = meander_cols[0]
            col_meander_total = meander_total_cols[0]
            col_cw = cw_cols[0]
            col_ccw = ccw_cols[0]
            
            try:
              df[col_dist] = pd.to_numeric(df[col_dist])
              df[col_velocity] = pd.to_numeric(df[col_velocity])
              df[col_moving] = pd.to_numeric(df[col_moving])
              df[col_not_moving] = pd.to_numeric(df[col_not_moving])
              df[col_turn_angle] = pd.to_numeric(df[col_turn_angle])
              df[col_ang_velocity] = pd.to_numeric(df[col_ang_velocity])
              df[col_meander] = pd.to_numeric(df[col_meander])
              df[col_meander_total] = pd.to_numeric(df[col_meander_total])
              df[col_cw] = pd.to_numeric(df[col_cw])
              df[col_ccw] = pd.to_numeric(df[col_ccw])

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

              # === ROTATION CW / CCW: среднее и отношение по сегментам ===
              cw_mean = df.groupby([col_well, "Segment"], sort=False)[col_cw].mean().unstack().add_prefix("CW_rotation_mean_")
              ccw_mean = df.groupby([col_well, "Segment"], sort=False)[col_ccw].mean().unstack().add_prefix("CCW_rotation_mean_")


              ratio_df = pd.DataFrame(index=cw_mean.index)
              for s in ["1", "2", "3", "4", "5"]:
                  cw_col = f"CW_rotation_mean_{s}"
                  ccw_col = f"CCW_rotation_mean_{s}"
                  ratio_col = f"Rotation_ratio_CW_CCW_{s}"
                  ratio_df[ratio_col] = cw_mean.get(cw_col, np.nan) / ccw_mean.get(ccw_col, np.nan).replace(0, np.nan)

              rotation_combined = pd.concat([cw_mean, ccw_mean, ratio_df], axis=1).reset_index()

              # === MEANDER TOTAL: сумма по модулю по всем сегментам ===
              meander_total_df = (
                  df.groupby(col_well, sort=False)[col_meander_total]
                  .apply(lambda x: x.abs().sum())
                  .reset_index()
                  .rename(columns={col_meander_total: "M_t_total_abs"})
              )

              # === MEANDER MEAN: сумма и среднее по модулю ===
              meander_df = (
                  df.groupby([col_well, "Segment"], sort=False)[col_meander]
                  .agg(['sum', lambda x: x.abs().mean()])
                  .rename(columns={"sum": "signed_sum", "<lambda_0>": "abs_mean"})
                  .reset_index()
              )

              meander_pivot_signed = meander_df.pivot(index=col_well, columns="Segment", values="signed_sum").add_prefix("M_m_sum_")
              meander_pivot_abs = meander_df.pivot(index=col_well, columns="Segment", values="abs_mean").add_prefix("M_m_mean_abs_")

              meander_combined = meander_pivot_abs.copy()
              for s in ["1", "2", "3", "4", "5"]:
                  meander_combined[f"M_m_sum_{s}"] = meander_pivot_signed[f"M_m_sum_{s}"]

              meander_combined = meander_combined.reset_index()

              # === ANGULAR VELOCITY: среднее по модулю и сумма ===
              angvel_df = (
                  df.groupby([col_well, "Segment"], sort=False)[col_ang_velocity]
                  .agg(['sum', lambda x: x.abs().mean()])
                  .rename(columns={"sum": "signed_sum", "<lambda_0>": "abs_mean"})
                  .reset_index()
              )

              angvel_pivot_signed = angvel_df.pivot(index=col_well, columns="Segment", values="signed_sum").add_prefix("A_v_sum_")
              angvel_pivot_abs = angvel_df.pivot(index=col_well, columns="Segment", values="abs_mean").add_prefix("A_v_mean_abs_")

              # Собрать в нужном порядке
              angvel_combined = angvel_pivot_abs.copy()
              for s in ["1", "2", "3", "4", "5"]:
                  angvel_combined[f"A_v_sum_{s}"] = angvel_pivot_signed[f"A_v_sum_{s}"]

              # Финализируем
              angvel_combined = angvel_combined.reset_index()

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
              # === Извлечение мета-информации по каждой лунке (первая строка на каждую well_id)
              meta_df = df[[col_well, "Test/control", "Compound", "Concentration"]].drop_duplicates(subset=col_well)
              
              merged_df = result_dist.merge(vel_mean_df, on=col_well, how="left")
              merged_df = merged_df.merge(vel_var_df, on=col_well, how="left")
              merged_df = merged_df.merge(moving_ratio_df, on=col_well, how="left")
              merged_df = merged_df.merge(angle_combined, on=col_well, how="left")
              merged_df = merged_df.merge(angvel_combined, on=col_well, how="left")
              merged_df = merged_df.merge(meander_combined, on=col_well, how="left")
              merged_df = merged_df.merge(meander_total_df, on=col_well, how="left")
              merged_df = merged_df.merge(rotation_combined, on=col_well, how="left")

              merged_df = merged_df.merge(meta_df, on=col_well, how="left")

              # Приводим Concentration к числовому дальше нужно будет
              if "Concentration" in merged_df.columns:
                  merged_df["Concentration"] = pd.to_numeric(merged_df["Concentration"], errors="coerce") #здесь тихая замена не страшна, т.к. в данном случае не должно быть None или чего-то еще

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
              meta_cols = ["Test/control", "Compound", "Concentration"]
              base_cols = [col_well] + meta_cols + ["Distance Total"]

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
              
              angvel_cols = []
              for s in ["1", "2", "3", "4", "5"]:
                  angvel_cols += [f"A_v_mean_abs_{s}", f"A_v_sum_{s}"]

              meander_cols = []
              for s in ["1", "2", "3", "4", "5"]:
                  meander_cols += [f"M_m_mean_abs_{s}", f"M_m_sum_{s}"]

              rotation_cols = []
              for s in ["1", "2", "3", "4", "5"]:
                  rotation_cols += [
                      f"CW_rotation_mean_{s}",
                      f"CCW_rotation_mean_{s}",
                      f"Rotation_ratio_CW_CCW_{s}",
                  ]

              # Объединение всех колонок
              final_cols = base_cols + dist_seg_cols + vel_cols + ratio_cols + moving_ratio_cols + angle_cols + angvel_cols + meander_cols
              final_cols += ["M_t_total_abs"]
              final_cols += rotation_cols
              merged_df = merged_df[[c for c in final_cols if c in merged_df.columns]]

              # Порядок лунок
              merged_df[col_well] = pd.Categorical(merged_df[col_well], categories=original_well_order, ordered=True)
              merged_df = merged_df.sort_values(col_well).reset_index(drop=True)

              # === ВСТАВИТЬ ЗДЕСЬ === Дополняем название препарата у контроля
              merged_df = fill_control_compound(merged_df)

              st.dataframe(merged_df, use_container_width=True)

              # Скачивание
              towrite = io.BytesIO()
              with pd.ExcelWriter(towrite, engine="openpyxl") as writer:
                  merged_df.to_excel(writer, index=False, sheet_name="Calculations_inside_well")
              towrite.seek(0)
              st.download_button(
                  label="⬇️ Скачать таблицу расчётов",
                  data=towrite,
                  file_name=f"{key_selected_calc}_summary.xlsx",
                  mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                  key="dl_summary_calc"
              )

              # Путь к справке
              help_path = pathlib.Path("docs/metrics_explained.txt")

              if help_path.exists():
                  with st.expander("ℹ️ Справка по метрикам", expanded=False):
                      with open(help_path, "r", encoding="utf-8") as f:
                          lines = f.readlines()

                      # Парсинг в словарь: ключ — метрика, значение — описание
                      help_dict = {}
                      for line in lines:
                          if "—" in line:
                              key, desc = line.split("—", 1)
                              help_dict[key.strip()] = desc.strip()

                      # Поле поиска
                      query = st.text_input("🔍 Введите название или часть описания метрики:")

                      if query:
                          # Поиск по ключу или описанию
                          results = {
                              k: v for k, v in help_dict.items()
                              if query.lower() in k.lower() or query.lower() in v.lower()
                          }
                          if results:
                              st.markdown("### 🔎 Результаты поиска:")
                              for k, v in results.items():
                                  st.markdown(f"**`{k}`** — {v}")
                          else:
                              st.info("❗ Ничего не найдено.")
                      else:
                          st.markdown("### 📋 Все метрики:")
                          for k, v in help_dict.items():
                              st.markdown(f"**`{k}`** — {v}")
              else:
                  st.warning("Файл справки по метрикам не найден: `docs/metrics_explained.txt`")

              st.markdown("### 📊 Агрегация по группам лунок")

              # Определяем ключевые колонки
              meta_cols = ["Test/control", "Compound", "Concentration"]
              well_col = merged_df.columns[0]  # обычно это well_id
              data_cols = [col for col in merged_df.columns if col not in [well_col] + meta_cols]

              # Добавим колонку со списком лунок
              well_df = (
                  merged_df
                  .groupby(meta_cols, sort=False, dropna=False)[well_col]
                  .apply(lambda x: sort_wells(x.dropna().astype(str)))
                  .reset_index()
                  .rename(columns={well_col: "Wells"})
              )

              # Считаем среднее и SD
              agg_df = (
                  merged_df
                  .groupby(meta_cols, sort=False, dropna=False)[data_cols]
                  .agg(['mean', 'std'])
                  .reset_index()
              )

              # Уплощаем MultiIndex колонок
              agg_df.columns = [
                  f"{col[0]} ({col[1]})" if col[1] else col[0]
                  for col in agg_df.columns.values
              ]

              # Объединяем с таблицей лунок
              agg_df = pd.merge(agg_df, well_df, on=meta_cols, how="left")

              # Переставим "Wells" в начало
              cols = ["Wells"] + [c for c in agg_df.columns if c != "Wells"]
              agg_df = agg_df[cols]

              # Показываем
              st.dataframe(agg_df, use_container_width=True)

              # Скачивание
              towrite = io.BytesIO()
              with pd.ExcelWriter(towrite, engine="openpyxl") as writer:
                  agg_df.to_excel(writer, index=False, sheet_name="Aggregated by group")
              towrite.seek(0)
              st.download_button(
                  label="⬇️ Скачать агрегированные данные по группам",
                  data=towrite,
                  file_name=f"{key_selected_calc}_group_summary.xlsx",
                  mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                  key="dl_group_summary"
              )
            
              ##################FC

              # Инициализация session state для fc_mode, если его еще нет
              if "fc_mode_neuro" not in st.session_state:
                  st.session_state["fc_mode_neuro"] = 'log₂(B/A)'  # Значение по умолчанию

              fc_mode = st.radio(
                    "Режим расчета Fold Change",
                    ('ratio (B/A)', 'difference ((B-A)/A)', 'log₂(B/A)'),
                    index=['ratio (B/A)', 'difference ((B-A)/A)', 'log₂(B/A)'].index(st.session_state["fc_mode_neuro"]),
                    key="key_radio_fc_mode_neuro"
                )

              if fc_mode != st.session_state["fc_mode_neuro"]:
                 st.session_state["fc_mode_neuro"] = fc_mode


              mode = 'ratio' if fc_mode == 'ratio (B/A)' else 'difference' if fc_mode == 'difference ((B-A)/A)' else 'log2_ratio'

              fc_df,warnings_fc = calculate_fold_change_with_pvalues(merged_df, mode=mode)

              st.session_state['fc_df_neuro'] = fc_df
              st.session_state['warnings_fc_neuro'] = warnings_fc

              if 'fc_df_neuro' in st.session_state:
                  fc_df = st.session_state['fc_df_neuro']
                  warnings_fc = st.session_state['warnings_fc_neuro']
                  
                  st.subheader(f"Результаты расчета Fold Change ({st.session_state['fc_mode_neuro']})")
                  st.dataframe(fc_df)

                  with st.sidebar:

                      if len(warnings_fc) !=0:

                          selected_warnings_fc = st.selectbox(
                              f"Отчет об ошибках в расчётах ({len(warnings_fc)} всего)",
                              warnings_fc,
                              index=0,
                              key="warnings_fc_select_neuro"
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
                      file_name=f"fold_change_results_{st.session_state['fc_mode_neuro']}.xlsx",
                      mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                      key = "Скачать полные результаты Fold Change_neuro"

                  )

                  # Добавляем Volcano Plot только для режима log₂(B/A)
                  if st.session_state["fc_mode_neuro"] == 'log₂(B/A)':
                      st.subheader("Volcano Plot")
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
                      available_drugs = fc_df['Compound'].unique()
                      available_concentrations= fc_df['Concentration'].unique()


                      if 'show_vulcano_plot_neuro' not in st.session_state:
                          st.session_state['show_vulcano_plot_neuro'] = False
                      
                      # Создаем отдельную форму для Volcano Plot
                      with st.form("volcano_form_neuro"):
                          # Независимый выбор препаратов
                          volcano_drugs = st.multiselect(
                              "Выберите препараты для Volcano Plot",
                              available_drugs,
                              default=available_drugs[:min(5, len(available_drugs))],  # Первые 5 по умолчанию
                              key="volcano_drugs_neuro"
                          )

                          volcano_conc = st.selectbox(
                              "Выберите концентрацию для Volcano Plot",
                              available_concentrations,
                              index=1,  # Первая тестовая
                              key="volcano_conc_neuro"
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
                                      help="Уровни значимости: *** - 0.1%, ** - 1%, * - 5%, ns - не значимо",
                                      key= "Порог p-value (статистическая значимость)_neuro"
                                  )
                                  
                                  log2fc_threshold = st.radio(
                                      "Порог log2FC (кратность изменения)",
                                      options=[1.0, 0.58],
                                      index=0,  # 1.0 по умолчанию
                                      format_func=lambda x: f"{x} ({'2-кратное' if x == 1.0 else '1.5-кратное'} изменение)",
                                      horizontal=True,
                                      key= "Порог log2FC (кратность изменения)_neuro"
                                  )
                              
                              with col2:
                                  # Настройки отображения
                                  show_legend = st.checkbox(
                                      "Показывать легенду", 
                                      value=True,
                                      key='volcano_show_legend_neuro'
                                  )
                                  
                                  # Раздельные настройки цветов линий
                                  hline_color = st.color_picker(
                                      "Цвет горизонтальной линии (p-value)",
                                      value='#FF0000',
                                      key='volcano_hline_color_neuro'
                                  )
                                  
                                  vline_color = st.color_picker(
                                      "Цвет вертикальных линий (log2FC)",
                                      value='#808080',
                                      key='volcano_vline_color_neuro'
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
                                                  key=f"volcano_color_{drug}_neuro"
                                              )
                          
                          submitted_volcano = st.form_submit_button("Построить Volcano Plot") #здесь не нужен ключ

                          if submitted_volcano:
                             st.session_state['show_vulcano_plot_neuro'] = True

                      if st.session_state['show_vulcano_plot_neuro'] and volcano_drugs and volcano_conc:
                          significant_df = plot_volcano(
                              fc_df, 
                              volcano_drugs,
                              volcano_conc, 
                              p_value_threshold=p_value_threshold,
                              log2fc_threshold=log2fc_threshold,
                              custom_colors=volcano_colors,
                              show_legend=show_legend,
                              hline_color=hline_color,
                              vline_color=vline_color
                          )
                          
                          if significant_df is not None and not significant_df.empty:
                              st.subheader("Значимые метрики")
                              st.write(f"Найдено {len(significant_df)} значимых метрик (p < {p_value_threshold}, |log2FC| > {log2fc_threshold})")
                              st.dataframe(significant_df)
                              
                              # Кнопка скачивания значимых метаболитов
                              output_significant = io.BytesIO()
                              with pd.ExcelWriter(output_significant, engine='openpyxl') as writer:
                                  significant_df.to_excel(writer, index=False, sheet_name='Significant_Metrics')
                              output_significant.seek(0)
                              
                              st.download_button(
                                  label="Скачать таблицу значимых метрик",
                                  data=output_significant,
                                  file_name=f"significant_metrics_p{p_value_threshold}_fc{log2fc_threshold}.xlsx",
                                  mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                                  key = "Скачать таблицу значимых метрик_neuro"
                              )
                          elif significant_df is not None and significant_df.empty:
                              st.warning("Нет значимых метрик при заданных параметрах.")


            except ValueError as e:
               st.error(f"❌ В данных есть необработанные значения ('-'). Перейдите во вкладку «⚙️ Анализ данных» и обработайте их.")
