import streamlit as st
import pandas as pd
from io import BytesIO
import matplotlib.pyplot as plt
import numpy as np
from utils.activity_cardio import *
from utils.hr_cardio import *


def cardiotoxicity_app():

    tab_activity, tab_hr = st.tabs([
        "📊 Активность (%)",
        "❤️ ЧСС"
    ])

    with tab_activity:

        st.title("📊 Активность (%)")

        uploaded_file = st.file_uploader(
            "Загрузите Excel файл (Raw Data)",
            type=["xlsx", "xls"]
        )

        if uploaded_file:
            xls = pd.ExcelFile(uploaded_file)

            combined_df = pd.DataFrame()
            name_counters = {}

            for sheet_name in xls.sheet_names:
                df = xls.parse(sheet_name, header=0)

                parsed = extract_time_percent(df)
                if parsed is None:
                    st.warning(f"⚠️ Лист '{sheet_name}' пропущен (нет Time (s) / %)")
                    continue

                base_name = extract_base_name(sheet_name)

                name_counters.setdefault(base_name, 0)
                name_counters[base_name] += 1

                col_suffix = f"{base_name}_{name_counters[base_name]}"

                parsed = parsed.rename(
                    columns={
                        "Time (s)": f"Time_{col_suffix}",
                        "%": f"Percent_{col_suffix}"
                    }
                )

                combined_df = pd.concat([combined_df, parsed], axis=1)

            st.subheader("📋 Исходная таблица")
            st.dataframe(combined_df, use_container_width=True)

            # ---- выгрузка ----
            buffer = BytesIO()
            with pd.ExcelWriter(buffer) as writer:
                combined_df.to_excel(writer, index=False, sheet_name="Result")

            st.download_button(
                label="⬇️ Скачать результат (Excel)",
                data=buffer.getvalue(),
                file_name="parsed_result.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

            st.subheader("🎛️ Выбор субъектов")

            available_subjects = get_available_subjects(combined_df)

            selected_subjects = st.multiselect(
                label="Показать на графике:",
                options=available_subjects,
                default=[],
                help="Выберите субъекты для отображения графиков"
            )

            st.subheader("📈 Графики активности")

            chart_type = st.radio(
                "Тип графика",
                ["Matplotlib (статичный)", "Plotly (интерактивный)"],
                index=0,
                horizontal=True,
                key="activity_chart_type",
            )

            title_activity = st.text_input(
                "Заголовок графика",
                value="Compound",
                key="activity_title",
            )

            colA, colB = st.columns(2)
            with colA:
                control_name = st.text_input(
                    "Название Control на графике",
                    value="Control",
                    key="activity_control_name",
                    help="Например: Control / Контроль"
                )
            with colB:
                subject_word = st.text_input(
                    "Название Subject на графике",
                    value="Subject",
                    key="activity_subject_word",
                    help="Например: Subject / Fish / Larva"
                )

            colx, coly = st.columns(2)
            with colx:
                x_label = st.text_input("X-label", value="Time (s)", key="activity_xlabel")
            with coly:
                y_label = st.text_input("Y-label", value="Activity (%)", key="activity_ylabel")

            if selected_subjects:
                if chart_type.startswith("Matplotlib"):
                    fig = plot_activity_curves(
                        combined_df,
                        selected_subjects,
                        title=title_activity,
                        xlabel=x_label,
                        ylabel=y_label,
                        control_name=control_name,
                        subject_word=subject_word,
                    )
                    st.pyplot(fig, use_container_width=True)
                else:
                    fig_p = plot_activity_curves_plotly(
                        combined_df,
                        selected_subjects,
                        title=title_activity,
                        xlabel=x_label,
                        ylabel=y_label,
                        control_name=control_name,
                        subject_word=subject_word,
                    )
                    st.plotly_chart(fig_p, use_container_width=True)
            else:
                st.info("Выбери субъектов в селекторе выше для построения графика.")
            
    with tab_hr:
        st.title("❤️ ЧСС")

        # ---------- Загрузка файла ----------
        hr_file = st.file_uploader(
            label="Загрузите Excel файл (Trial Statistics)",
            type=["xlsx", "xls"],
            key="hr_uploader"
        )

        if hr_file is not None:
            # ---------- Сбор базовой таблицы ----------
            try:
                hr_df = build_hr_table(hr_file)
            except Exception as e:
                st.error(f"Ошибка при обработке файла ЧСС: {e}")
                st.stop()

            # ---------- Таблица 1: исходные данные ----------
            st.subheader("📋 Исходная таблица ЧСС (по субъектам)")
            st.dataframe(hr_df, use_container_width=True)

            raw_excel = dataframe_to_excel_bytes(
                hr_df,
                sheet_name="HR_raw"
            )

            st.download_button(
                label="⬇️ Скачать исходную таблицу ЧСС (Excel)",
                data=raw_excel,
                file_name="heart_rate_raw.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

            # ---------- Агрегация ----------
            agg_df = aggregate_hr_by_compound(hr_df)
            summary_df = format_hr_summary_table(agg_df)

            # ---------- Таблица 2: агрегированная ----------
            st.subheader("📊 Агрегация ЧСС по Compound (N / Mean / SE)")
            st.dataframe(summary_df, use_container_width=True)

            summary_excel = dataframe_to_excel_bytes(
                summary_df,
                sheet_name="HR_summary"
            )

            st.download_button(
                label="⬇️ Скачать агрегированную таблицу ЧСС (Excel)",
                data=summary_excel,
                file_name="heart_rate_summary.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
            
            ####Статистика
            
            st.subheader("Статистический анализ")

            # ===== alpha selectors =====
            ALPHA_OPTIONS = [0.05, 0.01, 0.001]

            # ---------- Нормальность распределения (Shapiro–Wilk) ----------
            st.subheader("🧪 Проверка нормальности по группам (Shapiro–Wilk)")
            
            
            alpha_norm = st.selectbox(
                "α (нормальность, Shapiro)",
                ALPHA_OPTIONS,
                index=0,
                key="alpha_norm",
            )

            normality_df = build_normality_report_table(
                hr_df,
                metrics=["BPS Mean", "BPM Mean"],
                alpha=alpha_norm,
            )

            st.dataframe(normality_df, use_container_width=True)

            with st.expander("ⓘ Проверка нормальности (Shapiro–Wilk) — как интерпретировать"):
                st.markdown("""
                **Что проверяется**  
                Для каждой метрики (*BPS Mean*, *BPM Mean*) тест Shapiro–Wilk выполняется **отдельно в каждой группе Compound**.

                **Как читать таблицу**
                - `p_value > α` → нет оснований отвергнуть нормальность распределения  
                - `p_value ≤ α` → распределение статистически отличается от нормального  
                - `Normal (p>α)` — логический итог для конкретной группы  
                - `Note` показывает, был ли тест применён (например, при `n < 3` тест не проводится)

                **Важно**
                - При малых выборках отсутствие значимости **не доказывает** нормальность.
                - Для выбора ANOVA **все группы** должны быть нормальными.
                """)

            normality_excel = dataframe_to_excel_bytes(
                normality_df,
                sheet_name="Normality"
            )

            st.download_button(
                label="⬇️ Скачать таблицу нормальности (Excel)",
                data=normality_excel,
                file_name="hr_normality_shapiro.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
            
            # ---------- Гомогенность дисперсий (Levene) ----------
            st.subheader("🧪 Проверка гомогенности дисперсий (Levene)")

            alpha_levene = st.selectbox(
                "α (дисперсии, Levene)",
                ALPHA_OPTIONS,
                index=0,
                key="alpha_levene",
            )

            levene_df = build_levene_report_table(
                hr_df,
                metrics=["BPS Mean", "BPM Mean"],
                alpha=alpha_levene,
            )

            st.dataframe(levene_df, use_container_width=True)

            with st.expander("ⓘ Проверка равенства дисперсий (Levene) — как интерпретировать"):
                st.markdown("""
            **Что проверяется**  
            Тест Levene оценивает, равны ли дисперсии между всеми группами *Compound* для данной метрики.

            **Как читать таблицу**
            - `p_value > α` → дисперсии статистически не различаются  
            - `Equal_variance = True` → условие гомогенности дисперсий выполнено  
            - Используется вариант с `center = median`, устойчивый к выбросам

            **Важно**
            - Для применения **ANOVA** требуется равенство дисперсий  
            - При нарушении этого условия выбирается **Kruskal–Wallis**
            """)

            levene_excel = dataframe_to_excel_bytes(levene_df, sheet_name="Levene")
            st.download_button(
                label="⬇️ Скачать таблицу Levene (Excel)",
                data=levene_excel,
                file_name="hr_levene.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

            # ---------- Выбор теста + глобальный тест (единая таблица) ----------
            st.subheader("📌 Выбор теста + глобальный тест (ANOVA / Kruskal–Wallis)")


            alpha_global = st.selectbox(
                "α (глобальный тест)",
                ALPHA_OPTIONS,
                index=0,
                key="alpha_global",
            )
        
            combined_global_df = build_decision_and_global_table(
                hr_df=hr_df,
                normality_df=normality_df,
                levene_df=levene_df,
                metrics=["BPS Mean", "BPM Mean"],
                alpha=alpha_global,
            )

            st.dataframe(combined_global_df, use_container_width=True)

            with st.expander("ⓘ Глобальный тест (ANOVA / Kruskal–Wallis) — что это значит"):
                st.markdown("""
            **Как выбирается тест**
            - Если **все группы нормальны** и **дисперсии равны** → **ANOVA**
            - В остальных случаях → **Kruskal–Wallis**

            **Как читать таблицу**
            - `Statistic`:
            - ANOVA → F-статистика  
            - Kruskal–Wallis → H-статистика
            - `p_value < α` → есть статистически значимые различия **хотя бы между двумя группами**
            - `p_value ≥ α` → различий не обнаружено

            **Важно**
            - Глобальный тест **не показывает**, какие именно группы различаются  
            - Он только отвечает на вопрос: *есть ли различия вообще*
            """)

            combined_excel = dataframe_to_excel_bytes(combined_global_df, sheet_name="Decision+Global")
            st.download_button(
                label="⬇️ Скачать (выбор теста + глобальный тест) (Excel)",
                data=combined_excel,
                file_name="hr_decision_global.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            )

            # ---------- Post-hoc (Tukey / Dunn-Holm) ----------
            st.subheader("🔍 Post-hoc анализ (vs Control)")

            alpha_posthoc = st.selectbox(
                "α (post-hoc)",
                ALPHA_OPTIONS,
                index=0,
                key="alpha_posthoc",
            )

            posthoc_summary_df, posthoc_details = build_posthoc_reports(
                hr_df=hr_df,
                decision_global_df=combined_global_df,
                metrics=["BPS Mean", "BPM Mean"],
                alpha_posthoc=alpha_posthoc,
                alpha_global=alpha_global,
                control_label="Control",
                only_vs_control=True,
                dunn_p_adjust="holm",
            )

            st.dataframe(posthoc_summary_df, use_container_width=True)

            with st.expander("ⓘ Post-hoc анализ (Tukey / Dunn) — интерпретация"):
                st.markdown("""
            **Когда выполняется**
            - Только если глобальный тест значим (`p_value < α`)
            - Если глобальный тест незначим — post-hoc **не проводится** (это корректно)

            **Какие тесты используются**
            - После **ANOVA** → *Tukey HSD*
            - После **Kruskal–Wallis** → *Dunn с коррекцией Holm*

            ---

            ## Если выбран **Tukey HSD** (после ANOVA)

            **Что сравнивается**
            - Каждая пара групп (*Group 1* vs *Group 2*).  
            - Только сравнения конкретной дозы **vs Control**.

            **Колонки таблицы при Tukey**

            - `meandiff`  
            Разность средних: **Mean(Group 1) − Mean(Group 2)**.  
            Пример: `meandiff = -0.83` значит, что среднее в Control **меньше**, чем в Group 2, на 0.83 (в единицах метрики).  
            > Важно: знак показывает направление изменения.

            - `p_value_adj`  
            p-value **с учётом множественных сравнений**.  
            Это именно то значение, по которому судим о значимости.

            - `lower`, `upper`  
            Доверительный интервал (95%) для `meandiff`.  
            Если интервал **не включает 0**, различие статистически значимо.

            - `Reject_H0`  
            Итоговый флаг значимости (обычно эквивалентно `p_value_adj < α`).  
            `True` → нулевая гипотеза отвергается, изменение значимо.

            ---

            ## Если выбран **Dunn (Holm)** (после Kruskal–Wallis)

            **Что сравнивается**
            - Попарные сравнения между группами на основе рангов (не средних).

            **Колонки Dunn**
            - `Group 1`, `Group 2` — сравниваемые группы (`Group 1 = Control`)
            - `p_value_adj` — p-value после коррекции Holm (только для сравнений vs Control)
            - `Reject_H0` — `p_value_adj < α` (изменение значимо)

            **Важно**
            - В Dunn-тесте нет `meandiff` и доверительного интервала: он даёт только статистическую значимость (через p-values)
            - Коррекция Holm контролирует семейную ошибку первого рода при множественных сравнениях
            """)


            posthoc_summary_excel = dataframe_to_excel_bytes(posthoc_summary_df, sheet_name="Posthoc_summary")
            st.download_button(
                label="⬇️ Скачать сводку post-hoc (Excel)",
                data=posthoc_summary_excel,
                file_name="hr_posthoc_summary.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            )

            # Детальные таблицы
            if posthoc_details:
                for metric, ph_df in posthoc_details.items():
                    st.markdown(f"**{metric}**")
                    st.dataframe(ph_df, use_container_width=True)

                    ph_excel = dataframe_to_excel_bytes(ph_df, sheet_name="Posthoc")
                    st.download_button(
                        label=f"⬇️ Скачать post-hoc для {metric} (Excel)",
                        data=ph_excel,
                        file_name=f"hr_posthoc_{metric.replace(' ', '_')}.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    )

            ##############Графики
            def fig_to_png_bytes(fig) -> bytes:
                buf = BytesIO()
                fig.savefig(buf, format="png", dpi=600, bbox_inches="tight")
                buf.seek(0)
                return buf.getvalue()

            st.subheader("📊 Boxplots (BPS / BPM)")

            title = st.text_input("Заголовок", value="24 hpf", key="title")


            col1, col2, col3 = st.columns(3)
            with col1:
                ylabel_bps = st.text_input("Y-label BPS", value="Heart rate (bps)", key="bps_ylabel")
            with col2:
                ylabel_bpm = st.text_input("Y-label BPM", value="Heart rate (bpm)", key="bpm_ylabel")
            with col3:
                xlabel_common = st.text_input("X-label", value="", key="hr_xlabel", help = "Оставь пустым для автоподстановки названия соединения из таблицы или впиши свой вариант")
            
            control_tick = st.text_input("Подпись тика Control", value="Control", key="hr_control_tick")


            show_points = st.checkbox("Показывать точки (individual values)", value=True, key="hr_show_points")


            def p_to_sig(p: float | None) -> str:
                if p is None or (isinstance(p, float) and not np.isfinite(p)):
                    return "ns"
                if p < 0.001:
                    return "***"
                if p < 0.01:
                    return "**"
                if p < 0.05:
                    return "*"
                return "ns"
            
            def build_signif_map_for_metric(
                hr_df: pd.DataFrame,
                metric: str,
                posthoc_details: dict[str, pd.DataFrame],
                control_label: str = "Control",
            ) -> dict[str, str]:
                """
                Возвращает dict[full_group_name -> '***'/'**'/'*'/'ns'] для всех НЕ-control групп.
                """
                # default: ns для всех доз
                groups = sorted(hr_df["Compound"].astype(str).unique().tolist())
                out = {g: "ns" for g in groups if g.strip().lower() != control_label.strip().lower()}

                ph = posthoc_details.get(metric)
                if ph is None or ph.empty:
                    return out  # post-hoc не считался -> всё ns

                # ожидаем таблицу вида: Group 1 (Control), Group 2 (dose), p_value_adj
                if "Group 2" not in ph.columns or "p_value_adj" not in ph.columns:
                    return out

                for _, row in ph.iterrows():
                    g2 = str(row["Group 2"])
                    p = row["p_value_adj"]
                    p = float(p) if pd.notna(p) else np.nan
                    if g2.strip().lower() == control_label.strip().lower():
                        continue
                    out[g2] = p_to_sig(p)

                return out

            # --- BPS ---
            signif_bps = build_signif_map_for_metric(
                hr_df=hr_df,
                metric="BPS Mean",
                posthoc_details=posthoc_details,
                control_label="Control",
            )

            fig_bps = build_hr_boxplot(
                hr_df=hr_df,
                metric="BPS Mean",
                title=title,
                ylabel=ylabel_bps,
                xlabel=xlabel_common,
                show_points=show_points,
                control_tick=control_tick,
                signif_by_group=signif_bps,
            )
            st.pyplot(fig_bps, use_container_width=True)

            png_bps = fig_to_png_bytes(fig_bps)
            st.download_button(
                "⬇️ Скачать BPS boxplot (PNG)",
                data=png_bps,
                file_name="boxplot_bps.png",
                mime="image/png",
            )

            # --- BPM ---
            signif_bpm = build_signif_map_for_metric(
                hr_df=hr_df,
                metric="BPM Mean",
                posthoc_details=posthoc_details,
                control_label="Control",
            )

            fig_bpm = build_hr_boxplot(
                hr_df=hr_df,
                metric="BPM Mean",
                title=title,
                ylabel=ylabel_bpm,
                xlabel=xlabel_common,
                show_points=show_points,
                control_tick=control_tick,
                signif_by_group=signif_bpm,
            )
            st.pyplot(fig_bpm, use_container_width=True)

            png_bpm = fig_to_png_bytes(fig_bpm)
            st.download_button(
                "⬇️ Скачать BPM boxplot (PNG)",
                data=png_bpm,
                file_name="boxplot_bpm.png",
                mime="image/png",
            )