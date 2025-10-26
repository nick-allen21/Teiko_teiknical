"""
Drug analysis toolkit and dashboard.

Implements the plan described in README:
  - Part 1: Create SQLite schema and load `data/cell-count.csv`
  - Part 2: Build summary table of cell population relative frequencies per sample
  - Part 3: Compare responders vs non-responders (melanoma, PBMC, miraclib) with stats and boxplots
  - Part 4: Subset analysis for baseline melanoma PBMC samples on miraclib

Notes:
  - The dataset includes melanoma and carcinoma; analyses in Part 3 focus on melanoma only.
  - Run as a Streamlit app: `streamlit run drug_analysis.py`
  - Or use CLI: `python drug_analysis.py --help`
"""

from __future__ import annotations

import argparse
import os
import sqlite3
from dataclasses import dataclass
from typing import Iterable, List, Tuple
import plotly.express as px 
import numpy as np
import pandas as pd
import streamlit as st 



DATABASE_DEFAULT_PATH = "teiko.sqlite"
CSV_DEFAULT_PATH = os.path.join("data", "cell-count.csv")


POPULATION_COLUMNS: List[str] = [
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]


def get_db_connection(db_path: str = DATABASE_DEFAULT_PATH) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON;")
    return conn


def initialize_database(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS samples (
            sample TEXT PRIMARY KEY,
            project TEXT NOT NULL,
            subject TEXT NOT NULL,
            condition TEXT NOT NULL,
            age INTEGER,
            sex TEXT,
            treatment TEXT,
            response TEXT,
            sample_type TEXT NOT NULL,
            time_from_treatment_start INTEGER NOT NULL
        );

        CREATE TABLE IF NOT EXISTS cell_counts (
            sample TEXT NOT NULL,
            population TEXT NOT NULL,
            count INTEGER NOT NULL,
            PRIMARY KEY (sample, population),
            FOREIGN KEY (sample) REFERENCES samples(sample) ON DELETE CASCADE
        );

        CREATE INDEX IF NOT EXISTS idx_samples_subject ON samples(subject);
        CREATE INDEX IF NOT EXISTS idx_samples_filters ON samples(condition, treatment, sample_type, time_from_treatment_start);
        CREATE INDEX IF NOT EXISTS idx_cell_counts_population ON cell_counts(population);
        """
    )
    conn.commit()


def _standardize_response(value: str | float | None) -> str | None:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    text = str(value).strip().lower()
    if text in {"yes", "y", "true", "respond", "responder"}:
        return "yes"
    if text in {"no", "n", "false", "non-respond", "nonresponder", "non-responder"}:
        return "no"
    return text or None


def load_csv_to_db(csv_path: str, conn: sqlite3.Connection) -> None:
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV not found at {csv_path}")

    df = pd.read_csv(csv_path)

    required_columns = {
        "project",
        "subject",
        "condition",
        "age",
        "sex",
        "treatment",
        "response",
        "sample",
        "sample_type",
        "time_from_treatment_start",
    } | set(POPULATION_COLUMNS)

    missing = required_columns - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing columns: {sorted(missing)}")

    # Standardize response to yes/no/None
    df["response"] = df["response"].apply(_standardize_response)

    # Upsert samples
    samples_df = df[
        [
            "sample",
            "project",
            "subject",
            "condition",
            "age",
            "sex",
            "treatment",
            "response",
            "sample_type",
            "time_from_treatment_start",
        ]
    ].drop_duplicates(subset=["sample"]).copy()

    samples_records: List[Tuple] = list(
        samples_df.itertuples(index=False, name=None)
    )

    conn.executemany(
        """
        INSERT OR REPLACE INTO samples (
            sample, project, subject, condition, age, sex, treatment, response, sample_type, time_from_treatment_start
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        samples_records,
    )

    # Prepare long-form cell_counts
    long_rows: List[Tuple[str, str, int]] = []
    for _, row in df.iterrows():
        sample_id = row["sample"]
        for population in POPULATION_COLUMNS:
            count_value = int(row[population])
            long_rows.append((sample_id, population, count_value))

    conn.executemany(
        """
        INSERT OR REPLACE INTO cell_counts (sample, population, count)
        VALUES (?, ?, ?)
        """,
        long_rows,
    )
    conn.commit()



def fetch_counts_with_metadata(conn: sqlite3.Connection) -> pd.DataFrame:
    query = """
        SELECT s.sample,
               s.project,
               s.subject,
               s.condition,
               s.age,
               s.sex,
               s.treatment,
               s.response,
               s.sample_type,
               s.time_from_treatment_start,
               c.population,
               c.count
        FROM cell_counts c
        JOIN samples s ON s.sample = c.sample
    """
    return pd.read_sql_query(query, conn)


def compute_summary_table(conn: sqlite3.Connection) -> pd.DataFrame:
    df = fetch_counts_with_metadata(conn)
    if df.empty:
        return df

    totals = df.groupby("sample")["count"].sum().rename("total_count").reset_index()
    merged = df.merge(totals, on="sample", how="left")
    merged["percentage"] = merged["count"] / merged["total_count"].replace(0, np.nan) * 100.0

    # Ensure ordering of populations is stable
    merged["population"] = pd.Categorical(merged["population"], categories=POPULATION_COLUMNS, ordered=True)

    summary_cols = [
        "sample",
        "total_count",
        "population",
        "count",
        "percentage",
        # include metadata useful for downstream filtering/visuals
        "project",
        "subject",
        "condition",
        "age",
        "sex",
        "treatment",
        "response",
        "sample_type",
        "time_from_treatment_start",
    ]
    return merged[summary_cols].sort_values(["sample", "population"]).reset_index(drop=True)


@dataclass
class TestResult:
    population: str
    n_yes: int
    n_no: int
    median_yes: float
    median_no: float
    effect_size: float
    p_value: float
    p_adj: float
    significant: bool


def _cliffs_delta(x: Iterable[float], y: Iterable[float]) -> float:
    # Cliff's delta: probability that a value in x is greater than one in y minus the reverse
    x_arr = np.asarray(list(x), dtype=float)
    y_arr = np.asarray(list(y), dtype=float)
    n_x = len(x_arr)
    n_y = len(y_arr)
    if n_x == 0 or n_y == 0:
        return np.nan
    gt = 0
    lt = 0
    for xv in x_arr:
        gt += np.sum(xv > y_arr)
        lt += np.sum(xv < y_arr)
    return float(gt - lt) / float(n_x * n_y)


def _bh_adjust(p_values: List[float]) -> List[float]:
    m = len(p_values)
    if m == 0:
        return []
    order = np.argsort(p_values)
    ordered_p = np.array(p_values)[order]
    adjusted = np.empty(m, dtype=float)
    # Compute raw adjusted values
    raw = ordered_p * m / (np.arange(m) + 1)
    # Enforce monotonicity from the end
    adjusted[-1] = min(raw[-1], 1.0)
    for i in range(m - 2, -1, -1):
        adjusted[i] = min(raw[i], adjusted[i + 1])
    # Revert to original order
    result = np.empty(m, dtype=float)
    result[order] = adjusted
    return result.tolist()


def analyze_responders_vs_nonresponders(summary_df: pd.DataFrame, alpha: float = 0.05) -> Tuple[pd.DataFrame, List["object"]]:
    # Filter per spec: melanoma, PBMC, miraclib; compare yes vs no responders
    df = summary_df.copy()
    mask = (
        (df["condition"].str.lower() == "melanoma")
        & (df["sample_type"].str.upper() == "PBMC")
        & (df["treatment"].str.lower() == "miraclib")
        & (df["response"].isin(["yes", "no"]))
    )
    df = df.loc[mask]

    # Prepare figures 
    figures = []
   
    fig = px.box(
        df,
        x="population",
        y="percentage",
        color="response",
        points="all",
        hover_data=["sample", "subject", "time_from_treatment_start"],
        category_orders={"population": POPULATION_COLUMNS, "response": ["no", "yes"]},
        title="PBMC melanoma miraclib: Relative frequency by response",
    )
    fig.update_layout(legend_title_text="Responder")
    figures.append(fig)
   

    # Statistical tests per population: Mann-Whitney U
    from scipy.stats import mannwhitneyu

    results: List[TestResult] = []
    grouped = df.groupby("population")
    pvals: List[float] = []
    interim: List[Tuple[str, float, int, int, float, float]] = []  # pop, p, n_yes, n_no, med_yes, med_no
    for pop, g in grouped:
        yes_vals = g.loc[g["response"] == "yes", "percentage"].dropna()
        no_vals = g.loc[g["response"] == "no", "percentage"].dropna()
        n_yes = int(yes_vals.shape[0])
        n_no = int(no_vals.shape[0])
        if n_yes == 0 or n_no == 0:
            p = np.nan
        else:
            # two-sided test; use continuity False for exact small-sample if possible
            try:
                stat = mannwhitneyu(yes_vals, no_vals, alternative="two-sided")
                p = float(stat.pvalue)
            except Exception:
                p = np.nan
        med_yes = float(np.nanmedian(yes_vals)) if n_yes else np.nan
        med_no = float(np.nanmedian(no_vals)) if n_no else np.nan
        pvals.append(p if not np.isnan(p) else 1.0)
        interim.append((pop, p, n_yes, n_no, med_yes, med_no))

    padj = _bh_adjust(pvals)
    for (pop, p, n_yes, n_no, med_yes, med_no), q in zip(interim, padj):
        yes_vals = df.loc[(df["population"] == pop) & (df["response"] == "yes"), "percentage"].dropna()
        no_vals = df.loc[(df["population"] == pop) & (df["response"] == "no"), "percentage"].dropna()
        eff = _cliffs_delta(yes_vals, no_vals)
        results.append(
            TestResult(
                population=str(pop),
                n_yes=n_yes,
                n_no=n_no,
                median_yes=med_yes,
                median_no=med_no,
                effect_size=float(eff) if eff is not None else np.nan,
                p_value=float(p) if not np.isnan(p) else np.nan,
                p_adj=float(q),
                significant=bool(q < alpha) if not np.isnan(q) else False,
            )
        )

    results_df = pd.DataFrame([r.__dict__ for r in results]).sort_values("p_adj")
    return results_df, figures


def subset_baseline_melanoma_pbmc_miraclib(conn: sqlite3.Connection) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # Base subset table of samples meeting criteria
    base_query = """
        SELECT * FROM samples
        WHERE condition = 'melanoma'
          AND sample_type = 'PBMC'
          AND treatment = 'miraclib'
          AND time_from_treatment_start = 0
    """
    subset_df = pd.read_sql_query(base_query, conn)

    by_project = (
        subset_df.groupby("project")["sample"].nunique().rename("num_samples").reset_index()
    )

    # responder/non-responder counts at subject level
    subj_first = subset_df.dropna(subset=["subject"]).drop_duplicates(subset=["subject"]).copy()
    by_response = subj_first.groupby("response")["subject"].nunique().rename("num_subjects").reset_index()

    by_sex = subj_first.groupby("sex")["subject"].nunique().rename("num_subjects").reset_index()

    return subset_df, by_project, by_response, by_sex


def run_streamlit_app():

    st.set_page_config(page_title="Drug Analysis - Teiko", layout="wide")
    st.title("Drug Analysis Dashboard")
    st.caption(
        "Analyzing immune cell populations across samples. Dataset includes melanoma and carcinoma; key analysis focuses on melanoma responders vs non-responders on miraclib."
    )

    with st.sidebar:
        st.header("Data Setup")
        csv_path = st.text_input("CSV path", value=CSV_DEFAULT_PATH)
        db_path = st.text_input("SQLite DB path", value=DATABASE_DEFAULT_PATH)
        col1, col2 = st.columns(2)
        with col1:
            init_clicked = st.button("Initialize DB")
        with col2:
            load_clicked = st.button("Load CSV → DB")

    conn = get_db_connection(db_path)

    if init_clicked:
        initialize_database(conn)
        st.success("Database initialized.")

    if load_clicked:
        initialize_database(conn)
        load_csv_to_db(csv_path, conn)
        st.success("CSV loaded into database.")

    st.subheader("Part 2: Summary Table")
    if st.button("Compute Summary"):
        summary_df = compute_summary_table(conn)
        if summary_df.empty:
            st.warning("No data available. Initialize DB and load CSV first.")
        else:
            st.dataframe(summary_df.head(50))
            st.download_button(
                label="Download full summary CSV",
                data=summary_df.to_csv(index=False).encode("utf-8"),
                file_name="summary_table.csv",
                mime="text/csv",
            )

            st.subheader("Part 3: Responders vs Non-Responders (Melanoma, PBMC, Miraclib)")
            results_df, figures = analyze_responders_vs_nonresponders(summary_df)
            if figures:
                for fig in figures:
                    st.plotly_chart(fig, use_container_width=True)
            st.dataframe(results_df)

            sig = results_df.loc[results_df["significant"]]
            if not sig.empty:
                st.success(
                    f"Significant differences (BH FDR<0.05): {', '.join(sig['population'].tolist())}"
                )
            else:
                st.info("No significant differences at FDR < 0.05.")

            st.subheader("Part 4: Baseline Melanoma PBMC Samples on Miraclib")
            subset_df, by_project, by_response, by_sex = subset_baseline_melanoma_pbmc_miraclib(conn)
            st.markdown("Samples meeting criteria:")
            st.dataframe(subset_df)

            st.markdown("Counts by project:")
            st.dataframe(by_project)

            st.markdown("Subject counts by response:")
            st.dataframe(by_response)

            st.markdown("Subject counts by sex:")
            st.dataframe(by_sex)


def run_cli(csv_path: str, db_path: str) -> None:
    print("Initializing database and loading CSV...")
    conn = get_db_connection(db_path)
    initialize_database(conn)
    load_csv_to_db(csv_path, conn)
    print("Loaded.")

    print("\nComputing summary table...")
    summary_df = compute_summary_table(conn)
    print(summary_df.head(10).to_string(index=False))

    print("\nAnalyzing responders vs non-responders (melanoma, PBMC, miraclib)...")
    results_df, _ = analyze_responders_vs_nonresponders(summary_df)
    print(results_df.to_string(index=False))

    print("\nSubset analysis: baseline melanoma PBMC samples on miraclib")
    subset_df, by_project, by_response, by_sex = subset_baseline_melanoma_pbmc_miraclib(conn)
    print("\nSamples meeting criteria:")
    print(subset_df.to_string(index=False))
    print("\nCounts by project:")
    print(by_project.to_string(index=False))
    print("\nSubject counts by response:")
    print(by_response.to_string(index=False))
    print("\nSubject counts by sex:")
    print(by_sex.to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description="Drug analysis toolkit and dashboard")
    parser.add_argument("--csv", dest="csv_path", default=CSV_DEFAULT_PATH, help="Path to cell-count.csv")
    parser.add_argument("--db", dest="db_path", default=DATABASE_DEFAULT_PATH, help="Path to SQLite DB file")
    subparsers = parser.add_subparsers(dest="command")

    subparsers.add_parser("init-db", help="Initialize the database schema")
    subparsers.add_parser("load", help="Load CSV into the database (initializes if needed)")
    subparsers.add_parser("summary", help="Compute and print the summary table head")
    subparsers.add_parser("analyze", help="Run responder vs non-responder analysis")
    subparsers.add_parser("subset", help="Run baseline subset analysis (Part 4)")
    subparsers.add_parser("run-all", help="Run init, load, summary, analyze, subset in sequence")

    args = parser.parse_args()

    if args.command == "init-db":
        conn = get_db_connection(args.db_path)
        initialize_database(conn)
        print("Database initialized at:", args.db_path)
        return
    if args.command == "load":
        conn = get_db_connection(args.db_path)
        initialize_database(conn)
        load_csv_to_db(args.csv_path, conn)
        print("CSV loaded into DB.")
        return
    if args.command == "summary":
        conn = get_db_connection(args.db_path)
        df = compute_summary_table(conn)
        if df.empty:
            print("No data. Initialize DB and load CSV first.")
        else:
            print(df.head(50).to_string(index=False))
        return
    if args.command == "analyze":
        conn = get_db_connection(args.db_path)
        df = compute_summary_table(conn)
        results, _ = analyze_responders_vs_nonresponders(df)
        print(results.to_string(index=False))
        return
    if args.command == "subset":
        conn = get_db_connection(args.db_path)
        subset_df, by_project, by_response, by_sex = subset_baseline_melanoma_pbmc_miraclib(conn)
        print("\nSamples meeting criteria:")
        print(subset_df.to_string(index=False))
        print("\nCounts by project:")
        print(by_project.to_string(index=False))
        print("\nSubject counts by response:")
        print(by_response.to_string(index=False))
        print("\nSubject counts by sex:")
        print(by_sex.to_string(index=False))
        return
    if args.command == "run-all" or args.command is None:
        # default: run the full pipeline in CLI mode
        run_cli(args.csv_path, args.db_path)


if __name__ == "__main__":
    main()


