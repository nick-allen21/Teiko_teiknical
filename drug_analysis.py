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

