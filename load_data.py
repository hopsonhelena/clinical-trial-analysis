"""
load_data.py
Initialze SQLite data base and load cell-count.csv data (Part 1)

Schema design: 

Run:
    python load_data.py
"""

import sqlite3
import pandas as pd

DB_PATH = "clinical_trial.db"
CSV_PATH = "cell-count.csv"

def main():
    print("Initializing database...")

    # Connect to DB (creates file if not exists)
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    conn.execute("PRAGMA foreign_keys = ON;")

    try:
        # Drop existing tables 
        cursor.execute("DROP TABLE IF EXISTS cell_counts")
        cursor.execute("DROP TABLE IF EXISTS samples")
        cursor.execute("DROP TABLE IF EXISTS subjects")

        # Create tables

        # Subjects (patient info - one row per subjet/patient)
        cursor.execute("""
        CREATE TABLE subjects (
            subject TEXT PRIMARY KEY,
            age INTEGER,
            sex TEXT
            );         
        """)

        # Samples (experiment/sample info - one row per sample, linked to subject)
        cursor.execute("""
        CREATE TABLE samples (
            sample TEXT PRIMARY KEY,
            subject TEXT,
            project TEXT,
            sample_type TEXT,
            condition TEXT,
            treatment TEXT,
            response TEXT,
            time_from_treatment_start INTEGER,
            FOREIGN KEY (subject) REFERENCES subjects(subject) 
        );
        """)

        # Cell counts (long format, one row per population per sample)
        cursor.execute("""
        CREATE TABLE cell_counts (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sample TEXT,
            population TEXT,
            count INTEGER,
            FOREIGN KEY (sample) REFERENCES samples(sample) 
        );
        """)

        # Load CSV
        df = pd.read_csv(CSV_PATH)

        # Insert subjects into subjects table - what if a subject doesn't have information for treatmetnt response etc? We can insert NULL values for those fields, but we should ensure that the subject_id is still unique and serves as the primary key.
        subjects_df = df[["subject", "age", "sex"]].drop_duplicates()
        subjects_df.to_sql("subjects", conn, if_exists="append", index=False)

        # Insert samples into samples table
        samples_df = df[[
                "sample", "subject", "project", "condition",
                "treatment", "response", "sample_type",
                "time_from_treatment_start"
        ]]

        samples_df = samples_df.drop_duplicates()  # drop duplicates if any just in case

        samples_df.to_sql("samples", conn, if_exists="append", index=False)

        # Insert cell counts in long format into cell_counts table
        cell_cols = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

        df_long = df.melt(id_vars=["sample"], 
                          value_vars=cell_cols, 
                          var_name="population", 
                          value_name="count")

        df_long.to_sql("cell_counts", conn, if_exists="append", index=False)

        print("Database created and data loaded successfully!")

    finally:
        conn.close()


if __name__ == "__main__":
    main()
