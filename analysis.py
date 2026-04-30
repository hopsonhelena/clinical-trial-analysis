"""
analysis.py 

Analysis for Bob's immune cell project (parts 2-4).

Outputs: 
- output/pt2_cell_freq.csv
Summary table with relative frequencies of each cell population with samples as rows.

- output/pt3_responder_yn_stats.csv
Statistical comparison of differences in cell population frequencies.

- output/pt3_responder_yn_boxplot.png
Boxplots showing cell population frequencies by responder status.

- output/pt4_filtered_samples.csv
Filtered sample list (melanoma, PBMC, time=0, miraclib treatment) used for part 4 analysis.

- output/pt4_samples_per_project.csv
Number of samples per project.

- output/pt4_response_subject.csv
Number of responder vs non-responder subjects.

- output/pt4_sex_subject.csv
Number of male vs female subjects.

Run:
    python analysis.py
"""

import sqlite3
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

DB_PATH = "clinical_trial.db"
OUTPUT_DIR = "output"

# Check if output directory exists, if not, create it
def ensure_output_dir():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)


# Part 2: Initial Analysis - Data Overview

def create_frequency_table(conn):
    """
    For each sample, calculate the total number of cells by summing the counts across all five populations.
    Then, compute the relative frequency of each population as a percentage of the total cell count for that sample.
    
    Returns a DataFrame with columns: 
        sample, total_count, population, count, percentage
    """

    query = """
    SELECT 
        sample, 
        population,
        count,
        SUM(count) OVER (PARTITION BY sample) AS total_count
    FROM cell_counts
    """

    df = pd.read_sql_query(query, conn)

    #  Calculate total cell count for each sample 
    df["percentage"] = ((df["count"] / df["total_count"]) * 100).round(2)
    return df

# Part 3: Statistical Analysis 

def responders_vs_nonresponders_stats(conn, freq_df):
    """
    Compare differences in cell population frequencies between responder and non-responder melanoma patients
    receiving miraclib using PBMC samples only. 

    Returns a dictionary containing:
    - stats: statistical summary
    - data: merged dataset for plotting
    """

    query = """
    SELECT 
        sample,
        response
    FROM samples
    WHERE condition = 'melanoma' 
    AND treatment = 'miraclib' 
    AND sample_type = 'PBMC'
    """

    metadata_pt3 = pd.read_sql_query(query, conn)

    # Merge metadata with frequency table to get response
    merged = metadata_pt3.merge(freq_df, on = "sample")

    # Map population names for display
    population_map = {
    "b_cell": "B cell",
    "cd8_t_cell": "CD8 T cell",
    "cd4_t_cell": "CD4 T cell",
    "nk_cell": "NK cell",
    "monocyte": "Monocyte"
    }

    merged["population"] = merged["population"].replace(population_map)

    # Summarize statistics with t-test for each population 
    results = []
    populations = merged['population'].unique()

    for pop in populations:
        pop_data = merged[merged['population'] == pop]
        responders = pop_data[pop_data['response'] == 'yes']['percentage']
        non_responders = pop_data[pop_data['response'] == 'no']['percentage']

        # Get mean and mean difference for responders and non-responders
        mean_responders= responders.mean()
        mean_nonresponders = non_responders.mean()
        mean_diff = mean_responders - mean_nonresponders

        # Use Welch t-test if both groups have more than 30 samples, otherwise use Mann-Whitney U test
        if len(responders) > 30 and len(non_responders) > 30:
            t_stat, p_value = stats.ttest_ind(responders, non_responders, equal_var=False)
            test_used = "Welch's t-test"

        else:
            t_stat, p_value = stats.mannwhitneyu(responders, non_responders, alternative='two-sided')
            test_used = "Mann-Whitney U test"

    
        results.append({
            "population": pop,
            "p_value": p_value,
            "mean_responders": mean_responders,
            "mean_nonresponders": mean_nonresponders,
            "mean_diff": mean_diff,
            "statistic": t_stat,
            "test": test_used
        })

    stats_df = pd.DataFrame(results).sort_values("p_value").round(3)

    # Save stats table to CSV
    stats_df.to_csv(os.path.join(OUTPUT_DIR, "pt3_responder_yn_stats.csv"), index=False)
    
    # Print the stats table
    print("\nResponder vs Non-responder Immune Cell Frequencies\n")
    print(stats_df.to_string(index=False))

    # Report which cell populations have significant differences in relative frequency between responders and non-responders 
    significant = stats_df[stats_df["p_value"] < 0.05]
    print("\nSignificant differences:\n")

    if significant.empty:
        print("No statistically significant differences found (p >= 0.05 for all populations).")
    else:
        print(significant.to_string(index=False))

    return {
        "stats": stats_df,
        "data": merged
    }


# Boxplot 
def plot_boxplot(analysis):
    """
    Boxplot of immune cell frequencies for responders vs non-responders for each population.
    Saves the plot as a PNG file in the output directory.
    """
    data = analysis["data"]
    stats_df = analysis["stats"]

    # plot
    g = sns.catplot(
        data=data,
        x="response",
        y="percentage",
        col="population",
        kind="box",
        col_wrap=3,
        height=4,
        boxprops={'alpha': 0.4}, 
        showfliers=False
     )

    # Set y limit 0-100 
    for ax in g.axes.flat:
        ax.set_ylim(0, 100)

    # Loop through each population subplot
    populations = data['population'].unique()

    for ax, pop in zip(g.axes.flat, populations):
        
        pop_data = data[data['population'] == pop]

        # Add points, randomly sampled to reduce clutter
        sampled = pop_data.sample(n=80, random_state=1)

        sns.stripplot(
            data=sampled,
            x="response",
            y="percentage",
            ax=ax,
            color="black",
            alpha=0.3,
            size=4,
            jitter=0.25
        )

        # Annotate p-value from stats_df
        subset = stats_df[stats_df["population"] == pop]
        if not subset.empty:
            p_value = subset["p_value"].values[0]
            p_text = f"p = {p_value:.3f}"

        ax.text(
            0.95, 0.95,
            p_text,
            transform=ax.transAxes,
            ha='right',
            va='top',
            fontsize=10
        )

    # Axis labels and title
    g.set_axis_labels("Response", "Relative Frequency (%)")
    g.set_titles("{col_name}")
    g.fig.suptitle("Immune Cell Frequencies in \nMelanoma Patients receiving Miraclib (PBMC samples)", y=1.05)

    # Save plot
    plt.savefig(os.path.join(OUTPUT_DIR, "pt3_responder_yn_boxplot.png"), bbox_inches="tight")
    print(f"\nBoxplot saved to {OUTPUT_DIR}/pt3_responder_yn_boxplot.png")   



# Part 4: Data Subset Analysis
"""
Filter the database to allow Bob to identify all melanoma PBMC samples at baseline (time_from_treatment_start is 0) 
from patients who have been treated with miraclib.

Returns answers to Bob's questions:
- How many samples from each project
- How many subjects were responders/non-responders
- How many subjects were males/females
"""

# Base filter: melanoma, miraclib, time=0, PBMC samples
def base_filter(conn):
    query = """
    SELECT *
    FROM samples
    WHERE condition = 'melanoma'
    AND treatment = 'miraclib'
    AND sample_type = 'PBMC'
    AND time_from_treatment_start = 0
    """

    df = pd.read_sql_query(query, conn)

    output_path = os.path.join(OUTPUT_DIR, "pt4_filtered_samples.csv")
    df.to_csv(output_path, index=False)

    print("\nAll melanoma PBMC samples at baseline from patients who have been treated with miraclib "
          f"Saved to {output_path}")
    
    return df

# Among these samples, extend the query to determine:

# How many samples from each project
def count_samples_per_project(conn):
    query = """
    SELECT 
        project,
        COUNT(sample) AS sample_count
    FROM samples 
    WHERE condition = 'melanoma'
    AND treatment = 'miraclib'
    AND sample_type = 'PBMC'
    AND time_from_treatment_start = 0
    GROUP BY project
    """
    df = pd.read_sql_query(query, conn)

    output_path = os.path.join(OUTPUT_DIR, "pt4_samples_per_project.csv")
    df.to_csv(output_path, index=False)

    print("\nNumber of samples from each project:")
    print(df.to_string(index=False))
    print(f"\nSaved to {output_path}")

    return df

# How many subjects were responders/non-responders
def count_responders_nonresponders(conn):
    query = """
    SELECT 
        response,
        COUNT(DISTINCT subject) AS subject_count
    FROM samples 
    WHERE condition = 'melanoma'
    AND treatment = 'miraclib'
    AND sample_type = 'PBMC'
    AND time_from_treatment_start = 0
    GROUP BY response
    """
    df = pd.read_sql_query(query, conn)

    output_path = os.path.join(OUTPUT_DIR, "pt4_response_subject.csv")
    df.to_csv(output_path, index=False)

    print(f"\nNumber of subjects who were responders vs non-responders:")
    print(df.to_string(index=False))
    print(f"\nSaved to {output_path}")

    return df

# How many subjects were males/females
def count_males_females(conn):
    query = """
    SELECT 
        sub.sex,
        COUNT(DISTINCT sub.subject) AS subject_count
    FROM subjects sub
    JOIN samples s ON s.subject = sub.subject
    WHERE s.condition = 'melanoma'
    AND s.treatment = 'miraclib'
    AND s.sample_type = 'PBMC'
    AND s.time_from_treatment_start = 0
    GROUP BY sub.sex
    """
    df = pd.read_sql_query(query, conn)

    output_path = os.path.join(OUTPUT_DIR, "pt4_sex_subject.csv")
    df.to_csv(output_path, index=False)

    print("\nNumber of subjects who were males vs females:")
    print(df.to_string(index=False))
    print(f"\nSaved to {output_path}")
    
    return df


# Considering Melanoma males, what is the average number of B cells for responders at time=0?
# Use two decimals 
def average_b_cells_males_responders(conn):
    query = """
    SELECT AVG(c.count) AS avg_b_cells
    FROM samples s
    JOIN subjects sub ON s.subject = sub.subject
    JOIN cell_counts c ON s.sample = c.sample
    WHERE s.condition = 'melanoma'
    AND s.time_from_treatment_start = 0
    AND sub.sex = 'M'
    AND s.response = 'yes'
    AND c.population = 'b_cell'
    """
    df = pd.read_sql_query(query, conn)
    
    avg_b_cells = df['avg_b_cells'].iloc[0]
    
    if pd.isna(avg_b_cells):
        print("\nNo matching samples found for male melanoma responders at baseline.")
    else:
        print(f"\nAverage number of B cells for male melanoma responders at baseline: {avg_b_cells:.2f}")
    
    return avg_b_cells


# Main 

def main():
    
    # ensure output directory exists 
    ensure_output_dir()

    # connect to database
    conn = sqlite3.connect(DB_PATH)

    try:

        # Part 2
        print()
        print("\nRunning Part 2: Creating frequency table...")
        freq_df = create_frequency_table(conn)
        freq_df.to_csv(os.path.join(OUTPUT_DIR, "pt2_cell_freq.csv"), index=False)
        print(f"\nSummary table with relative frequencies saved to {OUTPUT_DIR}/pt2_cell_freq.csv")

        # Part 3
        print()
        print("\nRunning Part 3: Statistical analysis of responders vs non-responders...")
        analysis = responders_vs_nonresponders_stats(conn, freq_df)
        plot_boxplot(analysis)

        # Part 4
        print()
        print("\nRunning Part 4: Data subset analysis...")
        base_filter(conn)
        count_samples_per_project(conn)
        count_responders_nonresponders(conn)
        count_males_females(conn)
        average_b_cells_males_responders(conn)

        print("\nPipeline completed successfully!")

    finally:
        conn.close()

if __name__ == "__main__":
    main()