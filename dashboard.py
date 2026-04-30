import streamlit as st
import pandas as pd 
import sqlite3
import seaborn as sns
import matplotlib.pyplot as plt
import os
from analysis import create_frequency_table
from scipy import stats


st.title("Bob's Immune Cell Analysis Dashboard")

DB_PATH = "clinical_trial.db"
FREQ_PATH = "output/celltype_frequency_table.csv"

@st.cache_resource
def get_connection():
    return sqlite3.connect(DB_PATH)


# Part 2 - Function to load data 
@st.cache_data
def load_frequency_data():
    if os.path.exists(FREQ_PATH):
        # Read from local csv if exists (avoid re-running analysis)
        return pd.read_csv(FREQ_PATH)
    else:
        # If not, run the analysis to generate it
        conn = get_connection()
        df = create_frequency_table(conn)
        df.to_csv(FREQ_PATH, index=False)
        return df
    

# Part 3 - Function to load data 
@st.cache_data
def load_analysis_data(condition, treatment, sample_type):
    conn = get_connection()

    # join sample and count tables and
    query = f"""
    SELECT 
        s.sample,
        s.response,
        s.condition,
        s.treatment,
        c.population,
        c.count,
        SUM(c.count) OVER (PARTITION BY s.sample) AS total_count
    FROM samples s
    JOIN cell_counts c ON s.sample = c.sample
    WHERE s.condition = ?
    AND s.treatment = ?
    AND s.sample_type = ?
    """

    df = pd.read_sql(query, conn, params=[condition, treatment, sample_type])

    # Compute percentages
    df["percentage"] = ((df["count"] / df["total_count"]) * 100).round(2)

    # Clean population labels
    population_map = {
        "b_cell": "B cell",
        "cd8_t_cell": "CD8 T cell",
        "cd4_t_cell": "CD4 T cell",
        "nk_cell": "NK cell",
        "monocyte": "Monocyte"
    }
    df["population"] = df["population"].replace(population_map)

    return df

# Part 3 - Function to calculate statistics table 
@st.cache_data
def compute_stats(df):
    results = []

    populations = df["population"].unique()

    for pop in populations:
        pop_data = df[df["population"] == pop]

        responders = pop_data[pop_data["response"] == "yes"]["percentage"]
        non_responders = pop_data[pop_data["response"] == "no"]["percentage"]

        # skip empty groups
        if len(responders) == 0 or len(non_responders) == 0:
            continue

        mean_r = responders.mean()
        mean_nr = non_responders.mean()
        mean_diff = mean_r - mean_nr

        if len(responders) > 30 and len(non_responders) > 30:
            stat, p = stats.ttest_ind(responders, non_responders, equal_var=False)
            test = "Welch t-test"
        else:
            stat, p = stats.mannwhitneyu(responders, non_responders)
            test = "Mann-Whitney U"

        results.append({
            "population": pop,
            "p_value": round(p, 3),
            "mean_responders": round(mean_r, 2),
            "mean_nonresponders": round(mean_nr, 2),
            "mean_diff": round(mean_diff, 2),
            "test": test
        })

    return pd.DataFrame(results).sort_values("p_value")

# Part 3 - Function to add filters to boxplots 
@st.cache_data
def load_filter_options():
    conn = get_connection()
    conditions = pd.read_sql("SELECT DISTINCT condition FROM samples", conn)
    treatments = pd.read_sql("SELECT DISTINCT treatment FROM samples", conn)
    sample_types = pd.read_sql("SELECT DISTINCT sample_type FROM samples", conn)
    return conditions, treatments, sample_types


tab2, tab3, tab4 = st.tabs([
    "Part 2: Immune cell frequency", 
    "Part 3: Responder Analysis", 
    "Part 4: Subset"
    ])


### Part 2 tab ###
with tab2:
    st.header("Part 2: Summary of Immune Cell Frequencies")

    df = load_frequency_data()

    # Get population options
    populations = sorted(df["population"].unique())

    # Multiselect (same style as Part 3)
    selected_pops = st.multiselect(
        "Select immune cell populations",
        options=populations,
        default=populations
    )

    # Guard: nothing selected
    if len(selected_pops) == 0:
        st.warning("Please select at least one immune cell population.")
        st.stop()

    # Filter data
    filtered_df = df[df["population"].isin(selected_pops)]
    filtered_df = filtered_df.sort_values("population")

    # Show filtered table
    st.dataframe(filtered_df, hide_index=True)
    
    # Option to show full dataset
    with st.expander("Show full dataset"):
        st.dataframe(df, hide_index=True)


### Part 3 tab ###
with tab3:
    st.header("Part 3: Responder vs Non-Responder Analysis")

    # Add filters 
    conditions, treatments, sample_types = load_filter_options()

    condition = st.selectbox("Condition", conditions["condition"])
    treatment = st.selectbox("Treatment", treatments["treatment"])
    sample_type = st.selectbox("Sample Type", sample_types["sample_type"])


    # Load filtered data
    df = load_analysis_data(condition, treatment, sample_type)

    # Add population selector 
    populations = sorted(df["population"].unique())

    selected_pops = st.multiselect(
        "Select immune cell populations",
        options=populations,
        default=populations   # show all by default
    )

    # Guard: healthy and nothing selected 
    if condition.lower() == "healthy":
        st.error("Responder analysis is not applicable for healthy condition.")
        st.stop()

    if len(selected_pops) == 0:
        st.warning("Please select at least one immune cell population.")
        st.stop()


    # Filter data for populations 
    df = df[df["population"].isin(selected_pops)]
    
    # Dynamic layout
    col_wrap = min(5, len(selected_pops)) if selected_pops else 1

    # Plot
    g = sns.catplot(
        data=df,
        x="response",
        y="percentage",
        col="population",
        col_order=selected_pops,
        kind="box",
        col_wrap=col_wrap,
        height = 4,
        sharey=True,
        showfliers=False
    )

    plt.xticks(rotation=45)

    # Add points
    for ax, pop in zip(g.axes.flat, df["population"].unique()):
        pop_data = df[df["population"] == pop]

        sns.stripplot(
            data=pop_data,
            x="response",
            y="percentage",
            ax=ax,
            alpha=0.3,
            size=3,
            jitter=0.2
        )

    # Format all axes 
    for ax in g.axes.flat:
        ax.set_ylim(0, 100)
        ax.set_xlabel("Response")
        ax.set_ylabel("Relative Frequency (%)")

    # Render
    st.pyplot(g.fig)

    # Add stats table 
    stats_df = compute_stats(df)

    # Display table
    st.subheader("Statistical Summary")

    st.dataframe(stats_df, hide_index=True)

    # Filter significant rows 
    sig_only = st.checkbox("Show only significant (p < 0.05)")

    if sig_only:
        stats_df = stats_df[stats_df["p_value"] < 0.05]

    st.dataframe(stats_df, hide_index=True)