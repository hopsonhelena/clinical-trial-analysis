"""
load_dashboard.py

Streamlit dashboard for immune cell analysis.

Run:
    python -m streamlit run load_dashboard.py
"""

import io
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sqlite3
import streamlit as st
from analysis import create_frequency_table
from scipy import stats


# --- Constants ---

DB_PATH = "clinical_trial.db"
FREQ_PATH = "output/celltype_frequency_table.csv"

POPULATION_LABELS = {
    "b_cell": "B cell",
    "cd8_t_cell": "CD8 T cell",
    "cd4_t_cell": "CD4 T cell",
    "nk_cell": "NK cell",
    "monocyte": "Monocyte",
}

GROUP_BY_OPTIONS = [
    "sex",
    "response",
    "condition",
    "treatment",
    "time_from_treatment_start",
    "project",
    "age",
]


# --- Database helpers ---


def get_connection():
    return sqlite3.connect(DB_PATH)


@st.cache_data
def get_age_bounds():
    conn = get_connection()
    result = pd.read_sql(
        "SELECT MIN(age) AS min_age, MAX(age) AS max_age FROM subjects", conn
    )
    return int(result["min_age"][0]), int(result["max_age"][0])


# --- Database loaders ---


@st.cache_data
def load_pt2_data() -> pd.DataFrame:
    """Load or generate the cell-type frequency table."""
    if os.path.exists(FREQ_PATH):
        return pd.read_csv(FREQ_PATH)
    conn = get_connection()
    df = create_frequency_table(conn)
    df.to_csv(FREQ_PATH, index=False)
    return df


@st.cache_data
def load_pt3_data(condition: str, treatment: str, sample_type: str) -> pd.DataFrame:
    """Load per-sample cell counts with computed percentages for Part 3."""
    query = """
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
        WHERE s.condition   = ?
          AND s.treatment   = ?
          AND s.sample_type = ?
    
    """
    conn = get_connection()
    df = pd.read_sql(query, conn, params=[condition, treatment, sample_type])
    df["percentage"] = ((df["count"] / df["total_count"]) * 100).round(2)
    df["population"] = df["population"].replace(POPULATION_LABELS)
    return df


@st.cache_data
def load_pt4_data(
    condition=None,
    treatment=None,
    sample_type=None,
    time_from_treatment_start=None,
    project=None,
    response=None,
    sex=None,
    age_range=None,
) -> pd.DataFrame:
    """Load subject/sample data with optional filters for Part 4."""
    query = """
        SELECT
            s.sample,
            s.subject,
            s.project,
            s.sample_type,
            s.condition,
            s.treatment,
            s.response,
            s.time_from_treatment_start,
            sub.age,
            sub.sex
        FROM samples s
        JOIN subjects sub ON s.subject = sub.subject
        WHERE 1=1
    """
    params = []

    list_filters = [
        ("s.condition", condition),
        ("s.treatment", treatment),
        ("s.sample_type", sample_type),
        ("sub.sex", sex),
        ("s.response", response),
        ("s.time_from_treatment_start", time_from_treatment_start),
        ("s.project", project),
    ]
    for col, values in list_filters:
        if values:
            placeholders = ",".join(["?"] * len(values))
            query += f" AND {col} IN ({placeholders})"
            params += list(values)

    if age_range:
        query += " AND sub.age BETWEEN ? AND ?"
        params += list(age_range)

    conn = get_connection()
    return pd.read_sql(query, conn, params=params)


@st.cache_data
def load_pt3_filter_options():
    """Return distinct condition, treatment, and sample_type values for Part 3."""
    conn = get_connection()
    conditions = pd.read_sql("SELECT DISTINCT condition   FROM samples", conn)
    treatments = pd.read_sql("SELECT DISTINCT treatment   FROM samples", conn)
    sample_types = pd.read_sql("SELECT DISTINCT sample_type FROM samples", conn)
    return conditions, treatments, sample_types


@st.cache_data
def load_pt4_filter_options():
    """Return distinct filter values for Part 4."""
    conn = get_connection()
    conditions = pd.read_sql(
        "SELECT DISTINCT condition                 FROM samples", conn
    )
    treatments = pd.read_sql(
        "SELECT DISTINCT treatment                 FROM samples", conn
    )
    sample_types = pd.read_sql(
        "SELECT DISTINCT sample_type               FROM samples", conn
    )
    projects = pd.read_sql(
        "SELECT DISTINCT project                   FROM samples", conn
    )
    responses = pd.read_sql(
        "SELECT DISTINCT response                  FROM samples", conn
    )
    times = pd.read_sql("SELECT DISTINCT time_from_treatment_start FROM samples", conn)
    sexes = pd.read_sql("SELECT DISTINCT sex                       FROM subjects", conn)
    return conditions, treatments, sample_types, projects, sexes, responses, times


# --- Analysis helpers ---


@st.cache_data
def compute_pt3_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Run per-population responder vs. non-responder significance tests."""
    results = []
    for pop in df["population"].unique():
        pop_data = df[df["population"] == pop]
        responders = pop_data[pop_data["response"] == "yes"]["percentage"]
        non_responders = pop_data[pop_data["response"] == "no"]["percentage"]

        if responders.empty or non_responders.empty:
            continue

        if len(responders) > 30 and len(non_responders) > 30:
            stat, p = stats.ttest_ind(responders, non_responders, equal_var=False)
            test = "Welch t-test"
        else:
            stat, p = stats.mannwhitneyu(responders, non_responders)
            test = "Mann-Whitney U"

        results.append(
            {
                "population": pop,
                "p_value": round(p, 3),
                "mean_responders": round(responders.mean(), 2),
                "mean_nonresponders": round(non_responders.mean(), 2),
                "mean_diff": round(responders.mean() - non_responders.mean(), 2),
                "test": test,
            }
        )

    return pd.DataFrame(results).sort_values("p_value")


# --- Page layout ---

st.set_page_config(layout="wide")
st.title("Bob's Immune Cell Analysis Dashboard")

# Tabs
tab2, tab3, tab4 = st.tabs(
    ["Part 2: Immune cell frequency", "Part 3: Responder Analysis", "Part 4: Subset"]
)


# --- Tab2 - Immune cell frequencies ---

with tab2:
    st.header("Part 2: Summary of Immune Cell Frequencies")

    df2 = load_pt2_data()
    if df2.empty:
        st.warning("No frequency data available.")
        st.stop()

    populations = sorted(df2["population"].unique())
    selected_pops = st.multiselect(
        "Select immune cell populations",
        options=populations,
        default=populations,
    )

    if not selected_pops:
        st.warning("Please select at least one immune cell population.")
        st.stop()

    filtered_df2 = df2[df2["population"].isin(selected_pops)].sort_values("population")
    st.dataframe(filtered_df2, hide_index=True)

    with st.expander("Show full dataset"):
        st.dataframe(df2, hide_index=True)


# --- Tab3 - Responder vs Non-Responder Analysis ---

with tab3:
    st.header("Part 3: Responder vs Non-Responder Analysis")

    conditions_opts, treatments_opts, sample_types_opts = load_pt3_filter_options()

    cond3 = st.selectbox("Condition", conditions_opts["condition"])
    treat3 = st.selectbox("Treatment", treatments_opts["treatment"])
    samp_type3 = st.selectbox("Sample Type", sample_types_opts["sample_type"])

    df3 = load_pt3_data(cond3, treat3, samp_type3)
    if df3.empty:
        st.warning("No data available for this selection.")
        st.stop()

    if cond3.lower() == "healthy":
        st.error("Responder analysis is not applicable for healthy condition.")
        st.stop()

    populations3 = sorted(df3["population"].unique())
    selected_pops3 = st.multiselect(
        "Select immune cell populations",
        options=populations3,
        default=populations3,
    )

    if not selected_pops3:
        st.warning("Please select at least one immune cell population.")
        st.stop()

    df3 = df3[df3["population"].isin(selected_pops3)]
    col_wrap = min(5, len(selected_pops3))

    # Box plots
    g = sns.catplot(
        data=df3,
        x="response",
        y="percentage",
        col="population",
        col_order=selected_pops3,
        kind="box",
        col_wrap=col_wrap,
        height=6,
        aspect=1.2,
        sharey=True,
        showfliers=False,
    )

    for ax, pop in zip(g.axes.flat, selected_pops3):
        sns.stripplot(
            data=df3[df3["population"] == pop],
            x="response",
            y="percentage",
            ax=ax,
            alpha=0.3,
            size=3,
            jitter=0.2,
        )

    for ax in g.axes.flat:
        ax.set_ylim(0, 100)
        ax.set_xlabel("Response")
        ax.set_ylabel("Relative Frequency (%)")

    g.fig.set_size_inches(18, 6.5)
    g.fig.set_dpi(200)
    g.fig.tight_layout()

    st.pyplot(g.fig, use_container_width=True)

    # Download options
    FORMATS = {
        "PNG": ("image/png", "png"),
        "PDF": ("application/pdf", "pdf"),
    }
    fmt = st.radio("Download format", list(FORMATS.keys()), horizontal=True)
    mime, ext = FORMATS[fmt]

    buf = io.BytesIO()
    g.fig.savefig(
        buf,
        format=ext,
        bbox_inches="tight",
        **({"dpi": 150} if ext == "png" else {}),
    )
    buf.seek(0)
    st.download_button(
        label=f"Download plot as {fmt}",
        data=buf,
        file_name=f"responder_analysis_{cond3}_{treat3}_{samp_type3}.{ext}",
        mime=mime,
    )

    # Stats table
    st.subheader("Statistical Summary")
    stats_df = compute_pt3_stats(df3)
    if st.checkbox("Show only significant (p < 0.05)"):
        stats_df = stats_df[stats_df["p_value"] < 0.05]
    st.dataframe(stats_df, hide_index=True)


# --- Tab4 - Subset Analysis ---

with tab4:
    st.header("Part 4: Subset Tables")

    (
        cond_opts4,
        treat_opts4,
        samp_opts4,
        proj_opts4,
        sex_opts4,
        resp_opts4,
        time_opts4,
    ) = load_pt4_filter_options()

    min_age, max_age = get_age_bounds()

    col_a, col_b = st.columns(2)
    with col_a:
        condition4 = st.multiselect("Condition", cond_opts4["condition"].tolist())
        treatment4 = st.multiselect("Treatment", treat_opts4["treatment"].tolist())
        sample_type4 = st.multiselect("Sample Type", samp_opts4["sample_type"].tolist())
        response4 = st.multiselect("Response", resp_opts4["response"].tolist())
    with col_b:
        project4 = st.multiselect("Project", proj_opts4["project"].tolist())
        sex4 = st.multiselect("Sex", sex_opts4["sex"].tolist())
        time4 = st.multiselect(
            "Time from treatment start",
            sorted(time_opts4["time_from_treatment_start"].tolist()),
        )
        age_range4 = st.slider(
            "Age range",
            min_value=min_age,
            max_value=max_age,
            value=(min_age, max_age),
        )

    df4 = load_pt4_data(
        condition=condition4 or None,
        treatment=treatment4 or None,
        sample_type=sample_type4 or None,
        time_from_treatment_start=time4 or None,
        project=project4 or None,
        response=response4 or None,
        sex=sex4 or None,
        age_range=age_range4,
    )

    if df4.empty:
        st.warning("No data matches your filters.")
        st.stop()

    st.dataframe(df4)

    group_cols = st.multiselect("Group by", GROUP_BY_OPTIONS)
    if group_cols:
        grouped = df4.groupby(group_cols).size().reset_index(name="count")
        st.dataframe(grouped)
