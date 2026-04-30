---
title: "Bob's immune cell analysis"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# How to run 
make setup
make pipeline
make dashboard

# Database schema 
The database is normalized into separate subjects (patient level data), samples (biological sample data), and cell_counts (experimental data) tables to reduce redundancy and improve scalability for larger datasets. 

# Code structure 

load_data.py initializes database schema 

analysis.py contains the statistics and plots for specific questions (parts 2-4). 
Welch’s t-test was used when both groups had sufficient sample size (>30), otherwise a Mann–Whitney U test was applied

dashboard.py is the interactive dashboard. 

## Key Findings
CD4 T cell frequencies were significantly different between responders and non-responders (p < 0.05), however the absolute difference in means was small (~0.6%) and may not be biologically meaningful. 

# Dashboard link 
To access the dashboard: 
# clinical-trial-analysis
