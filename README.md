# Clinical Trial Immune Cell Analysis & Dashboard

## Overview

This project analyzes clinical trial data and immune cell populations for Bob.

---

## How to Run

Run the following commands:

* `make setup`: install dependencies
* `make pipeline`: run full analysis (data loading + outputs)
* `make dashboard`: launch interactive dashboard

---

## How to Access Dashboard (Codespaces)

1. In Codespaces, look at the bottom panel    
2. Click “Ports” tab   
3. You should see something like:  
**Port** 8501   
4. Under **Forwarded Address**, click the globe button  

---

## Database Schema

The SQLite database is normalized into three tables:

* **subjects**: patient-level demographic data (sex, age)
* **samples**: sample-level information (timepoint, treatment, etc.)
* **cell_counts**: immune cell population counts in long format (multiple rows per sample, one row per cell type)

### Rationale

* Reduces redundancy

  * subject level attributes (sex, age) are stored once in **subjects** 
  * sample metadata with one row per sample are stored in **samples**
  * measurements are stored in **cell_counts**

This avoids repeating subject information across multiple samples, reducing storage space and preventing inconsistencies.

* Improves scalability

  * New samples from existing subjects can be added to **samples** without duplicating or altering **subject** data
  * Additional cell populations or measurements can be added to **cell_counts** without schema changes
  * New metadata can be added to samples without affecting other tables 

---

## Code Structure
The project is organized into 3 main components:

* **load_data.py**
  Initializes the database schema, loads data from `cell-count.csv`.

* **analysis.py**
  Performs core analysis and generates output files for parts 2-4. 
  * Computes relative frequencies of immune cell populations
  * Performs statistical comparisons (responders vs non-responders)
  * Generates group summary tables
  
  Functions are written in a modular way to improve readability, maintainability, and enable reuse.  

* **dashboard.py**
  Provides a Streamlit-based dashboard for exploratory data analysis, while also enabling interactive reproduction of the outputs from `analysis.py`. 
  * Allows filtering, grouping, and visualization of results
  * One tab per analysis section (3 tabs for Parts 2-4)

  Design choices: 
  * Part 2 data is loaded from a precomputed CSV if available, or recomputed using analysis.py, since the analysis does not depend on user-driven filtering. 
  * Parts 3 and 4 query the database directly to support flexible, user-driven filtering.
  * Queries are executed in SQL to avoid loading unnecessary data into memory.
  * Querying the database directly ensures that results always reflect the most up-to-date data 
  * The dashboard is designed to reproduce the outputs from `analysis.py` interactively, while also enabling additional exploratory analyses
 
