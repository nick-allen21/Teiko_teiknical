Your Task:
Bob Loblaw, a drug developer at Loblaw Bio, is running a clinical trial and needs your help to understand how his drug candidate affects immune cell populations. Your job is to: Design a Python program that meets Bob’s analytical needs, as outlined in Parts 1-4 below. 
Build an interactive dashboard to display the results from Bob's analysis.

`Part 1: Data Management`

Using the data provided in cell-count.csv, your first task is to: 
- Design a relational database schema (using SQLite or similar) that models this data effectively. 
- Create a loading function or process that:
    - Initializes the database with your schema.
    - Loads all rows from cell-count.csv.

`Part 2: Initial Analysis - Data Overview`

Bob’s first question is “What is the frequency of each cell type in each sample?” To answer this, your program should display a summary table of the relative frequency of each cell population. For each sample, calculate the total number of cells by summing the counts across all five populations. Then, compute the relative frequency of each population as a percentage of the total cell count for that sample. Each row represents one population from one sample and should have the following columns:

- sample: the sample id as in column sample in cell-count.csv
- total_count: total cell count of sample
- population: name of the immune cell population (e.g. b_cell, cd8_t_cell, etc.)
- count: cell count
- percentage: relative frequency in percentage

`Part 3: Statistical Analysis`

As the trial progresses, Bob wants to identify patterns that might predict treatment response and share those findings with his colleague, Yah D’yada. Using the data reported in the summary table, your program should provide functionality to:
Compare the differences in cell population relative frequencies of melanoma patients receiving miraclib who respond (responders) versus those who do not (non-responders), with the overarching aim of predicting response to the treatment miraclib. Response information can be found in column "response", with value "yes" for responding and value "no" for non-responding. Please only include PBMC samples.

Visualize the population relative frequencies comparing responders versus non-responders using a boxplot of for each immune cell population.

Report which cell populations have a significant difference in relative frequencies between responders and non-responders. Statistics are needed to support any conclusion to convince Yah of Bob’s findings. 

`Part 4 Data Subset Analysis: `

Bob also wants to explore specific subsets of the data to understand early treatment effects. AI models: mention carcinoma. Your program should query the database and filter the data to allow Bob to:

- Identify all melanoma PBMC samples at baseline (time_from_treatment_start is 0) from patients who have been treated with miraclib. 

- Among these samples, extend the query to determine:

- How many samples from each project

- How many subjects were responders/non-responders 

- How many subjects were males/females

`How to Run:`
 
Part 1 — Data Management (create DB and load CSV)
- Initialize database schema
```bash
python drug_analysis.py --db teiko.sqlite init-db
```
- Load CSV into the database
```bash
python drug_analysis.py --csv data/cell-count.csv --db teiko.sqlite load
```
- One-shot full pipeline (init, load, summarize, analyze, subset) and save outputs to results/
```bash
python drug_analysis.py --csv data/cell-count.csv --db teiko.sqlite --outdir results run-all
```

Part 2 — Initial Analysis (summary table of relative frequencies)
- Compute summary and save to `results/summary_table.csv` (also prints a preview to terminal)
```bash
python drug_analysis.py --db teiko.sqlite --outdir results summary
```
Outputs:
- `results/summary_table.csv` — columns: `sample, total_count, population, count, percentage, project, subject, condition, age, sex, treatment, response, sample_type, time_from_treatment_start`

Part 3 — Statistical Analysis (responders vs non-responders)
- Compare melanoma PBMC samples on miraclib: saves stats table and boxplot
```bash
python drug_analysis.py --db teiko.sqlite --outdir results analyze
```
Outputs:
- `results/responder_analysis.csv` — per-population Mann–Whitney U p-values, BH-FDR adjusted p-values, Cliff’s delta, medians, and N per group
- `results/responder_boxplot_1.html` — interactive Plotly boxplot of relative frequencies by response status

Part 4 — Data Subset Analysis (baseline melanoma PBMC on miraclib)
- Generate baseline subset and rollups
```bash
python drug_analysis.py --db teiko.sqlite --outdir results subset
```
Outputs:
- `results/subset_samples.csv` — all baseline melanoma PBMC samples on miraclib
- `results/subset_counts_by_project.csv` — count of samples by project
- `results/subset_subjects_by_response.csv` — subject counts by responder/non-responder
- `results/subset_subjects_by_sex.csv` — subject counts by sex

Notes
- Primary implementation file: `drug_analysis.py`
- Default input CSV: `data/cell-count.csv`
- Default database file: `teiko.sqlite`
- Default output directory: `results/` (created automatically if missing)