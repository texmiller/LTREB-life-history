# LTREB-life-history

This repository includes code to analyze the effects of fungal endophyte symbiosis on the life histories of grass hosts, as described in the manuscript by Scherick et al. ("Fungal symbionts promote slower and steadier life histories of grass hosts"). The repository consists of three main folders for data storage and preparation ("data prep"), data analysis, modeling, and visualization ("analysis"), and storage of manuscript figures and tables ("manuscript"). The following scripts reproduce the analysis:

**data prep/LTREB_data_QAQC.R** – This script reads in raw data, implements all data wrangling, cleaning, and QA/QC checks, and outputs a tidy data file (ltreb_allspp_qaqc.csv) that is used for all subsequent analysis.

**analysis/age_endo_specific_vital_rates.R** – This script reads in the derived tidy data, applies age grouping (lumping advanced ages into an "old" age class), and fits hierarchical statistical models with calls to rstan (stan models are contained within the subfolder **analysis/Stan**). Fitted models are used to create figures for age-specific survival, age-specific fertility, seedling recruitment, and age of first reproduction. Finally, the script collects parameter values into demographic transition matrices and computes fitness and life history metrics, which are stored in lifehistorypost.csv.

**analysis/life_history_results.R** – This script reads posterior draws of the life history trait values, implements PCA, and analyzes and visualizes trait and fitness differences between S+ and S-.
