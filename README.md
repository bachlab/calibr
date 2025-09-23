# calibr: A Package for Experimental calibration using the Meta-Analysis of Bayes Factors

📌 **Overview**  
calibr is an R package that provides tools for experimental calibration, enabling researchers to estimate and update measurement method rankings across various datasets. 
We propose a principled approach via combining **Bayes Factors** for the object model, induced by a certain measurement, against the null model, which is simply an intercept (the worst possible measurement method). 
The package is particularly useful for combining evidence from multiple studies.

The infrastructure is a verbatim implementation of the methods in *Nikolakopoulos, S., & Ntzoufras, I. (2021). Meta analysis of Bayes factors. arXiv preprint arXiv:2103.13236*.
The logic underlying the choice of Bayes factors to indicate retrodictive validity is described in *Mancinelli, F., & Bach, D. (2025). Experiment-based calibration: inference and decision-making*.

---

## 🔧 **Installation**
To install calibr from GitLab, use the following command in R:

```r
# Install necessary dependencies
install.packages("devtools")  

# Install from GitLab
devtools::install_git("https://gitlab.com/bachlab/methods/calibr.git")
```

## 📊  Use case: computing a calibration summary
In this example, three datasets from different studies are analyzed in a between-subjects setting: the binary standard score (0/1) forms two groups. 
Per combination of study × standard score × measurement, we simulate 20 observations. “Measure_A” is better than “Measure_B”.

# 1. Data preparation

```r
library(dplyr)
set.seed(123)

# Parameters
n_studies <- 3
n_per_group <- 20
datasets <- paste0("Study", 1:n_studies)
measurements <- c("Measure_A", "Measure_B")
standard_scores <- c(0, 1)

# Simulate
df <- expand.grid(
  dataset = datasets,
  standard_score = standard_scores,
  measurement_name = measurements
) %>%
  rowwise() %>%
  do(data.frame(
    dataset = rep(.$dataset, n_per_group),
    standard_score = rep(.$standard_score, n_per_group),
    measurement_name = rep(.$measurement_name, n_per_group),
    measurement_value =
      if (.$standard_score == 0) {
        if (.$measurement_name == "Measure_A") rnorm(n_per_group, mean = 45, sd = 30) else rnorm(n_per_group, mean = 40, sd = 30)
      } else {
        if (.$measurement_name == "Measure_A") rnorm(n_per_group, mean = 70, sd = 30) else rnorm(n_per_group, mean = 55, sd = 30)
      }
  )) %>%
  ungroup()
```

# 2. Compute a calibration summary

```r
library(calibr)

result <- calibration_summary(
  data = df,
  dataset_col = "dataset",
  measurement_name_col = "measurement_name",
  measurement_values_col = "measurement_value",
  standard_scores_col = "standard_score"
)

# Main summary table (ordered by log BF)
print(result$summary_table)
# Matrix of pairwise log BF differences
print(result$log_bayes_factor_diffs)
# Short statement about the best measurement
cat(result$best_measurement, "\n")
# Subject overlap (if subject IDs were provided)
print(result$subject_overlap)
```

## 🙋 **FAQ**  

**I don’t understand the logic of subject overlap. When can I trust the computations?**  
`calibr` provides two functions to infer measurement rankings from data: `meta_bf_from_data` and `strict_meta_bf_from_data`.  

- **`meta_bf_from_data`** compares Bayes factors assuming data homogeneity. Since not all measurements are always available across datasets, the Bayes factors may not be based on the exact same data.  
- **`strict_meta_bf_from_data`** restricts comparisons to datasets where all relevant pairs of measurements are available.  

The `subject_overlap` metric in `meta_bf_from_data` quantifies how reliable each comparison is. With complete overlap, Bayes factors are directly comparable.  

---

**My BF ranking and Cohen’s d ranking don’t perfectly match. Why?**  
Bayes factors (BFs) and Cohen’s *d* are not in one-to-one correspondence. Discrepancies arise because each is combined across datasets differently. The weighting schemes emphasize different dataset features, which can lead to mismatches.  

---

**What if my standard score isn’t binary?**  
If your standard score is not binary, Cohen’s *d* is not meaningful. You will still receive all other outputs.  
