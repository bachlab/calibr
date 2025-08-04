# CalibR: A Package for Experimental Calibration using the Meta-Analysis of Bayes Factors

📌 **Overview**  
CalibR is an R package that provides tools for experimental calibration, enabling researchers to estimate and update measurement method rankings across various datasets. 
This is achieved via combining **Bayes Factors** for the object model, induced by a certain measurement, against the null model, which is simply an intercept (the worst possible measurement method). 
The package is particularly useful for combining evidence from multiple studies.

---

## 🔧 **Installation**
To install CalibR from GitLab, use the following command in R:

```r
# Install necessary dependencies
install.packages("devtools")  

# Install from GitLab
devtools::install_git("https://gitlab.com/bachlab/methods/calibr.git")
```

## 📊  Use case: computing a calibration summary
In this example, three datasets from three different studies are analysed. The standard score is binary (0/1); we adopt a between-subject setting, where
different standard scores correspond to different subject groups. Per study, standard score value, and  measure, we have a sample of 20 values. 
These numbers can be changed in the code. The measure 'Measure_A' will be better than 'Measure_B', because Measure_A (Measure_B) values come from a 
gaussian of mean 45 (40) when standard score is 0, and 70 (55) when standard score is 1. Standard deviations are fixed at 30.

# 1. Data preparation

```r
# Load library dplyr
library(dplyr)

# Set seed for reproducibility
set.seed(123)

# Define parameters
n_studies <- 3  # Number of studies
n_standard_scores <- 2  # Binary standard score: 0 or 1
n_per_group <- 20  # Number of samples per study-standard_score-measurement tuple

# Define dataset, measurement, and standard score combinations
datasets <- paste0("Study", 1:n_studies)
measurements <- c("Measure_A", "Measure_B")
standard_scores <- c(0, 1)  # Binary standard scores

# Generate example dataset
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
    measurement_value = if (.$standard_score == 0) { 
      if (.$measurement_name == "Measure_A") {
        rnorm(n_per_group, mean = 45, sd = 30)  # Lower for standard_score = 0
      } else {
        rnorm(n_per_group, mean = 40, sd = 30)
      }
    } else {  # standard_score == 1
      if (.$measurement_name == "Measure_A") {
        rnorm(n_per_group, mean = 70, sd = 30)  # Higher for standard_score = 1
      } else {
        rnorm(n_per_group, mean = 55, sd = 30)
      }
    }
  )) %>%
  ungroup()
```
# 2. Compute a Calibration summary

```r
# Load CalibR
library(CalibR)

# Compute the calibration summary
result <- calibration_summary(data = df, 
                              dataset_col = "dataset", 
                              measurement_name_col = "measurement_name",
                              measurement_values_col = "measurement_value", 
                              standard_scores_col = "standard_score"
                              )

# Print summary table
print(result$summary_table)

# Print Bayes Factor ratios
print(result$log_bayes_factor_diffs)

# Print best measurement statement
print(result$best_measurement)

print(result$subject_overlap)
```

