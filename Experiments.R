# - - - - - - - - - - - - - - - - - - - - - Overall vs. BF-obtained retrod.

# Load necessary library
library(dplyr)
library(stringr)

# Set seed for reproducibility
set.seed(123)

# Function to generate a dataset with correlated binary and real random vars
generate_dataset_binary <- function(dataset_name, n, rho, set_id) {

  # Calculate the threshold c
  c <- qnorm(1 - p)

  # Calculate the required correlation r
  phi_c <- dnorm(c)
  r <- rho * sqrt(p * (1 - p)) / phi_c

  # Create the covariance matrix
  Sigma <- matrix(c(1, r, r, 1), nrow = 2)

  # Generate bivariate normal data
  library(MASS)
  data <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  X <- data[, 1]
  Y <- data[, 2]

  # Generate the binary variable B
  B <- as.numeric(Y > c)

  # Return as dataframe
  data.frame(
    set_id = set_id,  # Identify which batch this dataset belongs to
    dataset = dataset_name,
    standardscore = as.vector(B),
    measurement = as.vector(X)
  )
}

# Function to generate a dataset with correlated real random vars
generate_dataset_continuous <- function(dataset_name, n, rho, set_id) {
  # Generate X (continuous)
  X <- scale(rnorm(n))  # Standardized for easier correlation control

  # Generate Y (continuous) based on X with specified correlation
  Y <- rho * X + sqrt(1 - rho^2) * rnorm(n)  # Continuous outcome

  # Return as dataframe
  data.frame(
    set_id = set_id,
    dataset = dataset_name,
    standardscore = as.vector(Y),
    measurement = scale(X)  # Continuous outcome
  )
}

t_to_correlation <- function(t, n) {
  if (n <= 2) {
    stop("Sample size must be greater than 2.")
  }
  r <- t / sqrt(n - 2 + t^2)
  return(r)
}

# Initialize an empty dataframe
all_datasets <- data.frame()

# Loop to generate 100 sets, each containing 3 datasets with varying correlations
for (i in 1:100) {

  df_measurement_u <- bind_rows(
    generate_dataset_binary(paste0("Dataset_A_", i), 100, runif(1, 0.4, 0.6), i),  # Weak correlation
    generate_dataset_binary(paste0("Dataset_B_", i), 100, runif(1, 0.4, 0.6), i),  # Moderate correlation
  )

  df_measurement_u$measurement_name = "Measurement_U"

  df_measurement_v <- bind_rows(
    generate_dataset_binary(paste0("Dataset_A_", i), 100, runif(1, 0.3, 0.5), i),  # Weak correlation
    generate_dataset_binary(paste0("Dataset_B_", i), 100, runif(1, 0.3, 0.5), i),  # Moderate correlation
  )

  df_measurement_v$measurement_name = "Measurement_V"

  # Append to the master dataframe
  all_datasets <- bind_rows(all_datasets, df_measurement_u, df_measurement_v)
}

# Print first few rows
head(all_datasets)

# Initialize storage lists
meta_correlations <- list()  # Stores correlations derived from Bayes Factor
overall_correlations <- list()  # Stores overall X-Y correlation per set_id

# Unique set IDs
set_ids <- unique(all_datasets$set_id)

# Loop over each set_id
for (i in set_ids) {
  # Subset data for the current set_id
  subset_data <- subset(all_datasets, set_id == i & measurement_name == "Measurement_V")

  # (1) Compute Meta-analytic Bayes Factor
  mbf <- meta_bf_from_data(subset_data, 'dataset', 'measurement_name', 'measurement', 'standardscore')

  # (2) Convert Bayes Factor to t-statistic
  mt <- bf2tstat(mbf$combined_bayes_factor, nrow(subset_data))

  # (3) Convert t-statistic to correlation
  meta_correlation <- t_to_correlation(mt, nrow(subset_data))

  # Store in list
  meta_correlations[[as.character(i)]] <- meta_correlation

  # Compute overall correlation between X (standardscore) and Y (measurement)
  overall_correlation <- cor(subset_data$standardscore, subset_data$measurement)

  # Store in list
  overall_correlations[[as.character(i)]] <- overall_correlation
}

# Convert lists to dataframes for easy inspection
meta_correlation_df <- data.frame(set_id = as.numeric(names(meta_correlations)),
                                  meta_correlation = unlist(meta_correlations))

overall_correlation_df <- data.frame(set_id = as.numeric(names(overall_correlations)),
                                     overall_correlation = unlist(overall_correlations))

# Merge both results into a single dataframe
final_results <- merge(meta_correlation_df, overall_correlation_df, by = "set_id")

# Print first few results
head(final_results)

plot(final_results$meta_correlation, final_results$overall_correlation)
hist(final_results$meta_correlation - final_results$overall_correlation)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - Sourav.

library(dplyr)
library(tidyr)
library(devtools)

sourav_output <- read.csv("~/Downloads/share_federico.csv")
sourav_output <- sourav_output %>% mutate(CS_at_USm_diff = CSpUSm - CSmUSm)
sourav_output_filtered <- sourav_output %>% filter(message == "Method ran successfully")

sourav_output_filtered <- tidyr::separate(sourav_output_filtered, dataset_name, into = c("dataset", "subject_ID"), sep = ", ", remove = FALSE)

sourav_output_filtered <- sourav_output_filtered %>% pivot_longer(
  cols = ends_with("USm"),
  names_to = "standard_score",
  values_to = "measurement_value"
)

sourav_output_filtered$standard_score <- ifelse(sourav_output_filtered$standard_score == "CSpUSm", 1, 0)

load_all()

cso <- calibration_summary(sourav_output_filtered,
                           dataset_col = "dataset",
                           measurement_name_col = "method_name",
                           measurement_values_col = "measurement_value",
                           standard_scores_col = "standard_score",
                           subject_id_col = "subject_ID"
)

mco <- meta_cohensd(sourav_output_filtered,
                    subject_id_col = "subject_ID",
                    dataset_col = "dataset",
                    measurement_name_col = "method_name",
                    measurement_values_col = "measurement_value",
                    standard_scores_col = "standard_score"
)

library(ggplot2)
library(dplyr)

# Extract method number and sort by decreasing Log(BF)
cso$summary_table <- cso$summary_table %>%
  mutate(
    method_number = as.numeric(gsub("method_", "", measurement)),  # Extract method number
    log_sample_size = log(sample_size)  # Compute log of sample size for dot size
  ) %>%
  arrange(desc(retrodiction_score))  # Sort by decreasing Log(BF)

# Create custom labels: Hide every other label
method_labels <- cso$summary_table$method_number

# Create ggplot
ggplot(cso$summary_table, aes(x = factor(method_number, levels = method_number),
                              y = retrodiction_score,
                              size = log_sample_size)) +  # Dot size by log(sample_size)
  geom_point(color = "blue", alpha = 0.4) +  # Points with transparency
  scale_x_discrete(labels = method_labels) +  # Apply the modified labels
  labs(
    x = "Method Number",
    y = "Meta retrodiction_score",
    size = "Log Sample Size"
  ) +
  theme_minimal(base_size = 14) +  # Elegant minimal theme
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),  # Tilted for readability
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"  # Keep legend readable
  )

# - - - - - - - - - - - - - - - - - - - - - - - - - - Generic example.

generate_study <- function(study_id, n_per_group, within = TRUE) {
  df <- expand.grid(
    study_id = study_id,
    subject_id = 1:n_per_group,
    standard_score = c(0, 1),
    measurement_type = c("A", "B")
  )
  if(!within){
    df$subject_id <- 1:nrow(df)
  }
  df$measurement <- with(df, rnorm(nrow(df),
                                   mean = ifelse(measurement_type == "A",
                                                 ifelse(standard_score == 1, 10, 5),
                                                 ifelse(standard_score == 1, 7, 6)), sd = 5))
  return(df)
}

all_studies <- do.call(rbind, list(
  generate_study("Study1", 30, within = TRUE),
  generate_study("Study2", 30, within = TRUE),
  generate_study("Study3", 30, within = FALSE)
))

result <- calibration_summary(
  data = all_studies, dataset_col = "study_id", standard_scores_col = "standard_score",
  measurement_name_col = "measurement_type", measurement_values_col = "measurement",
  subject_id_col = "subject_id"
)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - Olivier.

library(dplyr)
library(tidyr)
library(devtools)
load_all()

olivier_data <- read.csv("~/Desktop/SCRopt_data (1).csv")

cso_olivier <- calibration_summary(data = olivier_data,
                                   dataset_col = "database",
                                   measurement_name_col = "method",
                                   measurement_values_col = "FC_norm",
                                   standard_scores_col = "standard_score",
                                   subject_id_col = "ppid")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - CombiMeasures.

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(devtools)
load_all()

# import data.
data_orig = read.csv("../Calibration/estimates_multimodal_0921.csv");

# eliminate SPS measurements.
data_orig$SPS = c();

data_orig = subset(data_orig,
                   (HPR > min(HPR, na.rm = T) | is.na(HPR)) &
                     (SCR > min(SCR, na.rm = T) & SCR < max(SCR, na.rm = T) | is.na(SCR)))

data_orig.s =
  data_orig %>%
  group_by(Dataset) %>%
  mutate(
    SCR = SCR/sd(SCR,na.rm = T),
    PSR = PSR/sd(PSR,na.rm = T),
    HPR = HPR/sd(HPR,na.rm = T),
    RAR = RAR/sd(RAR,na.rm = T)
  )

data_orig.s[data_orig.s$Dataset == "FER01",]$RAR = NaN;

# Step 1: Assign subject_id as Dataset + row number within Dataset
data_orig.s <- data_orig.s %>%
  group_by(Dataset) %>%
  mutate(subject_id = paste0(Dataset, "_", row_number())) %>%
  ungroup()

# Tidy up, eliminate the peoples with more NAs than the minimum for each dataset.
data_clean <- data_orig.s %>%
  mutate(
    na_count = rowSums(across(c(PSR, SCR, RAR, HPR), is.na))
  ) %>%
  group_by(Dataset) %>%
  mutate(
    min_na_count = min(na_count)
  ) %>%
  filter(na_count == min_na_count) %>%
  ungroup()

# Step 2: Pivot longer and add both standard scores
estimates_long <- data_clean %>%
  pivot_longer(cols = c(PSR, SCR, RAR, HPR), names_to = "measurement", values_to = "value") %>%
  mutate(standard_score = 1) %>%
  bind_rows(
    data_clean %>%
      pivot_longer(cols = c(PSR, SCR, RAR, HPR), names_to = "measurement", values_to = "value") %>%
      mutate(value = 0, standard_score = 0)
  ) %>%
  arrange(subject_id, measurement, standard_score)

# Step 3: Tidy columns
estimates_long <- estimates_long %>%
  select(subject_id, Dataset, standard_score, measurement, value, everything())

estimates_long <- estimates_long %>%
  group_by(subject_id, Dataset, measurement) %>%
  filter(!any(is.na(value))) %>%
  ungroup()

# Step 1

step_1_result <- calibration_summary(
  data = estimates_long %>% filter(complete.cases(.)) %>% filter(Dataset %in% c("SC4B")),
  dataset_col = "Dataset",
  measurement_name_col = 'measurement',
  measurement_values_col = 'value',
  standard_scores_col = 'standard_score',
  subject_id_col = 'subject_id',
  power = 0.8,
  sig.level = 0.05,
  alternative = "two.sided",
  relevant_difference = 1
)

# Step 2

step_2_result <- calibration_summary(
  data = estimates_long %>% filter(complete.cases(.)) %>% filter(Dataset %in% c("PubFe", "SC4B")) ,
  dataset_col = "Dataset",
  measurement_name_col = 'measurement',
  measurement_values_col = 'value',
  standard_scores_col = 'standard_score',
  subject_id_col = 'subject_id',
  power = 0.8,
  sig.level = 0.05,
  alternative = "two.sided"
)

# Step 3

step_3_result <- calibration_summary(
  data = estimates_long %>% filter(complete.cases(.)) %>% filter(Dataset %in% c("VC7B", "SC4B", "PubFe")) ,
  dataset_col = "Dataset",
  measurement_name_col = 'measurement',
  measurement_values_col = 'value',
  standard_scores_col = 'standard_score',
  subject_id_col = 'subject_id',
  power = 0.8,
  sig.level = 0.05,
  alternative = "two.sided"
)

# Step 4

step_4_result  <- calibration_summary(
  data = estimates_long %>% filter(complete.cases(.)) %>% filter(Dataset %in% c("VC7B", "SC4B", "PubFe", "DoxMemP")) ,
  dataset_col = "Dataset",
  measurement_name_col = 'measurement',
  measurement_values_col = 'value',
  standard_scores_col = 'standard_score',
  subject_id_col = 'subject_id',
  power = 0.8,
  sig.level = 0.05,
  alternative = "two.sided"
)

# Step 5

step_5_result  <- calibration_summary(
  data = estimates_long %>% filter(complete.cases(.)) %>% filter(Dataset %in% c("VC7B", "SC4B", "PubFe", "DoxMemP", "FR")) ,
  dataset_col = "Dataset",
  measurement_name_col = 'measurement',
  measurement_values_col = 'value',
  standard_scores_col = 'standard_score',
  subject_id_col = 'subject_id',
  power = 0.8,
  sig.level = 0.05,
  alternative = "two.sided"
)

# Step 6

step_6_result  <- calibration_summary(
  data = estimates_long %>% filter(complete.cases(.)) %>% filter(Dataset %in% c("VC7B", "SC4B", "PubFe", "DoxMemP", "FR", "TC")) ,
  dataset_col = "Dataset",
  measurement_name_col = 'measurement',
  measurement_values_col = 'value',
  standard_scores_col = 'standard_score',
  subject_id_col = 'subject_id',
  power = 0.8,
  sig.level = 0.05,
  alternative = "two.sided"
)

# Step 7

step_7_result  <- calibration_summary(
  data = estimates_long %>% filter(complete.cases(.)) %>% filter(Dataset %in% c("VC7B", "SC4B", "PubFe", "DoxMemP", "FR", "TC", "FSS6B")) ,
  dataset_col = "Dataset",
  measurement_name_col = 'measurement',
  measurement_values_col = 'value',
  standard_scores_col = 'standard_score',
  subject_id_col = 'subject_id',
  power = 0.8,
  sig.level = 0.05,
  alternative = "two.sided"
)

# Step 8

step_8_result  <- calibration_summary(
  data = estimates_long %>% filter(complete.cases(.)) %>% filter(Dataset %in% c("VC7B", "SC4B", "PubFe", "DoxMemP", "FR", "TC", "FSS6B", "FER02")) ,
  dataset_col = "Dataset",
  measurement_name_col = 'measurement',
  measurement_values_col = 'value',
  standard_scores_col = 'standard_score',
  subject_id_col = 'subject_id',
  power = 0.8,
  sig.level = 0.05,
  alternative = "two.sided"
)

# Step 9

step_9_result  <- calibration_summary(
  data = estimates_long %>% filter(complete.cases(.)) %>% filter(Dataset %in% c("VC7B", "SC4B", "PubFe", "DoxMemP", "FR", "TC", "FSS6B", "FER02","FER01")) ,
  dataset_col = "Dataset",
  measurement_name_col = 'measurement',
  measurement_values_col = 'value',
  standard_scores_col = 'standard_score',
  subject_id_col = 'subject_id',
  power = 0.8,
  sig.level = 0.05,
  alternative = "two.sided",
  relevant_difference = 1
)


# - - - - - - - - - - - - The lone logBF values - - - - - - - - - - - - - - - -

# Combine summary tables from each step, with a "step" indicator
# List all step results in order
step_results_list <- mget(paste0("step_", 1:9, "_result"))

# Extract, add step, and bind together
step_results_combined <- bind_rows(
  Map(
    function(res, n) {
      res$summary_table %>%
        select(measurement, log_combined_bayes_factor, retrodiction_score, meta_hedge_g, sample_size) %>%
        mutate(step = paste("Step", n))
    },
    step_results_list, 1:9
  )
) %>%
  select(step, measurement, log_combined_bayes_factor, retrodiction_score, meta_hedge_g, sample_size)

step_results_combined <- step_results_combined %>%
  mutate(
    step_num = as.integer(stringr::str_extract(step, "\\d+")),
    step = case_when(
      step_num == 1 ~ "1",
      step_num == 2 ~ "1,2",
      step_num > 2 ~ paste0("1..", step_num)
    )
  ) %>%
  select(-step_num)

step_results_combined$retrodiction_score <- abs(step_results_combined$retrodiction_score)
step_results_combined$meta_hedge_g <- abs(step_results_combined$meta_hedge_g)


# - - - - - - - - - - - - All the plots - - - - - - - - - - - - - - - -

plot_abs_lbf <- step_results_combined %>%
  ggplot(aes(x = step,
             y = log_combined_bayes_factor,
             group = measurement,
             color = measurement)) +
  geom_line(aes(alpha = sqrt(sample_size), linewidth = sqrt(sample_size)), show.legend = TRUE) +
  # geom_point(aes(alpha = sqrt(sample_size), size = sqrt(sample_size)),
  #            shape = 21, fill = "white", stroke = 1.2, show.legend = TRUE) +
  scale_alpha(range = c(0.5, 1), guide = "none") +
  scale_linewidth(range = c(0.8, 2.5), guide = "none") +
  scale_size(range = c(2, 6), guide = "none") +
  scale_color_brewer(palette = "Set1", name = "Measurement") +
  labs(
    title = "Log(BF) ",
    subtitle = "Line alpha/width reflect overall sample size",
    x = "Datasets integrated",
    y = "Log(BF)",
    color = "Measurement"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

plot_abs_retr <- step_results_combined %>%
  ggplot(aes(x = step,
             y = retrodiction_score,
             group = measurement,
             color = measurement)) +
  geom_line(aes(alpha = sqrt(sample_size), linewidth = sqrt(sample_size)), show.legend = TRUE) +
  # geom_point(aes(alpha = sqrt(sample_size), size = sqrt(sample_size)),
  #            shape = 21, fill = "white", stroke = 1.2, show.legend = TRUE) +
  scale_alpha(range = c(0.5, 1), guide = "none") +
  scale_linewidth(range = c(0.8, 2.5), guide = "none") +
  scale_size(range = c(2, 6), guide = "none") +
  scale_color_brewer(palette = "Set1", name = "Measurement") +
  labs(
    title = "Retrodiction Score ",
    subtitle = "Line alpha/width reflects overall sample size",
    x = "Datasets integrated",
    y = "Pearson r",
    color = "Measurement"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

plot_abs_hedg <- step_results_combined %>%
  ggplot(aes(x = step,
             y = meta_hedge_g,
             group = measurement,
             color = measurement)) +
  geom_line(aes(alpha = sqrt(sample_size), linewidth = sqrt(sample_size)), show.legend = TRUE) +
  # geom_point(aes(alpha = sqrt(sample_size), size = sqrt(sample_size)),
  #            shape = 21, fill = "white", stroke = 1.2, show.legend = TRUE) +
  scale_alpha(range = c(0.5, 1), guide = "none") +
  scale_linewidth(range = c(0.8, 2.5), guide = "none") +
  scale_size(range = c(2, 6), guide = "none") +
  scale_color_brewer(palette = "Set1", name = "Measurement") +
  labs(
    title = "Meta Hedge's g ",
    subtitle = "Line alpha/width reflects overall sample size",
    x = "Datasets integrated",
    y = "Hedge's g",
    color = "Measurement"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )


# - - - - - - - - - - - - The differences - - - - - - - - - - - - - - - -

extract_pairwise_df <- function(step_result, step_name) {
  diffs <- as.data.frame(step_result$log_bayes_factor_diffs)
  overlap <- as.data.frame(step_result$subject_overlap)
  summary_table <- as.data.frame(step_result$summary_table)

  if (ncol(diffs) <= 1) return(data.frame())

  diffs$measurement1 <- rownames(diffs)
  overlap$measurement1 <- rownames(overlap)

  # Pivot longer for pairwise matrices
  diffs_long <- diffs %>%
    tidyr::pivot_longer(-measurement1, names_to = "measurement2", values_to = "log_bayes_factor_diff")

  overlap_long <- overlap %>%
    tidyr::pivot_longer(-measurement1, names_to = "measurement2", values_to = "subject_overlap")

  # Get sample sizes per measurement from summary_table
  sample_sizes <- summary_table %>%
    dplyr::select(measurement, sample_size = sample_size)

  # Join, create pair, and merge in sample sizes
  pairwise_df <- diffs_long %>%
    dplyr::left_join(overlap_long, by = c("measurement1", "measurement2")) %>%
    dplyr::filter(measurement1 != measurement2) %>%
    # Sort names for unique unordered pairs
    dplyr::mutate(
      measurement_low = pmin(as.character(measurement1), as.character(measurement2)),
      measurement_high = pmax(as.character(measurement1), as.character(measurement2)),
      measurements_involved = paste(measurement_low, measurement_high, sep = "-"),
      step = step_name
    ) %>%
    # Add sample sizes for each measurement in the pair
    dplyr::left_join(sample_sizes, by = c("measurement1" = "measurement")) %>%
    dplyr::rename(sample_size_1 = sample_size) %>%
    dplyr::left_join(sample_sizes, by = c("measurement2" = "measurement")) %>%
    dplyr::rename(sample_size_2 = sample_size) %>%
    dplyr::select(step, measurements_involved, log_bayes_factor_diff, subject_overlap, sample_size_1, sample_size_2) %>%
    dplyr::distinct(measurements_involved, .keep_all = TRUE) # Only one per pair

  return(pairwise_df)
}

all_pairs <- list()
for (i in 1:9) {
  step_name <- paste0("Step ", i)
  res <- get(paste0("step_", i, "_result"))
  out <- tryCatch(extract_pairwise_df(res, step_name), error = function(e) data.frame())
  all_pairs[[i]] <- out
}

all_pairs_df <- bind_rows(all_pairs)

all_pairs_df <- all_pairs_df %>%
  mutate(
    step_num = as.integer(stringr::str_extract(step, "\\d+")),
    step = case_when(
      step_num == 1 ~ "1",
      step_num == 2 ~ "1,2",
      step_num > 2 ~ paste0("1..", step_num)
    )
  ) %>%
  select(-step_num)

plot_diff <- all_pairs_df %>%
  filter(str_starts(measurements_involved, "HPR")) %>%
  ggplot(aes(x = step,
             y = log_bayes_factor_diff,
             group = measurements_involved,
             color = measurements_involved
             )) +
  geom_line(aes(alpha = subject_overlap, linewidth = subject_overlap), show.legend = TRUE) +
  # geom_point(aes(alpha = sqrt(subject_overlap), size = sqrt(subject_overlap)),
  #            shape = 21, fill = "white", stroke = 1.2, show.legend = TRUE) +
  scale_alpha(range = c(0.5, 1), guide = "none") +
  scale_linewidth(range = c(0.8, 2.5), guide = "none") +
  scale_size(range = c(2, 6), guide = "none") +
  scale_color_grey(start = 0.1, end = 0.9, name = "Comparison Pair") +
  labs(
    title = "Log(BF) Difference (HPR pairs)",
    subtitle = "Line width reflects subject overlap",
    x = "Datasets integrated",
    y = "Log(BF) Difference",
    color = "Comparison Pair"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

# - - - - - - - - - - - - Combined plots - - - - - - - - - - - - - - - -

library(patchwork)

da_plot <- plot_abs_lbf + plot_abs_retr + plot_diff

ggsave("../_Temporary/ensemble_img_calibration.png", plot = da_plot, width = 9, height = 2.5)


# - - - - - - - - - - - - Checks - - - - - - - - - - - - - - - - - - - -

