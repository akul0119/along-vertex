# Along-Vertex Analysis for DTI Data: TBI vs PTE Classification

################################################################################
# Libraries and Setup
################################################################################
library(dplyr)
library(ggplot2)
library(glmnet)
library(pROC)
library(caret)

# Set paths
base_path <- "~/Desktop/bundles_along_vertex/datasets"
roc_plot_path <- "~/Desktop/bundles_along_vertex/roc"
model_save_path <- "~/Desktop/bundles_along_vertex/models"
plot_save_path <- "~/Desktop/bundles_along_vertex/plots"

demo_file <- "~/Desktop/tract_bundles_p3/patient_demo.csv"

# Load demo data
demo <- read.csv(demo_file, header = TRUE)

# File paths for DTI data
file_paths <- c(
  "bundles.along.vertex.dti.map_FA_mean.csv",
  "bundles.along.vertex.dti.map_MD_mean.csv",
  "bundles.along.vertex.dti.map_RD_mean.csv"
)

################################################################################
# Functions
################################################################################

# Function to preprocess data
preprocess_data <- function(data_path, demo, outcome_data) {
  data <- read.csv(data_path, header = TRUE)
  merged_data <- merge(outcome_data, data, by = "subject")
  merged_data <- merge(demo, merged_data, by = "subject")

  processed_data <- merged_data %>%
    filter(late == 0 | late == 1, subject != 22) %>%
    mutate(late = as.numeric(late))

  # Remove NA values and impute missing values
  processed_data <- processed_data %>%
    select(where(~ mean(!is.na(.)) > 0.9)) %>%
    mutate_if(is.numeric, ~ifelse(is.na(.), mean(., na.rm = TRUE), .))

  return(processed_data)
}

# Function to perform linear modeling
perform_linear_modeling <- function(data, variables, model_save_path, file_name) {
  significant_results <- data.frame(
    Model = character(),
    Variable = character(),
    Estimate = numeric(),
    Std.Error = numeric(),
    t_value = numeric(),
    p_value = numeric(),
    Overall_p_value = numeric(),
    stringsAsFactors = FALSE
  )

  for (variable in variables) {
    formula <- as.formula(paste(variable, "~ late + age + sex"))
    model <- lm(formula, data = data)
    summary_model <- summary(model)
    overall_p_value <- pf(summary_model$fstatistic[1], summary_model$fstatistic[2], summary_model$fstatistic[3], lower.tail = FALSE)

    if (summary_model$coefficients["late", "Pr(>|t|)"] < 0.05) {
      significant_results <- rbind(significant_results, data.frame(
        Model = variable,
        Variable = "late",
        Estimate = summary_model$coefficients["late", "Estimate"],
        Std.Error = summary_model$coefficients["late", "Std. Error"],
        t_value = summary_model$coefficients["late", "t value"],
        p_value = summary_model$coefficients["late", "Pr(>|t|)"],
        Overall_p_value = overall_p_value
      ))
    }
  }

  save_file <- file.path(model_save_path, paste0("significant_results_", gsub(".csv", "", file_name), ".csv"))
  write.csv(significant_results, save_file, row.names = FALSE)

  return(significant_results)
}

# Function to plot significant variables
plot_significant_variables <- function(data, significant_results, save_path, file_name) {
  data$late <- as.factor(data$late)

  for (i in seq_len(nrow(significant_results))) {
    variable <- significant_results$Model[i]
    p_value <- significant_results$p_value[i]

    plot <- ggplot(data, aes(x = late, y = .data[[variable]])) +
      geom_boxplot() +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(title = paste("Comparison of", variable, "by late"),
           subtitle = paste("P-value:", format(p_value, scientific = TRUE))) +
      theme_minimal()

    plot_file <- file.path(save_path, paste0("plot_", gsub(".csv", "", file_name), "_", variable, ".png"))
    ggsave(plot_file, plot, width = 8, height = 6)
  }
}

# Function to perform Elastic Net Regression
perform_enr <- function(data, model_save_path, roc_plot_path, file_name) {
  metrics <- names(data)[!(names(data) %in% c("subject", "site", "age", "sex", "pid", "acute", "late"))]
  data$late <- factor(data$late, levels = c(0, 1), labels = c("TBI", "PTS"))

  weights <- ifelse(data$late == "PTS", 1 / table(data$late)["PTS"], 1 / table(data$late)["TBI"])

  X <- as.matrix(data[, metrics])
  X <- scale(X, center = TRUE, scale = apply(X, 2, function(x) ifelse(sd(x) > 0, sd(x), 1)))
  y <- as.numeric(data$late) - 1

  cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = 0.5, grouped = FALSE, weights = weights)

  best_lambda <- cv_fit$lambda.min
  best_coefficients <- coef(cv_fit, s = "lambda.min")

  probabilities <- predict(cv_fit, newx = X, s = "lambda.min", type = "response")
  predictions <- ifelse(probabilities > 0.5, 1, 0)

  roc_obj <- roc(y, probabilities)
  auc_value <- auc(roc_obj)

  save_file <- file.path(model_save_path, paste0("late_enr_", gsub(".csv", "", file_name), ".txt"))
  sink(save_file)
  cat("Elastic Net Regression Results\n")
  cat(sprintf("Best Lambda: %f\n", best_lambda))
  print(roc_obj)
  sink()

  roc_plot_file <- file.path(roc_plot_path, paste0("late_ROC_", gsub(".csv", "", file_name), ".png"))
  png(roc_plot_file)
  plot(roc_obj, main = paste("ROC Curve -", gsub(".csv", "", file_name)))
  dev.off()
}

################################################################################
# Main Script
################################################################################
outcome_data <- read.csv("~/Desktop/outcome_data.csv")
roc_list <- list()

for (file_name in file_paths) {
  data_path <- file.path(base_path, file_name)
  processed_data <- preprocess_data(data_path, demo, outcome_data)

  Vars <- names(processed_data)[!(names(processed_data) %in% c("subject", "acute", "late", "site", "age", "sex", "pid"))]

  significant_results <- perform_linear_modeling(processed_data, Vars, model_save_path, file_name)
  plot_significant_variables(processed_data, significant_results, plot_save_path, file_name)

  perform_enr(processed_data, model_save_path, roc_plot_path, file_name)
}
