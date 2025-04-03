# --------------------------- 1. Load Libraries ---------------------------
library(readr)
library(ggplot2)
library(tidyr)
library(MASS)  # For stepwise model selection
library(car)  # For VIF calculation
library(stargazer)  # For model summary table
library(corrplot)
library(xtable)  # Load the package
library(dplyr)
library(knitr)
library(kableExtra)
library(RColorBrewer)


# --------------------------- 2. Data Import and Variable Definitions ---------------------------
# Mortgage Yield dataset from University of Florida
data_mort <- read.csv("http://users.stat.ufl.edu/~winner/data/myield.csv")
str(data_mort)
row.names(data_mort)

# --------------------------- 3. Exploratory Data Analysis (EDA) ---------------------------

# 3.1 Univariate Numerical Analysis: Summary Statistics
summary(data_mort)

# 3.2 Univariate Graphical Analysis: Histograms
par(mfrow = c(2, 3))  # 2 rows, 3 columns layout
for(var in c("X1", "X2", "X3", "X4", "X5", "X6")){
  hist(data_mort[[var]], main = paste("Histogram of", var),
       xlab = var, col = "lightblue", border = "white")
}

#  3.3 Pairwise simple correlations

# The numeric strength of association is computed as for any pair of variables with
# a correlation coecient such as Pearson’s. Since these only consider two variables
# at a time, they are called simple coecients.

pairs( ~ X1 + X2 + X3 + X4 + X5 +X6, data=data_mort)

# 3.4 Scatter Plots with Trend Lines for comparison with dependent variable 
par(mfrow = c(2, 3))
for (var in c("X1", "X2", "X3", "X4", "X5", "X6")) {
  plot(data_mort[[var]], data_mort$mortYld, 
       main = paste("Mortage Yield vs", var),
       xlab = var, ylab = "mortYld",
       col = "blue", pch = 16)  
  lines(lowess(data_mort[[var]], data_mort$mortYld), col = "red", lwd = 2)
}

#3.3.4.2 Multivariate correlation 

# = relationship independent of other variables 
# = correlation between two variables after correcting for all others
# better than normal correlation because, when we see high corr it is still 
# open question if it is because of a 3rd variable 

# how is it done?
# isolating the unique relationship between two variables by regressing 
# each on all other variables and correlating the resulting residuals

# Partial Correlation Function
p.cor <- function(x){
  inv <- solve(var(x))
  sdi <- diag(1/sqrt(diag(inv)))
  p.cor.mat <- -(sdi %*% inv %*% sdi)
  diag(p.cor.mat) <- 1
  rownames(p.cor.mat) <- colnames(p.cor.mat) <- colnames(x)
  return(p.cor.mat)
}

# Compute Partial Correlations
cat("\nPartial Correlations:\n")
p.cor(data_mort[c(2:8)])

# Compute Partial Correlation Matrix
partial_cor_matrix <- p.cor(data_mort[, c("mortYld", "X1", "X2", "X3", "X4", "X5", "X6")])
par(mfrow = c(1, 1))
# Plot the Partial Correlation Matrix
corrplot(partial_cor_matrix, method = "color", type = "upper",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         addCoef.col = "black", # Show values
         tl.col = "black", tl.srt = 45) # Rotate labels for readability

# Compute Variance Inflation Factor (VIF)
vif_values <- vif(lm(mortYld ~ X1 + X2 + X3 + X4 + X5 + X6, data = data_mort))


cat("\nVIF Values:\n")
print(vif_values)

barplot(vif_values, col = "skyblue", main = "Variance Inflation Factor (VIF)")


# Explore individual relationships as done in the paper - > how much of the variance
# is explanied by each variable 

r2_values <- c()
variables <- c("X1", "X2", "X3", "X4", "X5", "X6")
for (var in variables) {
  model <- lm(mortYld ~ data_mort[[var]], data = data_mort)
  r2_values <- c(r2_values, summary(model)$r.squared)
}
# summary table
r2_table <- data.frame(Variable = variables, R2 = r2_values)
print(r2_table)

barplot(r2_table$R2,
        names.arg = r2_table$Variable,
        col = "skyblue",
        main = "R² Values for Univariate Models",
        xlab = "Predictor Variables",
        ylab = expression(R^2))
# by here we see that X1, X2,X3 explain the most, should we start from there?

# --------------------------- 4. Model Fitting ---------------------------
# Mathematical formulation of the full linear model:
#   mortYld = β0 + β1*X1 + β2*X2 + β3*X3 + β4*X4 + β5*X5 + β6*X6 + ε


# 4.1 Full Linear Model
full_model <- lm(mortYld ~ X1 + X2 + X3 + X4 + X5 + X6, data = data_mort)

# 4.2 Stepwise Model Selection to determine best subset of predictors
#we want Akaike’s Information Criterion (AIC) to be as lower as possible so R tries iteratively
#to remove variables until there is no bigger improvement 
selected_model <- step(full_model, data = data_mort, direction = 'both')

#From our prior analysis, we can drop X4 due to its correlation to X3, a VIF>5, and low R^2, meaning it
# has low explanatory power

reduced_model <- lm(mortYld ~ X1 + X2 + X3 + X5 + X6, data = data_mort)

# 4.2 Stepwise Model Selection to determine best subset of predictors
#we want Akaike’s Information Criterion (AIC) to be as lower as possible so R tries iteratively
#to remove variables until there is no bigger improvement 
selected_reduced_model <- step(reduced_model, data = data_mort, direction = 'both')

#----------
# Compute Variance Inflation Factor (VIF)
vif_values_reduced <- vif(selected_reduced_model)
cat("\nVIF Values:\n")
print(vif_values)

barplot(vif_values_reduced, col = "skyblue", main = "Variance Inflation Factor (VIF) for reduced model")



#Plot residuals against each predictor
par(mfrow = c(2, 2))
for (var in c("X1", "X2", "X3")) {
  plot(data_mort[[var]], residuals(selected_reduced_model), 
       main = paste("Residuals vs", var),
       xlab = var, ylab = "Residuals")
  lines(lowess(data_mort[[var]], residuals(selected_reduced_model)), col = "red")
}

#Partial regression plots
library(car)
avPlots(selected_reduced_model)  # Look for nonlinear patterns in partial residuals

#Implementing polynomial terms
poly_model <- lm(mortYld ~ X1 + I(X1^2) + X2 + I(X2^2) + X3 + I(X3^2), data = data_mort)
poly_model <- lm(mortYld ~ I(X1^2) + I(X2^2) + I(X3^2), data = data_mort)

#check_model(poly_model, check = c("all"))

#comparing performance
summary(selected_reduced_model)
summary(poly_model)

#Implementing log model on outcome
log_model_outcome <- log(mortYld) ~ X1 + X2 + X3
log.model.lm <- lm(log_model_outcome, data = data_mort)
summary(log.model.lm)

#log_model <- update(formula(selected_reduced_model), log(.) ~ .)
#log.model.lm <- lm(log_model, data = data_mort)
#summary(log.model.lm)


#Log model on predictors
delta = 1e-16
log_model1 <- mortYld ~ log(X1) + log(X2+delta) + log(X3)
log.model.lm1 <- lm(log_model1, data = data_mort)
summary(log.model.lm1)

# ----------- Comparing models --------------
library(performance)

# List of models
models <- list(
  "Reduced" = selected_reduced_model,
  "Polynomial" = poly_model,
  "Log-Response" = log.model.lm,
  "Log-Predictors" = log.model.lm1
)

# Create performance comparison table
model_metrics <- compare_performance(models, metrics = c("AIC", "BIC", "R2", "Adj_R2", "RMSE"))
print(model_metrics)


# Melt metrics for plotting
plot_data <- model_metrics %>%
  pivot_longer(cols = c(AIC, BIC, RMSE), names_to = "Metric", values_to = "Value")

# Plot AIC values separately
p1 <- ggplot(plot_data[plot_data$Metric == "AIC",], aes(x = Name, y = Value, fill = Name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Value, 2)), vjust = -0.5) +
  labs(title = "AIC Values for All Models", y = "AIC Value", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Spectral")

# Plot BIC values separately
p2 <- ggplot(plot_data[plot_data$Metric == "BIC",], aes(x = Name, y = Value, fill = Name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Value, 2)), vjust = -0.5) +
  labs(title = "BIC Values for All Models", y = "BIC Value", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Spectral")

# Plot AIC values separately
p3 <- ggplot(plot_data[plot_data$Metric == "RMSE",], aes(x = Name, y = Value, fill = Name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Value, 2)), vjust = -0.5) +
  labs(title = "RMSE Values for All Models", y = "RMSE", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Spectral")


# Calculate adjusted R² values
adj_r2 <- c(
  summary(selected_reduced_model)$adj.r.squared,
  summary(poly_model)$adj.r.squared,
  summary(log.model.lm)$adj.r.squared,
  summary(log.model.lm1)$adj.r.squared
)




# Create data frame for adjusted R² values
adj_r2_df <- data.frame(
  Model = c("Reduced", "Polynomial", "Log-Response", "Log-Predictors"),
  Adj_R2 = adj_r2
)

# Merge AIC and adjusted R² results
#comparison_df <- merge(aic_df, adj_r2_df, by = "Model", all = TRUE)

# Plot Adjusted R^2
p4 <- ggplot(adj_r2_df, aes(x = Model, y = Adj_R2, fill = Model)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Adj_R2, 2)), vjust = -0.5) +
  labs(title = expression(paste("Adjusted ",R^2," for All Models")), y = expression(paste("Adjusted",R^2)), x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette = "Spectral") 

library(patchwork)
p1 
p2
p3
p4

par(mfrow = c(2, 2))
plot(selected_reduced_model, which = 1, main = "Reduced Model")
plot(poly_model, which = 1, main = "Polynomial Model")
plot(log.model.lm, which = 1, main = "Log-Response Model")
plot(log.model.lm1, which = 1, main = "Log-Predictors Model")

# Generate predictions
data_mort$pred_reduced <- predict(selected_reduced_model)
data_mort$pred_poly <- predict(poly_model)
data_mort$pred_logresp <- exp(predict(log.model.lm))  # Reverse log-transform for Y
data_mort$pred_logpred <- predict(log.model.lm1)

data_long <- data_mort %>%
  pivot_longer(
    cols = starts_with("pred_"),
    names_to = "Model",
    values_to = "Predicted"
  ) %>%
  mutate(
    Model = case_when(
      Model == "pred_reduced" ~ "Reduced",
      Model == "pred_poly" ~ "Polynomial",
      Model == "pred_logresp" ~ "Log-Response",
      Model == "pred_logpred" ~ "Log-Predictors"
    )
  )
ggplot(data_long, aes(x = mortYld, y = Predicted, color = Model)) +
  geom_point(alpha = 1) +  # Add points
  #geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +  # Linear trend lines
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Actual vs. Predicted Mortgage Yields with Trend Lines",
    x = "Actual",
    y = "Predicted"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Spectral")
 
#-------------------------------


check_model(selected_reduced_model, check="all") #stepwise reduced model
check_model(poly_model, check = "all")  # Polynomial model
check_model(log.model.lm, check = "all") #Log-response model
check_model(log.model.lm1, check = "all")  # Log-predictors model
check_model(log.model.lm1, check = "qq")  # Log-predictors model
