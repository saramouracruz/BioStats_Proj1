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
# is explained by each variable 

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
selected_model <- step(full_model, data = data_mort)

#as mentioned in the paper we can see that there is a more or less log relationship between
#independent variable and Xs

# 4.3 Log-Transformed Model (addressing potential non-linear relationships)
delta = 1e-8
log_model <- lm(mortYld ~ log(X1+delta) + log(X2+delta) + X3 + X4 + log(X5+delta) + log(X6+delta), data = data_mort)

# 4.4 Stepwise Model Selection to determine best subset of predictors !!!looks like the best metrics? cooks distance is fucked 
selected_model_log <- step(log_model, data = data_mort)
summary(selected_model_log)
plot(selected_model_log, which = 1:4)

#. 4.5 Try interaction model


# --------------------------- 5. Model Assessment ---------------------------
# Model Assumptions:
# 1. Errors have mean 0.
# 2. Errors are homoscedastic (constant variance).
# 3. Errors are uncorrelated.
# 4. Errors are normally distributed.
#
# Diagnostic Plots for the Selected Model:

par(mfrow = c(2, 3))
#plot(selected_model, which = 1:6, main = "Diagnostics: Selected Model")
plot(selected_model, which = 1:6)



# QQ Plot of Residuals
tmp <- qqnorm(residuals(selected_model), pch=20, main = "QQ Plot of Residuals")
qqline(residuals(selected_model), col = "red")
# Compute the differences between the QQ plot points and the line
diff <- (tmp$x - tmp$y)
### label the residuals that deviate too far from the line and by that try to understand 
### why are they problematic 
text(tmp$x, tmp$y, ifelse((abs(diff) > 0.9), names(diff), ""), pos=2)
rm(tmp,diff)

# --------------------------------------------------
# Observed vs. Fitted Plot for the Selected Mortgage Yield Model
# --------------------------------------------------
plot(data_mort$mortYl ~ fitted(selected_model), 
     pch = 20,
     # Optionally adjust x and y limits based on your data range
     xlim = c(min(fitted(selected_model)), max(fitted(selected_model))),
     ylim = c(min(data_mort$mortYld), max(data_mort$mortYld)),
     xlab = "Fitted Mortgage Yield",
     ylab = "Observed Mortgage Yield",
     main = "Observed vs. Fitted Mortgage Yield")
abline(0, 1, col = "red", lwd = 2)
grid(col = "black")
par(mfrow = c(1,1))


#plot(selected_model$fitted.values, selected_model$residuals,
     ##xlab = "Fitted Values", ylab = "Residuals",
     #pch = 16, col = "darkgreen")
#abline(h = 0, col = "red", lwd = 2)


# --------------------------- 6. Final Model and Estimated Equation ---------------------------

# Assuming the selected_model is final, its mathematical form can be written as:
#   Ȳ = β0̂ + β1̂ * X1 + β3̂ * X3 + β4̂ * X4
# (Adjust the formula when we decide what the fuck we are choosing)
# Here, Ȳ denotes the predicted Mortgage Yield (%)

final_coef <- round(coef(selected_model), 2)
cat("Final Model Equation:\n")
cat("Mortgage Yield (%) = ",
    paste0(final_coef[1], " + ", 
           ifelse("X1" %in% names(final_coef), paste0(final_coef["X1"], "*X1"), ""),
           ifelse("X3" %in% names(final_coef), paste0(" + ", final_coef["X3"], "*X3"), ""),
           ifelse("X4" %in% names(final_coef), paste0(" + ", final_coef["X4"], "*X4"), ""), "\n"))

# --------------------------- 7. Summary Table of Models ---------------------------

# summary table to compare all model variations
model_list <- list(
  "Full Model" = full_model,
  "Selected Model" = selected_model,
  "Log Model" = log_model,
  "Log improved" = selected_model_log
)

model_summary <- data.frame(
  Model = names(model_list),
  p_value = sapply(model_list, function(mod) {
    fstat <- summary(mod)$fstatistic
    pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  }),
  Adj_R2 = sapply(model_list, function(mod) summary(mod)$adj.r.squared),
  AIC = sapply(model_list, AIC),
  BIC = sapply(model_list, BIC)
)

print(model_summary)

plot(selected_model_log, which = 1:4) #-> doesn't look good at all.....

#ok so we started with full models 

# Compute Variance Inflation Factor (VIF)
vif(lm(mortYld ~ X1 + X2 + X3 + X4 + X5 + X6, data = data_mort))

#according to this criteria variables X3 (and X4) is highly correlated wit other!!!
#re-run the model starting with this variable out 

vif(lm(mortYld ~ X1 + X2 + X4 + X5 + X6, data = data_mort))

#take X3 and see what happens
AIC(lm(mortYld ~ X1 + X2 + X4 + X5 + X6, data = data_mort)) #-20

AIC(step(lm(mortYld ~ X1 + X2 + X4 + X5 + X6, data = data_mort))) #-24

model_final <- lm(mortYld ~ X1 + X2 + X5, data = data_mort)
summary(model_final)
coef(model_final)

#take X4 and see what happens
AIC(lm(mortYld ~ X1 + X2 + X3 + X5 + X6, data = data_mort)) #-22

AIC(step(lm(mortYld ~ X1 + X2 + X3 + X5 + X6, data = data_mort))) #-25

model_final <- lm(mortYld ~ X1 + X2 + X3, data = data_mort)
summary(model_final)
coef(model_final)



#blah i dont fully get it why X3 is more correlated but take it out does not help


#---------------------------------------------------------------
# ------ How much of the total variability of the predictand 
# ------     is explained by each of the models? 
#---------------------------------------------------------------

r2_values <- round(sapply(c("X1", "X2", "X3", "X4", "X5", "X6"), function(var) {
  model <- lm(mortYld ~ data_mort[[var]], data = data_mort)
  summary(model)$r.squared
}), 3)  # Round to 3 decimal places

r2_row <- c("R² (from simple regression)", "", r2_values)

# Create table with state names, mortgage yield, and explanatory variables
table_data <- data_mort %>%
  select(smsa, mortYld, X1, X2, X3, X4, X5, X6) %>%
  rename(
    `Mortgage Yield (%)` = mortYld,
    `Avg Loan/Mortgage Ratio` = X1,
    `Distance from Boston (miles)` = X2,
    `Savings per New Unit Built` = X3,
    `Savings per Capita` = X4,
    `Population Increase 1950-1960 (%)` = X5,
    `% First Mortgage from Inter-regional Banks` = X6
  )


# Convert to dataframe for table formatting
final_table <- rbind(table_data, r2_row)

# Print table with formatting
kable(final_table, format = "html", caption = "Mortgage Yield and Explanatory Variables") %>%
  kable_styling(full_width = FALSE) %>%
  row_spec(nrow(final_table), bold = TRUE, background = "#f0f0f0")  # Highlight R² row


# Convert to dataframe for table formatting
final_table <- rbind(table_data, r2_row)
colnames(final_table) <- c(
  "State",
  "\\rotatebox{90}{Mortgage Yield (\\%)}",
  "\\rotatebox{90}{Avg Loan/Mortgage Ratio}",
  "\\rotatebox{90}{Distance from Boston (miles)}",
  "\\rotatebox{90}{Savings per New Unit Built}",
  "\\rotatebox{90}{Savings per Capita}",
  "\\rotatebox{90}{Population Increase 1950-1960 (\\%)}",
  "\\rotatebox{90}{\\% First Mortgage from Inter-regional Banks}"
)

print(xtable(final_table, caption = "Mortgage Yield and Explanatory Variables\\newline (Simple Regression Results)"), 
      type = "latex", 
      caption.placement = "top", 
      include.rownames = FALSE,
      sanitize.colnames.function = identity)
