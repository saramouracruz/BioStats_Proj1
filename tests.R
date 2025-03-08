# You should predict Mortgage Yield in % (mortYld) as a function of:
# Average Loan/Mortgage Ratio  High values ==> Low Down Payments   (X1)
# Distance from Boston in miles  (X2)
# Savings per new unit built   (X3)
# Savings per capita   (X4)
# Pop inc 1950-1960 in %   (X5)
# Percent of first mortgage from inter-regional banks  (X6)


# import data 
library(readr)
library(ggplot2)
library(tidyr)

data_mort <- read.csv("myield.csv")

#kable(data_mort, caption = "Mortgage Yield Data")

str(data_mort)

#visual inspection--------------------------------------------------------------------

#layout to show multiple plots
par(mfrow = c(2, 3))  # 2 rows, 3 columns 

# Loop through each independent variable
for (var in c("X1", "X2", "X3", "X4", "X5", "X6")) {
  plot(data_mort[[var]], data_mort$mortYld, 
       main = paste("mortYld vs", var),
       xlab = var, ylab = "mortYld",
       col = "blue", pch = 16)  
  
  # Add a trend line
  lines(lowess(data_mort[[var]], data_mort$mortYld), col = "red", lwd = 2)
}

#as mentioned in the paper we can see that there is a more or less log relantionship between
#independent variable and Xs

#partial correlations (function given by prof) ----------------------------------------
# = relationship independent of other variables 
# = correlation between two variables after correcting for all others
# better than normal correlation because, when we see high corr it is still 
# open question if it is because of a 3rd variable 

p.cor <- function(x){
  inv <- solve(var(x))
  sdi <- diag(1/sqrt(diag(inv)))
  p.cor.mat <- -(sdi %*% inv %*% sdi)
  diag(p.cor.mat) <- 1
  rownames(p.cor.mat) <- colnames(p.cor.mat) <- colnames(x)
  return(p.cor.mat) }

#adapted function to see corr>threshold

p.cor.threshold <- function(x, threshold = 0.5) {
  inv <- solve(var(x))  
  sdi <- diag(1/sqrt(diag(inv)))  
  p.cor.mat <- -(sdi %*% inv %*% sdi)  
  diag(p.cor.mat) <- 1  
  rownames(p.cor.mat) <- colnames(p.cor.mat) <- colnames(x)  
  
  # Get upper triangle (avoid duplicates)
  p.cor.mat[lower.tri(p.cor.mat, diag = TRUE)] <- NA  
  
  # Find pairs with correlation > threshold
  significant_pairs <- which(abs(p.cor.mat) > threshold, arr.ind = TRUE)
  
  # Print results
  if (nrow(significant_pairs) > 0) {
    for (i in 1:nrow(significant_pairs)) {
      row <- rownames(p.cor.mat)[significant_pairs[i, 1]]
      col <- colnames(p.cor.mat)[significant_pairs[i, 2]]
      value <- p.cor.mat[significant_pairs[i, 1], significant_pairs[i, 2]]
      cat(sprintf("Correlation between %s and %s: %.3f\n", row, col, value))
    }
  } else {
    cat("No correlations above the threshold found.\n")
  }
  
  return(p.cor.mat)  # Return full matrix 
}

p.cor.threshold(data_mort[c(2:8)])

#conclusions: 
# related to the predictant (mortageyield) = all variables haev corr!=0 so might be useful for predictions
#the variables with higher partial corr are to mortageyield are X1 and X3
# between predictors: we see at least 3 pairs high correlated which might be useful to know afterwards 
# X2,X3
# X2,X4
# X3,X4

# building models (linear) ----------------------------------------------------------------------
#we can use step in order to let R fit the model, and subtract the variables one by one 

# let stepwise pick the best from a full model

mortage.formula <- mortYld ~ (X1 + X2 + X3 + X4 + X5 + X6)
mortage.lm <- lm(mortage.formula, data = data_mort)

a <- step(lm(mortage.formula, data = data_mort))

#how to interpret?
#we want Akaike’s Information Criterion (AIC) to be as lower as possible so R tries iteratively
#to remove variables until there is no bigger improvement 
# best model = Step:  AIC=-82.81
# mortYld ~ X1 + X3 + X4 

best.formula <- mortYld ~ (X1 + X3 + X4)
best.lm <- lm(best.formula, data = data_mort)

layout(matrix(1:6,ncol=3))
#plot full model
plot(mortage.lm, which = c(1,2,3,4,5,6))
summary(mortage.lm)
coef(mortage.lm)

#plot reduced model
plot(best.lm, which = c(1,2,3,4,5,6))
summary(best.lm)
coef(best.lm)
#better R^2 and p-value


#explore better the plots but looks better: cooks distance and residuals

#let's try new stuff - paper uses log!

log_mortage.formula <- mortYld ~ (log(X1) + log(X3) + log(X4))
a <- step(lm(log_mortage.formula, data = data_mort))
log.lm <- lm(log_mortage.formula, data = data_mort)
summary(log.lm)
plot(log.lm, which = c(1,2,3,4,5,6))

log_mortage.formula <- mortYld ~ (log(X1) + log(X2+0.00000001) + log(X3) + log(X5) + log(X6))
#ok here acctually gives the same variables that the paper identifies which is gutten
a <- step(lm(log_mortage.formula, data = data_mort))

log_mortage_corrected.formula <- mortYld ~ (log(X1) + log(X2+0.00000001)+ log(X5))
log_corrected.lm <- lm(log_mortage_corrected.formula, data = data_mort)
summary(log_corrected.lm)
plot(log.lm, which = c(1,2,3,4,5,6))



