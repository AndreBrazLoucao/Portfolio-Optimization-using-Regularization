#---Install relevant packages---
#install.packages("quadprog")
#install.packages("quantmod")
#install.packages("MASS")
#install.packages("corrplot")
#install.packages("Hmisc")
#install.packages("emdbook")
#install.packages("dplyr")
#install.packages("pracma")
#install.packages("optiSolve")
#---Load Packages---
library(optiSolve)
library(dplyr)
library(plyr)
library(quantmod)
library(MASS)
library(Hmisc)
library(emdbook)

#List of required packages
packages <- c(
  "optiSolve",
  "dplyr",
  "plyr",
  "quantmod",
  "MASS",
  "Hmisc",
  "emdbook")

# Install and load missing packages
packages_New <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(packages_New)) install.packages(packages_New)
lapply(packages, require, character.only=TRUE)

####################################################################################################
#The functions contained in R Script "RegularizedPortfolioFunctions" are needed to run this script
####################################################################################################

#Setting a seed for reproducable results
set.seed(1654359)

#Set working directory
setwd("H:/Thesis/R Project - Portfolio Optimization")

#--------------------------------------------------------------
#Read data
#--------------------------------------------------------------

#48 Industry Portfolios
indPF <- read.delim(file = "48_Industry_Portfolios_Daily.txt", header = T, dec = ".", sep = "", row.names=NULL)

#Risk Free Rate
RiskFreeRate <-  read.delim(file = "FF_3Factors.txt", header = T, dec = ".", sep = "", row.names=NULL)[-c(2,3,4)]

#--------------------------------------------------------------
#Data Preparation
#--------------------------------------------------------------
#---48 Industry Portfolios---
#Convert dates
colnames(indPF)[1] <- "Date"

indPF$Date <- as.Date(indPF$Date, "%Y%m%d")

#Remove irrelevant rows
indPF <- indPF[-c(0:45011),]

#Convert Returns to numeric
indPF[,2:49] <- sapply(indPF[,2:49],as.numeric)

#Convert to time series object
indPF <- xts(indPF[,2:49], order.by = as.Date(indPF[,1]))

#Split into months for rolling window
splitPF <-  split(indPF, format(index(indPF), "%Y-%m"))


#---Risk Free Rate---
#Convert dates
colnames(RiskFreeRate)[1] <- "Date"

RiskFreeRate$Date <- as.Date(RiskFreeRate$Date, "%Y%m%d")

#Remove irrelevant rows
RiskFreeRate <- RiskFreeRate[RiskFreeRate[["Date"]] >= "2001-07-02", ]

#Split into months for rolling window
splitRF <-  split(RiskFreeRate, format(index(indPF), "%Y-%m"))

#--------------------------------------------------------------
#----Preliminary Data Analysis----
#--------------------------------------------------------------

#Check for NA values
which(is.na.data.frame(indPF))

#Check for values -99.99 and -999 (values that are used by Fama and French to indicate missing values)
which(indPF == -99.99)

which(indPF == -999)



#--------------------------------------------------------------
#----Benchmark Equal Weights Portfolio----
#--------------------------------------------------------------

#Remove values that other models use for sample covariance matrix estimation 
trunc_indPF <- indPF[-c(1:247)]

#Split into months
splitEW <-  split(trunc_indPF, format(index(trunc_indPF), "%Y-%m"))

#Portfolio Returns with equal Weights
returnsEW <- rowSums(1/48 * trunc_indPF)
#Convert return to time series object
returnsEW <- xts(returnsEW, order.by = index(trunc_indPF) )


#Create weight vector of equal weights
EqualWeights <- as.vector( rep(1/48, 48) )

#Pre allocate vector for monthly variance values
MonthlyRisk <- vector("double", length(splitPF) )

#Compute equal weights portfolio variance over whole dataset
for (i in seq_along(splitPF)) {
  MonthlyRisk[[i]] <- t(EqualWeights) %*% cov(splitPF[[i]]) %*% EqualWeights
}

#Average variance of portfolio over whole period
EWRisk <- mean(MonthlyRisk)
print(EWRisk)

#Pre allocate Sharpe ratio vector
sharpeEW <- vector("double", length(splitEW) )

#Compute Sharpe ratio for all iterations of rolling window
for (i in seq_along(splitEW)) {
  sharpeEW[[i]] <- mean( rowSums(splitEW[[i]] * EqualWeights ) -  splitRF[[i]][1,2])  / sd(rowSums(splitEW[[i]] * EqualWeights)) * sqrt(252)
}

#Take mean of all Sharpe ratios
EWsharpe <- mean(sharpeEW)
print(EWsharpe)
#--------------------------------------------------------------
#----Minimum variance portfolio----
#--------------------------------------------------------------

#Pre-allocate data frame and vectors with values for evaluation, weights and returns
Evaluation <- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Validation_Risk", "Test_Risk", "Sharpe"))))

weights <- vector(mode = "list", length = ( length(splitPF)-13 ))

returns <- vector(mode = "list", length = (length(splitPF)-13))

#Loop over full data set for each rolling window iteration
for (i in 1:227) {
  #Covariance matrix estimation time frame
  covMat_end <- i+11
  
  #Validation Month is equivalent to test month due to missing cross-validation
  val_month <- i+12
  
  #Theoretical test month if cross-validation would be used
  test_month <- i+13
  
  #Unlist returns to compute covariance in next line
  ret <-  ldply(splitPF[i:covMat_end], data.frame, .id = NULL)
  
  #Fit model
  fit <- MinimumVariance(Cov = cov(ret))

  #Record evaluation criteria, weights and returns
  Validation_Risk <- t(fit$Weights) %*% cov(splitPF[[val_month]]) %*% fit$Weights
  
  Test_Risk <- t(fit$Weights) %*% cov(splitPF[[test_month]]) %*% fit$Weights
  
  sharpe <- mean( rowSums(splitPF[[test_month]] * fit$Weights) -  splitRF[[i]][1,2])  / sd( rowSums(splitPF[[test_month]] * fit$Weights)) * sqrt(252)
  
  weights[[i]] <- fit$Weights
  
  returns[[i]] <- t(t(splitPF[[test_month]]) * fit$Weights)
  
  #Add row of i'th portfolios evaluation criteria to evaluation data frame
  Evaluation[nrow(Evaluation) + 1,] = c(Validation_Risk, Test_Risk, sharpe)
}

MinimumVariance <- list(
  Evaluation = Evaluation,
  Weights = weights,
  Returns = returns
)

#--------------------------------------------------------------
#----No Short Sale Portfolio----
#--------------------------------------------------------------

#Pre-allocate data frame with values for evaluation, weights and returns
weights <- vector(mode = "list", length =  (length(splitPF)-13) )

Evaluation <- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Validation_Risk", "Test_Risk", "Sharpe"))))

returns <- vector(mode = "list", length = (length(splitPF)-13))

#Loop over full data set for each rolling window iteration
for (i in 1:227) {
  #Define sample covariance, validation and test set
  covMat_end <- i+11
  val_month <- i+12
  test_month <- i+13
  
  #Compute covariance matrix of assets
  cov <-  cov(do.call("rbind", splitPF[i:covMat_end]))
  
  #Fit model
  fit <- ConstraintLasso(Cov = cov, lambda = 1)
  
  #Record evaluation criteria, weights and returns
  Validation_Risk <- t(fit$Weights) %*% cov(splitPF[[val_month]]) %*% fit$Weights
  
  Test_Risk <- t(fit$Weights) %*% cov(splitPF[[test_month]]) %*% fit$Weights
  
  sharpe <- mean( rowSums(splitPF[[test_month]] * fit$Weights) - splitRF[[i]][1,2]) / sd( rowSums(splitPF[[test_month]] * fit$Weights)) * sqrt(252)
  
  weights[[i]] <- fit$Weights
  
  returns[[i]] <- t(t(splitPF[[test_month]]) * fit$Weights)
  
  #Add row of i'th portfolios evaluation criteria to evaluation data frame
  Evaluation[nrow(Evaluation) + 1,] = c(Validation_Risk, Test_Risk, sharpe)
}

NoShort <- list(
  Evaluation = Evaluation,
  Weights = weights,
  Returns = returns
)

#--------------------------------------------------------------
#----Lasso----
#--------------------------------------------------------------
#--Weight computation for whole data set--
#Pre-allocate evaluation variables
risks <- vector(mode = "list", length = ( length(splitPF)-13 ))
lambdas <- vector(mode = "list", length = ( length(splitPF)-13 ))
weights <- vector(mode = "list", length = ( length(splitPF)-13 ))
risks_test <- vector(mode = "list", length = ( length(splitPF)-13 ))
returns <- vector(mode = "list", length = (length(splitPF)-13))

#Data frame with evaluation criteria 
Evaluation <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("Lambda", "Validation_Risk", "Test_Risk", "Sharpe"))))

#Loop over every iteration of rolling window
for (i in 1:227) {
  #Define sample covariance, validation and test set
  covMat_end <- i+11
  val_month <- i+12
  test_month <- i+13
  
  #Fit all models using cross-validation
  cv_fit <- lasso_cv(Input = splitPF[i:covMat_end], Target = splitPF[[val_month]])
  
  #Evaluation Metrics
  weights[[i]] <- cv_fit$WeightsOpt
 
  risks_test[[i]] <- t(cv_fit$WeightsOpt) %*% cov(splitPF[[test_month]]) %*% cv_fit$WeightsOpt 
  
  sharpe <- mean( rowSums(weights[[i]] * splitPF[[test_month]])-  splitRF[[i]][1,2]) / sd( rowSums(weights[[i]] * splitPF[[test_month]])) * sqrt(252)
  
  returns[[i]] <- t(t(splitPF[[test_month]]) * cv_fit$Weights)
  
  Evaluation[nrow(Evaluation) + 1,] = c(cv_fit$LambdaOpt[[1]], cv_fit$LambdaOpt[[2]], risks_test[[i]], sharpe)
}

Lasso <- list(
  Evaluation = Evaluation,
  Weights = weights,
  Returns = returns
)

#--------------------------------------------------------------
#----Ridge----
#--------------------------------------------------------------
#--Weight Computation for whole data set--
#The procedure is identical to the Lasso Portfolio with the exception of using function ridge_cv() instead of lasso_cv()

#Pre-allocate evaluation objects
risks <- vector(mode = "list", length = ( length(splitPF)-13 ))
lambdas <- vector(mode = "list", length = ( length(splitPF)-13 ))
weights <- vector(mode = "list", length = ( length(splitPF)-13 ))
risks_test <- vector(mode = "list", length = ( length(splitPF)-13 ))
returns <- vector(mode = "list", length = (length(splitPF)-13))

Evaluation <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("Lambda", "Validation_Risk", "Test_Risk", "Sharpe"))))

#Loop over full data set
for (i in 1:227) {
  
  covMat_end <- i+11
  val_month <- i+12
  test_month <- i+13
  
  #Fit all models using cross-validation
  cv_fit <- ridge_cv(Input = splitPF[i:covMat_end], Target = splitPF[[val_month]])
  
  #Evaluation Metrics
  weights[[i]] <- cv_fit$WeightsOpt
  
  risks_test[[i]] <- t(cv_fit$WeightsOpt) %*% cov(splitPF[[test_month]]) %*% cv_fit$WeightsOpt 
  
  sharpe <- mean( rowSums(weights[[i]] * splitPF[[test_month]]) -  splitRF[[i]][1,2]) / sd( rowSums(weights[[i]] * splitPF[[test_month]]))* sqrt(252)
  
  returns[[i]] <- t(t(splitPF[[test_month]]) * cv_fit$Weights)
  
  Evaluation[nrow(Evaluation) + 1,] = c(cv_fit$LambdaOpt[[1]], cv_fit$LambdaOpt[[2]], risks_test[[i]], sharpe)

}


#--Verifying Budget Constraint--

#List all sums of the weight vectors for each held portfolio
budget_constr <- unlist(lapply(weights, sum))

print("All values need to equal 1 for Budget constraint to hold:")
print(budget_constr)

#--Verifying Ridge Constraint--
#-Check for percentage deviation of ridge constraint parameter from actual sum of squared weights-

#Compute squared weights
SquaredWeights <- lapply(weights, function(weights) weights^2)

#Define Vector to track percentage deviations
ConstraintDeviation <- vector( length = ( length(splitPF)-13 ))

#Compute deviation of actual sum of squares versus passed lambda value for all months
for (i in 1:227) {
  ConstraintDeviation[i] <- abs(((sum(SquaredWeights[[i]]) - Evaluation[i,1]) / Evaluation[i,1]) * 100)
}

#Maximum Deviation over whole period
max_dev <- max(ConstraintDeviation)

#Mean Deviation over whole period
mean_dev <- mean(ConstraintDeviation)

print("Maximum deviation from ridge constraint:")
print(mean_dev)

#Mean Deviation over whole period
mean_dev <- mean(ConstraintDeviation)

print("Mean deviation from ridge constraint:")
print(mean_dev)

#Add Approximated Lambda Value to Evaluation data frame
Evaluation$LamdbaApprox <- unlist(lapply(SquaredWeights, sum))

Ridge <- list(
  Evaluation = Evaluation,
  Weights = weights,
  Returns = returns
)


#--------------------------------------------------------------
#----Partially Egalitarian Lasso Step 2: Simple Average----
#--------------------------------------------------------------
#--Weight Computation for whole data set--
#The procedure is identical to the Lasso and Ridge portfolio, with peLassoAV_cv

risks <- vector(mode = "list", length = ( length(splitPF)-13 ))
lambdas <- vector(mode = "list", length = ( length(splitPF)-13 ))
weights <- vector(mode = "list", length = ( length(splitPF)-13 ))
risks_test <- vector(mode = "list", length = ( length(splitPF)-13 ))
returns <- vector(mode = "list", length = (length(splitPF)-13))

Evaluation <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("Lambda", "Validation_Risk", "Test_Risk", "Sharpe"))))

for (i in 1:227) {
  
  covMat_end <- i+11
  val_month <- i+12
  test_month <- i+13
  
  cv_fit <- peLassoAV_cv(Input = splitPF[i:covMat_end], Target = splitPF[[val_month]])
  
  #Evaluation Metrics
  weights[[i]] <- cv_fit$WeightsOpt
  
  risks_test[[i]] <- t(cv_fit$WeightsOpt) %*% cov(splitPF[[test_month]]) %*% cv_fit$WeightsOpt 
  
  sharpe <- mean( rowSums(weights[[i]] * splitPF[[test_month]]) -  splitRF[[i]][1,2]) / sd( rowSums(weights[[i]] * splitPF[[test_month]])) * sqrt(252)
  
  returns[[i]] <- t(t(splitPF[[test_month]]) * cv_fit$Weights)
  
  Evaluation[nrow(Evaluation) + 1,] = c(cv_fit$LambdaOpt[[1]], cv_fit$LambdaOpt[[2]], risks_test[[i]], sharpe)
}

peLassoAverage <- list(
  Evaluation = Evaluation,
  Weights = weights,
  Returns = returns
)

#--------------------------------------------------------------
#----Partially egalitarian lasso (Ridge)----
#--------------------------------------------------------------
#--Weight Computation for whole Dataset--
#Pre-allocate evaluation objects
risks <- vector(mode = "list", length = ( length(splitPF)-13 ))
lambdas <- vector(mode = "list", length = ( length(splitPF)-13 ))
weights <- vector(mode = "list", length = ( length(splitPF)-13 ))
risks_test <- vector(mode = "list", length = ( length(splitPF)-13 ))
returns <- vector(mode = "list", length = (length(splitPF)-13))

Evaluation <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("Lambda (Lasso)","Lambda (Ridge)", "Validation_Risk", "Test_Risk", "Sharpe", "NonZero Weights"))))

#Loop over all rolling window iterations
for (i in 1:227) {
  #Define sample covariance, validation and test set
  covMat_end <- i+11
  val_month <- i+12
  test_month <- i+13
  
  #First step: Cross-validation using the lasso
  cv_fit <- lasso_cv(Input = splitPF[i:covMat_end], Target = splitPF[[val_month]])
  
  #Store index of weights that are not set to zero
  NonZeroIndex <- which(cv_fit$WeightsOpt != 0)
  
  #Compute covariance of non zero assets
  CovNoNull <- cov(do.call("rbind", splitPF[i:covMat_end]))[c(NonZeroIndex),c(NonZeroIndex)]
  
  #Store number of weight that are not zero
  K <- ncol(CovNoNull)
  
  #Second Step: Cross-validation using the ridge and the weights that were not set to one in the second step
  cv_fitRidge <- ridge_cv(Input = CovNoNull, Target = splitPF[[val_month]][,c(NonZeroIndex)], CovInput = T)
  
  #Evaluation Metrics
  weights[[i]] <- cv_fitRidge$WeightsOpt
  
  risks_test[[i]] <- t(cv_fitRidge$WeightsOpt) %*% cov(splitPF[[test_month]][,c(NonZeroIndex)]) %*% cv_fitRidge$WeightsOpt 
  
  sharpe <- mean( rowSums(weights[[i]] * splitPF[[test_month]][,c(NonZeroIndex)])-  splitRF[[i]][1,2]) / sd( rowSums(weights[[i]] * splitPF[[test_month]][,c(NonZeroIndex)])) * sqrt(252)
  
  returns[[i]] <- t(t(splitPF[[test_month]]) * cv_fit$Weights)
  
  Evaluation[nrow(Evaluation) + 1,] = c(cv_fit$LambdaOpt[[1]], cv_fitRidge$LambdaOpt[[1]], cv_fitRidge$LambdaOpt[[2]], risks_test[[i]], sharpe, K)
}

#--Verifying Ridge Constraint--
#-Check for percentage deviation of ridge constraint parameter from actual sum of squared weights-

#Compute squared weights
SquaredWeights <- lapply(weights, function(weights) weights^2)

#Define Vector to track percentage deviations
ConstraintDeviation <- vector( length = ( length(splitPF)-13 ))

#Compute deviation of actual sum of squares versus passed lambda value for all months
for (i in 1:227) {
  ConstraintDeviation[i] <- abs(((sum(SquaredWeights[[i]]) - Evaluation[i,2]) / Evaluation[i,2]) * 100)
}

#--Verifying Budget Constraint--

#List all sums of the weight vectors for each held portfolio
budget_constr <- unlist(lapply(weights, sum))

print("All values need to equal 1 for Budget constraint to hold:")
print(budget_constr)

#Maximum Deviation over whole period
max_dev <- max(ConstraintDeviation)

print("Maximum deviation from ridge constraint:")
print(mean_dev)

#Mean Deviation over whole period
mean_dev <- mean(ConstraintDeviation)

print("Mean deviation from ridge constraint:")
print(mean_dev)

#Add Approximated Lambda Value to Evaluation data frame
Evaluation$LamdbaApprox <- unlist(lapply(SquaredWeights, sum))

peLassoRidge <- list(
  Evaluation = Evaluation,
  Weights = weights,
  Returns = returns
)


#------------------------------------------------------------------------------------------
#Lasso Portfolio with fixed Lambda Values
#Lambda = 1.5
Lasso15 <- FixedLasso(lambda = 1.5)

#Lambda = 2
Lasso2 <- FixedLasso(lambda = 2)

#Lambda = 3
Lasso3 <- FixedLasso(lambda = 3)



#########################
#Remove all but irrelevant objects
rm(list= ls()[!(ls() %in% c('indPF',
                            'Lasso','Lasso15','Lasso2','Lasso3',
                            'MinimumVariance','NoShort','returnsEW','Ridge','RiskFreeRate','splitEW','splitRF','splitPF',
                            'EWRisk',
                            'ConstraintDeviation',
                            'peLassoAverage',
                            'peLassoRidge'
                            ))])
#------------------------------------------------------------------------------------------------------------
#----Evaluation----
#------------------------------------------------------------------------------------------------------------
#--Average of Evaluation Criteria--
#NoShort (Lambda = 1)
print("Average of Evaluation Criteria (Short Sale Constrained Portfolio):")
colMeans(NoShort$Evaluation)

#Lambda = 1.5
print("Average of Evaluation Criteria (Lasso (Lambda = 1.5)):")
colMeans(Lasso15$Evaluation)

#Lambda = 2
print("Average of Evaluation Criteria (Lasso (Lambda = 2)):")
colMeans(Lasso2$Evaluation)

#Lambda = 3
print("Average of Evaluation Criteria (Lasso (Lambda = 1.5))):")
colMeans(Lasso3$Evaluation)

#MinimumVariance
print("Average of Evaluation Criteria (Minimum Variance Portfolio):")
colMeans(MinimumVariance$Evaluation)

#Lasso
print("Average of Evaluation Criteria (Lasso (cross-validated)):")
colMeans(Lasso$Evaluation)

#Ridge
print("Average of Evaluation Criteria (Ridge):")
colMeans(Ridge$Evaluation)

#peLasso Average
print("Average of Evaluation Criteria (peLasso (Average)):")
colMeans(peLassoAverage$Evaluation)

#peLasso Ridge
print("Average of Evaluation Criteria (peLaso (Ridge)):")
colMeans(peLassoRidge$Evaluation)

#Non-Zero and Short Positions
positions <- data.frame(c("Minimum Variance",
                          "Short Sale Constrained",
                          "Lasso (Lambda = 1.5)",
                          "Lasso (Lambda = 2)",
                          "Lasso (Lambda = 3)",
                          "Lasso (cross-validated)",
                          "Ridge",
                          "peLasso (Average)",
                          "peLasso (Ridge)"),
                        c(
                          (sum(as.numeric(lapply(MinimumVariance$Weights, function(weights) sum(weights != 0)))) / 227),
                          (sum(as.numeric(lapply(NoShort$Weights, function(weights) sum(weights != 0)))) / 227),
                          (sum(as.numeric(lapply(Lasso15$Weights, function(weights) sum(weights != 0)))) / 227),
                          (sum(as.numeric(lapply(Lasso2$Weights, function(weights) sum(weights != 0)))) / 227),
                          (sum(as.numeric(lapply(Lasso3$Weights, function(weights) sum(weights != 0)))) / 227),
                          (sum(as.numeric(lapply(Lasso$Weights, function(weights) sum(weights != 0)))) / 227),
                          (sum(as.numeric(lapply(Ridge$Weights, function(weights) sum(weights != 0)))) / 227),
                          (sum(as.numeric(lapply(peLassoAverage$Weights, function(weights) sum(weights != 0)))) / 227),
                          (sum(as.numeric(lapply(peLassoRidge$Weights, function(weights) sum(weights != 0)))) / 227)),
                        c(
                          (sum(as.numeric(lapply(MinimumVariance$Weights, function(weights) sum(weights < 0)))) / 227),
                          #Short sale constrained portfolio
                          0,
                          (sum(as.numeric(lapply(Lasso15$Weights, function(weights) sum(weights < 0)))) / 227),
                          (sum(as.numeric(lapply(Lasso2$Weights, function(weights) sum(weights < 0)))) / 227),
                          (sum(as.numeric(lapply(Lasso3$Weights, function(weights) sum(weights < 0)))) / 227),
                          (sum(as.numeric(lapply(Lasso$Weights, function(weights) sum(weights < 0)))) / 227),
                          (sum(as.numeric(lapply(Ridge$Weights, function(weights) sum(weights < 0)))) / 227),
                          (sum(as.numeric(lapply(peLassoAverage$Weights, function(weights) sum(weights < 0)))) / 227),
                          (sum(as.numeric(lapply(peLassoRidge$Weights, function(weights) sum(weights < 0)))) / 227)
                        ))


print(positions)

positions <- data.frame( Portfolio = c("Minimum Variance",
                                       "Short Sale Constrained",
                                       "Lasso (Lambda = 1.5)",
                                       "Lasso (Lambda = 2)",
                                       "Lasso (Lambda = 3)",
                                       "Lasso (cross-validated)",
                                       "Ridge",
                                       "peLasso (Average)",
                                       "peLasso (Ridge)"),
                         ActivePositions = c((sum(as.numeric(lapply(MinimumVariance$Weights, function(weights) sum(weights != 0)))) / 227),
                                             (sum(as.numeric(lapply(NoShort$Weights, function(weights) sum(weights != 0)))) / 227),
                                             (sum(as.numeric(lapply(Lasso15$Weights, function(weights) sum(weights != 0)))) / 227),
                                             (sum(as.numeric(lapply(Lasso2$Weights, function(weights) sum(weights != 0)))) / 227),
                                             (sum(as.numeric(lapply(Lasso3$Weights, function(weights) sum(weights != 0)))) / 227),
                                             (sum(as.numeric(lapply(Lasso$Weights, function(weights) sum(weights != 0)))) / 227),
                                             (sum(as.numeric(lapply(Ridge$Weights, function(weights) sum(weights != 0)))) / 227),
                                             (sum(as.numeric(lapply(peLassoAverage$Weights, function(weights) sum(weights != 0)))) / 227),
                                             (sum(as.numeric(lapply(peLassoRidge$Weights, function(weights) sum(weights != 0)))) / 227)),
                         ShortSales = c((sum(as.numeric(lapply(MinimumVariance$Weights, function(weights) sum(weights < 0)))) / 227),
                                        #Short sale constrained portfolio
                                        0,
                                        (sum(as.numeric(lapply(Lasso15$Weights, function(weights) sum(weights < 0)))) / 227),
                                        (sum(as.numeric(lapply(Lasso2$Weights, function(weights) sum(weights < 0)))) / 227),
                                        (sum(as.numeric(lapply(Lasso3$Weights, function(weights) sum(weights < 0)))) / 227),
                                        (sum(as.numeric(lapply(Lasso$Weights, function(weights) sum(weights < 0)))) / 227),
                                        (sum(as.numeric(lapply(Ridge$Weights, function(weights) sum(weights < 0)))) / 227),
                                        (sum(as.numeric(lapply(peLassoAverage$Weights, function(weights) sum(weights < 0)))) / 227),
                                        (sum(as.numeric(lapply(peLassoRidge$Weights, function(weights) sum(weights < 0)))) / 227)
                         ))

print(positions)
