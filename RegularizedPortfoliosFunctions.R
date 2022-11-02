#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Author:             Andre Braz Loucao
# Short Description:  Functions for regularized portfolio optimization
#----------------------------------------------------------------------
#----------------------------------------------------------------------

# Clear Workspace
####rm(list=ls())

# Which Packages are required?
packages <- c(
  "optiSolve",
  "Matrix"
)

# Installing and loading missing packages
packages_New <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(packages_New)) install.packages(packages_New)
lapply(packages, require, character.only=TRUE)

#----------------------------------------------------------------------
#----------------------------------------------------------------------

#--------------------------------------------------------------
#----Function to round and preserve the sum----
#--------------------------------------------------------------


round_preserve_sum <- function(x, digits = 0) {
  up <- 10 ^ digits
  x <- x * up
  z <- floor(x)
  indices <- tail(order(x-z), round(sum(x)) - sum(z))
  z[indices] <- z[indices] + 1
  z / up
}




#--------------------------------------------------------------
#----Function for calculating the weights for the Constraint Lasso----
#--------------------------------------------------------------

# Cov -- Covariance Matrix N x N
# lambda -- lasso Parameter (default is 1)

ConstraintLasso <- function(Cov, lambda = 1){
  # Number of Assets 
  N <- ncol(Cov)
  
  # There are sometimes Problems if both Unity and Lasso Constraint are 1
  # Therefore, I use a if condition. In the first case for lambda = 1
  # the lasso constraint is equal to the unity constraint and therefore not considered.
  
  
  
  # lambda == 1
  ##################################
  
  if(lambda == 1){
    # Objective Function
    Q  <- optiSolve::quadfun(Cov)
    
    # All Variables greater zero
    lb <- optiSolve::lbcon(val = 0, id = 1:N)
    
    # Linear Constraint
    A <- optiSolve::lincon(
      A   = matrix(rep(1,N), nrow = 1, byrow = T),
      dir = c("=="),
      val = c(1),
      name = c("Unity")
    )
    
    # Define the Optimization Problem
    OptimizationProblem <- optiSolve::cop(f = Q, max = FALSE, lb = lb, lc = A)
    # Solve the Optimization Problem
    result <- optiSolve::solvecop(OptimizationProblem, quiet = TRUE)
    
    # Round Sums preserving sum (4 Digits)
    weights <- round_preserve_sum(
      x      = result$x[1:N],
      digits = 4
    )
  } else {
    
    # lambda > 1
    ##################################
    
    # Transforming the Matrix for the lasso constraint
    Cov <- cbind(rbind(Cov, -Cov), rbind(-Cov, Cov))
    
    # Objective Function
    Q  <- optiSolve::quadfun(Cov)
    
    # All Variables greater zero
    lb <- optiSolve::lbcon(val = 0, id = 1:(2*N))
    
    # Linear Constraint
    A <- optiSolve::lincon(
      A   = matrix(c(rep(1,N), rep(-1,N),rep(1,2*N)), nrow = 2, byrow = T),
      dir = c("==", "<="),
      val = c(1, lambda),
      name = c("Unity", "Lasso")
    )
    
    # Define the Optimization Problem
    OptimizationProblem <- optiSolve::cop(f = Q, max = FALSE, lb = lb, lc = A)
    # Solve the Optimization Problem
    result <- optiSolve::solvecop(OptimizationProblem, quiet = TRUE)
    
    # Round Sums preserving sum (4 Digits)
    weights <- round_preserve_sum(
      x      = result$x[1:N] - result$x[(N+1):(2*N)],
      digits = 4
    )
  }
  
  
  # Objective Function
  obj <- as.numeric(t(weights) %*% Cov[1:N,1:N] %*% weights)
  
  return(list(
    Weights            = weights,
    LassoParameter     = lambda,
    ObjectiveFunction  = obj,
    SumWeights         = sum(weights),
    SumAbsoluteWeights = sum(abs(weights))
  ))
}

#--------------------------------------------------------------
#----Function for calculating the weights for the constrained ridge----
#--------------------------------------------------------------

# Cov -- Covariance Matrix N x N
# lambda -- Ridge Parameter (default is 1)

ConstraintRidge <- function(Cov, lambda = 1){
  # Number of Assets 
  N <- ncol(Cov)
  
  
  
  # Objective Function
  Q  <- optiSolve::quadfun(Cov)
  
  # Linear Constraint (budget constraint)
  A <- optiSolve::lincon(
    A   = matrix(rep(1,N), nrow = 1, byrow = T),
    dir = c("=="),
    val = c(1),
    name = c("Unity")
  )
  
  # Quadratic Constraint (sum of squared weights)
  Aquad <- optiSolve::quadcon(
    Q = diag(N),
    dir = c("<="),
    val = c(lambda),
    name = c("Ridge")
    
  )
  
  
  # Define the Optimization Problem
  OptimizationProblem <- optiSolve::cop(f = Q, max = FALSE, lc = A, qc = Aquad)
  # Solve the Optimization Problem
  
  result <- optiSolve::solvecop(OptimizationProblem, quiet = TRUE)
  
  # Round Sums preserving sum (4 Digits)
  weights <- round_preserve_sum(
    x      = result$x[1:N],
    digits = 4
  )
  
  
  # Objective Function
  obj <- as.numeric(t(weights) %*% Cov[1:N,1:N] %*% weights)
  
  return(list(
    Weights            = weights,
    RidgeParameter     = lambda,
    ObjectiveFunction  = obj,
    SumWeights         = sum(weights),
    SumSquaredWeights = sum((weights)^2)
  ))
}

#--------------------------------------------------------------
#----Function for calculating the weights for the minimum variance portfolio----
#--------------------------------------------------------------

# Cov -- Covariance Matrix N x N
# lambda -- lasso Parameter (default is 1)

MinimumVariance <- function(Cov){
  # Number of Assets 
  N <- ncol(Cov)
  
  
  
  # Objective Function
  Q  <- optiSolve::quadfun(Cov)
  
  # Linear Constraint
  A <- optiSolve::lincon(
    A   = matrix(rep(1,N), nrow = 1, byrow = T),
    dir = c("=="),
    val = c(1),
    name = c("Unity")
  )
  
  
  # Define the Optimization Problem
  OptimizationProblem <- optiSolve::cop(f = Q, max = FALSE, lc = A)
  # Solve the Optimization Problem
  
  result <- optiSolve::solvecop(OptimizationProblem, quiet = TRUE)
  
  # Round Sums preserving sum (4 Digits)
  weights <- round_preserve_sum(
    x      = result$x[1:N],
    digits = 4
  )
  
  
  # Objective Function
  obj <- as.numeric(t(weights) %*% Cov[1:N,1:N] %*% weights)
  
  return(list(
    Weights            = weights,
    ObjectiveFunction  = obj,
    SumWeights         = sum(weights)
  ))
}


#--------------------------------------------------------------
#----Function for calculating the weights for the partially egalitarian lasso (Average)----
#--------------------------------------------------------------

# Cov -- Covariance Matrix N x N
# lambda -- lasso Parameter (default is 1)

peLasso_Average <- function(Cov, lambda = 1){
  # Number of Assets 
  N <- ncol(Cov)
  
  # There are sometimes Problems if both Unity and Lasso Constraint are 1
  # Therefore, I use a if condition. In the first case for lambda = 1
  # the lasso constraint is equal to the unity constraint and therefore not considered.
  
  
  
  # lambda == 1
  ##################################
  
  if(lambda == 1){
    # Objective Function
    Q  <- optiSolve::quadfun(Cov)
    
    # All Variables greater zero
    lb <- optiSolve::lbcon(val = 0, id = 1:N)
    
    # Linear Constraint
    A <- optiSolve::lincon(
      A   = matrix(rep(1,N), nrow = 1, byrow = T),
      dir = c("=="),
      val = c(1),
      name = c("Unity")
    )
    
    # Define the Optimization Problem
    OptimizationProblem <- optiSolve::cop(f = Q, max = FALSE, lb = lb, lc = A)
    # Solve the Optimization Problem
    result <- optiSolve::solvecop(OptimizationProblem, quiet = TRUE)
    
    # Round Sums preserving sum (4 Digits)
    weights <- round_preserve_sum(
      x      = result$x[1:N],
      digits = 4
    )
    
    #Averaging weights by taking the simple average
    #Number of nonzero weights
    K <- sum(weights != 0)
    #Set surviving weights to 1/K
    weights[which(weights !=0)] <- 1/K
    
  } else {
    
    # lambda > 1
    ##################################
    
    # Transforming the Matrix for the lasso constraint
    Cov <- cbind(rbind(Cov, -Cov), rbind(-Cov, Cov))
    
    # Objective Function
    Q  <- optiSolve::quadfun(Cov)
    
    # All Variables greater zero
    lb <- optiSolve::lbcon(val = 0, id = 1:(2*N))
    
    # Linear Constraint
    A <- optiSolve::lincon(
      A   = matrix(c(rep(1,N), rep(-1,N),rep(1,2*N)), nrow = 2, byrow = T),
      dir = c("==", "<="),
      val = c(1, lambda),
      name = c("Unity", "Lasso")
    )
    
    # Define the Optimization Problem
    OptimizationProblem <- optiSolve::cop(f = Q, max = FALSE, lb = lb, lc = A)
    # Solve the Optimization Problem
    result <- optiSolve::solvecop(OptimizationProblem, quiet = TRUE)
    
    # Round Sums preserving sum (4 Digits)
    weights <- round_preserve_sum(
      x      = result$x[1:N] - result$x[(N+1):(2*N)],
      digits = 4
    )
    
    #Number of nonzero weights
    K <- sum(weights != 0)
    #Set surviving weights to 1/K
    NonZero <- which(weights !=0) 
    weights[NonZero] <- 1/K
  }
  
  
  # Objective Function
  obj <- as.numeric(t(weights) %*% Cov[1:N,1:N] %*% weights)
  
  return(list(
    Weights            = weights,
    LassoParameter     = lambda,
    ObjectiveFunction  = obj,
    SumWeights         = sum(weights),
    SumAbsoluteWeights = sum(abs(weights))
  ))
}

#--------------------------------------------------------------
#----Function for Lasso Portfolio with fixed Lambda----
#--------------------------------------------------------------
#Wheight Computation for whole Dataset
#Define Constraints
FixedLasso <- function(lambda = 1){
  
  #Pre allocate weight vector
  Weights <- vector(mode = "list", length = ( length(splitPF)-13 ))
  
  #Create dataframe for evaluation metrics
  Evaluation <- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("Validation_Risk", "Test_Risk", "Sharpe"))))
  
  #Loop over full dataset
  for (i in 1:227) {
    #Define covariance matrix estimation time frame
    covMat_end <- i+11
    #Validation Set
    val_month <- i+12
    #Test set
    test_month <- i+13
    
    #Sample Covariance Matrix
    cov <-  cov(do.call("rbind", splitPF[i:covMat_end]))
    
    #Fit model
    fit <- ConstraintLasso(Cov = cov, lambda = lambda)
    
    #Assign evaluation metrics
    Validation_Risk <- t(fit$Weights) %*% cov(splitPF[[val_month]]) %*% fit$Weights
    
    Test_Risk <- t(fit$Weights) %*% cov(splitPF[[test_month]]) %*% fit$Weights
    
    sharpe <- mean( rowSums(splitPF[[test_month]] * fit$Weights) - splitRF[[i]][1,2]) / sd( rowSums(splitPF[[test_month]] * fit$Weights)) * sqrt(252)
    
    Weights[[i]] <- fit$Weights
    
    returns[[i]] <- t(t(splitPF[[test_month]]) * fit$Weights)
    
    #Add row of i'th portfolio evaluation
    Evaluation[nrow(Evaluation) + 1,] = c(Validation_Risk, Test_Risk, sharpe)
    
  }
  return(list(
    Evaluation = Evaluation,
    Weights = Weights,
    Returns = returns
  ))
}

#--------------------------------------------------------------
#----Function for Lasso cross-validation----
#--------------------------------------------------------------
lasso_cv <- function(Input, Target) {
  
  #Transform passed Input to covariance
  cov <- cov(do.call("rbind", Input))
  
  #Find maximum lambda value
  lambda_max <- ConstraintLasso(Cov = cov, lambda = 100)$SumAbsoluteWeights
  
  #Compute grid of lambda values 
  grid = lseq(from = 1, to = lambda_max, length.out = 100)
  
  #Pre-allocate cross-validation metrics and weight vector
  CrossVal <<- data.frame(matrix(ncol=2,nrow=100, dimnames=list(NULL, c("Lambda", "Risk"))))
  weights <- vector(mode = "list", length = 100)
  
  #Fit model for all lambda values in grid
  for (i in seq_along(grid)) {
    
    #Vector of weights for every lambda value
    weights[[i]] <- ConstraintLasso(Cov = cov, lambda = grid[i])$Weights
    
    #Create data frame with values
    #Lambda
    CrossVal$Lambda[i] <<- grid[i]
    #Risk on Validation Set
    CrossVal$Risk[i] <<- t(weights[[i]]) %*% cov(Target) %*% weights[[i]]
  } 
  cv_result <<- list(Values = CrossVal, LambdaOpt = CrossVal[which.min(CrossVal$Risk),], WeightsOpt = weights[[which.min(CrossVal$Risk)]])
}


#--------------------------------------------------------------
#----Function for ridge cross-validation----
#--------------------------------------------------------------
ridge_cv <- function(Input, Target, CovInput = F) {
  
  #If/else statement distinguishes between input of covariance matrix vs. input of returns, due to peLasso (Ridge) implementation
  #Input as returns
  if(CovInput == F){
    #Transform returns in Input to covariance matrix
    cov <- cov(do.call("rbind", Input))
    
    #Find maximum lambda value
    lambda_max <- ConstraintRidge(Cov = cov, lambda = 100)$SumSquaredWeights
    
    #Compute grid of lambda values 
    grid = lseq(from =  1.1/ncol(cov), to = lambda_max, length.out = 100)
    
    #Pre-allocate cross-validation metrics and weight vector
    CrossVal <<- data.frame(matrix(ncol=2,nrow=100, dimnames=list(NULL, c("Lambda", "Risk"))))
    weights <- vector(mode = "list", length = 100)
    
    #Fit model for all values in lambda grid
    for (i in seq_along(grid)) {
      
      #Vector of weights for every lambda value
      weights[[i]] <- ConstraintRidge(Cov = cov, lambda = grid[i])$Weights
      
      #Create dataframe with evaluation for all lambdas
      #Lambda
      CrossVal$Lambda[i] <<- grid[i]
      #Risk on Validation Set
      CrossVal$Risk[i] <<- t(weights[[i]]) %*% cov(Target) %*% weights[[i]]
    } 
    #Store Lambdas, Risks, and optimal weights
    cv_result <<- list(Values = CrossVal, LambdaOpt = CrossVal[which.min(CrossVal$Risk),], WeightsOpt = weights[[which.min(CrossVal$Risk)]])
  
    
    #Computation of optimal lambda for peLasso (Ridge)  
  } else {
    #Declare input as covariance matrix for implementation of second step in peLasso
    cov <- Input
    
    #Find lambda max
    lambda_max <- ConstraintRidge(Cov = cov, lambda = 100)$SumSquaredWeights
    
    #Compute grid of lambda values 
    grid = lseq(from =  1.1/ncol(cov), to = lambda_max, length.out = 100)

    #Pre allocate cross-validation metrics and weight vector    
    CrossVal <<- data.frame(matrix(ncol=2,nrow=100, dimnames=list(NULL, c("Lambda", "Risk"))))
    weights <- vector(mode = "list", length = 100)
    
    #Fit model for all values in lambda grid
    for (i in seq_along(grid)) {
      
      #Vector of weights for every lambda value
      weights[[i]] <- ConstraintRidge(Cov = cov, lambda = grid[i])$Weights
      
      #Create dataframe with values
      #Lambda
      CrossVal$Lambda[i] <<- grid[i]
      #Risk on Validation Set
      CrossVal$Risk[i] <<- t(weights[[i]]) %*% cov(Target) %*% weights[[i]]
    } 
    
    
    cv_result <<- list(Values = CrossVal,
                       LambdaOpt = CrossVal[which.min(CrossVal$Risk),],
                       WeightsOpt = weights[[which.min(CrossVal$Risk)]]
    )
  }
}

#--------------------------------------------------------------
#----Function for partially egalitarian lasso (Average) cross-validation----
#--------------------------------------------------------------
peLassoAV_cv <- function(Input, Target) {
  #Transform returns in Input to covariance matrix
  cov <- cov(do.call("rbind", Input))
  
  #Find lambda max
  lambda_max <- ConstraintLasso(Cov = cov, lambda = 100)$SumAbsoluteWeights
  
  #Compute grid of lambda values 
  grid = lseq(from = 1, to = lambda_max, length.out = 100)
  
  CrossVal <<- data.frame(matrix(ncol=2,nrow=100, dimnames=list(NULL, c("Lambda", "Risk"))))
  weights <- vector(mode = "list", length = 100)
  
  #Find Solutions for all lambda values
  for (i in seq_along(grid)) {
    
    weights[[i]] <- peLasso_Average(Cov = cov, lambda = grid[i])$Weights
    
    #Create dataframe with values
    #Lambda
    CrossVal$Lambda[i] <<- grid[i]
    #Risk on Validation Set
    CrossVal$Risk[i] <<- t(weights[[i]]) %*% cov(Target) %*% weights[[i]]
  } 
  
  cv_result <<- list(Values = CrossVal, LambdaOpt = CrossVal[which.min(CrossVal$Risk),], WeightsOpt = weights[[which.min(CrossVal$Risk)]])
}