# This script runs the simulation for the talk Bayesian sequential designs Wed 29 July 2020
# 
# /* 
# ----------------------------- General stuff ---------------------------
# */
# Setting seed
set.seed(911225)

# Set library path
.libPaths('/home/aq01/R/x86_64-pc-linux-gnu-library/3.5')

# Libraries
library(rslurm)
library(BayesFactor)
library(assortedRFunctions)

# Job parameters
n_nodes       <- 1
cpus_per_node <- 6
nIter         <- 10000

# /* 
# ----------------------------- Traditional design ---------------------------
# */
# Function
helperfunction <- function(n, d, side){
  results   <- list()
  
  data      <- rnorm(n, d, 1)
  if(side == 2){
    # Two sided
    bf        <- reportBF(ttestBF(data), 4)
  } else {
    # One sided
    bf        <- reportBF(ttestBF(data, nullInterval = c(-Inf, 0))[2], 4)
  }
  
  results$side <- side
  results$d    <- d
  results$n    <- n
  results$bf   <- bf
  return(results)
}

# Setting parameters 
n           <- seq(2, 1000, 10)
d0          <- 0.0
d1          <- 0.5

params      <- data.frame(n = c(rep(n, nIter), rep(n, nIter)),
                          d = c(rep(d0, nIter*length(n)), rep(d1, nIter*length(n))))

params      <- rbind(params, params)
params$side <- c(rep(2, nrow(params)/2), rep(1, nrow(params)/2))

# Submittig jobs
sjob1 <- slurm_apply(helperfunction, params, jobname = 'traditionalDesign',
                     nodes = n_nodes, cpus_per_node = cpus_per_node, submit = FALSE)

# /* 
# ----------------------------- Sequential design without limit ---------------------------
# */
# Function
helperfunction <- function(minN, d, crit1, crit2, batchSize, side){
  bf        <- c()
  results   <- list()
  
  # Create minium sample and calculate BF
  i    <- 1
  n    <- as.numeric(minN)
  data <- rnorm(n, d, 1)
  if(side == 2){
    # Two sided
    bf        <- reportBF(ttestBF(data), 4)
  } else {
    # One sided
    bf        <- reportBF(ttestBF(data, nullInterval = c(-Inf, 0))[2], 4)
  }
  
  
  # Within simulation loop
  while(bf[length(bf)] < crit1 & bf[length(bf)] > crit2){
    n         <- n + batchSize
    data      <- c(data, rnorm(batchSize, d, 1))
    if(side == 2){
      # Two sided
      bf[i + 1] <- reportBF(ttestBF(data), 4)
    } else {
      # One sided
      bf[i + 1] <- reportBF(ttestBF(data, nullInterval = c(-Inf, 0))[2], 4)
    }
    
    i         <- i + 1
  }
  
  
  # Return results
  results$side      <- side
  results$d         <- d
  results$n         <- n
  results$bf        <- bf
  results$crit1     <- crit1
  results$crit2     <- crit2
  results$batchSize <- batchSize
  return(results)
}

# Parameters
params      <- data.frame(minN      = rep(10, nIter*2),
                          d         = c(rep(d0, nIter), rep(d1, nIter)),
                          crit1     = rep(10, nIter*2),
                          crit2     = rep(1/6, nIter*2),
                          batchSize = rep(5, nIter*2))

params      <- rbind(params, params)
params$side <- c(rep(2, nrow(params)/2), rep(1, nrow(params)/2))

# Create job
sjob1 <- slurm_apply(helperfunction, params, jobname = 'SequentialDesignWithoutLimit',
                     nodes = n_nodes, cpus_per_node = cpus_per_node, submit = FALSE)

# /* 
# ----------------------------- Sequential design with limit ---------------------------
# */
# Function
helperfunction <- function(minN, d, crit1, crit2, batchSize, limit){
  bf        <- c()
  results   <- list()
  
  # Create minium sample and calculate BF
  i    <- 1
  n    <- as.numeric(minN)
  data <- rnorm(n, d, 1)
  bf   <- reportBF(ttestBF(data), 4)
  
  
  # Within simulation loop
  while(bf[length(bf)] < crit1 & bf[length(bf)] > crit2 & n < limit){
    n         <- n + batchSize
    data      <- c(data, rnorm(batchSize, d, 1))
    bf[i + 1] <- reportBF(ttestBF(data), 4)
    i         <- i + 1
  }
  
  
  # Return results
  results$d         <- d
  results$n         <- n
  results$bf        <- bf
  results$crit1     <- crit1
  results$crit2     <- crit2
  results$batchSize <- batchSize
  results$limit     <- limit
  return(results)
}

#minN, d, crit1, crit2, batchSize, side
# Parameters
params      <- data.frame(minN      = rep(10, nIter*2),
                          d         = c(rep(d0, nIter), rep(d1, nIter)),
                          crit1     = rep(10, nIter*2),
                          crit2     = rep(1/6, nIter*2),
                          batchSize = rep(5, nIter*2),
                          limit     = rep(100, nIter*2))

# Create job
sjob1 <- slurm_apply(helperfunction, params, jobname = 'SequentialDesignWithtLimit',
                     nodes = n_nodes, cpus_per_node = cpus_per_node, submit = FALSE)

# /* 
# ----------------------------- Sequential design with limit and with low minimum size ---------------------------
# */
# Parameters
params      <- data.frame(minN      = rep(2, nIter*2),
                          d         = c(rep(d0, nIter), rep(d1, nIter)),
                          crit1     = rep(10, nIter*2),
                          crit2     = rep(1/6, nIter*2),
                          batchSize = rep(5, nIter*2),
                          limit     = rep(100, nIter*2))

# Create job
sjob1 <- slurm_apply(helperfunction, params, jobname = 'SequentialDesignWithtLimit_lowMinN',
                     nodes = n_nodes, cpus_per_node = cpus_per_node, submit = FALSE)


# /* 
# ----------------------------- Sequential design with limit and checking every time ---------------------------
# */
# Parameters
params      <- data.frame(minN      = rep(10, nIter*2),
                          d         = c(rep(d0, nIter), rep(d1, nIter)),
                          crit1     = rep(10, nIter*2),
                          crit2     = rep(1/6, nIter*2),
                          batchSize = rep(1, nIter*2),
                          limit     = rep(100, nIter*2))

# Create job
sjob1 <- slurm_apply(helperfunction, params, jobname = 'SequentialDesignWithtLimit_batcheSize1',
                     nodes = n_nodes, cpus_per_node = cpus_per_node, submit = FALSE)