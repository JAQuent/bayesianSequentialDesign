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
  bf_1        <- reportBF(ttestBF(data, nullInterval = c(-Inf, 0))[2], 4)
  bf_2        <- reportBF(ttestBF(data), 4)
  
  results <- data.frame(d = d,
                        n = n,
                        bf_1 = bf_1,
                        bf_2 = bf_2)

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

# Submitting jobs
sjob1 <- slurm_apply(helperfunction, params, jobname = 'traditionalDesign',
                     nodes = n_nodes, cpus_per_node = cpus_per_node, submit = FALSE)

# /* 
# ----------------------------- Sequential design without limit ---------------------------
# */
# Function
helperfunction <- function(minN, d, crit1, crit2, batchSize, side){
  # Create minimum sample and calculate BF
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
  
  # Create DF
  results <- data.frame(index = 1,
                        d = d, 
                        n = n,
                        bf = bf,
                        crit1 = crit1,
                        crit2 = crit2,
                        batchSize = batchSize)
  
  # Within simulation loop
  while(bf < crit1 & bf > crit2){
    n         <- n + batchSize
    data      <- c(data, rnorm(batchSize, d, 1))
    if(side == 2){
      # Two sided
      bf <- reportBF(ttestBF(data), 4)
    } else {
      # One sided
      bf <- reportBF(ttestBF(data, nullInterval = c(-Inf, 0))[2], 4)
    }
    
    i         <- i + 1
    # Bind to DF
    results <- rbind(results,  data.frame(index = i,
                                          d = d, 
                                          n = n,
                                          bf = bf,
                                          crit1 = crit1,
                                          crit2 = crit2,
                                          batchSize = batchSize))
  }
  
  # Return results
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
# Note this function can also be run with no limit by setting limit = Inf.
# Function
helperfunction <- function(minN, d, crit1, crit2, batchSize, limit){
  # Create minimum sample and calculate BF
  i    <- 1
  n    <- as.numeric(minN)
  data <- rnorm(n, d, 1)
  bf   <- reportBF(ttestBF(data), 4)
  
  # Create DF
  results <- data.frame(index = 1,
                        d = d, 
                        n = n,
                        bf = bf,
                        crit1 = crit1,
                        crit2 = crit2,
                        batchSize = batchSize,
                        limit = limit)
  
  
  # Within simulation loop
  while(bf < crit1 & bf> crit2 & n < limit){
    n         <- n + batchSize
    data      <- c(data, rnorm(batchSize, d, 1))
    bf        <- reportBF(ttestBF(data), 4)
    i         <- i + 1
    
    # Bind to DF
    results <- rbind(results,  data.frame(index = i,
                                          d = d, 
                                          n = n,
                                          bf = bf,
                                          crit1 = crit1,
                                          crit2 = crit2,
                                          batchSize = batchSize,
                                          limit = limit))
  }
  
  
  # Return results
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