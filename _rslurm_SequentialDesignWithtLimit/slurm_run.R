library(base, quietly = TRUE)
library(methods, quietly = TRUE)
library(datasets, quietly = TRUE)
library(utils, quietly = TRUE)
library(grDevices, quietly = TRUE)
library(graphics, quietly = TRUE)
library(stats, quietly = TRUE)
library(usethis, quietly = TRUE)
library(devtools, quietly = TRUE)
library(rslurm, quietly = TRUE)
library(coda, quietly = TRUE)
library(Matrix, quietly = TRUE)
library(BayesFactor, quietly = TRUE)
library(assortedRFunctions, quietly = TRUE)


.rslurm_func <- readRDS('f.RDS')
.rslurm_params <- readRDS('params.RDS')
.rslurm_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
.rslurm_istart <- .rslurm_id * 20000 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 20000, nrow(.rslurm_params))
.rslurm_result <- do.call(parallel::mcmapply, c(
    FUN = .rslurm_func,
    .rslurm_params[.rslurm_istart:.rslurm_iend, , drop = FALSE],
    mc.cores = 6,
    mc.preschedule = TRUE,
    SIMPLIFY = FALSE))

saveRDS(.rslurm_result, file = paste0('results_', .rslurm_id, '.RDS'))
