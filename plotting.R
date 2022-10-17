### 
# get packages
library(ggplot2)
library(bayestestR)

### 
# function to read logs
read_logs <- function(refs, nReps, base_dir, mod_dir = NULL, cut) {
  # create data frame list
  fbd_logs <- list()
  
  # for each ref
  for (ref in refs) {
    # for each rep
    for (rep in 1:nReps) {
      # ref dir
      ref_dir <- paste0(base_dir, "/", ref)
      
      # get rep directory
      rep_dir <- paste0(ref_dir, "/rep_", rep)
      
      # joint log file
      joint_log <- if (is.null(mod_dir)) paste0(rep_dir, "/ssefbd.log") else
        paste0(rep_dir, "/", mod_dir, "/ssefbd.log")
      
      # read log
      fbd_logs[[rep]] <- read.delim(joint_log)[, -1:-cut]
    }
  }
  
  # return logs
  return(fbd_logs)
}

###
# read logs

##
# pure FBD

# refs - only first one
refs <- "1_base"

# number of reps ran for FBD
nReps <- 500

# directory where all output for a given rep is 
base_dir <- "/Users/petrucci/Documents/research/ssefbd_evol22/power/analysis/output/2_med_psi"

# mod-fbd
mod_dir <- "fbd"

# get logs (cut 6 - organizational cols, likelihoods, and div)
fbd_logs <- read_logs(refs, nReps, base_dir, mod_dir, 6)


### 
# making plots

# set up expected values
expec <- data.frame(lambda = 0.1,
                    mu = 0.03,
                    psi = 0.05)

# list to hold data frame of means and 95% CI
means <- list()

# one to hold the mode
modes <- list()

# and a data frame to hold the coverage
coverages <- list()

# get col names from each log
cnames <- colnames(fbd_logs[[1]])

# set coverage colnames to that
colnames(coverages) <- cnames

# start our new colnames
cnamesMeans <- c()

# for each of the previous colnames
for (cname in cnames) {
  # get a low95, mean, and up95
  cnamesMeans <- c(cnamesMeans, paste0(cname, "low95"), paste0(cname, "mean"),
                   paste0(cname, "high95"))
}

# for each ref
for (ref in refs) {
  # initialize data frames
  means[[ref]] <- data.frame(matrix(ncol = 9, nrow = 0))

  modes[[ref]] <- coverages[[ref]] <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(modes[[ref]]) <- colnames(coverages[[ref]]) <- cnames
  
  # for each log
  for (i in 1:length(fbd_logs)) {
    # initialize vectors of results
    resMeans <- c()
    resModes <- c()
    
    # get this log
    log <- fbd_logs[[i]]
    
    # for each column in log
    for (col in 1:ncol(log)) {
      # 97.5% and 2.5% quantiles
      quants <- quantile(log[, col], c(0.025, 0.975))
      
      # append to res with mean
      resMeans <- c(resMeans, quants[1], mean(log[, col]), quants[2])
      
      # get a mode estimate
      mode <- map_estimate(log[, col])
      
      # append to resModes
      resModes <- c(resModes, mode)
    }
    
    # add res to data frame
    means[[ref]] <- rbind(means[[ref]], resMeans)
    modes[[ref]] <- rbind(modes[[ref]], resModes)
    
    # auxiliar dfs
    low <- means[[ref]][, c(1, 4, 7)]
    high <- means[[ref]][, c(3, 6, 9)]
    
    # calculate coverage
    coverages[[ref]] <- as.data.frame(lapply(1:ncol(log), function(x) 
      sum(low[, x] < expec[x] & high[, x] > expec[x])))
    colnames(coverages[[ref]]) <- colnames(modes[[ref]]) <- cnames
  }
  
  # colnames
  colnames(means[[ref]]) <- cnamesMeans
}

# set up 6 plots
par(mfrow = c(3, 1))

# for each ref
for (comb in 1:length(refs)) {
  # get ref
  ref <- refs[comb]
  
  # get means data frame
  meansRef <- means[[ref]]
  
  # for each parameter
  for (i in 1:ncol(fbd_logs[[1]])) {
    # get parameter name
    name <- colnames(fbd_logs[[1]])[i]
    
    # get vectors of means and 95 CI
    low95 <- meansRef[, 3*(i - 1) + 1]
    mean <- meansRef[, 3*(i - 1) + 2]
    high95 <- meansRef[, 3*i]
    
    # plot the means
    plot(mean, ylab = "rate (events/lineage/my)", xlab = "rep",
         main = paste0("Mean and 95% CI of ", name, " for ", ref),
         ylim = c(0, max(high95) + max(high95)/5))
    
    # add lines for 95% CI
    lines(low95, col = "RED")
    lines(high95, col = "RED")
    
    # add expected line
    abline(h = expec[comb, i], col = "green")
  }
}

# get div
modes[[1]]$div <- modes[[1]]$lambda - modes[[1]]$mu

# ggplot for the modes
ggplot(modes[[1]], aes(x = mu, y = div)) +
  geom_point() + geom_hline(yintercept = 0.07) + geom_vline(xintercept = 0.03)
