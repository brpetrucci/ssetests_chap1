#####################################
##            Chapter 1            ##
##   BiSSE + FBD Simulation Study  ##
##            Plotting             ##
##    Bruno do Rosario Petrucci    ##
#####################################

###
## reading logs

# create list of logs
logs <- list()

# create a list of ESS values
ess <- list()

# and create a dataframe to hold comb values for each ref
refs_df <- data.frame(matrix(nrow = 0, ncol = 4))

# number of parameter combinations
n_refs <- 76

# number of reps
n_reps <- 100

# vectors of combination names
psiRefs <- c("1_low_psi", "2_med_psi", "3_high_psi")
parRefs <- c("1_base", "2_high_lambda1", "3_high_mu0", "4_low_q10")
modelRefs <- c("fbd", "bisse", "both")
traitRefs <- c("real", "low", "mid", "high")

# base directory
base_dir <- "/Users/petrucci/Documents/research/ssetests_chap1/"

# for each combination
for (i in 1:n_refs) {
  # source the ref in question
  source(paste0(base_dir, "analysis/bisse/refs/refs_", i, ".Rev"))
  
  # if model is 1 (fbd), trait is not necessary, same for 2 (bisse) and psi
  if (modelComb == 1) {
    traitComb <- 1
  } else if (modelComb == 2) {
    psiComb <- 1
  }
  
  # psi reference
  psiRef <- psiRefs[psiComb]
  
  # parameter reference
  parRef <- parRefs[parComb]
  
  # model reference
  modelRef <- modelRefs[modelComb]
  
  # trait reference
  traitRef <- traitRefs[traitComb]
  
  # append to refs_df
  refs_df <- rbind(refs_df, c(psiComb, parComb, modelComb, traitComb))
  
  # if model is 1, no traitRef field is necessary
  if (modelComb == 1) {
    # get the directory where the reps are
    target_dir <- paste0(base_dir, "output/bisse/", psiRef, "/", parRef, "/",
                         modelRef, "/")
  } else {
    # get the directory where the reps are
    target_dir <- paste0(base_dir, "output/bisse/", psiRef, "/", parRef, "/",
                         modelRef, "/", traitRef, "/")
  }
  
  # get each rep
  reps <- lapply(1:n_reps, function(x) {
    table <- read.delim(paste0(target_dir, "rep_", x, "/bisse.log"), 
                        header = TRUE, sep = "\t");
    table[, 5:ncol(table)]
    })
  
  # get the ess vales
  ess_reps <- lapply(1:n_reps, function(x) effectiveSize(reps[[x]]))
  col_names <- names(ess_reps[[1]])
  ess_reps <- do.call(rbind.data.frame, ess_reps)
  colnames(ess_reps) <- col_names
  
  # put it in the list of logs
  logs[[i]] <- reps
  ess[[i]] <- ess_reps
}

# name refs
colnames(refs_df) <- c("psiComb", "parComb", "modelComb", "traitComb")

###
## testing accuracy - FBD

# get log indices for the FBD reps
fbd_refs <- refs_df[which(refs_df[, 3] == 1), ]

# add true values to the data frame
fbd_refs$lambda <- 0.1
fbd_refs$lambda2 <- c(0.1, 0.2)[(fbd_refs$parComb == 2) + 1]
fbd_refs$mu <- 0.03
fbd_refs$mu2 <- c(0.03, 0.06)[(fbd_refs$parComb == 3) + 1]
fbd_refs$psi <- c(0.01, 0.05, 0.1)[fbd_refs$psiComb]

# make data frames for 95% CI and median, and one for mean
fbd_low <- fbd_med <- fbd_high <- fbd_mean <- 
  data.frame(matrix(nrow = 0, ncol = 4))

# iterate through logs
for (i in as.numeric(rownames(fbd_refs))) {
  # and reps
  for (j in 1:n_reps) {
    # get log
    log <- logs[[i]][[j]]
    
    # get quantiles
    quants <- do.call(rbind.data.frame, lapply(c(0.025, 0.5, 0.975), 
                                                    function(x)
      unlist(lapply(1:ncol(log), function(y) quantile(log[, y], x)))))
    quants <- cbind(rep(i, 3), quants)
    colnames(quants) <- c("comb", "lambda", "mu", "psi")
    
    # fill data frames
    fbd_low <- rbind(fbd_low, quants[1, ])
    fbd_med <- rbind(fbd_med, quants[2, ])
    fbd_high <- rbind(fbd_high, quants[3, ])
    fbd_mean <- rbind(fbd_mean, c(i, colMeans(log)))
  }
}

# name fbd_mean
colnames(fbd_mean) <- colnames(fbd_low)
