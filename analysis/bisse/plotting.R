#####################################
##            Chapter 1            ##
##   BiSSE + FBD Simulation Study  ##
##            Plotting             ##
##    Bruno do Rosario Petrucci    ##
#####################################

### 
## packages

# ggplot2
library(ggplot2)

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
fbd_refs$comb <- as.numeric(rownames(fbd_refs))

# add true values to the data frame
fbd_refs$lambda1 <- 0.1
fbd_refs$lambda2 <- c(0.1, 0.2)[(fbd_refs$parComb == 2) + 1]
fbd_refs$mu1 <- 0.03
fbd_refs$mu2 <- c(0.03, 0.06)[(fbd_refs$parComb == 3) + 1]
fbd_refs$psi <- c(0.01, 0.05, 0.1)[fbd_refs$psiComb]

# make data frames for 95% CI and median, and one for mean
fbd_low <- fbd_med <- fbd_high <- fbd_mean <- 
  data.frame(matrix(nrow = 0, ncol = 4))
fbd_cov <- data.frame(matrix(nrow = 0, ncol = 5))

# iterate through logs
for (i in 1:nrow(fbd_refs)) {
  # get ref
  ref <- fbd_refs$comb[i]
  
  # and reps
  for (j in 1:n_reps) {
    # get log
    log <- logs[[ref]][[j]]
    
    # apply burnout
    log <- log[(nrow(log)/5):nrow(log), ]
    
    # get quantiles
    quants <- do.call(rbind.data.frame, lapply(c(0.025, 0.5, 0.975), 
                                                    function(x)
      unlist(lapply(1:ncol(log), function(y) quantile(log[, y], x)))))
    quants <- cbind(rep(ref, 3), quants)
    colnames(quants) <- c("comb", "lambda", "mu", "psi")
    
    # fill data frames
    fbd_low <- rbind(fbd_low, quants[1, ])
    fbd_med <- rbind(fbd_med, quants[2, ])
    fbd_high <- rbind(fbd_high, quants[3, ])
    fbd_mean <- rbind(fbd_mean, c(ref, colMeans(log)))
  }
  
  # low and high dfs for this combination
  fbd_low_i <- fbd_low[fbd_low[, 1] == ref, ]
  fbd_high_i <- fbd_high[fbd_high[, 1] == ref, ]
  
  # coverage
  fbd_cov <- rbind(fbd_cov, 
                   c(sum(fbd_low_i[, 2] < fbd_refs$lambda1[i] & 
                     fbd_high_i[, 2] > fbd_refs$lambda1[i]),
                   sum(fbd_low_i[, 2] < fbd_refs$lambda2[i] &
                     fbd_high_i[, 2] > fbd_refs$lambda2[i]),
                   sum(fbd_low_i[, 3] < fbd_refs$mu1[i] &
                     fbd_high_i[, 3] > fbd_refs$mu1[i]),
                   sum(fbd_low_i[, 3] < fbd_refs$mu2[i] &
                     fbd_high_i[, 3] > fbd_refs$mu2[i]),
                   sum(fbd_low_i[, 4] < fbd_refs$psi[i] &
                     fbd_high_i[, 4] > fbd_refs$psi[i])))
  colnames(fbd_cov) <- c("lambda1", "lambda2", "mu1", "mu2", "psi")
}

# normalize and name fbd_cov
fbd_cov <- fbd_cov / n_reps
rownames(fbd_cov) <- rownames(fbd_refs)

# name fbd_mean
colnames(fbd_mean) <- colnames(fbd_low)

# mean boxplot
fbd_mean_bplot_lambda <- ggplot(fbd_mean, aes(factor(comb), lambda)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = fbd_refs, aes(x = factor(comb), y = lambda1), 
             col = "red", size = 2) + 
  geom_point(data = fbd_refs, aes(x = factor(comb), y = lambda2), 
             col = "red", size = 2) + 
  labs(x = "Parameter combination") +
  theme_bw()
fbd_mean_bplot_mu <- ggplot(fbd_mean, aes(factor(comb), mu)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = fbd_refs, aes(x = factor(comb), y = mu1), 
             col = "red", size = 2) + 
  geom_point(data = fbd_refs, aes(x = factor(comb), y = mu2), 
             col = "red", size = 2) + 
  ylim(0, 0.1) +
  labs(x = "Parameter combination") +
  theme_bw()
fbd_mean_bplot_psi <- 
  ggplot(fbd_mean, aes(reorder(factor(comb), psi, FUN = mean), psi)) + 
  geom_boxplot(outlier.shape = NA) + geom_point(data = fbd_refs,
                              aes(x = reorder(factor(comb), psi), y = psi),
                              col = "red", size = 2) + 
  labs(x = "Parameter combination") + 
  theme_bw()

###
## testing accuracy - BiSSE (real trait only)

# get log indices for the BiSSE reps with real traits
bisse_refs <- refs_df[which(refs_df[, 3] == 2 & refs_df[, 4] == 1), ]
bisse_refs$comb <- as.numeric(rownames(bisse_refs))

# add true values to the data frame
bisse_refs$lambda1 <- 0.1
bisse_refs$lambda2 <- c(0.1, 0.2)[(bisse_refs$parComb == 2) + 1]
bisse_refs$mu1 <- 0.03
bisse_refs$mu2 <- c(0.03, 0.06)[(bisse_refs$parComb == 3) + 1]
bisse_refs$q01 <- 0.01
bisse_refs$q10 <- c(0.01, 0.005)[(bisse_refs$parComb == 4) + 1]

# make data frames for 95% CI and median, one for mean, and one for coverage
bisse_low <- bisse_med <- bisse_high <- bisse_mean <- bisse_cov <-
  data.frame(matrix(nrow = 0, ncol = 6))

# iterate through logs
for (i in 1:nrow(bisse_refs)) {
  # get ref
  ref <- bisse_refs$comb[i]
  
  # and reps
  for (j in 1:n_reps) {
    # get log
    log <- logs[[ref]][[j]]
    
    # apply burnout and take out pi
    log <- log[(nrow(log)/5):nrow(log), -c(5, 6)]
    
    # get quantiles
    quants <- do.call(rbind.data.frame, lapply(c(0.025, 0.5, 0.975), 
                                               function(x)
      unlist(lapply(1:ncol(log), function(y) quantile(log[, y], x)))))
    quants <- cbind(rep(ref, 3), quants)
    colnames(quants) <- c("comb", "lambda1", "lambda2", 
                          "mu1", "mu2", "q01", "q10")
    
    # fill data frames
    bisse_low <- rbind(bisse_low, quants[1, ])
    bisse_med <- rbind(bisse_med, quants[2, ])
    bisse_high <- rbind(bisse_high, quants[3, ])
    bisse_mean <- rbind(bisse_mean, c(ref, colMeans(log)))
  }
  
  # low and high dfs for this combination
  bisse_low_i <- bisse_low[bisse_low[, 1] == ref, ]
  bisse_high_i <- bisse_high[bisse_high[, 1] == ref, ]
  
  # coverage
  bisse_cov <- rbind(bisse_cov, 
                   c(sum(bisse_low_i[, 2] < bisse_refs$lambda1[i] & 
                           bisse_high_i[, 2] > bisse_refs$lambda1[i]),
                     sum(bisse_low_i[, 3] < bisse_refs$lambda2[i] &
                           bisse_high_i[, 3] > bisse_refs$lambda2[i]),
                     sum(bisse_low_i[, 4] < bisse_refs$mu1[i] &
                           bisse_high_i[, 4] > bisse_refs$mu1[i]),
                     sum(bisse_low_i[, 5] < bisse_refs$mu2[i] &
                           bisse_high_i[, 5] > bisse_refs$mu2[i]),
                     sum(bisse_low_i[, 6] < bisse_refs$q01[i] &
                           bisse_high_i[, 6] > bisse_refs$q01[i]),
                     sum(bisse_low_i[, 7] < bisse_refs$q10[i] &
                           bisse_high_i[, 7] > bisse_refs$q10[i])))
  colnames(bisse_cov) <- c("lambda1", "lambda2", "mu1", "mu2", "q01", "q10")
}

# normalize and name bisse_cov
bisse_cov <- bisse_cov / n_reps
rownames(bisse_cov) <- rownames(bisse_refs)

# name bisse_mean
colnames(bisse_mean) <- colnames(bisse_low)

# mean boxplot
bisse_mean_bplot_lambda1 <- ggplot(bisse_mean, aes(factor(comb), lambda1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = bisse_refs, aes(x = factor(comb), y = lambda1), 
             col = "red", size = 2) +
  labs(x = "Parameter combination") +
  theme_bw()
bisse_mean_bplot_lambda2 <- ggplot(bisse_mean, aes(factor(comb), lambda2)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = bisse_refs, aes(x = factor(comb), y = lambda2), 
             col = "red", size = 2) +
  labs(x = "Parameter combination") +
  theme_bw()
bisse_mean_bplot_mu1 <- ggplot(bisse_mean, aes(factor(comb), mu1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = bisse_refs, aes(x = factor(comb), y = mu1), 
             col = "red", size = 2) + 
  ylim(0, 0.05) +
  labs(x = "Parameter combination") +
  theme_bw()
bisse_mean_bplot_mu2 <- ggplot(bisse_mean, aes(factor(comb), mu2)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = bisse_refs, aes(x = factor(comb), y = mu2), 
             col = "red", size = 2) + 
  ylim(0, 0.061) +
  labs(x = "Parameter combination") +
  theme_bw()
bisse_mean_bplot_q01 <- ggplot(bisse_mean, aes(factor(comb), q01)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = bisse_refs, aes(x = factor(comb), y = q01), 
             col = "red", size = 2) + 
  labs(x = "Parameter combination") +
  theme_bw()
bisse_mean_bplot_q10 <- ggplot(bisse_mean, aes(factor(comb), q10)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = bisse_refs, aes(x = factor(comb), y = q10), 
             col = "red", size = 2) + 
  labs(x = "Parameter combination") +
  theme_bw()

# plots equivalent to Maddison 2007 - means
bisse_comp_lambda <- ggplot(bisse_mean, aes(lambda1, lambda2, 
                                            color = factor(comb == 11))) +
  geom_point() +
  theme_bw()
bisse_comp_mu <- ggplot(bisse_mean, aes(mu1, mu2, 
                                            color = factor(comb == 15))) +
  geom_point() +
  theme_bw()
bisse_comp_q <- ggplot(bisse_mean, aes(q01, q10, 
                                            color = factor(comb == 19))) +
  geom_point() +
  theme_bw()

###
## testing accuracy - BiSSE + FBD