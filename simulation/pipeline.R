#######################################
##            Chapter 1              ##
##   HiSSE + FBD Simulation Study    ##
##     Bruno do Rosario Petrucci     ##
#######################################

###
# goal: replicate a portion of Beaulieu & O'Meara 2016 
# HiSSE and Maddison et al 2007 BiSSE tests using
# trees including fossil samples as well (FBD)

### 
# packages required

# devtools 
suppressMessages(library(devtools))

# paleobuddy - need to remember to be in development branch
suppressMessages(library(paleobuddy))

# APE (for saving trees)
library(ape)

# readr (for write_tsv)
library(readr)

# load paleobuddy
load_all()

###
# simulation pipeline settings:
# Number of traits: 1 real, 3 neutral (q from 0.01 to 1)
# Trait evolution: Bi/HiSSE
# states: c(0A, 0B, 1A, 1B)
# X0: 0A
# Transition rates: qjk = q01 = 0.01, and q10 = 0.005
# mu: mujk = mu0A = 0.03; mu0B = 0.06
# lambda: lambdajk = lambda0A = 0.1; lambda1B = 0.2
# psi: 0.01, 0.05, 0.1
# N: 500

###
# set up the baseline simulation settings

# initial number of species
n0 <- 1

# expected number of extant species
N <- 500

# base (null model) rates
# equivalent to BiSSE pipeline
# 0 = 0A, 1 = 1A, 2 = 0B, 3 = 1B
lambda0 <- lambda1 <- lambda2 <- lambda3 <- 0.1

mu0 <- mu1 <- mu2 <- mu3 <- 0.03

q01 <- q10 <- q23 <- q32 <- q02 <- q20 <- q13 <- q31 <- 0.01

q03 <- q30 <- q12 <- q21 <- 0

psi <- 0.05

# create variations we want to test -- also from BiSSE
lambda1Var <- lambda3Var <- 0.2

mu0Var <- mu2Var <- 0.06

q10Var <- q32Var <- 0.005

lambda3Var2 <- 0.2

mu3Var2 <- 0.06

###
# create directories for saving data

# make a smarter dir.create function
smart.dir.create <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir)
}

###
# consider parameter changes and run simulations

# collect baseline parameters in a data
# frame that will later become the key

key <- data.frame(ref = "base",
                  lambda0 = lambda0,
                  lambda1 = lambda1,
                  lambda2 = lambda2,
                  lambda3 = lambda3,
                  mu0 = mu0,
                  mu1 = mu1, 
                  mu2 = mu2,
                  mu3 = mu3,
                  q01 = q01,
                  q10 = q10,
                  q02 = q02,
                  q20 = q20,
                  q03 = q03,
                  q30 = q30,
                  q12 = q12,
                  q21 = q21,
                  q13 = q13,
                  q31 = q31,
                  q23 = q23,
                  q32 = q32,
                  psi = psi)

# make a vector to hold the baseline parameters
base <- key[1, ]

# reference vector for each parameter change
refs <- c("high_lambda1", "high_mu0", "low_q10", "high_lambda1B", "high_mu0B")

# the corresponding column numbers to change in the base
change <- list(c(3, 5), c(6, 8), c(11, 21), 5, 8)

# the new values
new <- list(c(lambda1Var, lambda3Var),  c(mu0Var, mu2Var), c(q10Var, q32Var),
            lambda3Var2, mu3Var2)

# copy the baseline a bunch of times to fill the key
key <- rbind(key, key[rep(1, length(change)), ])
rownames(key) <- 1:nrow(key)

# populate it with the values in aux
for (i in 2:nrow(key)) {
  key$ref[i] <- refs[i - 1]
  key[i, change[[i - 1]]] <- new[[i - 1]]
}

###
# auxiliary functions

# correcting root edge after dropping some tips
correct_root_edge <- function(phy, sim, fossils) {
  if (max(node.depth.edgelength(phy)) + phy$root.edge != max(sim$TS)) {
    # if the MRCA is a speciation event
    if (max(node.depth.edgelength(phy)) %in% sim$TS) {
      # get species that was born
      spBorn <- which(sim$TS == max(node.depth.edgelength(phy)))
      
      # get its parent
      spPar <- sim$PAR[spBorn]
      
      # get fossils for parent species
      fossilsSpPar <- fossils[fossils$Species == paste0("t", spPar), ]
      
      # check if spPar has any fossils older than spBorn
      if (nrow(fossilsSpPar) > 0 && 
          any(fossilsSpPar$SampT > sim$TS[spBorn])) {
        # root edge will be the difference between 
        # the speciation time and oldest fossil
        phy$root.edge <- sim$TS[spPar] - max(fossilsSpPar$SampT)
      } else {
        # if not, root edge is just the time between speciations
        phy$root.edge <- sim$TS[spPar] - sim$TS[spBorn]
      }
    } else {
      # if not, it is a fossil, find which species'
      spFossil <- fossils$Species[which(fossils$SampT == 
                                          max(node.depth.edgelength(phy)))]
      spFossil <- as.numeric(gsub('.([0-9]+)*','\\1', spFossil))
      
      # root edge then is the difference between speciation and fossilization
      phy$root.edge <- sim$TS[spFossil] - 
        max(node.depth.edgelength(phy))
    }
  }
  
  return(phy)
}

# create simulation function for one rep
simulate_rep <- function(lambda0, lambda1, lambda2, lambda3,
                         mu0, mu1, mu2, mu3, 
                         q01, q10, q02, q20, q03, q30,
                         q12, q21, q13, q31, q23, q32, N, psi) {
  ### set up parameters
  # initial number of species
  n0 <- 1
  
  # speciation
  lambda <- c(lambda0, lambda1, lambda2, lambda3)
  
  # extinction
  mu <- c(mu0, mu1, mu2, mu3)
  
  # transition rate matrix
  Q <- list(matrix(c(0, q10, q20, q30, q01, 0, q21, q31, 
                     q02, q12, 0, q32, q03, q13, q23, 0), 4, 4),
            matrix(c(0, 0.01, 0.01, 0, 0.01, 0, 0, 0.01, 
                     0.01, 0, 0, 0.01, 0, 0.01, 0.01, 0), 4, 4),
            matrix(c(0, 0.1, 0.1, 0, 0.1, 0, 0, 0.1, 
                     0.1, 0, 0, 0.1, 0, 0.1, 0.1, 0), 4, 4),
            matrix(c(0, 1, 1, 0, 1, 0, 0, 1, 
                     1, 0, 0, 1, 0, 1, 1, 0), 4, 4))
  
  # initialize check
  bounds <- FALSE
  
  # counter to make sure we don't keep insisting on bad parameter values
  counter <- 1
  ratioSampled <- c()
  ratioTrait <- c()
  
  while (!bounds) {
    # run simulation
    sim <- bd.sim.traits(n0, lambda, mu, N = N, 
                         nTraits = 4, nStates = 2, nHidden = 2, Q = Q)
    
    # run again until we have some extinct species
    while (sum(!sim$SIM$EXTANT) == 0) {
      sim <- bd.sim.traits(n0, lambda, mu, N = N, 
                           nTraits = 4, nStates = 2, nHidden = 2, Q = Q)
      
    }
    
    # traits
    traits <- sim$TRAITS
    
    # and get the sim object
    sim <- sim$SIM
    
    # get tMax
    tMax <- max(sim$TS)
    
    # sample data frame
    sample <- data.frame()
    while (nrow(sample) == 0) {
      sample <- suppressMessages(sample.clade(sim, psi, max(sim$TS)))
    }
    
    # get complete tree
    tree <- make.phylo(sim)
    
    # create an ultrametric tree
    treeUltra <- drop.fossil(tree)
    
    # get an FBD tree
    treeFBD <- make.phylo(sim, sample, returnTrueExt = FALSE)
    
    # correct root edge if necessary
    treeFBD <- correct_root_edge(treeFBD, sim, sample)
    
    # get origin
    originFBD <- max(node.depth.edgelength(treeFBD)) + treeFBD$root.edge
    
    # make the trait lists
    # for the ultrametric tree
    traitsSummaryUltra <- traits.summary(sim, traits, selection = "extant")
    
    treeUltraReal <- traitsSummaryUltra$trait1
    
    treeUltraNeutralLow <- traitsSummaryUltra$trait2
    
    treeUltraNeutralMid <- traitsSummaryUltra$trait3
    
    treeUltraNeutralHigh <- traitsSummaryUltra$trait4
    
    # and for the fbd tree
    traitsSummaryFBD <- traits.summary(sim, traits, 
                                       sample, selection = "sampled")
    
    treeFBDReal <- traitsSummaryFBD$trait1
    
    treeFBDNeutralLow <- traitsSummaryUltra$trait2
    
    treeFBDNeutralMid <- traitsSummaryUltra$trait3
    
    treeFBDNeutralHigh <- traitsSummaryUltra$trait4
    
    # check the minimum percentage of the rare trait
    ratioTrait <- min(sum(treeUltraReal)/sum(sim$EXTANT),
                      sum(treeFBDReal)/(length(sim$TS) + 
                                          nrow(sample) - sum(!sim$EXTANT)))
    
    # check that both traits are represented at least 10%
    bounds <- (ratioTrait >= 0.1) &&
      (ratioTrait <= 0.9)
    
    # increase counter
    counter <- counter + 1
    
    # if counter is higher than 10, maybe rethink the parameters
    if (counter > 10) {
      print(ratioTrait)
      stop("Hard to find replicate within bounds")
    }
  }
  
  # make a list to return
  res <- list(SIM = sim, SAMPLE = sample,
              TREE = tree, ULTRATREE = treeUltra, FBDTREE = treeFBD,
              ORIGIN = originFBD,
              REALULTRA = treeUltraReal, REALFBD = treeFBDReal,
              NEUTRALLOWULTRA = treeUltraNeutralLow, 
              NEUTRALLOWFBD = treeFBDNeutralLow,
              NEUTRALMIDULTRA = treeUltraNeutralMid, 
              NEUTRALMIDFBD = treeFBDNeutralMid,
              NEUTRALHIGHULTRA = treeUltraNeutralHigh, 
              NEUTRALHIGHFBD = treeFBDNeutralHigh)
  
  return(res)
}

# function to save simulation lists
save_sims <- function(nReps, simReps, targetDir) {
  # extract simulations
  simList <- lapply(1:nReps, function(x) simReps[[x]]$SIM)
  
  # trait data 
  realUltra <- lapply(1:nReps, function(x) simReps[[x]]$REALULTRA)
  realFBD <- lapply(1:nReps, function(x) simReps[[x]]$REALFBD)
  neutralLowUltra <- lapply(1:nReps, function(x) simReps[[x]]$NEUTRALLOWULTRA)
  neutralLowFBD <- lapply(1:nReps, function(x) simReps[[x]]$NEUTRALLOWFBD)
  neutralMidUltra <- lapply(1:nReps, function(x) simReps[[x]]$NEUTRALMIDULTRA)
  neutralMidFBD <- lapply(1:nReps, function(x) simReps[[x]]$NEUTRALMIDFBD)
  neutralHighUltra <- lapply(1:nReps, function(x) simReps[[x]]$NEUTRALMIDULTRA)
  neutralHighFBD <- lapply(1:nReps, function(x) simReps[[x]]$NEUTRALMIDFBD)
  
  # fossil records
  fossilsList <- lapply(1:nReps, function(x) simReps[[x]]$SAMPLE)
  
  # trees
  treeList <-  lapply(1:nReps, function(x) simReps[[x]]$TREE)
  treeUltraList <- lapply(1:nReps, function(x) simReps[[x]]$ULTRATREE)
  treeFBDList <- lapply(1:nReps, function(x) simReps[[x]]$FBDTREE)
  
  # origin times
  origin <- lapply(1:nReps, function(x) simReps[[x]]$ORIGIN)
  
  # save simulation RData
  save(simList, file = paste0(targetDir, "sim_list.RData"))
  
  # simulations directory
  simsDir <- paste0(targetDir, "sims/")
  smart.dir.create(simsDir)
  
  # trait lists directories
  traitsDir <- paste0(targetDir, "traits/")
  smart.dir.create(traitsDir)
  
  traitsUltraDir <- paste0(traitsDir, "ultrametric/")
  smart.dir.create(traitsUltraDir)
  
  neutralUltraDir <- paste0(traitsUltraDir, "neutral/")
  smart.dir.create(neutralUltraDir)
  
  neutralUltraLowDir <- paste0(neutralUltraDir, "low/")
  smart.dir.create(neutralUltraLowDir)
  
  neutralUltraMidDir <- paste0(neutralUltraDir, "mid/")
  smart.dir.create(neutralUltraMidDir)
  
  neutralUltraHighDir <- paste0(neutralUltraDir, "high/")
  smart.dir.create(neutralUltraHighDir)
  
  traitsFBDDir <- paste0(traitsDir, "fbd/")
  smart.dir.create(traitsFBDDir)
  
  neutralFBDDir <- paste0(traitsFBDDir, "neutral/")
  smart.dir.create(neutralFBDDir)
  
  neutralFBDLowDir <- paste0(neutralFBDDir, "low/")
  smart.dir.create(neutralFBDLowDir)
  
  neutralFBDMidDir <- paste0(neutralFBDDir, "mid/")
  smart.dir.create(neutralFBDMidDir)
  
  neutralFBDHighDir <- paste0(neutralFBDDir, "high/")
  smart.dir.create(neutralFBDHighDir)
  
  # fossils directory
  fossilsDir <- paste0(targetDir, "fossils/")
  smart.dir.create(fossilsDir)
  
  # and tree directories
  treeDir <- paste0(targetDir, "trees/")
  smart.dir.create(treeDir)
  
  treeUltraDir <- paste0(treeDir, "ultrametric/")
  smart.dir.create(treeUltraDir)
  
  treeFBDDir <- paste0(treeDir, "fbd/")
  smart.dir.create(treeFBDDir)
  
  # print simulations to file
  invisible(lapply(1:nReps, function(x)
    capture.output(print(simList[[x]]),
                   file = paste0(simsDir, "sim_", x, ".txt")))) 
  invisible(lapply(1:nReps, function(x)
    capture.output(lapply(simList[[x]], function(x) {
      cat("\n")
      x}),
      file = paste0(simsDir, "sim_", x, ".txt"),
      append = TRUE))) 
  
  # save trait data as as .nex
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(realUltra[[x]], 
                     file = paste0(traitsUltraDir, "real_ultra_", x, ".nex"),
                     format = "standard")))
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(neutralLowUltra[[x]], 
                     file = paste0(neutralUltraLowDir, 
                                   "neutral_low_ultra_", x, ".nex"),
                     format = "standard")))
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(neutralMidUltra[[x]], 
                     file = paste0(neutralUltraMidDir, 
                                   "neutral_mid_ultra_", x, ".nex"),
                     format = "standard")))
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(neutralHighUltra[[x]], 
                     file = paste0(neutralUltraHighDir, 
                                   "neutral_high_ultra_", x, ".nex"),
                     format = "standard")))
  
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(realFBD[[x]], 
                     file = paste0(traitsFBDDir, "real_fbd_", x, ".nex"),
                     format = "standard")))
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(neutralLowFBD[[x]], 
                     file = paste0(neutralFBDLowDir, 
                                   "neutral_low_fbd_", x, ".nex"),
                     format = "standard")))
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(neutralMidFBD[[x]], 
                     file = paste0(neutralFBDMidDir, 
                                   "neutral_mid_fbd_", x, ".nex"),
                     format = "standard")))
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(neutralHighFBD[[x]], 
                     file = paste0(neutralFBDHighDir, 
                                   "neutral_high_fbd_", x, ".nex"),
                     format = "standard")))
  
  # save fossil records as tsv
  invisible(lapply(1:nReps, function(x)
    write_tsv(fossilsList[[x]], 
              file = paste0(fossilsDir, "fossils_", x, ".tsv"))))
  
  # save ultrametric trees
  invisible(lapply(1:nReps, function(x)
    write.nexus(treeUltraList[[x]], 
                file = paste0(treeUltraDir, "tree_ultra_", x, ".nex"))))
  # and fbd trees
  invisible(lapply(1:nReps, function(x)
    write.nexus(treeFBDList[[x]], 
                file = paste0(treeFBDDir, "tree_fbd_", x, ".nex"))))
  # write origins to FBD directory
  invisible(write.table(origin, quote = FALSE, 
                        col.names = FALSE, row.names = FALSE,
                        file = paste0(treeFBDDir, "origins_fbd.tsv"),
                        sep = "\t"))
  
  # save complete trees as well
  invisible(lapply(1:nReps, function(x)
    write.nexus(treeList[[x]], 
                file = paste0(treeDir, "tree_", x, ".nex"))))
}

# create function to run simulations for a list of parameters
simulate <- function(seeds, nReps, comb, key, simDir, N, psi) {
  ### recover parameters from key
  # ref
  ref <- key$ref[comb]
  
  # speciation rates
  lambda0 <- key$lambda0[comb]
  lambda1 <- key$lambda1[comb]
  lambda2 <- key$lambda2[comb]
  lambda3 <- key$lambda3[comb]
  
  # extinction rates
  mu0 <- key$mu0[comb]
  mu1 <- key$mu1[comb]
  mu2 <- key$mu2[comb]
  mu3 <- key$mu3[comb]
  
  # transition rates
  q01 <- key$q01[comb]
  q10 <- key$q10[comb]
  q02 <- key$q02[comb]
  q20 <- key$q20[comb]
  q03 <- key$q03[comb]
  q30 <- key$q30[comb]
  q12 <- key$q12[comb]
  q21 <- key$q21[comb]
  q13 <- key$q13[comb]
  q31 <- key$q31[comb]
  q23 <- key$q23[comb]
  q32 <- key$q32[comb]
  
  # fossil sampling rates
  psi <- psi
  
  ## create directories
  # base directory for simulations, if it wasn't created
  smart.dir.create(simDir)
  
  # create directory for this combination of simulations
  combDir <- paste0(simDir, comb, "_", ref, "/")
  smart.dir.create(combDir)
  
  # run simulations
  simReps <- lapply(1:nReps, function(x) {
    print(paste0("psi: ", psi, " comb: ", ref, " rep: ", x, 
                 " seed: ", seeds[x]))
    set.seed(seeds[x])
    simulate_rep(lambda0, lambda1, lambda2, lambda3, 
                 mu0, mu1, mu2, mu3,
                 q01, q10, q02, q20, q03, q30,
                 q12, q21, q13, q31, q23, q32, N, psi)
  })
  
  # save them
  invisible(save_sims(nReps, simReps, combDir))
}

###
# number of reps to run each combination of parameters
nReps <- 100

# psi ref
psiRefs <- c("1_low_psi", "2_med_psi", "3_high_psi")

# psi vector
psiVar <- c(0.01, 0.05, 0.1)

# for each psi
for (i in 1:length(psiVar)) {
  # simulations directory
  simDir <- paste0("/Users/petrucci/Documents/research/ssetests_chap1/", 
                   "simulation/replicates/", psiRefs[i], "/")
  smart.dir.create(simDir)
  
  # run simulations for each combination of parameters
  for (comb in 1:nrow(key)) {
    # get a seed for each rep
    seeds <- runif(nReps, (comb - 1)*nReps, comb*nReps)
    
    # create directory for this combination of simulations
    combDir <- paste0(simDir, comb, "_", key$ref[comb], "/")
    smart.dir.create(combDir)
    
    # save seeds
    save(seeds, file = paste0(combDir, "seeds.RData"))
    
    # simulate the reps with this parameter combination
    simulate(seeds, nReps, comb, key, simDir, N, psiVar[i])
  }
}
