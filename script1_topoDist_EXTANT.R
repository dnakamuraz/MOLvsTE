# Script 1: Calculating topological distances

#################
# LOAD PACKAGES #
#################

library(ape)
library(dplyr)
library(ggridges)
library(gridExtra)
library(phytools)
library(RNODE)
library(tidyverse)
library(TreeDist)
library(viridis)

##############
# READ TREES #
##############

setwd("/Users/labanfibios/Desktop/Doutorado/Project/B2_TEvsMOL/GitHub/Trees_extant")

# List all files in the directory
files <- list.files()
# Extract unique prefixes (assumes filenames start with numbers)
prefixes <- unique(sub("^(\\d+).*", "\\1", files))

# MP MOL STRICT CONSENSUS
mp_mol_cons <- list()
# Iterate over each prefix
for (prefix in prefixes) {
  # Filter files belonging to the current prefix
  prefix_files <- files[grep(paste0("^", prefix), files)]
  # Apply the appropriate function based on the filename pattern
  if (any(grepl("a_strictConsensus", prefix_files))) {
    mp_mol_cons[[paste0("mp_mol_cons_", prefix)]] <- read.tree(prefix_files[grepl("a_strictConsensus", prefix_files)])
  }
}

# MP MOL MPTs
mp_mol_mpts <- list()
# Iterate over each prefix
for (prefix in prefixes) {
  # Filter files belonging to the current prefix
  prefix_files <- files[grep(paste0("^", prefix), files)]
  # Apply the appropriate function based on the filename pattern
  if (any(grepl("a_MOL_TNT_results.out", prefix_files))) {
    mp_mol_mpts[[prefix]][["mp_mol_mpts"]] <- ReadTntTree(
      prefix_files[grepl("a_MOL_TNT_results.out", prefix_files)]
    )
  }
}

# MP MOL BOOTSTRAP (header in the files were manually deleted)
mp_mol_bs <- list()
# Iterate over each prefix
for (prefix in prefixes) {
  # Filter files belonging to the current prefix
  prefix_files <- files[grep(paste0("^", prefix), files)]
  # Apply the appropriate function based on the filename pattern
  if (any(grepl("b_MOL_BS_TNT.nwk", prefix_files))) {
    mp_mol_bs[[paste0("mp_mol_bs", prefix)]] <- read.tree(prefix_files[grepl("b_MOL_BS_TNT.nwk", prefix_files)])
  }
}

# MP TE STRICT CONSENSUS
mp_te_cons <- list()
# Iterate over each prefix
for (prefix in prefixes) {
  # Filter files belonging to the current prefix
  prefix_files <- files[grep(paste0("^", prefix), files)]
  # Apply the appropriate function based on the filename pattern
  if (any(grepl("c_strictConsensus_TE_TNT_results.nwk", prefix_files))) {
    mp_te_cons[[paste0("mp_te_cons", prefix)]] <- read.tree(prefix_files[grepl("c_strictConsensus_TE_TNT_results.nwk", prefix_files)])
  }
}

# MP TE MPTs
mp_te_mpts <- list()
# Iterate over each prefix
for (prefix in prefixes) {
  # Filter files belonging to the current prefix
  prefix_files <- files[grep(paste0("^", prefix), files)]
  # Apply the appropriate function based on the filename pattern
  if (any(grepl("c_TE_TNT_results.out", prefix_files))) {
    mp_te_mpts[[prefix]][["mp_te_mpts"]] <- ReadTntTree(
      prefix_files[grepl("c_TE_TNT_results.out", prefix_files)]
    )
  }
}

# MP TE BOOTSTRAP (header in the files were manually deleted)
mp_te_bs <- list()
# Iterate over each prefix
for (prefix in prefixes) {
  # Filter files belonging to the current prefix
  prefix_files <- files[grep(paste0("^", prefix), files)]
  # Apply the appropriate function based on the filename pattern
  if (any(grepl("d_TE_BS_TNT.nwk", prefix_files))) {
    mp_te_bs[[paste0("mp_te_bs", prefix)]] <- read.tree(prefix_files[grepl("d_TE_BS_TNT.nwk", prefix_files)])
  }
}

# ML MOL (optimal tree)
ml_mol_best = list()
# Iterate over each prefix
for (prefix in prefixes) {
  # Filter files belonging to the current prefix
  prefix_files <- files[grep(paste0("^", prefix), files)]
  # Apply the appropriate function based on the filename pattern
  if (any(grepl("_MOL_IQTREE.treefile", prefix_files))) {
    ml_mol_best[[paste0("ml_mol_best_", prefix)]] <- read.tree(prefix_files[grepl("_MOL_IQTREE.treefile", prefix_files)])
  }
}

# ML MOL (BS)
ml_mol_bs = list()
# Iterate over each prefix
for (prefix in prefixes) {
  # Filter files belonging to the current prefix
  prefix_files <- files[grep(paste0("^", prefix), files)]
  # Apply the appropriate function based on the filename pattern
  if (any(grepl("_MOL_IQTREE.contree", prefix_files))) {
    ml_mol_bs[[paste0("ml_mol_bs_", prefix)]] <- read.tree(prefix_files[grepl("_MOL_IQTREE.contree", prefix_files)])
  }
}

# ML TE ASC (optimal tree)
ml_te_asc_best = list()
# Iterate over each prefix
for (prefix in prefixes) {
  # Filter files belonging to the current prefix
  prefix_files <- files[grep(paste0("^", prefix), files)]
  # Apply the appropriate function based on the filename pattern
  if (any(grepl("_TE_ASC_IQTREE.treefile", prefix_files))) {
    ml_te_asc_best[[paste0("ml_te_asc_best_", prefix)]] <- read.tree(prefix_files[grepl("_TE_ASC_IQTREE.treefile", prefix_files)])
  }
}

# ML TE ASC (BS)
ml_te_asc_bs = list()
# Iterate over each prefix
for (prefix in prefixes) {
  # Filter files belonging to the current prefix
  prefix_files <- files[grep(paste0("^", prefix), files)]
  # Apply the appropriate function based on the filename pattern
  if (any(grepl("_TE_ASC_IQTREE.contree", prefix_files))) {
    ml_te_asc_bs[[paste0("ml_te_asc_bs_", prefix)]] <- read.tree(prefix_files[grepl("_TE_ASC_IQTREE.contree", prefix_files)])
  }
}

# ML TE noASC (optimal tree)
ml_te_noasc_best = list()
# Iterate over each prefix
for (prefix in prefixes) {
  # Filter files belonging to the current prefix
  prefix_files <- files[grep(paste0("^", prefix), files)]
  # Apply the appropriate function based on the filename pattern
  if (any(grepl("_TE_ASC_IQTREE.treefile", prefix_files))) {
    ml_te_noasc_best[[paste0("ml_te_noasc_best_", prefix)]] <- read.tree(prefix_files[grepl("_TE_noASC_IQTREE.treefile", prefix_files)])
  }
}

# ML TE noASC (BS)
ml_te_noasc_bs = list()
# Iterate over each prefix
for (prefix in prefixes) {
  # Filter files belonging to the current prefix
  prefix_files <- files[grep(paste0("^", prefix), files)]
  # Apply the appropriate function based on the filename pattern
  if (any(grepl("_TE_noASC_IQTREE.contree", prefix_files))) {
    ml_te_noasc_bs[[paste0("ml_te_noasc_bs_", prefix)]] <- read.tree(prefix_files[grepl("_TE_noASC_IQTREE.contree", prefix_files)])
  }
}

# Check number of datasets
length(mp_mol_cons) # 57
length(mp_mol_mpts) # 57
length(mp_mol_bs) # 57
length(mp_te_cons) # 57
length(mp_te_mpts) # 57
length(mp_te_bs) # 57
length(ml_mol_best) # 57
length(ml_mol_bs) # 57
length(ml_te_asc_best) # 57
length(ml_te_asc_bs) # 57
length(ml_te_noasc_best) # 57
length(ml_te_noasc_bs) # 57

#########################
# TOPOLOGICAL DISTANCES #
#########################

# A) MP MOL vs MP TE
MPmol_MPte = vector("list", length(mp_mol_cons))
# Compute SPR
SPR <- vector("list", length(mp_mol_mpts)) # empty list
batch_size <- 5 # run analyses in batches of 5 to avoid crash
n <- length(mp_mol_mpts) # number of analyses
for (start in seq(1, n, by = batch_size)) {
  end <- min(start + batch_size - 1, n)
  cat("\nProcessing batch:", start, "to", end, "\n")
  for (i in start:end) {
    SPR[[i]] <- RNODE::multiSPR(
      mp_mol_mpts[[i]]$mp_mol_mpts,
      mp_te_mpts[[i]]$mp_te_mpts,
      normalization = TRUE,
      method = "all"
    )
    cat("Result for i =", i, "\n")
    print(SPR[[i]])
  }
  saveRDS(SPR, file = "SPR_partial.rds")
  cat("Saved progress up to index", end, "\n")
  gc()
}
SPR <- lapply(SPR, function(x) {
  # If the element is a single value
  if (length(x) == 1) {
    data.frame(
      minRF = x,
      maxRF = x,
      meanRF = x
    )
  } else {
    # Otherwise keep it as is
    x
  }
})
# Compute RF
RF <- vector("list", length(mp_mol_mpts)) # empty list
batch_size <- 5 # run analyses in batches of 5 to avoid crash
n <- length(mp_mol_mpts) # number of analyses
for (start in seq(1, n, by = batch_size)) {
  end <- min(start + batch_size - 1, n)
  cat("\nProcessing batch:", start, "to", end, "\n")
  for (i in start:end) {
    RF[[i]] <- RNODE::multiRF(
      mp_mol_mpts[[i]]$mp_mol_mpts,
      mp_te_mpts[[i]]$mp_te_mpts,
      normalization = TRUE,
      method = "all"
    )
    cat("Result for i =", i, "\n")
    print(RF[[i]])
  }
  saveRDS(RF, file = "SPR_partial.rds")
  cat("Saved progress up to index", end, "\n")
  gc()
}
RF <- lapply(RF, function(x) {
  # If the element is a single value
  if (length(x) == 1) {
    data.frame(
      minRF = x,
      maxRF = x,
      meanRF = x
    )
  } else {
    # Otherwise keep it as is
    x
  }
})
# Compute CID
CID <- vector("list", length(mp_mol_mpts)) # empty list
batch_size <- 5 # run analyses in batches of 5 to avoid crash
n <- length(mp_mol_mpts) # number of analyses
for (start in seq(1, n, by = batch_size)) {
  end <- min(start + batch_size - 1, n)
  cat("\nProcessing batch:", start, "to", end, "\n")
  for (i in start:end) {
    CID[[i]] <- RNODE::multiCID(
      mp_mol_mpts[[i]]$mp_mol_mpts,
      mp_te_mpts[[i]]$mp_te_mpts,
      normalization = TRUE,
      method = "all"
    )
    cat("Result for i =", i, "\n")
    print(CID[[i]])
  }
  saveRDS(CID, file = "SPR_partial.rds")
  cat("Saved progress up to index", end, "\n")
  gc()
}
CID <- lapply(CID, function(x) {
  # If the element is a single value
  if (length(x) == 1) {
    data.frame(
      minCID = x,
      maxCID = x,
      meanCID = x
    )
  } else {
    # Otherwise keep it as is
    x
  }
})
# Summarize in a dataframe
MPmol_MPte_df <- data.frame(
  index = seq_along(SPR),
  setNames(do.call(rbind, lapply(SPR, unlist)), c("minSPR","maxSPR","meanSPR")),
  setNames(do.call(rbind, lapply(RF,  unlist)), c("minRF","maxRF","meanRF")),
  setNames(do.call(rbind, lapply(CID, unlist)), c("minCID","maxCID","meanCID"))
)
# Add a new column with the category of comparison
MPmol_MPte_df$comparison = "MPmol_MPte"
MPmol_MPte_df
#write_csv(MPmol_MPte_df, "../MPmol_MPte_extant.csv")

# B) ML MOL vs ML TE ASC
MLmol_MLteAsc = vector("list", length(ml_mol_best))
# For each dataset, calculate metrics of topological distance
for (i in seq_along(ml_mol_best)) {
  MLmol_MLteAsc[[i]] = summaryTopologicalDist(ml_mol_best[[i]], ml_te_asc_best[[i]])
}
# Convert vector to a df
MLmol_MLteAsc_df <- do.call(rbind, lapply(MLmol_MLteAsc, as.data.frame))
# Compute SPR
SPR = vector("list", length(ml_mol_best)) # empty list
for (i in seq_along(ml_mol_best)) {SPR[[i]] = RNODE::multiSPR(ml_mol_best[[i]],
                                                  ml_te_asc_best[[i]],
                                                  normalization=T, method='minSPR')}
MLmol_MLteAsc_df$SPR = SPR # Append vector of SPR to df
# Add a new column with the category of comparison
MLmol_MLteAsc_df$comparison = "MLmol_MLteAsc"
# Add a new column with the dataset name
MLmol_MLteAsc_df$dataset = names(ml_mol_best)
MLmol_MLteAsc_df <- MLmol_MLteAsc_df[, c(ncol(MLmol_MLteAsc_df), 1:(ncol(MLmol_MLteAsc_df) - 1))]
MLmol_MLteAsc_df
#write_csv(MLmol_MLteAsc_df, "../MLmol_MLteAsc_df.csv")

# No. polytomies: MOL
ml_mol_polytomies = vector("list", length(ml_mol_best))
for (i in seq_along(ml_mol_best)) {
  ml_mol_polytomies[[i]] = RNODE:::howManyPolytomies(ml_mol_best[[i]])}
values <- unlist(ml_mol_polytomies)
mean_val <- mean(values, na.rm = TRUE)
min_val  <- min(values, na.rm = TRUE)
max_val  <- max(values, na.rm = TRUE)
sd_val   <- sd(values, na.rm = TRUE)
summary_stats <- data.frame(
  Mean = mean_val,
  Min  = min_val,
  Max  = max_val,
  SD   = sd_val
)
summary_stats
# No. polytomies: TE ASC
ml_te_asc_polytomies = vector("list", length(ml_te_asc_best))
for (i in seq_along(ml_te_asc_best)) {
  ml_te_asc_polytomies[[i]] = RNODE:::howManyPolytomies(ml_te_asc_best[[i]])}
values <- unlist(ml_te_asc_polytomies)
mean_val <- mean(values, na.rm = TRUE)
min_val  <- min(values, na.rm = TRUE)
max_val  <- max(values, na.rm = TRUE)
sd_val   <- sd(values, na.rm = TRUE)
summary_stats <- data.frame(
  Mean = mean_val,
  Min  = min_val,
  Max  = max_val,
  SD   = sd_val
)
summary_stats

# C) ML MOL vs ML TE noASC
MLmol_MLteNoAsc = vector("list", length(ml_mol_best))
# For each dataset, calculate metrics of topological distance
for (i in seq_along(ml_mol_best)) {
  MLmol_MLteNoAsc[[i]] = summaryTopologicalDist(ml_mol_best[[i]], ml_te_noasc_best[[i]])
}
# Convert vector to a df
MLmol_MLteNoAsc_df <- do.call(rbind, lapply(MLmol_MLteNoAsc, as.data.frame))
# Compute SPR
SPR = vector("list", length(ml_mol_best)) # empty list
for (i in seq_along(ml_mol_best)) {SPR[[i]] = RNODE::multiSPR(ml_mol_best[[i]],
                                                              ml_te_noasc_best[[i]],
                                                              normalization=T, method='minSPR')}
MLmol_MLteNoAsc_df$SPR = SPR # Append vector of SPR to df
# Add a new column with the category of comparison
MLmol_MLteNoAsc_df$comparison = "MLmol_MLteNoAsc"
# Add a new column with the dataset name
MLmol_MLteNoAsc_df$dataset = names(ml_mol_best)
MLmol_MLteNoAsc_df <- MLmol_MLteNoAsc_df[, c(ncol(MLmol_MLteNoAsc_df), 1:(ncol(MLmol_MLteNoAsc_df) - 1))]
View(MLmol_MLteNoAsc_df)
#write.csv(MLmol_MLteNoAsc_df, "../MLmol_MLteNoAsc_df.csv")
# No. polytomies: TE noASC
ml_te_noasc_polytomies = vector("list", length(ml_te_noasc_best))
for (i in seq_along(ml_te_noasc_best)) {
  ml_te_noasc_polytomies[[i]] = RNODE:::howManyPolytomies(ml_te_noasc_best[[i]])}
values <- unlist(ml_te_noasc_polytomies)
mean_val <- mean(values, na.rm = TRUE)
min_val  <- min(values, na.rm = TRUE)
max_val  <- max(values, na.rm = TRUE)
sd_val   <- sd(values, na.rm = TRUE)
summary_stats <- data.frame(
  Mean = mean_val,
  Min  = min_val,
  Max  = max_val,
  SD   = sd_val
)
summary_stats
