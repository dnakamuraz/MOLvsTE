# Script 0: Preparing Trees_extinct

#################
# Load packages #
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
# Read trees #
##############

setwd("/Users/labanfibios/Desktop/Doutorado/Project/B2_TEvsMOL/GitHub/Trees_extinct")

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
length(mp_mol_cons) # 22
length(mp_mol_mpts) # 22
length(mp_mol_bs) # 22
length(mp_te_cons) # 22
length(mp_te_mpts) # 22
length(mp_te_bs) # 22
length(ml_mol_best) # 22
length(ml_mol_bs) # 22
length(ml_te_asc_best) # 22
length(ml_te_asc_bs) # 22
length(ml_te_noasc_best) # 22
length(ml_te_noasc_bs) # 22

####################
# DATA EXPLORATION #
####################

# Check number of nodes
collect_multiple_lists_Nnode_df <- function(lists_named, base_names, from = 3, to = 58) {
  indices <- sprintf("%03d", seq(from, to))
  df <- data.frame(matrix(NA_integer_, nrow = length(indices), ncol = length(base_names)),
                   row.names = indices,
                   stringsAsFactors = FALSE)
  colnames(df) <- base_names
  for (i in seq_along(indices)) {
    idx_str <- indices[i]
    for (j in seq_along(lists_named)) {
      lst <- lists_named[[j]]
      base <- base_names[j]
      nm1 <- paste0(base, idx_str)       # e.g. ml_mol_best003
      nm2 <- paste0(base, "_", idx_str)  # e.g. ml_mol_best_003
      el <- lst[[nm1]]
      if (is.null(el)) el <- lst[[nm2]]
      if (!is.null(el) && !is.null(el$Nnode)) {
        df[i, j] <- as.integer(el$Nnode)
      }
    }
  }
  df <- df[rowSums(!is.na(df)) > 0, , drop = FALSE]
  return(df)
}
df_Nnode <- collect_multiple_lists_Nnode_df(
  lists_named = list(
    mp_mol_cons,
    mp_te_cons,
    ml_mol_best,
    ml_te_asc_best,
    ml_te_noasc_best
  ),
  base_names = c(
    "mp_mol_cons",
    "mp_te_cons",
    "ml_mol_best",
    "ml_te_asc_best",
    "ml_te_noasc_best"
  )
)
View(df_Nnode)

# Check mean BS
collect_bootstrap_means_df <- function(from = 3, to = 58) {
  indices <- seq(from, to)
  idx_strs <- sprintf("%03d", indices)
  df <- data.frame(
    mp_mol_bs       = rep(NA_real_, length(indices)),
    mp_te_bs        = rep(NA_real_, length(indices)),
    ml_mol_bs       = rep(NA_real_, length(indices)),
    ml_te_asc_bs    = rep(NA_real_, length(indices)),
    ml_te_noasc_bs  = rep(NA_real_, length(indices)),
    row.names = idx_strs,
    stringsAsFactors = FALSE
  )
  get_mean_node_label <- function(lst, base, idx) {
    name1 <- paste0(base, idx)       # e.g. ml_mol_bs003
    name2 <- paste0(base, "_", idx)  # e.g. ml_mol_bs_003
    el <- lst[[name1]]
    if (is.null(el)) el <- lst[[name2]]
    if (is.null(el) || is.null(el$node.label)) return(NA_real_)
    vals <- as.numeric(el$node.label)
    if (all(is.na(vals))) return(NA_real_)
    round(mean(na.omit(vals)), 1)
  }
  for (i in seq_along(indices)) {
    idx <- idx_strs[i]
    df$mp_mol_bs[i]      <- get_mean_node_label(mp_mol_bs,      "mp_mol_bs",      idx)
    df$mp_te_bs[i]       <- get_mean_node_label(mp_te_bs,       "mp_te_bs",       idx)
    df$ml_mol_bs[i]      <- get_mean_node_label(ml_mol_bs,      "ml_mol_bs",      idx)
    df$ml_te_asc_bs[i]   <- get_mean_node_label(ml_te_asc_bs,   "ml_te_asc_bs",   idx)
    df$ml_te_noasc_bs[i] <- get_mean_node_label(ml_te_noasc_bs, "ml_te_noasc_bs", idx)
  }
  df <- df[rowSums(!is.na(df)) > 0, , drop = FALSE]
  return(df)
}
df_bootstrap_means <- collect_bootstrap_means_df()
View(df_bootstrap_means)

# Check number of polytomies (necessary because I want to report the number of polytomies in the tree with fossils, before pruning)
# No. polytomies: MP MOL
mp_mol_polytomies = vector("list", length(mp_mol_cons))
for (i in seq_along(mp_mol_cons)) {
  mp_mol_polytomies[[i]] = RNODE:::howManyPolytomies(mp_mol_cons[[i]])}
values <- unlist(mp_mol_polytomies)
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
# No. polytomies: MP TE
mp_te_polytomies = vector("list", length(mp_te_cons))
for (i in seq_along(mp_te_cons)) {
  mp_te_polytomies[[i]] = RNODE:::howManyPolytomies(mp_te_cons[[i]])}
values <- unlist(mp_te_polytomies)
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
# No. polytomies: ML MOL
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
# No. polytomies: TE ASC (MKv)
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
# No. polytomies: TE noASC (MK)
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

#########################################
# PREPARING TE TREES BY PRUNING FOSSILS #
#########################################

# It is necessary because multiSPR does not prune automatically.

# MP TE: Identify leaves to be deleted (= fossils)
stopifnot(length(mp_mol_cons) == length(mp_te_cons))
mp_te_cons_pruned <- vector("list", length(mp_te_cons))
mp_te_cons_pruned_tips <- vector("list", length(mp_te_cons))
for (i in seq_along(mp_te_cons)) {
  tree_mol <- mp_mol_cons[[i]]
  tree_te  <- mp_te_cons[[i]]
  if (inherits(tree_mol, "phylo") && inherits(tree_te, "phylo")) {
    # Find tips unique to mp_te_cons
    unique_tips <- setdiff(tree_te$tip.label, tree_mol$tip.label)
    mp_te_cons_pruned_tips[[i]] <- unique_tips
    if (length(unique_tips) > 0) {
      mp_te_cons_pruned[[i]] <- drop.tip(tree_te, unique_tips)
    } else {
      mp_te_cons_pruned[[i]] <- tree_te
    }
  } else {
    mp_te_cons_pruned[[i]] <- NULL
    mp_te_cons_pruned_tips[[i]] <- NULL
  }
}
names(mp_te_cons_pruned_tips) <- names(mp_te_cons)
# ---- Prune leaves in TE MPTs ----
mp_te_mpts_pruned <- vector("list", length(mp_te_mpts))
for (X in seq_along(mp_te_mpts)) {
  inner_list <- mp_te_mpts[[X]]$mp_te_mpts
  tips_to_remove <- mp_te_cons_pruned_tips[[X]]
  new_inner <- vector("list", length(inner_list))
  for (Y in seq_along(inner_list)) {
    tree_obj <- inner_list[[Y]]
    if (inherits(tree_obj, "phylo") && !is.null(tips_to_remove) && length(tips_to_remove) > 0) {
      new_inner[[Y]] <- drop.tip(tree_obj, tips_to_remove)
    } else {
      new_inner[[Y]] <- tree_obj
    }
  }
  if (!is.null(names(inner_list))) names(new_inner) <- names(inner_list)
  class(new_inner) <- class(inner_list)
  mp_te_mpts_pruned[[X]] <- list(mp_te_mpts = new_inner)
}
names(mp_te_mpts_pruned) <- names(mp_te_mpts)
# ---- Now do the same for MP MOL ----
mp_mol_cons_pruned_tips <- vector("list", length(mp_mol_cons))
mp_mol_cons_pruned <- vector("list", length(mp_mol_cons))
for (i in seq_along(mp_mol_cons)) {
  tree_mol <- mp_mol_cons[[i]]
  tree_te  <- mp_te_cons[[i]]
  if (inherits(tree_mol, "phylo") && inherits(tree_te, "phylo")) {
    # Find tips unique to mp_mol_cons
    unique_tips <- setdiff(tree_mol$tip.label, tree_te$tip.label)
    mp_mol_cons_pruned_tips[[i]] <- unique_tips
    if (length(unique_tips) > 0) {
      mp_mol_cons_pruned[[i]] <- drop.tip(tree_mol, unique_tips)
    } else {
      mp_mol_cons_pruned[[i]] <- tree_mol
    }
  } else {
    mp_mol_cons_pruned[[i]] <- NULL
    mp_mol_cons_pruned_tips[[i]] <- NULL
  }
}
names(mp_mol_cons_pruned_tips) <- names(mp_mol_cons)
mp_mol_mpts_pruned <- vector("list", length(mp_mol_mpts))
for (X in seq_along(mp_mol_mpts)) {
  inner_list <- mp_mol_mpts[[X]]$mp_mol_mpts
  tips_to_remove <- mp_mol_cons_pruned_tips[[X]]
  new_inner <- vector("list", length(inner_list))
  for (Y in seq_along(inner_list)) {
    tree_obj <- inner_list[[Y]]
    if (inherits(tree_obj, "phylo") && !is.null(tips_to_remove) && length(tips_to_remove) > 0) {
      new_inner[[Y]] <- drop.tip(tree_obj, tips_to_remove)
    } else {
      new_inner[[Y]] <- tree_obj
    }
  }
  if (!is.null(names(inner_list))) names(new_inner) <- names(inner_list)
  class(new_inner) <- class(inner_list)
  mp_mol_mpts_pruned[[X]] <- list(mp_mol_mpts = new_inner)
}
names(mp_mol_mpts_pruned) <- names(mp_mol_mpts)
# ---- Outputs ----
mp_te_mpts_pruned
mp_mol_mpts_pruned

# ML ASC: Identify leaves to be deleted (= fossils)
stopifnot(length(ml_mol_best) == length(ml_te_asc_best))
ml_mol_best_pruned <- vector("list", length(ml_mol_best))
ml_te_asc_best_pruned <- vector("list", length(ml_te_asc_best))
pruned_tips_list <- vector("list", length(ml_mol_best))  # optional: store pruned tips per pair
for (i in seq_along(ml_mol_best)) {
  tree1 <- ml_mol_best[[i]]
  tree2 <- ml_te_asc_best[[i]]
  if (inherits(tree1, "phylo") && inherits(tree2, "phylo")) {
    shared_terminals <- intersect(tree1$tip.label, tree2$tip.label)
    # Identify which tips will be removed from each
    pruned_tips_tree1 <- setdiff(tree1$tip.label, shared_terminals)
    pruned_tips_tree2 <- setdiff(tree2$tip.label, shared_terminals)
    # Store info for inspection
    pruned_tips_list[[i]] <- list(
      from_tree1 = pruned_tips_tree1,
      from_tree2 = pruned_tips_tree2
    )
    # Prune
    ml_mol_best_pruned[[i]] <- if (length(pruned_tips_tree1) > 0) {
      drop.tip(tree1, pruned_tips_tree1)
    } else {
      tree1
    }
    ml_te_asc_best_pruned[[i]] <- if (length(pruned_tips_tree2) > 0) {
      drop.tip(tree2, pruned_tips_tree2)
    } else {
      tree2
    }
  } else {
    ml_mol_best_pruned[[i]] <- NULL
    ml_te_asc_best_pruned[[i]] <- NULL
    pruned_tips_list[[i]] <- NULL
  }
}
# Preserve names if they exist
names(ml_mol_best_pruned) <- names(ml_mol_best)
names(ml_te_asc_best_pruned) <- names(ml_te_asc_best)
names(pruned_tips_list) <- names(ml_mol_best)

# ML noASC: Identify leaves to be deleted (= fossils)
stopifnot(length(ml_mol_best) == length(ml_te_noasc_best))
ml_mol_best_pruned_noasc <- vector("list", length(ml_mol_best))
ml_te_noasc_best_pruned <- vector("list", length(ml_te_noasc_best))
pruned_tips_list_noasc <- vector("list", length(ml_mol_best))
for (i in seq_along(ml_mol_best)) {
  tree1 <- ml_mol_best[[i]]
  tree2 <- ml_te_noasc_best[[i]]
  if (inherits(tree1, "phylo") && inherits(tree2, "phylo")) {
    shared_terminals <- intersect(tree1$tip.label, tree2$tip.label)
    pruned_tips_tree1 <- setdiff(tree1$tip.label, shared_terminals)
    pruned_tips_tree2 <- setdiff(tree2$tip.label, shared_terminals)
    pruned_tips_list_noasc[[i]] <- list(
      from_tree1 = pruned_tips_tree1,
      from_tree2 = pruned_tips_tree2
    )
    ml_mol_best_pruned_noasc[[i]] <- if (length(pruned_tips_tree1) > 0) {
      drop.tip(tree1, pruned_tips_tree1)
    } else {
      tree1
    }
    ml_te_noasc_best_pruned[[i]] <- if (length(pruned_tips_tree2) > 0) {
      drop.tip(tree2, pruned_tips_tree2)
    } else {
      tree2
    }
  } else {
    ml_mol_best_pruned_noasc[[i]] <- NULL
    ml_te_noasc_best_pruned[[i]] <- NULL
    pruned_tips_list_noasc[[i]] <- NULL
  }
}
# Preserve names if they exist
names(ml_mol_best_pruned_noasc) <- names(ml_mol_best)
names(ml_te_noasc_best_pruned) <- names(ml_te_noasc_best)
names(pruned_tips_list_noasc) <- names(ml_mol_best)

#####################################
# TOPOLOGICAL DISTANCES EXCEPT ACID #
#####################################

# A) MP MOL vs MP TE
MPmol_MPte = vector("list", length(mp_mol_cons))
# For each dataset, calculate metrics of topological distance
for (i in seq_along(mp_mol_cons)) {
  MPmol_MPte[[i]] = summaryTopologicalDist(mp_mol_cons[[i]], mp_te_cons[[i]])
}
# Convert list to a df
MPmol_MPte_df <- do.call(rbind, lapply(MPmol_MPte, as.data.frame))
# Compute SPR (CHECK VALUES AND, IF NA IS PRESENT, CALCULATE MANUALLY)
SPR = vector("list", length(mp_mol_mpts)) # empty list
for (i in seq_along(mp_mol_mpts)) {SPR[[i]] = RNODE::multiSPR(mp_mol_mpts[[i]]$mp_mol_mpts,
                                                              mp_te_mpts_pruned[[i]]$mp_te_mpts, 
                                                              normalization=T, 
                                                              method = 'minSPR')}
MPmol_MPte_df$SPR = SPR # Append vector of SPR to df
# Add a new column with the category of comparison
MPmol_MPte_df$comparison = "MPmol_MPte"
# Add a new column with the dataset name
MPmol_MPte_df$dataset = names(mp_mol_cons)
MPmol_MPte_df <- MPmol_MPte_df[, c(ncol(MPmol_MPte_df), 1:(ncol(MPmol_MPte_df) - 1))]
View(MPmol_MPte_df)
#write_csv(MPmol_MPte_df, "../MPmol_MPte.csv")

# B) ML MOL vs ML TE ASC
MLmol_MLteAsc = vector("list", length(ml_mol_best))
# For each dataset, calculate metrics of topological distance
for (i in seq_along(ml_mol_best)) {
  MLmol_MLteAsc[[i]] = summaryTopologicalDist(ml_mol_best[[i]], ml_te_asc_best_pruned[[i]])
}
# Convert vector to a df
MLmol_MLteAsc_df <- do.call(rbind, lapply(MLmol_MLteAsc, as.data.frame))
# Compute SPR
SPR = vector("list", length(ml_mol_best)) # empty list
for (i in seq_along(ml_mol_best)) {SPR[[i]] = RNODE::multiSPR(ml_mol_best[[i]],
                                                              ml_te_asc_best_pruned[[i]],
                                                              normalization=T, method='meanSPR')}
MLmol_MLteAsc_df$SPR = SPR # Append vector of SPR to df
# Add a new column with the category of comparison
MLmol_MLteAsc_df$comparison = "MLmol_MLteAsc"
# Add a new column with the dataset name
MLmol_MLteAsc_df$dataset = names(ml_mol_best)
MLmol_MLteAsc_df <- MLmol_MLteAsc_df[, c(ncol(MLmol_MLteAsc_df), 1:(ncol(MLmol_MLteAsc_df) - 1))]
MLmol_MLteAsc_df
#write_csv(MLmol_MLteAsc_df, "../MLmol_MLteAsc_df.csv")

# C) ML MOL vs ML TE noASC
MLmol_MLteNoAsc = vector("list", length(ml_mol_best))
# For each dataset, calculate metrics of topological distance
for (i in seq_along(ml_mol_best)) {
  MLmol_MLteNoAsc[[i]] = summaryTopologicalDist(ml_mol_best[[i]], ml_te_noasc_best_pruned[[i]])
}
# Convert vector to a df
MLmol_MLteNoAsc_df <- do.call(rbind, lapply(MLmol_MLteNoAsc, as.data.frame))
# Compute SPR
SPR = vector("list", length(ml_mol_best)) # empty list
for (i in seq_along(ml_mol_best)) {SPR[[i]] = RNODE::multiSPR(ml_mol_best[[i]],
                                                              ml_te_noasc_best_pruned[[i]],
                                                              normalization=T, method='minSPR')}
MLmol_MLteNoAsc_df$SPR = SPR # Append vector of SPR to df
# Add a new column with the category of comparison
MLmol_MLteNoAsc_df$comparison = "MLmol_MLteNoAsc"
# Add a new column with the dataset name
MLmol_MLteNoAsc_df$dataset = names(ml_mol_best)
MLmol_MLteNoAsc_df <- MLmol_MLteNoAsc_df[, c(ncol(MLmol_MLteNoAsc_df), 1:(ncol(MLmol_MLteNoAsc_df) - 1))]
View(MLmol_MLteNoAsc_df)
df_flat <- MLmol_MLteNoAsc_df
df_flat[] <- lapply(df_flat, function(col) {
  if (is.list(col)) sapply(col, function(x) paste(x, collapse = ";")) else col
})
write.csv(df_flat, "../MLmol_MLteNoAsc_df.csv", row.names = FALSE)

################################
# PREPARING TREES FOR TNT ACID #
################################

# Set directory for pruned trees
output_dir <- "../Trees_extinct_ACID/pruned_trees"
log_file <- file.path(output_dir, "pruning_log.txt")
# Create output directory if it does not exist
if (!dir.exists(output_dir)) dir.create(output_dir)

# Function to prepare TE trees pruning the fossils to compute ACID
prune_tree_lists <- function(list1, list2, output_dir, log = TRUE) {
  # Check lengths
  if (length(list1) != length(list2)) {
    stop("Both lists must have the same length (same number of trees).")
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Prepare log collector if needed
  if (log) {
    log_lines <- c("Pruning Report", "================", "")
    log_file <- file.path(output_dir, "pruning_log.txt")
  }
  
  # Internal helper to write the tree with header and footer
  write_with_header <- function(tree, filepath) {
    con <- file(filepath, open = "w")
    writeLines(c(
      "#NEXUS",
      "begin trees ;",
      "tree tagged_tree = [&U]"
    ), con)
    write.tree(tree, file = con, append = TRUE)
    writeLines("end ;", con)
    close(con)
  }
  
  # Iterate over pairs of trees
  for (i in seq_along(list1)) {
    tree1 <- list1[[i]]
    tree2 <- list2[[i]]
    
    # Get tip labels
    tips1 <- tree1$tip.label
    tips2 <- tree2$tip.label
    
    # Find common and unique tips
    common_tips <- intersect(tips1, tips2)
    unique1 <- setdiff(tips1, tips2)
    unique2 <- setdiff(tips2, tips1)
    
    # Prune trees
    pruned1 <- drop.tip(tree1, setdiff(tips1, common_tips))
    pruned2 <- drop.tip(tree2, setdiff(tips2, common_tips))
    
    # File names ending with .nwk
    tree1_file <- file.path(output_dir, paste0("list1_pruned_", i, ".nwk"))
    tree2_file <- file.path(output_dir, paste0("list2_pruned_", i, ".nwk"))
    
    # Write pruned trees with header and footer
    write_with_header(pruned1, tree1_file)
    write_with_header(pruned2, tree2_file)
    
    # Append log info if requested
    if (log) {
      log_lines <- c(
        log_lines,
        paste0("Pair ", i, ":"),
        paste("  Number of common tips:", length(common_tips)),
        paste("  Unique tips in list1[[", i, "]] (", length(unique1), "):",
              if (length(unique1) > 0) paste(unique1, collapse = ", ") else "None"),
        paste("  Unique tips in list2[[", i, "]] (", length(unique2), "):",
              if (length(unique2) > 0) paste(unique2, collapse = ", ") else "None"),
        ""
      )
    }
  }
  # Write log file if needed
  if (log) {
    writeLines(log_lines, con = log_file)
    message("Log file written to: ", log_file)
  }
  message("Pruning complete. Pruned trees saved in: ", output_dir)
}

# MP
prune_tree_lists(mp_mol_cons, mp_te_cons, output_dir="../Trees_extinct_ACID/pruned_trees_MP/", log=T)

# ML ASC
prune_tree_lists(ml_mol_best, ml_te_asc_best, output_dir="../Trees_extinct_ACID/pruned_trees_ML_ASC/", log=T)

# ML noASC
prune_tree_lists(ml_mol_best, ml_te_noasc_best, output_dir="../Trees_extinct_ACID/pruned_trees_ML_noASC/", log=T)

