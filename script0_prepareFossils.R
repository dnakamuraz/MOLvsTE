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

############################
# PREPARING TREES FOR ACID #
############################

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

