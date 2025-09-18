# Question 4. How do bootstrap values of shared clades differ between MOL and TE trees?

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

setwd("/Users/labanfibios/Desktop/Doutorado/Project/B2_TEvsMOL/GitHub/Trees_extant/")

##############
# READ TREES #
##############

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

# MP MOL BOOTSTRAP
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

# MP TE BOOTSTRAP
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

# ML MOL
ml_mol = list()
# Iterate over each prefix
for (prefix in prefixes) {
  # Filter files belonging to the current prefix
  prefix_files <- files[grep(paste0("^", prefix), files)]
  # Apply the appropriate function based on the filename pattern
  if (any(grepl("_MOL_IQTREE.contree", prefix_files))) {
    ml_mol[[paste0("ml_mol_cons_", prefix)]] <- read.tree(prefix_files[grepl("_MOL_IQTREE.contree", prefix_files)])
  }
}

# ML TE ASC
ml_te_asc = list()
# Iterate over each prefix
for (prefix in prefixes) {
  # Filter files belonging to the current prefix
  prefix_files <- files[grep(paste0("^", prefix), files)]
  # Apply the appropriate function based on the filename pattern
  if (any(grepl("_TE_ASC_IQTREE.contree", prefix_files))) {
    ml_te_asc[[paste0("ml_te_asc", prefix)]] <- read.tree(prefix_files[grepl("_TE_ASC_IQTREE.contree", prefix_files)])
  }
}

# ML TE noASC
ml_te_noasc = list()
# Iterate over each prefix
for (prefix in prefixes) {
  # Filter files belonging to the current prefix
  prefix_files <- files[grep(paste0("^", prefix), files)]
  # Apply the appropriate function based on the filename pattern
  if (any(grepl("_TE_noASC_IQTREE.contree", prefix_files))) {
    ml_te_noasc[[paste0("ml_te_noasc", prefix)]] <- read.tree(prefix_files[grepl("_TE_noASC_IQTREE.contree", prefix_files)])
  }
}

# Check number of datasets
length(mp_mol_cons) # 57
length(mp_mol_bs) # 57
length(mp_te_cons) # 57
length(mp_te_bs) # 57
length(ml_mol) # 57
length(ml_te_asc) # 57
length(ml_te_noasc) # 57

########################################
# DATASET: BOOTSTRAP OF SHARED CLADES  #
########################################

# A) MP MOL vs MP TE
MPmol_MPte = vector("list", length(mp_mol_bs))
# For each dataset, save BS of shared clades
for (i in seq_along(mp_mol_bs)) {
  MPmol_MPte[[i]] = sharedNodes(mp_mol_bs[[i]], mp_te_bs[[i]], spearman=F)
}
# Use lapply to remove rows where support is empty (i.e. root)
MPmol_MPte <- lapply(MPmol_MPte, function(df) {
  df[!is.na(df$Support_Tree_1) & df$Support_Tree_1 != "", , drop = FALSE]
})
# Add Approach1 and Approach2 columns
MPmol_MPte <- lapply(MPmol_MPte, function(df) {
  df$Approach1 <- "MOL" # Add column Approach1
  df$Approach2 <- "TE"  # Add column Approach2
  df
})
# Add the Dataset column
MPmol_MPte <- lapply(seq_along(MPmol_MPte), function(i) {
  df <- MPmol_MPte[[i]]
  df$Dataset <- i # Add a new column with the dataset number
  df
})
# Add ID for paired purposes
MPmol_MPte <- lapply(MPmol_MPte, function(df) {
  df$ID <- paste0(df$Dataset, "_", seq_len(nrow(df)))  # Combine Dataset and unique row number
  return(df)  # Return the modified dataframe
})
'
# Combine Node_Tree_2 below Node_Tree_1, Support_Tree_2 below Support_Tree_1, and Approach2 below Approach1
MPmol_MPte <- lapply(MPmol_MPte, function(df) {
  data.frame(
    Node = c(df$Node_Tree_1, df$Node_Tree_2),
    Support = c(df$Support_Tree_1, df$Support_Tree_2),
    Approach = c(df$Approach1, df$Approach2)
  )
})
'
# Convert the ID column to a factor in each dataframe
MPmol_MPte <- lapply(MPmol_MPte, function(df) {
  df$ID <- as.factor(df$ID) # Convert the Node column to factor
  df
})
# Combine the list into a single dataframe using base R
MPmol_MPte <- do.call(rbind, MPmol_MPte)
# Add the Optimality column
MPmol_MPte$Optimality <- "MP"


# B) ML MOL vs ML TE ASC
MLmol_MLteAsc = vector("list", length(ml_mol))
# For each dataset, save BS of shared clades
for (i in seq_along(ml_mol)) {
  MLmol_MLteAsc[[i]] = sharedNodes(ml_mol[[i]], ml_te_asc[[i]], spearman=F)
}
# Use lapply to remove rows where support is empty (i.e. root)
MLmol_MLteAsc <- lapply(MLmol_MLteAsc, function(df) {
  df[!is.na(df$Support_Tree_1) & df$Support_Tree_1 != "", , drop = FALSE]
})
# Add Approach1 and Approach2 columns
MLmol_MLteAsc <- lapply(MLmol_MLteAsc, function(df) {
  df$Approach1 <- "MOL" # Add column Approach1
  df$Approach2 <- "TE"  # Add column Approach2
  df
})
# Add the Dataset column
MLmol_MLteAsc <- lapply(seq_along(MLmol_MLteAsc), function(i) {
  df <- MLmol_MLteAsc[[i]]
  df$Dataset <- i # Add a new column with the dataset number
  df
})
# Add ID for paired purposes
MLmol_MLteAsc <- lapply(MLmol_MLteAsc, function(df) {
  df$ID <- paste0(df$Dataset, "_", seq_len(nrow(df)))  # Combine Dataset and unique row number
  return(df)  # Return the modified dataframe
})
'
# Combine Node_Tree_2 below Node_Tree_1, Support_Tree_2 below Support_Tree_1, and Approach2 below Approach1
MLmol_MLteAsc <- lapply(MLmol_MLteAsc, function(df) {
  data.frame(
    Node = c(df$Node_Tree_1, df$Node_Tree_2),
    Support = c(df$Support_Tree_1, df$Support_Tree_2),
    Approach = c(df$Approach1, df$Approach2)
  )
})
'
# Convert the ID column to a factor in each dataframe
MLmol_MLteAsc <- lapply(MLmol_MLteAsc, function(df) {
  df$ID <- as.factor(df$ID) # Convert the ID column to factor
  df
})
# Combine the list into a single dataframe using base R
MLmol_MLteAsc <- do.call(rbind, MLmol_MLteAsc)
# Add the Optimality column
MLmol_MLteAsc$Optimality <- "ML-ASC"


# C) ML MOL vs ML TE noASC
MLmol_MLteNoAsc = vector("list", length(ml_mol))
# For each dataset, save BS of shared clades
for (i in seq_along(ml_mol)) {
  MLmol_MLteNoAsc[[i]] = sharedNodes(ml_mol[[i]], ml_te_noasc[[i]], spearman=F)
}
# Use lapply to remove rows where support is empty (i.e. root)
MLmol_MLteNoAsc <- lapply(MLmol_MLteNoAsc, function(df) {
  df[!is.na(df$Support_Tree_1) & df$Support_Tree_1 != "", , drop = FALSE]
})
# Add Approach1 and Approach2 columns
MLmol_MLteNoAsc <- lapply(MLmol_MLteNoAsc, function(df) {
  df$Approach1 <- "MOL" # Add column Approach1
  df$Approach2 <- "TE"  # Add column Approach2
  df
})
# Add the Dataset column
MLmol_MLteNoAsc <- lapply(seq_along(MLmol_MLteNoAsc), function(i) {
  df <- MLmol_MLteNoAsc[[i]]
  df$Dataset <- i # Add a new column with the dataset number
  df
})
# Add ID for paired purposes
MLmol_MLteNoAsc <- lapply(MLmol_MLteNoAsc, function(df) {
  df$ID <- paste0(df$Dataset, "_", seq_len(nrow(df)))  # Combine Dataset and unique row number
  return(df)  # Return the modified dataframe
})
'
# Combine Node_Tree_2 below Node_Tree_1, Support_Tree_2 below Support_Tree_1, and Approach2 below Approach1
MLmol_MLteNoAsc <- lapply(MLmol_MLteNoAsc, function(df) {
  data.frame(
    Node = c(df$Node_Tree_1, df$Node_Tree_2),
    Support = c(df$Support_Tree_1, df$Support_Tree_2),
    Approach = c(df$Approach1, df$Approach2)
  )
})
'
# Convert the ID column to a factor in each dataframe
MLmol_MLteNoAsc <- lapply(MLmol_MLteNoAsc, function(df) {
  df$ID <- as.factor(df$ID) # Convert the ID column to factor
  df
})
# Combine the list into a single dataframe using base R
MLmol_MLteNoAsc <- do.call(rbind, MLmol_MLteNoAsc)
# Add the Optimality column
MLmol_MLteNoAsc$Optimality <- "ML-noASC"

# Combine the dataframes
combined_df <- rbind(MPmol_MPte, MLmol_MLteAsc, MLmol_MLteNoAsc)
combined_df$ID = paste0(combined_df$ID, "_", combined_df$Optimality)
# Convert support to numeric vectors
combined_df$Support_Tree_1 <- as.numeric(as.character(combined_df$Support_Tree_1))
combined_df$Support_Tree_2 <- as.numeric(as.character(combined_df$Support_Tree_2))
# Remove <50 BS values (present in IQTREE)
combined_df <- combined_df[combined_df$Support_Tree_1 >= 50 &
                             combined_df$Support_Tree_2 >= 50, ]
# Write to CSV file
write.csv(combined_df, file = "Dataset_R_q4_09-25_shared.csv", row.names = FALSE)

###################
# BS UNIQUE NODES #
###################

# A) MP MOL vs MP TE
MPmol_MPte_unique = vector("list", length(mp_mol_bs))
# For each dataset, save BS of shared clades
for (i in seq_along(mp_mol_bs)) {
  MPmol_MPte_unique[[i]] = uniqueNodes(mp_mol_bs[[i]], mp_te_bs[[i]], composition=F)
}
# Convert nested lists to a flat dataframe
MP_unique_df <- bind_rows(
  lapply(seq_along(MPmol_MPte_unique), function(i) {
    mol_df <- MPmol_MPte_unique[[i]][[1]] %>%
      mutate(Analysis = i, Type = "MOL", Method = "MP")
    
    te_df <- MPmol_MPte_unique[[i]][[2]] %>%
      mutate(Analysis = i, Type = "TE", Method = "MP")
    
    bind_rows(mol_df, te_df)
  })
)
# Reorder columns
MP_unique_df <- MP_unique_df %>% select(Analysis, Type, Method, Node, Support)
View(MP_unique_df)

# B) ML MOL vs ML TE ASC
MLmol_MLteAsc_unique = vector("list", length(ml_mol))
# For each dataset, save BS of shared clades
for (i in seq_along(ml_mol)) {
  MLmol_MLteAsc_unique[[i]] = uniqueNodes(ml_mol[[i]], ml_te_asc[[i]], composition=F)
}
# Convert nested lists to a flat dataframe
ML_ASC_unique_df <- bind_rows(
  lapply(seq_along(MLmol_MLteAsc_unique), function(i) {
    mol_df <- MLmol_MLteAsc_unique[[i]][[1]] %>%
      mutate(Analysis = i, Type = "MOL", Method = "ML_ASC")
    te_df <- MLmol_MLteAsc_unique[[i]][[2]] %>%
      mutate(Analysis = i, Type = "TE", Method = "ML_ASC")
    bind_rows(mol_df, te_df)
  })
)
# Reorder columns
ML_ASC_unique_df <- ML_ASC_unique_df %>% select(Analysis, Type, Method, Node, Support)
# Remove <50% BS from IQ-TREE
ML_ASC_unique_df <- ML_ASC_unique_df %>%
  mutate(Support = as.numeric(Support)) %>%
  filter(Support >= 50)
# View
View(ML_ASC_unique_df) 

# C) ML MOL vs MP TE noASC
MLmol_MLteNoAsc_unique = vector("list", length(ml_mol))
# For each dataset, save BS of shared clades
for (i in seq_along(ml_mol)) {
  MLmol_MLteNoAsc_unique[[i]] = uniqueNodes(ml_mol[[i]], ml_te_noasc[[i]], composition=F)
}
# Convert nested lists to a flat dataframe
ML_noASC_unique_df <- bind_rows(
  lapply(seq_along(MLmol_MLteNoAsc_unique), function(i) {
    mol_df <- MLmol_MLteNoAsc_unique[[i]][[1]] %>%
      mutate(Analysis = i, Type = "MOL", Method = "ML_noASC")
    te_df <- MLmol_MLteNoAsc_unique[[i]][[2]] %>%
      mutate(Analysis = i, Type = "TE", Method = "ML_noASC")
    bind_rows(mol_df, te_df)
  })
)
# Reorder columns
ML_noASC_unique_df <- ML_noASC_unique_df %>% select(Analysis, Type, Method, Node, Support)
# Remove <50% BS from IQ-TREE
ML_noASC_unique_df <- ML_noASC_unique_df %>%
  mutate(Support = as.numeric(Support)) %>%
  filter(Support >= 50)
# View
View(ML_noASC_unique_df) 

# Ensure Support is numeric in all dataframes
MP_unique_df <- MP_unique_df %>% mutate(Support = as.numeric(Support))
ML_ASC_unique_df <- ML_ASC_unique_df %>% mutate(Support = as.numeric(Support))
ML_noASC_unique_df <- ML_noASC_unique_df %>% mutate(Support = as.numeric(Support))
# Combine dataframes
all_unique_df <- bind_rows(MP_unique_df, ML_ASC_unique_df, ML_noASC_unique_df)
View(all_unique_df)
# Save as CSV
write.csv(all_unique_df, "./Dataset_R_q4_09-25_unique.csv", row.names = FALSE)



'
########################
# GLMMs in Posit Cloud #
########################

# Install and load the package
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  install.packages("glmmTMB")
}
if (!requireNamespace("emmeans", quietly = TRUE)) {
  install.packages("emmeans")
}
if (!requireNamespace("cowplot", quietly = TRUE)) {
  install.packages("cowplot")
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}
library("dplyr")
library("gridExtra")
library("cowplot")
library("emmeans")
library("ggplot2")
library("glmmTMB")
library("tidyr")

# Load the data
setwd("/Users/labanfibios/Desktop/Doutorado/Project/B2_TEvsMOL/GitHub/")
df <- read.csv("Dataset_R_q4_09-25.csv")
df$Optimality <- factor(df$Optimality, levels = c("MP", "ML-MKv", "ML-MK"))
df$SupportDiff = df$Support_Tree_1 - df$Support_Tree_2
View(df)

##########
# 4.a MP #
##########

df_MP <- df[grepl("MP", df$Optimality), ]
View(df_MP)

# Exploratory data analysis
mean(df_MP$Support_Tree_1) # MOL
mean(df_MP$Support_Tree_2) # TE
sd(df_MP$Support_Tree_1) # MOL
sd(df_MP$Support_Tree_2) # TE
min(df_MP$Support_Tree_1) # MOL
min(df_MP$Support_Tree_2) # TE
max(df_MP$Support_Tree_1) # MOL
max(df_MP$Support_Tree_2) # TE
sum(df_MP$Support_Tree_1 == df_MP$Support_Tree_2) # number of clades with BS MOL = BS TE
sum(df_MP$Support_Tree_1 > df_MP$Support_Tree_2) # number of clades with BS MOL > BS TE
sum(df_MP$Support_Tree_1 < df_MP$Support_Tree_2) # number of clades with BS MOL < BS TE

hist(1/df_MP$Support_Tree_1, breaks = 40, col = "lightblue", main = "Left-Skewed Distribution of Bootstrap Values", xlab = "Support Values")
hist(1/df_MP$Support_Tree_2, breaks = 40, col = "lightblue", main = "Left-Skewed Distribution of Bootstrap Values", xlab = "Support Values")

# Shapiro-Wilk test for normality
shapiro.test(df_MP$Support_Tree_1) # P < 0.05, non-parametric
shapiro.test(df_MP$Support_Tree_2) # P < 0.05, non-parametric

# Perform the Wilcoxon signed-rank test for paired values
mp_wilcoxon <- wilcox.test(as.numeric(df_MP$Support_Tree_1),
                           as.numeric(df_MP$Support_Tree_2),
                           paired = TRUE)
mp_wilcoxon

# Dataset for GLMMs
df_MP_glmm <- data.frame(
  Node = c(df_MP$Node_Tree_1, df_MP$Node_Tree_2),
  Approach = c(df_MP$Approach1, df_MP$Approach2),
  Support = c(df_MP$Support_Tree_1, df_MP$Support_Tree_2),
  ID = c(df_MP$ID, df_MP$ID)
)
mean(df_MP_glmm$Support) # Mean support of shared clades in MP analyses
sd(df_MP_glmm$Support) # SD support of shared clades in MP analyses
min(df_MP_glmm$Support) # Minimum support of shared clades in MP analyses
max(df_MP_glmm$Support) # Maximum support of shared clades in MP analyses

# GLMM: Beta (inappropriate due to the presence of zero and one)
#mp_beta <- glmmTMB(Support/100 ~ Approach + (1 | ID),
#                    family = beta_family(link = "logit"),
#                    zi = ~1,
#                    data = df_MP_glmm)
#summary(mp_beta)

# GLMM: Binomial
df_MP_glmm$Successes <- df_MP_glmm$Support
df_MP_glmm$Failures <- 100 - df_MP_glmm$Support
mp_binomial <- glmmTMB(cbind(Successes, Failures) ~ Approach + (1 | ID),
                       family = binomial(link = "logit"),
                       zi = ~1,
                       data = df_MP_glmm)
summary(mp_binomial)

# GLMM: Gamma (BS transformed to right-skewed distribution)
mp_gamma <- glmmTMB(1/Support ~ Approach + (1 | ID),
                    zi = ~1,
                    family = Gamma(link = "inverse"),
                    data = df_MP_glmm)
summary(mp_gamma)

##############
# 4.b ML ASC #
##############

df_ML_ASC <- df[grepl("ML-ASC", df$Optimality), ]
View(df_ML_ASC)

# Exploratory data analysis
mean(df_ML_ASC$Support_Tree_1) # MOL
mean(df_ML_ASC$Support_Tree_2) # TE
sd(df_ML_ASC$Support_Tree_1) # MOL
sd(df_ML_ASC$Support_Tree_2) # TE
var(df_ML_ASC$Support_Tree_1) # MOL
var(df_ML_ASC$Support_Tree_2) # TE
sum(df_ML_ASC$Support_Tree_1 == df_ML_ASC$Support_Tree_2) # number of clades with BS MOL = BS TE
sum(df_ML_ASC$Support_Tree_1 > df_ML_ASC$Support_Tree_2) # number of clades with BS MOL > BS TE
sum(df_ML_ASC$Support_Tree_1 < df_ML_ASC$Support_Tree_2) # number of clades with BS MOL < BS TE

hist(df_ML_ASC$Support_Tree_1, breaks = 40, col = "lightblue", main = "Left-Skewed Distribution of Bootstrap Values", xlab = "Support Values")
hist(df_ML_ASC$Support_Tree_2, breaks = 40, col = "lightblue", main = "Left-Skewed Distribution of Bootstrap Values", xlab = "Support Values")

# Shapiro-Wilk test for normality
shapiro.test(df_ML_ASC$Support_Tree_1) # P < 0.05, non-parametric
shapiro.test(df_ML_ASC$Support_Tree_2) # P < 0.05, non-parametric

# Perform the Wilcoxon signed-rank test for paired values
ml_asc_wilcoxon <- wilcox.test(as.numeric(df_ML_ASC$Support_Tree_1),
                               as.numeric(df_ML_ASC$Support_Tree_2),
                               paired = TRUE)
ml_asc_wilcoxon

# Dataset for GLMMs
df_ML_ASC_glmm <- data.frame(
  Node = c(df_ML_ASC$Node_Tree_1, df_ML_ASC$Node_Tree_2),
  Approach = c(df_ML_ASC$Approach1, df_ML_ASC$Approach2),
  Support = c(df_ML_ASC$Support_Tree_1, df_ML_ASC$Support_Tree_2),
  ID = c(df_ML_ASC$ID, df_ML_ASC$ID)
)
mean(df_ML_ASC_glmm$Support) # Mean support of shared clades in MP analyses
sd(df_ML_ASC_glmm$Support) # SD support of shared clades in MP analyses
min(df_ML_ASC_glmm$Support) # Minimum support of shared clades in MP analyses
max(df_ML_ASC_glmm$Support) # Maximum support of shared clades in MP analyses

# GLMM: Binomial
df_ML_ASC_glmm$Successes <- df_ML_ASC_glmm$Support
df_ML_ASC_glmm$Failures <- 100 - df_ML_ASC_glmm$Support
ml_asc_binomial <- glmmTMB(cbind(Successes, Failures) ~ Approach + (1 | ID),
                           family = binomial(link = "logit"),
                           zi = ~1,
                           data = df_ML_ASC_glmm)
summary(ml_asc_binomial)

# GLMM: Gamma (BS transformed to right-skewed distribution)
ml_asc_gamma <- glmmTMB(1/Support ~ Approach + (1 | ID),
                        zi = ~1,
                        family = Gamma(link = "inverse"),
                        data = df_ML_ASC_glmm)
summary(ml_asc_gamma)

#################
# 4.c ML no ASC #
#################

df_ML_noASC <- df[grepl("ML-noASC", df$Optimality), ]
View(df_ML_noASC)

# Exploratory data analysis
mean(df_ML_noASC$Support_Tree_1) # MOL
mean(df_ML_noASC$Support_Tree_2) # TE
sd(df_ML_noASC$Support_Tree_1) # MOL
sd(df_ML_noASC$Support_Tree_2) # TE
var(df_ML_noASC$Support_Tree_1) # MOL
var(df_ML_noASC$Support_Tree_2) # TE
sum(df_ML_noASC$Support_Tree_1 == df_ML_noASC$Support_Tree_2) # number of clades with BS MOL = BS TE
sum(df_ML_noASC$Support_Tree_1 > df_ML_noASC$Support_Tree_2) # number of clades with BS MOL > BS TE
sum(df_ML_noASC$Support_Tree_1 < df_ML_noASC$Support_Tree_2) # number of clades with BS MOL < BS TE

hist(df_ML_noASC$Support_Tree_1, breaks = 40, col = "lightblue", main = "Left-Skewed Distribution of Bootstrap Values", xlab = "Support Values")
hist(df_ML_noASC$Support_Tree_2, breaks = 40, col = "lightblue", main = "Left-Skewed Distribution of Bootstrap Values", xlab = "Support Values")

# Shapiro-Wilk test for normality
shapiro.test(df_ML_noASC$Support_Tree_1) # P < 0.05, non-parametric
shapiro.test(df_ML_noASC$Support_Tree_2) # P < 0.05, non-parametric

# Perform the Wilcoxon signed-rank test for paired values
ml_noasc_wilcoxon <- wilcox.test(as.numeric(df_ML_noASC$Support_Tree_1),
                                 as.numeric(df_ML_noASC$Support_Tree_2),
                                 paired = TRUE)
ml_noasc_wilcoxon

# Dataset for GLMMs
df_ML_noASC_glmm <- data.frame(
  Node = c(df_ML_noASC$Node_Tree_1, df_ML_noASC$Node_Tree_2),
  Approach = c(df_ML_noASC$Approach1, df_ML_noASC$Approach2),
  Support = c(df_ML_noASC$Support_Tree_1, df_ML_noASC$Support_Tree_2),
  ID = c(df_ML_noASC$ID, df_ML_noASC$ID)
)
mean(df_ML_noASC_glmm$Support) # Mean support of shared clades in MP analyses
sd(df_ML_noASC_glmm$Support) # SD support of shared clades in MP analyses
min(df_ML_noASC_glmm$Support) # Minimum support of shared clades in MP analyses
max(df_ML_noASC_glmm$Support) # Maximum support of shared clades in MP analyses

# GLMM: Binomial
df_ML_noASC_glmm$Successes <- df_ML_noASC_glmm$Support
df_ML_noASC_glmm$Failures <- 100 - df_ML_noASC_glmm$Support
ml_noasc_binomial <- glmmTMB(cbind(Successes, Failures) ~ Approach + (1 | ID),
                             family = binomial(link = "logit"),
                             zi = ~1,
                             data = df_ML_noASC_glmm)
summary(ml_asc_binomial)

# GLMM: Gamma (BS transformed to right-skewed distribution)
ml_noasc_gamma <- glmmTMB(1/Support ~ Approach + (1 | ID),
                          zi = ~1,
                          family = Gamma(link = "inverse"),
                          data = df_ML_noASC_glmm)
summary(ml_noasc_gamma)

################
# 4.d PLOTTING #
################

# Create a new column with the optimality criterion
df_MP_glmm$comparison = "MPmol_MPte"
df_ML_ASC_glmm$comparison = "MLmol_MLteAsc"
df_ML_noASC_glmm$comparison = "MLmol_MLteNoAsc"

# Merge all comparisons into the same dataset
df_plot = bind_rows(df_MP_glmm, df_ML_ASC_glmm, df_ML_noASC_glmm)
View(df_plot)

# Set order of comparisons
df_plot$comparison <- factor(df_plot$comparison, levels = c("MPmol_MPte", "MLmol_MLteAsc", "MLmol_MLteNoAsc"))
df_plot = df_plot[!is.na(df_plot$comparison), ]

# Set colors
custom_palette <- c("MPmol_MPte" = "#339999", "MLmol_MLteAsc" = "#6f77b4", "MLmol_MLteNoAsc" = "#cccccc")
# Create custom legend labels
custom_labels <- c(
  "MPmol_MPte" = "MP-MOL vs MP-TE",
  "MLmol_MLteAsc" = "ML-MOL vs ML-TE (with ASC)",
  "MLmol_MLteNoAsc" = "ML-MOL vs ML-TE (without ASC)"
)

# Box Plot MP
plot_mp = ggplot(df_MP_glmm, aes(x = Approach, y = Support)) +
  geom_boxplot(outlier.shape = NA) +  # Create boxplot without default outliers
  geom_jitter(color = "#339999", alpha = .05, width = 0.2) +  # Add transparent dots
  labs(x = "", y = "\nBootstrap\n", tag="A") +
  theme_minimal() +
  #theme(axis.text.x = element_blank(), legend.position = "none") +  # Hide legend if not needed
  scale_color_manual(values = custom_palette, labels = custom_labels)
# Box Plot ML ASC
plot_ml_asc = ggplot(df_ML_ASC_glmm, aes(x = Approach, y = Support)) +
  geom_boxplot(outlier.shape = NA) +  # Create boxplot without default outliers
  geom_jitter(color = "#6f77b4", alpha = .05, width = 0.2) +  # Add transparent dots
  labs(x = "", y = "\nBootstrap\n", tag="B") +
  theme_minimal() +
  scale_y_continuous(limits = c(50, NA)) + # Define o mínimo para 50, o máximo será ajustado automaticamente
  #theme(axis.text.x = element_blank(), legend.position = "none") +  # Hide legend if not needed
  scale_color_manual(values = custom_palette, labels = custom_labels)
# Box Plot ML noASC
plot_ml_noasc = ggplot(df_ML_noASC_glmm, aes(x = Approach, y = Support)) +
  geom_boxplot(outlier.shape = NA) +  # Create boxplot without default outliers
  geom_jitter(color = "#cccccc", alpha = .08, width = 0.2) +  # Add transparent dots
  labs(x = "", y = "\nBootstrap\n", tag="C") +
  theme_minimal() +
  scale_y_continuous(limits = c(50, NA)) + # Define o mínimo para 50, o máximo será ajustado automaticamente
  #theme(axis.text.x = element_blank(), legend.position = "none") +  # Hide legend if not needed
  scale_color_manual(values = custom_palette, labels = custom_labels)

# Extract the legend as a separate ggplot object
legend <- ggplot(df_plot, aes(x = comparison, y = Support, color = comparison)) +
  geom_point() +
  scale_color_manual(values = custom_palette, labels = custom_labels) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  )+
  guides(color = guide_legend(override.aes = list(size = 2, alpha=.7)))  # Increase dot size in legend
# Convert the plot to a grob
plot_grob <- ggplotGrob(legend)
# Extract the legend from the grob
legend_grob <- gtable::gtable_filter(plot_grob, "guide-box")

# Plot 2x2
tiff("figS2_BS.tiff", width = 6, height = 6, units = "in", res = 300)
grid.arrange(arrangeGrob(plot_mp + theme(plot.margin = margin(1,1,1,.01)),
                         plot_ml_asc + theme(plot.margin = margin(1,1,1,.01)),
                         plot_ml_noasc + theme(plot.margin = margin(1,1,1,.01)),
                         nrow=1, ncol=3),
             legend_grob,
             nrow=2,
             heights= c(20,1))
dev.off()

################
# UNIQUE NODES #
################

# A) MP MOL vs MP TE
uniqueMPmol_MPte = vector("list", length(mp_mol_bs))
# For each dataset, calculate no. unique clades in MOL and TE
for (i in seq_along(mp_mol_bs)) {
  uniqueMPmol_MPte[[i]] = uniqueNodes(mp_mol_bs[[i]], mp_te_bs[[i]], composition=F)
}
# MP: Check the number of shared clades in all datasets
mp_shared_clades = c() # empty vector
for (i in seq_along(MPmol_MPte)){
  mp_shared_clades[i] = length(MPmol_MPte[[i]]$Node_Tree_1)
}
mp_shared_clades = sum(mp_shared_clades)
mp_shared_clades
# MP: Check the total number of unique clades in MOL trees
mp_unique_MOL_clades = c() # empty vector
for (i in seq_along(uniqueMPmol_MPte)){
  mp_unique_MOL_clades[i] = as.numeric(nrow(uniqueMPmol_MPte[[i]][[1]]))
}
mp_unique_MOL_clades = sum(mp_unique_MOL_clades)
mp_unique_MOL_clades
# MP: Check the total number of unique clades in TE trees
mp_unique_TE_clades = c() # empty vector
for (i in seq_along(uniqueMPmol_MPte)){
  mp_unique_TE_clades[i] = as.numeric(nrow(uniqueMPmol_MPte[[i]][[2]]))
}
mp_unique_TE_clades = sum(mp_unique_TE_clades)
mp_unique_TE_clades
# MP: Check the total number of clades
mp_total_clades = mp_shared_clades+mp_unique_MOL_clades+mp_unique_TE_clades
mp_total_clades
# MP: Check the proportion of shared/total
mp_shared_clades/mp_total_clades

# B) ML MOL vs MP TE ASC
uniqueMLmol_MLteAsc = vector("list", length(ml_mol))
# For each dataset, calculate no. unique clades
for (i in seq_along(ml_mol)) {
  uniqueMLmol_MLteAsc[[i]] = uniqueNodes(ml_mol[[i]], ml_te_asc[[i]], composition=F)
}
# ML: Check the number of shared clades in all datasets
MLmol_MLteAsc_list <- split(MLmol_MLteAsc, MLmol_MLteAsc$Dataset)
mlasc_shared_clades = c() # empty vector
for (i in seq_along(MLmol_MLteAsc_list)){
  mlasc_shared_clades[i] = nrow(MLmol_MLteAsc_list[[i]])
}
mlasc_shared_clades = sum(mlasc_shared_clades)
mlasc_shared_clades
# ML: Check the total number of unique clades in MOL trees
mlasc_unique_MOL_clades = c() # empty vector
for (i in seq_along(uniqueMLmol_MLteAsc)){
  mlasc_unique_MOL_clades[i] = as.numeric(nrow(uniqueMLmol_MLteAsc[[i]][[1]]))
}
mlasc_unique_MOL_clades = sum(mlasc_unique_MOL_clades)
mlasc_unique_MOL_clades
# ML: Check the total number of unique clades in TE trees
mlasc_unique_TE_clades = c() # empty vector
for (i in seq_along(uniqueMLmol_MLteAsc)){
  mlasc_unique_TE_clades[i] = as.numeric(nrow(uniqueMLmol_MLteAsc[[i]][[2]]))
}
mlasc_unique_TE_clades = sum(mlasc_unique_TE_clades)
mlasc_unique_TE_clades
# ML: Check the total number of clades
mlasc_total_clades = mlasc_shared_clades+mlasc_unique_MOL_clades+mlasc_unique_TE_clades
mlasc_total_clades
# ML: Check the proportion of shared/total
mlasc_shared_clades/mlasc_total_clades

# C) ML MOL vs MP TE noASC
uniqueMLmol_MLteNoAsc = vector("list", length(ml_mol))
# For each dataset, calculate no. unique clades
for (i in seq_along(ml_mol)) {
  uniqueMLmol_MLteNoAsc[[i]] = uniqueNodes(ml_mol[[i]], ml_te_noasc[[i]], composition=F)
}
# ML: Check the number of shared clades in all datasets
MLmol_MLteNoAsc_list <- split(MLmol_MLteNoAsc, MLmol_MLteNoAsc$Dataset)
mlnoasc_shared_clades = c() # empty vector
for (i in seq_along(MLmol_MLteNoAsc_list)){
  mlnoasc_shared_clades[i] = nrow(MLmol_MLteNoAsc_list[[i]])
}
mlnoasc_shared_clades = sum(mlnoasc_shared_clades)
mlnoasc_shared_clades
# ML: Check the total number of unique clades in MOL trees
mlnoasc_unique_MOL_clades = c() # empty vector
for (i in seq_along(uniqueMLmol_MLteNoAsc)){
  mlnoasc_unique_MOL_clades[i] = as.numeric(nrow(uniqueMLmol_MLteNoAsc[[i]][[1]]))
}
mlnoasc_unique_MOL_clades = sum(mlnoasc_unique_MOL_clades)
mlnoasc_unique_MOL_clades
# ML: Check the total number of unique clades in TE trees
mlnoasc_unique_TE_clades = c() # empty vector
for (i in seq_along(uniqueMLmol_MLteNoAsc)){
  mlnoasc_unique_TE_clades[i] = as.numeric(nrow(uniqueMLmol_MLteNoAsc[[i]][[2]]))
}
mlnoasc_unique_TE_clades = sum(mlnoasc_unique_TE_clades)
mlnoasc_unique_TE_clades
# ML: Check the total number of clades
mlnoasc_total_clades = mlnoasc_shared_clades+mlnoasc_unique_MOL_clades+mlasc_unique_TE_clades
mlnoasc_total_clades
# ML: Check the proportion of shared/total
mlnoasc_shared_clades/mlnoasc_total_clades
'