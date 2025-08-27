# Script 1: Calculating topological distances

setwd("/Users/labanfibios/Desktop/Doutorado/Project/B2_TEvsMOL/GitHub/")

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

##################
# EXAMPLE 1: 046 #
##################

# Trees from 046 (without fossils, MP, Fig. 2A, Table S3)
# Map BS to optimal trees
mp_046_mol = read.tree("Trees_extant/046a_strictConsensus_MOL_TNT_results.nwk")
mp_046_mol_bs = read.tree("Trees_extant/046b_MOL_BS_TNT.nwk")
RNODE::mapSupport(mp_046_mol, mp_046_mol_bs,  
                  plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(1.2,-0.5), tree.cex=.8, 
                  tree.output = "Figures/Fig2_046_OPTvsBS_MP_MOL.pdf")
mp_046_te = read.tree("Trees_extant/046c_strictConsensus_TE_TNT_results.nwk")
mp_046_te_bs = read.tree("Trees_extant/046d_TE_BS_TNT.nwk")
RNODE::mapSupport(mp_046_te, mp_046_te_bs,
                  plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(1.2,-0.5), tree.cex=.8, 
                  tree.output = "Figures/Fig2_046_OPTvsBS_MP_TE.pdf")
# Identify shared nodes and plot tangle trees
mp_shared = RNODE::sharedNodes(mp_046_mol_bs, mp_046_te_bs, dataframe=F, root="Morone_saxatilis",
                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=15, tree.fsize=1, tree.adj=c(-1.5,0.5), tree.cex=0.6, 
                   output.tree1="Figures/Fig2_046_MP_MOL_extant.pdf", output.tree2="Figures/Fig2_046_MP_TE_extant.pdf",
                   tanglegram = T, output.tangletree = "Figures/Fig2_046_MP_MOLvsTE_extant_tangle.pdf", tanglegram.lab.cex=0.4, tanglegram.margin=5,
                   write.pruned=T, write.pruned1.name="Figures/Fig2_046_MP_MOL_extant_pruned.nwk", write.pruned2.name="Figures/Fig2_046_MP_TE_extant_pruned.nwk")
# No. of nodes where BS MOL = BS TE
a = sum(mp_shared$Support_Tree_1 == mp_shared$Support_Tree_2, na.rm = TRUE) # number
a
round(a*100 / as.numeric(nrow(mp_shared)), digits=3) # proportion
# No. of nodes where BS MOL > BS TE
a = sum(mp_shared$Support_Tree_1 > mp_shared$Support_Tree_2, na.rm = TRUE) # number
a
round(a*100 / as.numeric(nrow(mp_shared)), digits=3) # proportion
# No. of nodes where BS MOL < BS TE
a = sum(mp_shared$Support_Tree_1 < mp_shared$Support_Tree_2, na.rm = TRUE) # number
a
round(a*100 / as.numeric(nrow(mp_shared)), digits=3) # proportion
# Identify and plot unique nodes
RNODE::uniqueNodes(mp_046_mol, mp_046_te,
                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(-1.5,0.5), tree.cex=2, output.tree="Figures/Fig2_046_MP_MOLvsTE_extant_unique.pdf")

# Trees from 046 (with fossils, MP, Fig. 2B, Table S3)
# Map BS to optimal trees
mp_046_mol_f = read.tree("Trees_extinct//046a_strictConsensus_MOL_TNT_results.nwk")
mp_046_mol_bs_f = read.tree("Trees_extinct/046b_MOL_BS_TNT.nwk")
RNODE::mapSupport(mp_046_mol_f, mp_046_mol_bs_f,  
                  plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(1.2,-0.5), tree.cex=.8, 
                  tree.output = "Figures/Fig2_046_OPTvsBS_MP_MOL_f.pdf")
mp_046_te_f = read.tree("Trees_extinct//046c_strictConsensus_TE_TNT_results.nwk")
mp_046_te_bs_f = read.tree("Trees_extinct/046d_TE_BS_TNT.nwk")
RNODE::mapSupport(mp_046_te_f, mp_046_te_bs_f,
                  plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(1.2,-0.5), tree.cex=.8, 
                  tree.output = "Figures/Fig2_046_OPTvsBS_MP_TE_f.pdf")
# Identify shared nodes and plot tangle trees
mp_f_shared = RNODE::sharedNodes(mp_046_mol_bs_f, mp_046_te_bs_f, dataframe=F, root="Morone_saxatilis",
                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=15, tree.fsize=1, tree.adj=c(-1.5,0.5), tree.cex=0.6, 
                   output.tree1="Figures/Fig2_046_MP_MOL_extinct.pdf", output.tree2="Figures/Fig2_046_MP_TE_extinct.pdf",
                   tanglegram = T, output.tangletree = "Figures/Fig2_046_MP_MOLvsTE_extinct_tangle.pdf", tanglegram.lab.cex=0.4, tanglegram.margin=5,
                   write.pruned=T, write.pruned1.name="Figures/Fig2_046_MP_MOL_extinct_pruned.nwk", write.pruned2.name="Figures/Fig2_046_MP_TE_extinct_pruned.nwk")
# No. of nodes where BS MOL = BS TE
a = sum(mp_f_shared$Support_Tree_1 == mp_f_shared$Support_Tree_2, na.rm = TRUE) # number
a
round(a*100 / as.numeric(nrow(mp_f_shared)), digits=3) # proportion
# No. of nodes where BS MOL > BS TE
a = sum(mp_f_shared$Support_Tree_1 > mp_f_shared$Support_Tree_2, na.rm = TRUE) # number
a
round(a*100 / as.numeric(nrow(mp_f_shared)), digits=3) # proportion
# No. of nodes where BS MOL < BS TE
a = sum(mp_f_shared$Support_Tree_1 < mp_f_shared$Support_Tree_2, na.rm = TRUE) # number
a
round(a*100 / as.numeric(nrow(mp_f_shared)), digits=3) # proportion
# Identify and plot unique nodes
RNODE::uniqueNodes(mp_046_mol_f, mp_046_te_f,
                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(-1.5,0.5), tree.cex=2, output.tree="Figures/Fig2_046_MP_MOLvsTE_extinct_unique.pdf")


# Trees from 046 (without fossils, ML MKV, Fig. 2C, Table S3)
# Map BS to optimal trees
ml_046_mol = read.tree("Trees_extant/046_MOL_IQTREE.treefile")
ml_046_mol_bs = read.tree("Trees_extant/046_MOL_IQTREE.contree")
RNODE::mapSupport(ml_046_mol, ml_046_mol_bs,  
                  plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(1.2,-0.5), tree.cex=.8, 
                  tree.output = "Figures/Fig2_046_OPTvsBS_ML_MOL.pdf")
ml_046_te = read.tree("Trees_extant/046_TE_ASC_IQTREE.treefile")
ml_046_te_bs = read.tree("Trees_extant/046_TE_ASC_IQTREE.contree")
RNODE::mapSupport(ml_046_te, ml_046_te_bs,
                  plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(1.2,-0.5), tree.cex=.8, 
                  tree.output = "Figures/Fig2_046_OPTvsBS_ML_TE.pdf")
# Identify shared nodes and plot tangle trees
ml_mkv_shared = RNODE::sharedNodes(ml_046_mol_bs, ml_046_te_bs, dataframe=F, root="Morone_saxatilis",
                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=15, tree.fsize=1, tree.adj=c(-1.5,0.5), tree.cex=0.6, 
                   output.tree1="Figures/Fig2_046_ML_MOL_extant.pdf", output.tree2="Figures/Fig2_046_ML_TE_extant.pdf",
                   tanglegram = T, output.tangletree = "Figures/Fig2_046_ML_MOLvsTE_extant_tangle.pdf", tanglegram.lab.cex=0.4, tanglegram.margin=5,
                   write.pruned=T, write.pruned1.name="Figures/Fig2_046_ML_MOL_extant_pruned.nwk", write.pruned2.name="Figures/Fig2_046_ML_TE_extant_pruned.nwk")
# No. of nodes where BS MOL = BS TE
a = sum(ml_mkv_shared$Support_Tree_1 == ml_mkv_shared$Support_Tree_2, na.rm = TRUE) # number
a
round(a*100 / as.numeric(nrow(ml_mkv_shared)), digits=3) # proportion
# No. of nodes where BS MOL > BS TE
a = sum(ml_mkv_shared$Support_Tree_1 > ml_mkv_shared$Support_Tree_2, na.rm = TRUE) # number
a
round(a*100 / as.numeric(nrow(ml_mkv_shared)), digits=3) # proportion
# No. of nodes where BS MOL < BS TE
a = sum(ml_mkv_shared$Support_Tree_1 < ml_mkv_shared$Support_Tree_2, na.rm = TRUE) # number
a
round(a*100 / as.numeric(nrow(ml_mkv_shared)), digits=3) # proportion# Identify and plot unique nodes
RNODE::uniqueNodes(ml_046_mol, ml_046_te,
                   plotTrees=T, node.numbers=F, tree.width=14, tree.height=17, tree.fsize=0.8, tree.adj=c(-1.5,0.5), tree.cex=2, output.tree="Figures/Fig2_046_ML_MOLvsTE_extant_unique.pdf")


# Trees from 046 (with fossils, ML MKV, Fig. 2D, Table S3)
# Map BS to optimal trees
ml_046_mol_f = read.tree("Trees_extinct//046_MOL_IQTREE.treefile")
ml_046_mol_bs_f = read.tree("Trees_extinct/046_MOL_IQTREE.contree")
RNODE::mapSupport(ml_046_mol_f, ml_046_mol_bs_f,  
                  plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(1.2,-0.5), tree.cex=.8, 
                  tree.output = "Figures/Fig2_046_OPTvsBS_ML_MOL_f.pdf")
ml_046_te_f = read.tree("Trees_extinct//046_TE_ASC_IQTREE.treefile")
ml_046_te_bs_f = read.tree("Trees_extinct/046_TE_ASC_IQTREE.contree")
RNODE::mapSupport(ml_046_te_f, ml_046_te_bs_f,
                  plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(1.2,-0.5), tree.cex=.8, 
                  tree.output = "Figures/Fig2_046_OPTvsBS_ML_TE_f.pdf")
# Identify shared nodes and plot tangle trees
ml_mkv_f_shared = RNODE::sharedNodes(ml_046_mol_bs_f, ml_046_te_bs_f, dataframe=F, root="Morone_saxatilis",
                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=15, tree.fsize=1, tree.adj=c(-1.5,0.5), tree.cex=0.6, 
                   output.tree1="Figures/Fig2_046_ML_MOL_extinct.pdf", output.tree2="Figures/Fig2_046_ML_TE_extinct.pdf",
                   tanglegram = T, output.tangletree = "Figures/Fig2_046_ML_MOLvsTE_extinct_tangle.pdf", tanglegram.lab.cex=0.4, tanglegram.margin=5,
                   write.pruned=T, write.pruned1.name="Figures/Fig2_046_ML_MOL_extinct_pruned.nwk", write.pruned2.name="Figures/Fig2_046_ML_TE_extinct_pruned.nwk")
# No. of nodes where BS MOL = BS TE
a = sum(ml_mkv_shared$Support_Tree_1 == ml_mkv_shared$Support_Tree_2, na.rm = TRUE) # number
a
round(a*100 / as.numeric(nrow(ml_mkv_shared)), digits=3) # proportion
# No. of nodes where BS MOL > BS TE
a = sum(ml_mkv_shared$Support_Tree_1 > ml_mkv_shared$Support_Tree_2, na.rm = TRUE) # number
a
round(a*100 / as.numeric(nrow(ml_mkv_shared)), digits=3) # proportion
# No. of nodes where BS MOL < BS TE
a = sum(ml_mkv_shared$Support_Tree_1 < ml_mkv_shared$Support_Tree_2, na.rm = TRUE) # number
a
round(a*100 / as.numeric(nrow(ml_mkv_shared)), digits=3) # proportion# Identify and plot unique nodes
# Identify and plot unique nodes
RNODE::uniqueNodes(ml_046_mol_f, ml_046_te_f,
                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(-1.5,0.5), tree.cex=2, output.tree="Figures/Fig2_046_ML_MOLvsTE_extinct_unique.pdf")

# Trees from 046 (without fossils, ML MK, Table S3)
# Map BS to optimal trees
ml_046_mol = read.tree("Trees_extant/046_MOL_IQTREE.treefile")
ml_046_mol_bs = read.tree("Trees_extant/046_MOL_IQTREE.contree")
RNODE::mapSupport(ml_046_mol, ml_046_mol_bs,  
                  plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(1.2,-0.5), tree.cex=.8, 
                  tree.output = "Figures/Fig2_046_OPTvsBS_ML_MOL.pdf")
ml_046_te_mk = read.tree("Trees_extant/046_TE_noASC_IQTREE.treefile")
ml_046_te_mk_bs = read.tree("Trees_extant/046_TE_noASC_IQTREE.contree")
RNODE::mapSupport(ml_046_te_mk, ml_046_te_mk_bs,
                  plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(1.2,-0.5), tree.cex=.8, 
                  tree.output = "Figures/Fig2_046_OPTvsBS_ML_TEmk.pdf")
# Identify shared nodes and plot tangle trees
ml_mk_shared = RNODE::sharedNodes(ml_046_mol_bs, ml_046_te_mk_bs, dataframe=F, root="Morone_saxatilis",
                                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=15, tree.fsize=1, tree.adj=c(-1.5,0.5), tree.cex=0.6, 
                                   output.tree1="Figures/Fig2_046_ML_MOL_extant.pdf", output.tree2="Figures/Fig2_046_ML_TE_extant.pdf",
                                   tanglegram = T, output.tangletree = "Figures/Fig2_046_ML_MOLvsTE_extant_tangle.pdf", tanglegram.lab.cex=0.4, tanglegram.margin=5,
                                   write.pruned=T, write.pruned1.name="Figures/Fig2_046_ML_MOL_extant_pruned.nwk", write.pruned2.name="Figures/Fig2_046_ML_TE_extant_pruned.nwk")
# No. of nodes where BS MOL = BS TE
sum(ml_mk_shared$Node_Tree_1 == ml_mk_shared$Node_Tree_2, na.rm = TRUE) # number
round(sum(ml_mk_shared$Node_Tree_1 == ml_mk_shared$Node_Tree_2, na.rm = TRUE)*100 / length(ml_mk_shared$Node_Tree_1), digits=3) # proportion
# No. of nodes where BS MOL > BS TE
sum(ml_mk_shared$Node_Tree_1 > ml_mk_shared$Node_Tree_2, na.rm = TRUE)
round(sum(ml_mk_shared$Node_Tree_1 > ml_mk_shared$Node_Tree_2, na.rm = TRUE)*100 / length(ml_mk_shared$Node_Tree_1), digits=3) # proportion
# No. of nodes where BS MOL < BS TE
sum(ml_mk_shared$Node_Tree_1 < ml_mk_shared$Node_Tree_2, na.rm = TRUE)
round(sum(ml_mk_shared$Node_Tree_1 < ml_mk_shared$Node_Tree_2, na.rm = TRUE)*100 / length(ml_mk_shared$Node_Tree_1), digits=3) # proportion
# Identify and plot unique nodes
RNODE::uniqueNodes(ml_046_mol, ml_046_te,
                   plotTrees=T, node.numbers=F, tree.width=14, tree.height=17, tree.fsize=0.8, tree.adj=c(-1.5,0.5), tree.cex=2, output.tree="Figures/Fig2_046_ML_MOLvsTE_extant_unique.pdf")


# Trees from 046 (with fossils, ML MK, Table S3)
# Map BS to optimal trees
ml_046_mol_f = read.tree("Trees_extinct//046_MOL_IQTREE.treefile")
ml_046_mol_bs_f = read.tree("Trees_extinct/046_MOL_IQTREE.contree")
RNODE::mapSupport(ml_046_mol_f, ml_046_mol_bs_f,  
                  plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(1.2,-0.5), tree.cex=.8, 
                  tree.output = "Figures/Fig2_046_OPTvsBS_ML_MOL_f.pdf")
ml_046_te_mk_f = read.tree("Trees_extinct//046_TE_noASC_IQTREE.treefile")
ml_046_te_mk_bs_f = read.tree("Trees_extinct/046_TE_noASC_IQTREE.contree")
RNODE::mapSupport(ml_046_te_mk_f, ml_046_te_mk_bs_f,
                  plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(1.2,-0.5), tree.cex=.8, 
                  tree.output = "Figures/Fig2_046_OPTvsBS_ML_TE_f.pdf")
# Identify shared nodes and plot tangle trees
ml_mk_f_shared = RNODE::sharedNodes(ml_046_mol_bs_f, ml_046_te_mk_bs_f, dataframe=F, root="Morone_saxatilis",
                                     plotTrees=T, node.numbers=F, tree.width=12, tree.height=15, tree.fsize=1, tree.adj=c(-1.5,0.5), tree.cex=0.6, 
                                     output.tree1="Figures/Fig2_046_ML_MOL_extinct.pdf", output.tree2="Figures/Fig2_046_ML_TE_extinct.pdf",
                                     tanglegram = T, output.tangletree = "Figures/Fig2_046_ML_MOLvsTE_extinct_tangle.pdf", tanglegram.lab.cex=0.4, tanglegram.margin=5,
                                     write.pruned=T, write.pruned1.name="Figures/Fig2_046_ML_MOL_extinct_pruned.nwk", write.pruned2.name="Figures/Fig2_046_ML_TE_extinct_pruned.nwk")
# No. of nodes where BS MOL = BS TE
sum(ml_mk_f_shared$Node_Tree_1 == ml_mk_f_shared$Node_Tree_2, na.rm = TRUE) # number
round(sum(ml_mk_f_shared$Node_Tree_1 == ml_mk_f_shared$Node_Tree_2, na.rm = TRUE)*100 / length(ml_mk_f_shared$Node_Tree_1), digits=3) # proportion
# No. of nodes where BS MOL > BS TE
sum(ml_mk_f_shared$Node_Tree_1 > ml_mk_f_shared$Node_Tree_2, na.rm = TRUE)
round(sum(ml_mk_f_shared$Node_Tree_1 > ml_mk_f_shared$Node_Tree_2, na.rm = TRUE)*100 / length(ml_mk_f_shared$Node_Tree_1), digits=3) # proportion
# No. of nodes where BS MOL < BS TE
sum(ml_mk_f_shared$Node_Tree_1 < ml_mk_f_shared$Node_Tree_2, na.rm = TRUE)
round(sum(ml_mk_f_shared$Node_Tree_1 < ml_mk_f_shared$Node_Tree_2, na.rm = TRUE)*100 / length(ml_mk_f_shared$Node_Tree_1), digits=3) # proportion
# Identify and plot unique nodes
RNODE::uniqueNodes(ml_046_mol_f, ml_046_te_f,
                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(-1.5,0.5), tree.cex=2, output.tree="Figures/Fig2_046_ML_MOLvsTE_extinct_unique.pdf")

############
# PLOTTING #
############

# Add Analysis column to each dataframe
mp_shared$Analysis       <- "mp_shared"
mp_f_shared$Analysis     <- "mp_f_shared"
ml_mkv_shared$Analysis   <- "ml_mkv_shared"
ml_mkv_f_shared$Analysis <- "ml_mkv_f_shared"
ml_mk_shared$Analysis    <- "ml_mk_shared"
ml_mk_f_shared$Analysis  <- "ml_mk_f_shared"

# Merge all into one dataframe
df <- rbind(
  mp_shared,
  mp_f_shared,
  ml_mkv_shared,
  ml_mkv_f_shared,
  ml_mk_shared,
  ml_mk_f_shared
)

library(dplyr)
library(tidyr)
library(ggplot2)

# Clean and reshape
df_clean <- df %>%
  filter(
    !is.na(Support_Tree_1),
    !is.na(Support_Tree_2),
    Support_Tree_1 != "",
    Support_Tree_2 != ""
  ) %>%
  mutate(
    Support_Tree_1 = as.numeric(Support_Tree_1),
    Support_Tree_2 = as.numeric(Support_Tree_2),
    Analysis = recode(
      Analysis,
      "ml_mk_f_shared"   = "ML MK (with fossils)",
      "ml_mk_shared"     = "ML MK (without fossils)",
      "ml_mkv_f_shared"  = "ML MKv (with fossils)",
      "ml_mkv_shared"    = "ML MKv (without fossils)",
      "mp_f_shared"      = "MP (with fossils)",
      "mp_shared"        = "MP (without fossils)"
    )
  )

# 2. Reshape to long format
df_long <- df_clean %>%
  pivot_longer(
    cols = c(Support_Tree_1, Support_Tree_2),
    names_to = "Tree",      # <- This becomes the new column name
    values_to = "Support"
  )

# Check structure
head(df_long)

# 3. Plot a normal boxplot
ggplot(df_long, aes(x = Tree, y = Support)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  facet_wrap(~Analysis) +
  theme_minimal()

# Density plots with fixed axes and custom legend
ggplot(df_long, aes(x = Support, fill = Tree, color = Tree)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~Analysis) +
  theme_minimal() +
  labs(title = "", x = "\nSupport", y = "Density\n") +
  coord_cartesian(xlim = c(50, 100), ylim = c(0, 1)) +
  scale_fill_manual(
    values = c("Support_Tree_1" = "blue", "Support_Tree_2" = "red"),
    labels = c("Support_Tree_1" = "MOL BS", "Support_Tree_2" = "TE BS")
  ) +
  scale_color_manual(
    values = c("Support_Tree_1" = "blue", "Support_Tree_2" = "red"),
    labels = c("Support_Tree_1" = "MOL BS", "Support_Tree_2" = "TE BS")
  )

# Histogram
jpeg("Figures/046_support.jpg", width=6, height=6, units="in",res=300)
ggplot(df_long, aes(x = Support, fill = Tree)) +
  geom_histogram(binwidth = 2, color = "black", alpha = 0.2, position = "identity") +
  facet_wrap(~Analysis) +
  scale_fill_manual(
    values = c("Support_Tree_1" = "blue", "Support_Tree_2" = "red"),
    labels = c("Support_Tree_1" = "MOL BS", "Support_Tree_2" = "TE BS")
  ) +
  theme_minimal() +
  labs(title = "", x = "\nSupport", y = "Count\n") +
  coord_cartesian(xlim = c(50, 100), ylim = c(0, 20))
dev.off() 

# Exploratory data
library(dplyr)
summary_stats <- df_long %>%
  group_by(Analysis, Tree) %>%
  summarise(
    mean_support = mean(Support, na.rm = TRUE),
    min_support  = min(Support, na.rm = TRUE),
    max_support  = max(Support, na.rm = TRUE),
    .groups = "drop"
  )
summary_stats

library(purrr)

# Function to run paired Wilcoxon test for each Analysis
paired_tests <- df_long %>%
  group_by(Analysis) %>%
  summarise(
    test = list(
      wilcox.test(
        Support ~ Tree,
        data = cur_data(),
        paired = TRUE,
        exact = FALSE
      )
    ),
    .groups = "drop"
  ) %>%
  mutate(
    W = map_dbl(test, ~ .$statistic),
    p_value = map_dbl(test, ~ .$p.value)
  ) %>%
  select(Analysis, W, p_value)
paired_tests
