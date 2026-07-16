############
# PREAMBLE #
############

# Set working directory
setwd("/Users/labanfibios/Desktop/Doutorado/Project/B2_TEvsMOL/GitHub/")

# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(RNODE)

#############
# LOAD DATA #
#############

# Load bootstrap trees
MOL_008_BS = read.tree("Trees_extant/008b_MOL_BS_TNT.nwk")
TE_008_BS = read.tree("Trees_extant/008d_TE_BS_TNT.nwk")
MOL_018_BS = read.tree("Trees_extant/018b_MOL_BS_TNT.nwk")
TE_018_BS = read.tree("Trees_extant/018d_TE_BS_TNT.nwk")
MOL_046_BS = read.tree("Trees_extant/046b_MOL_BS_TNT.nwk")
TE_046_BS = read.tree("Trees_extant/046d_TE_BS_TNT.nwk")
MOL_057_BS = read.tree("Trees_extant/057b_MOL_BS_TNT.nwk")
TE_057_BS = read.tree("Trees_extant/057d_TE_BS_TNT.nwk")
MOL_059_BS = read.tree("Trees_extant/059b_MOL_BS_TNT.nwk")
TE_059_BS = read.tree("Trees_extant/059d_TE_BS_TNT.nwk")


# Compute REP
setwd("/Users/labanfibios/Desktop/Doutorado/Project/B2_TEvsMOL/GitHub/GB_trees/")
print("The lenght of the worst binary tree computed in POY using negative weights: 
        - 008_MOL = 20592
        - 008_MORPH = 876
        - 008_TE = 21339
        - 018_MOL = 74934
        - 018_MORPH = 687
        - 018_TE =  74267
        - 046_MOL = 9532086
        - 046_MORPH = 1376
        - 046_TE = 9354305
        - 057_MOL = 409081
        - 057_MORPH = 1027
        - 057_TE = 409560
        - 059_MOL = 212953
        - 059_MORPH = 135
        - 059_TE = 211159")

MOL_008 = read.tree("008gb_ab_trees_MOL_TNT_results.nwk")[[2]]
MOL_008_REP = addREP(tree=MOL_008, best_tree_length = 15355, worst_tree_length = 20592)
TE_008 = read.tree("008gb_cd_trees_TE_TNT_results.nwk")[[2]]
TE_008_REP = addREP(tree=TE_008, best_tree_length = 15753, worst_tree_length = 21339)

MOL_018 = read.tree("018gb_ab_trees_MOL_TNT_results.nwk")[[2]]
MOL_018_REP = addREP(tree=MOL_018, best_tree_length = 9252, worst_tree_length = 74934)
TE_018 = read.tree("018gb_cd_trees_TE_TNT_results.nwk")[[2]]
TE_018_REP = addREP(tree=TE_018, best_tree_length = 9677, worst_tree_length = 74267)

MOL_046 = read.tree("046gb_ab_trees_MOL_TNT_results.nwk")[[2]]
MOL_046_REP = addREP(tree=MOL_046, best_tree_length = 435512, worst_tree_length = 9532086)
TE_046 = read.tree("046gb_cd_trees_TE_TNT_results.nwk")[[2]]
TE_046_REP = addREP(tree=TE_046, best_tree_length = 435773, worst_tree_length = 9354305)

MOL_057 = read.tree("057gb_ab_trees_MOL_TNT_results.nwk")[[2]]
MOL_057_REP = addREP(tree=MOL_057, best_tree_length = 94454, worst_tree_length = 409081)
TE_057 = read.tree("057gb_cd_trees_TE_TNT_results.nwk")[[2]]
TE_057_REP = addREP(tree=TE_057, best_tree_length = 94621, worst_tree_length = 409560)

MOL_059 = read.tree("059gb_ab_trees_MOL_TNT_results.nwk")[[2]]
MOL_059_REP = addREP(tree=MOL_059, best_tree_length = 89062, worst_tree_length = 212953)
TE_059 = read.tree("059gb_cd_trees_TE_TNT_results.nwk")[[2]]
TE_059_REP = addREP(tree=TE_059, best_tree_length=89155, worst_tree_length = 211159)

# Example of plotting REP
plot(ladderize(MOL_008_REP))
nodelabels(MOL_008_REP$node.label,
           frame = "none",
           cex = 0.7,
           adj = c(1.2, -0.2))

# Unify BS and REP data into a dataframe
df_MOL_008 = sharedNodes(MOL_008_BS, MOL_008_REP)
df_TE_008 = sharedNodes(TE_008_BS, TE_008_REP)
df_MOL_018 = sharedNodes(MOL_018_BS, MOL_018_REP)
df_TE_018 = sharedNodes(TE_018_BS, TE_018_REP)
df_MOL_046 = sharedNodes(MOL_046_BS, MOL_046_REP)
df_TE_046 = sharedNodes(TE_046_BS, TE_046_REP)
df_MOL_057 = sharedNodes(MOL_057_BS, MOL_057_REP)
df_TE_057 = sharedNodes(TE_057_BS, TE_057_REP)
df_MOL_059 = sharedNodes(MOL_059_BS, MOL_059_REP)
df_TE_059 = sharedNodes(TE_059_BS, TE_059_REP)
df <- bind_rows(
  mutate(df_MOL_008, dataset = "MOL", analysis = "008"),
  mutate(df_TE_008,  dataset = "TE",  analysis = "008"),
  mutate(df_MOL_018, dataset = "MOL", analysis = "018"),
  mutate(df_TE_018,  dataset = "TE",  analysis = "018"),
  mutate(df_MOL_046, dataset = "MOL", analysis = "046"),
  mutate(df_TE_046,  dataset = "TE",  analysis = "046"),
  mutate(df_MOL_057, dataset = "MOL", analysis = "057"),
  mutate(df_TE_057,  dataset = "TE",  analysis = "057"),
  mutate(df_MOL_059, dataset = "MOL", analysis = "059"),
  mutate(df_TE_059,  dataset = "TE",  analysis = "059")
)

# Delete empty rows or with question marks (trivial support values from TNT)
df <- df[rowSums(df == "?" | df == "", na.rm = TRUE) == 0, ]
df$Support_Tree_1 = as.numeric(df$Support_Tree_1)
df$Support_Tree_2 = as.numeric(df$Support_Tree_2)
df <- na.omit(df)

# Unify MOL and TE data into a second dataframe
df_BS_008 = sharedNodes(MOL_008_BS, TE_008_BS)
df_REP_008 = sharedNodes(MOL_008_REP, TE_008_REP)
df_BS_018 = sharedNodes(MOL_018_BS, TE_018_BS)
df_REP_018 = sharedNodes(MOL_018_REP, TE_018_REP)
df_BS_046 = sharedNodes(MOL_046_BS, TE_046_BS)
df_REP_046 = sharedNodes(MOL_046_REP, TE_046_REP)
df_BS_057 = sharedNodes(MOL_057_BS, TE_057_BS)
df_REP_057 = sharedNodes(MOL_057_REP, TE_057_REP)
df_BS_059 = sharedNodes(MOL_059_BS, TE_059_BS)
df_REP_059 = sharedNodes(MOL_059_REP, TE_059_REP)
df2 <- bind_rows(
  mutate(df_BS_008, support = "Bootstrap", analysis = "008"),
  mutate(df_REP_008,  support = "REP",  analysis = "008"),
  mutate(df_BS_018, support = "Bootstrap", analysis = "018"),
  mutate(df_REP_018,  support = "REP",  analysis = "018"),
  mutate(df_BS_046, support = "Bootstrap", analysis = "046"),
  mutate(df_REP_046,  support = "REP",  analysis = "046"),
  mutate(df_BS_057, support = "Bootstrap", analysis = "057"),
  mutate(df_REP_057,  support = "REP",  analysis = "057"),
  mutate(df_BS_059, support = "Bootstrap", analysis = "059"),
  mutate(df_REP_059,  support = "REP",  analysis = "059")
)
# Delete empty rows or with question marks (trivial support values from TNT)
df2 <- df2[rowSums(df2 == "?" | df2 == "", na.rm = TRUE) == 0, ]
df2$Support_Tree_1 = as.numeric(df2$Support_Tree_1)
df2$Support_Tree_2 = as.numeric(df2$Support_Tree_2)
df2 <- na.omit(df2)

##############
# QUESTION 4 #
##############

# Histograms
hist_bs_mol = ggplot(subset(df, dataset == "MOL"), aes(Support_Tree_1)) +
  geom_histogram(
    bins = 25,
    fill = "grey40",
    color = "white",
    linewidth = 0.2
  ) +
  labs(
    title = "MOL",
    x = "\nBootstrap",
    y = "Frequency\n"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "plain"), 
    plot.title = element_text(
      face = "plain",
      hjust = 0.5)
  )

hist_bs_te = ggplot(subset(df, dataset == "TE"), aes(Support_Tree_1)) +
  geom_histogram(
    bins = 25,
    fill = "black",
    color = "white",
    linewidth = 0.2
  ) +
  labs(
    title = "TE",
    x = "\nBootstrap",
    y = "Frequency\n"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "plain"), 
    plot.title = element_text(
      face = "plain",
      hjust = 0.5)
  )

hist_rep_mol = ggplot(subset(df, dataset == "MOL"), aes(Support_Tree_2)) +
  geom_histogram(
    bins = 25,
    fill = "grey40",
    color = "white",
    linewidth = 0.2
  ) +
  labs(
    title = "MOL",
    x = "\nREP",
    y = "Frequency\n"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "plain"), 
    plot.title = element_text(
      face = "plain",
      hjust = 0.5)
  )

hist_rep_te = ggplot(subset(df, dataset == "TE"), aes(Support_Tree_2)) +
  geom_histogram(
    bins = 25,
    fill = "black",
    color = "white",
    linewidth = 0.2
  ) +
  labs(
    title = "TE",
    x = "\nREP",
    y = "Frequency\n"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "plain"), 
    plot.title = element_text(
      face = "plain",
      hjust = 0.5)
  )
combined_plot <-
  (hist_rep_mol + hist_rep_te) /
  (hist_bs_mol + hist_bs_te) +
  plot_annotation(tag_levels = "a")
combined_plot
ggsave(
  "../Figures/figBS_REP_histograms.jpg",
  plot = combined_plot,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

# Correlation between BS and REP
BSvsREP = ggplot(df, aes(
  x = Support_Tree_1,
  y = Support_Tree_2,
  color = dataset,
)) +
  geom_point(size = 3.5, alpha = 0.8) +
  geom_smooth(
    aes(group = dataset),
    method = "lm",
    formula = y ~ poly(x, 2),
    se = T,
    linewidth = 1
  ) +
  scale_color_manual(values = c(MOL = "#1B9E77", TE = "#D95F02")) +
  labs(
    x = "Bootstrap support (%)",
    y = "REP support (%)",
    color = "Dataset") +
  theme_classic(base_size = 14)
ggsave(
  "../Figures/figBSvsREP.jpg",
  plot = BSvsREP,
  width = 6,
  height = 4,
  units = "in",
  dpi = 600
)

# Correlation between MOL and TE
# Bootstrap spearman
df_bs <- subset(df2, support == "Bootstrap")
cor_bs <- cor.test(df_bs$Support_Tree_1, df_bs$Support_Tree_2,
                   method = "spearman")
cor_bs
# REP spearman
df_rep <- subset(df2, support == "REP")
cor_rep <- cor.test(df_rep$Support_Tree_1, df_rep$Support_Tree_2,
                    method = "spearman")
cor_rep
# Plot
bs <- ggplot(subset(df2, support == "Bootstrap"),
             aes(x = Support_Tree_1, y = Support_Tree_2)) +
  geom_point(size = 3.5, alpha = 0.8, color = "#4E79A7") +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 1,
    color = "#4E79A7",
    fill = "#4E79A7",
    alpha = 0.2
  ) +
  annotate(
    "text",
    x = 71, y = 103,
    hjust = 1.05, vjust = 1.5,
    label = sprintf("Spearman's \u03c1 = %.2f\nP = %.3g",
                    cor_bs$estimate,
                    cor_bs$p.value),
    size = 3
  ) +
  labs(
    x = "\nMOL support (BS)",
    y = "TE support (BS)\n"
  ) 
rep <- ggplot(subset(df2, support == "REP"),
              aes(x = Support_Tree_1, y = Support_Tree_2)) +
  geom_point(size = 3.5, alpha = 0.8, color = "#E07A5F") +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 1,
    color = "#E07A5F",
    fill = "#E07A5F",
    alpha = 0.2
  ) +
  annotate(
    "text",
    x = 0.015, y = 0.04,
    hjust = 1.05, vjust = 1.5,
    label = sprintf("Spearman's \u03c1 = %.2f\nP = %.3g",
                    cor_rep$estimate,
                    cor_rep$p.value),
    size = 3
  ) +
  labs(
    x = "\nMOL support (REP)",
    y = "TE support (REP)\n"
  ) 
combined_plot2 <-
  (bs + rep) +
  plot_annotation(tag_levels = "a")
combined_plot2
ggsave(
  "../Figures/figBS_REP_correlationsMOLvsTE.jpg",
  plot = combined_plot2,
  width = 8,
  height = 4,
  units = "in",
  dpi = 600
)

# Perform the Wilcoxon signed-rank test for paired values in Bootstrap
df_bootstrap <- subset(df2, support == "Bootstrap")
wilcox.test(as.numeric(df_bootstrap$Support_Tree_1),
            as.numeric(df_bootstrap$Support_Tree_2),
            paired = TRUE)

# Perform the Wilcoxon signed-rank test for paired values in REP
df_rep <- subset(df2, support == "REP")
wilcox.test(as.numeric(df_rep$Support_Tree_1),
            as.numeric(df_rep$Support_Tree_2),
            paired = TRUE)

#################################################
# Question 5: Support MOL vs Clade occurence TE #
#################################################


