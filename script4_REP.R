############
# PREAMBLE #
############

# Set working directory
setwd("/Users/labanfibios/Desktop/Doutorado/Project/B2_TEvsMOL/GitHub/")

# Load packages
library(dplyr)
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

# Unify data into a dataframe
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

# Delete empty rows or with question marks
df <- df[rowSums(df == "?" | df == "", na.rm = TRUE) == 0, ]

########################
# STATISTICAL ANALYSES #
########################

# Correlation between BS and REP

plot(df$Support_Tree_1, df$Support_Tree_2)

# Question 4: BS MOL vs BS TE 

# Question 5

