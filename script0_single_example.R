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

# Trees from 046 (without fossils, MP)
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

RNODE::sharedNodes(mp_046_mol, mp_046_te, dataframe=F, root="Morone_saxatilis",
                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=15, tree.fsize=1, tree.adj=c(-1.5,0.5), tree.cex=0.6, 
                   output.tree1="Figures/Fig2_046_MP_MOL_extant.pdf", output.tree2="Figures/Fig2_046_MP_TE_extant.pdf",
                   tanglegram = T, output.tangletree = "Figures/Fig2_046_MP_MOLvsTE_extant_tangle.pdf", tanglegram.lab.cex=0.4, tanglegram.margin=5,
                   write.pruned=T, write.pruned1.name="Figures/Fig2_046_MP_MOL_extant_pruned.nwk", write.pruned2.name="Figures/Fig2_046_MP_TE_extant_pruned.nwk")
RNODE::uniqueNodes(mp_046_mol, mp_046_te,
                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(-1.5,0.5), tree.cex=2, output.tree="Figures/Fig2_046_MP_MOLvsTE_extant_unique.pdf")

# Trees from 046 (without fossils, ML)
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

RNODE::sharedNodes(ml_046_mol, ml_046_te, dataframe=F, root="Morone_saxatilis",
                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=15, tree.fsize=1, tree.adj=c(-1.5,0.5), tree.cex=0.6, 
                   output.tree1="Figures/Fig2_046_ML_MOL_extant.pdf", output.tree2="Figures/Fig2_046_ML_TE_extant.pdf",
                   tanglegram = T, output.tangletree = "Figures/Fig2_046_ML_MOLvsTE_extant_tangle.pdf", tanglegram.lab.cex=0.4, tanglegram.margin=5,
                   write.pruned=T, write.pruned1.name="Figures/Fig2_046_ML_MOL_extant_pruned.nwk", write.pruned2.name="Figures/Fig2_046_ML_TE_extant_pruned.nwk")
RNODE::uniqueNodes(ml_046_mol, ml_046_te,
                   plotTrees=T, node.numbers=F, tree.width=14, tree.height=17, tree.fsize=0.8, tree.adj=c(-1.5,0.5), tree.cex=2, output.tree="Figures/Fig2_046_ML_MOLvsTE_extant_unique.pdf")

# Trees from 046 (with fossils, MP)
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

RNODE::sharedNodes(mp_046_mol_f, mp_046_te_f, dataframe=F, root="Morone_saxatilis",
                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=15, tree.fsize=1, tree.adj=c(-1.5,0.5), tree.cex=0.6, 
                   output.tree1="Figures/Fig2_046_MP_MOL_extinct.pdf", output.tree2="Figures/Fig2_046_MP_TE_extinct.pdf",
                   tanglegram = T, output.tangletree = "Figures/Fig2_046_MP_MOLvsTE_extinct_tangle.pdf", tanglegram.lab.cex=0.4, tanglegram.margin=5,
                   write.pruned=T, write.pruned1.name="Figures/Fig2_046_MP_MOL_extinct_pruned.nwk", write.pruned2.name="Figures/Fig2_046_MP_TE_extinct_pruned.nwk")
RNODE::uniqueNodes(mp_046_mol_f, mp_046_te_f,
                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(-1.5,0.5), tree.cex=2, output.tree="Figures/Fig2_046_MP_MOLvsTE_extinct_unique.pdf")


# Trees from 046 (with fossils, ML)
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

RNODE::sharedNodes(ml_046_mol_f, ml_046_te_f, dataframe=F, root="Morone_saxatilis",
                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=15, tree.fsize=1, tree.adj=c(-1.5,0.5), tree.cex=0.6, 
                   output.tree1="Figures/Fig2_046_ML_MOL_extinct.pdf", output.tree2="Figures/Fig2_046_ML_TE_extinct.pdf",
                   tanglegram = T, output.tangletree = "Figures/Fig2_046_ML_MOLvsTE_extinct_tangle.pdf", tanglegram.lab.cex=0.4, tanglegram.margin=5,
                   write.pruned=T, write.pruned1.name="Figures/Fig2_046_ML_MOL_extinct_pruned.nwk", write.pruned2.name="Figures/Fig2_046_ML_TE_extinct_pruned.nwk")
RNODE::uniqueNodes(ml_046_mol_f, ml_046_te_f,
                   plotTrees=T, node.numbers=F, tree.width=12, tree.height=17, tree.fsize=0.8, tree.adj=c(-1.5,0.5), tree.cex=2, output.tree="Figures/Fig2_046_ML_MOLvsTE_extinct_unique.pdf")
