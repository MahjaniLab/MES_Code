# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Calculates % of genes are expressed above baseline per tissue

args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
`%ni%` <- Negate(`%in%`)

# MES genes input:
genes = read.delim(file = "passed_genes.txt", header = FALSE)
MES_genes = genes$V1

# Tissue names (removed spaces)
tissues = read.delim(file = "all_tissues.txt", header = FALSE)
nsamp = length(MES_genes)

# From GTEx, gene expression per gene
EXP = read.delim(file = "rna_tissue_consensus.tsv", header = TRUE)
all_genes = unique(EXP$Genename)

for (region in tissues$V1){    
  #region = "cerebralcortex"
  #region = "adrenalgland"

  All_genes_expression = as.numeric(EXP$nTPM[(EXP$Tissue == region)])
  MES_genes_expression = as.numeric(EXP$nTPM[(EXP$Tissue == region) & (EXP$Genename %in% MES_genes)])

  cat(region, "\t", 
	sum(All_genes_expression>1)/length(All_genes_expression), "\t",
	sum(MES_genes_expression>1)/length(MES_genes_expression), "\t", 
        (sum(MES_genes_expression>1)/length(MES_genes_expression)) - (sum(All_genes_expression>1)/length(All_genes_expression)), "\t")

  ## Samples and calculates p-values
   expressions_samps = vector(length = 1000)
   # 1000 iterations of sampling
   for (i in 1:1000){
     samp_genes = sample(x=all_genes, size=nsamp, replace=FALSE) 
     samp_genes_expression = as.numeric(EXP$nTPM[(EXP$Tissue == region) & (EXP$Genename %in% samp_genes)])
     expressions_samps[i] = sum(samp_genes_expression>1)/length(samp_genes_expression);
   }

   p_value <- pnorm(sum(MES_genes_expression>1)/length(MES_genes_expression), mean = mean(expressions_samps), sd = sd(expressions_samps))
   # P-value is like a percentile, so abs 1-value gives a more proper "p-value"
   cat(abs(1 - p_value), "\t",  mean(expressions_samps), "\t",  sd(expressions_samps), "\n")

}

