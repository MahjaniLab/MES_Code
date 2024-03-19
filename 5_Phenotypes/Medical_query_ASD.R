# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Testing comorbidities for MES gene carriers with ASD

options(tidyverse.quiet = TRUE)
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

#MES genes
GENES = read.delim(file = "MES_genes.txt", header = FALSE)

# SPARK metadata
FAM = read.delim(file = "SPARK.iWES_v1.mastertable.2022_02.tsv", header = TRUE)
unaffected = c(FAM$spid[FAM$asd == 1])
affected = c(FAM$spid[FAM$asd == 2])

# SPARK medical background
PHENO = read.csv(file = "basic_medical_screening-2022-06-03.csv", header = TRUE)
NAMES = PHENO$subject_sp_id

QUES = read.delim(file = "basic_med_questions.txt", header = FALSE)

samp_w_one_MES = vector()
for (gene in 1:length(GENES$V1)){
  infile = paste0("intermediates/", GENES$V1[gene], ".", GENES$V2[gene], ".cap")
  INFILE = read.delim(infile, header = FALSE)
  samp_w_one_MES = c(samp_w_one_MES, INFILE$V2)
}
# Removes samples with at least one MES gene
affected_no_MES = affected[!(affected %in% samp_w_one_MES)]

for (gene in 1:length(GENES$V1)){
  infile = paste0("intermediates/", GENES$V1[gene], ".", GENES$V2[gene], ".cap")
  INFILE = read.delim(infile, header = FALSE)

  ASD_W_GENE = INFILE$V2[(INFILE$V2 %in% unaffected) & (INFILE$V1 == GENES$V1[gene])]
  ASD_WO_GENE = affected_no_MES[!(affected_no_MES %in% ASD_W_GENE)]

  for (column in QUES$V1){
    PHENO_sub = as.vector(unlist(PHENO %>% select(all_of(column))))
    matrix_chr = matrix(nrow = 2, ncol = 2)

    matrix_chr[1,1] = sum(PHENO_sub[NAMES %in% ASD_W_GENE] == 1, na.rm = TRUE)
    matrix_chr[1,2] = length(PHENO_sub[NAMES %in% ASD_W_GENE]) - sum(PHENO_sub[NAMES %in% ASD_W_GENE] == 1, na.rm = TRUE)
    matrix_chr[2,1] = sum(PHENO_sub[NAMES %in% ASD_WO_GENE] == 1, na.rm = TRUE)
    matrix_chr[2,2] = length(PHENO_sub[NAMES %in% ASD_WO_GENE]) - sum(PHENO_sub[NAMES %in% ASD_WO_GENE] == 1, na.rm = TRUE)

    A <- chisq.test(matrix_chr, simulate.p.value = TRUE)

    cat(GENES$V1[gene], "\t", column, "\t", A$p.value, "\t", 
	matrix_chr[1,1], "\t", matrix_chr[1,2], "\t", matrix_chr[2,1], "\t", matrix_chr[2,2], "\n", sep = "")
  }

 ASD_W_GENE = INFILE$V2[(INFILE$V2 %in% affected) & (INFILE$V1 == GENES$V2[gene])]
 ASD_WO_GENE = affected_no_MES[!(affected_no_MES %in% ASD_W_GENE)]

  for (column in QUES$V1){
    PHENO_sub = as.vector(unlist(PHENO %>% select(all_of(column))))
    matrix_chr = matrix(nrow = 2, ncol = 2)

  matrix_chr[1,1] = sum(PHENO_sub[NAMES %in% ASD_W_GENE] == 1, na.rm = TRUE)
    matrix_chr[1,2] = length(PHENO_sub[NAMES %in% ASD_W_GENE]) - sum(PHENO_sub[NAMES %in% ASD_W_GENE] == 1, na.rm = TRUE)
    matrix_chr[2,1] = sum(PHENO_sub[NAMES %in% ASD_WO_GENE] == 1, na.rm = TRUE)
    matrix_chr[2,2] = length(PHENO_sub[NAMES %in% ASD_WO_GENE]) - sum(PHENO_sub[NAMES %in% ASD_WO_GENE] == 1, na.rm = TRUE)

    A <- chisq.test(matrix_chr, simulate.p.value = TRUE)

    cat(GENES$V2[gene], "\t", column, "\t", A$p.value, "\t",
        matrix_chr[1,1], "\t", matrix_chr[1,2], "\t", matrix_chr[2,1], "\t", matrix_chr[2,2], "\n", sep = "")
  }

}
