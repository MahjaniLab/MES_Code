# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Removing variants that are Mendelian errors.

library(stringr)
args = commandArgs(trailingOnly=TRUE)
chr = args[1]

# SPARK families (a variation on a ped file). Header:
# family_ID  father  mother  siblings        sibling_sexes   sibling_asd
FAM <- read.delim(file = "families.txt", header = TRUE)
system(paste0("rm ", chr, ".mendelian.report.txt"))

for (f in FAM$family){

  # Individual level variant information grouped by family. Header:
  #CHROM   POS     REF     ALT     ChildX       ChildY       ParentA       ParentB
  #chr1    942451  T       C       1/1     1/1     1/1     1/1

  input_file <- paste0(chr, "/", f, ".var")
  data <- read.delim(input_file, header = TRUE, stringsAsFactors = FALSE)

  # The last two columns are the parents
  child_cols <- (ncol(data) - 2):1  

  mendelian_errors <- data.frame()
  non_errors <- data.frame()

  for (i in 1:nrow(data)) {
    children_genotypes <- data[i, child_cols]
    parent1_genotype <- data[i, ncol(data) - 1]
    parent2_genotype <- data[i, ncol(data)]
  
    if (any(children_genotypes == "1/1") && parent1_genotype == "0/0" && parent2_genotype == "0/0") {
      mendelian_errors <- rbind(mendelian_errors, data[i, ])
    } else if (any(children_genotypes == "0/1") && parent1_genotype == "0/0" && parent2_genotype == "0/0") {
      mendelian_errors <- rbind(mendelian_errors, data[i, ])
    } else if (any(children_genotypes == "0/0") && parent1_genotype == "1/1" && parent2_genotype == "1/1") {
      mendelian_errors <- rbind(mendelian_errors, data[i, ])
    } else if (any(children_genotypes == "0/1") && parent1_genotype == "1/1" && parent2_genotype == "1/1") {
      mendelian_errors <- rbind(mendelian_errors, data[i, ])
    } else if (any(children_genotypes == "1/1") && parent1_genotype == "0/0" && parent2_genotype == "0/1") {
      mendelian_errors <- rbind(mendelian_errors, data[i, ])
    } else if (any(children_genotypes == "1/1") && parent1_genotype == "0/1" && parent2_genotype == "0/0") {
      mendelian_errors <- rbind(mendelian_errors, data[i, ])
    } else if (any(children_genotypes == "0/0") && parent1_genotype == "1/1" && parent2_genotype == "0/1") {
      mendelian_errors <- rbind(mendelian_errors, data[i, ])
    } else if (any(children_genotypes == "0/0") && parent1_genotype == "0/1" && parent2_genotype == "1/1") {
      mendelian_errors <- rbind(mendelian_errors, data[i, ])
    } else {
    # Store rows without Mendelian errors in the non_errors data frame
    non_errors <- rbind(non_errors, data[i, ])
    }
  }

  mendelian_errors_file <- paste0(chr, "/", f, ".menderrors.var")
  write.table(mendelian_errors, mendelian_errors_file, sep = "\t", row.names = FALSE, quote = FALSE)

  # Write the non-errors to a separate file
  non_errors_file <- paste0(chr, "/", f, ".nonmenderrors.var")
  write.table(non_errors, non_errors_file, sep = "\t", row.names = FALSE, quote = FALSE)

  cat(f, "\t", chr, "\t", nrow(mendelian_errors)/nrow(data), "\n", sep = "", file = paste0(chr, ".mendelian.report.txt"), append = TRUE)
}
