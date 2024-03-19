# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Parsing a variant file and reorganizing information into families.

library(data.table)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

# Arg is family (Father_Mother) id
system(paste0("rm intermediates/", args[1], "_valid_probands.txt"))

# Output header
cat("Variant", "Gene", "Consequence", "Family", "FatherID", "MotherID", 
  "Paternal_genotype", "Maternal_genotype",
  "N_offspring", "Offspring_IDs", 
  "N_probands", "Probands", "N_siblings", "Siblings", 
  "N_probands_w_variant", "Probands_w_variant", "Probands_w_variant_genotypes", 
  "N_probands_wo_variant", "Probands_wo_variant", "Probands_wo_variant_genotypes",
  "N_siblings_w_variant", "Siblings_w_variant", "Siblings_w_variant_genotypes", 
  "N_siblings_wo_variant", "Siblings_wo_variant", "Siblings_wo_variant_genotypes",
  "\n", sep = "\t", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))

# SPARK families with all members WES
FAM<- read.delim(file = "families.txt", header = TRUE)

for (f in 1:length(FAM$family)){ 
  if(FAM$family[f] == args[1]){
    mother = FAM$mother[f]
    father = FAM$father[f]
    family = FAM$family[f]
    siblings = unlist(str_split(FAM$siblings[f], ","))
    sibling_status = unlist(str_split(FAM$sibling_asd[f], ","))

    # Variants filtered to not have Mendelian errors
    IN = read.delim(file = paste0("family, ".nonmenderrors.var"), header = TRUE)

    for (var in 1:length(IN$CHROM)){
      variant = paste0(gsub("chr", "", IN$CHROM[var]), ":", IN$POS[var], ":", IN$REF[var], ":", IN$ALT[var])    

      cat(variant, "\t", ".", "\t", ".", "\t", family, "\t", father, "\t", mother, "\t", sep = "", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))
      cat(IN[var,colnames(IN)==father], "\t", IN[var,colnames(IN)==mother], "\t", length(siblings), "\t", sep = "", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))    
      cat(siblings, sep = ",", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))
      cat("\t", length(siblings[sibling_status == 2]), "\t", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))
      cat(siblings[sibling_status == 2], sep = ",", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))
      cat("\t", length(siblings[sibling_status == 1]), "\t", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))
      cat(siblings[sibling_status == 1], sep = ",", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))

      probands = siblings[sibling_status == 2]
      non_probands = siblings[sibling_status == 1]

      m = match(probands, colnames(IN))
      proband_geno = unlist(IN[var,m])
      proband_names = colnames(IN)[m]

      proband_geno_het = (proband_geno == "0/1") | (proband_geno == "1/1")

      # N_probands_w_variant
      cat("\t", length(proband_geno[proband_geno_het]), "\t", sep = "", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))
  
      # Probands_w_variant
      cat(proband_names[proband_geno_het], sep = ",", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))
      cat("\t", sep = ",", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))
  
      # Probands_w_variant_genotypes
      cat(proband_geno[proband_geno_het], sep = ",", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))

      proband_geno_hom = (proband_geno == "0/0")

      # N_probands_wo_variant
      cat("\t", length(proband_geno[proband_geno_hom]), "\t", sep = "", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))

      # Probands_wo_variant
      cat(proband_names[proband_geno_hom], sep = ",", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))
      cat("\t", sep = ",", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))

      # Probands_wo_variant_genotypes
      cat(proband_geno[proband_geno_hom], sep = ",", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))

      m = match(non_probands, colnames(IN))
      non_proband_geno = unlist(IN[var,m])
      non_proband_names = colnames(IN)[m]

      non_proband_geno_het = (non_proband_geno == "0/1") | (non_proband_geno == "1/1")

      # N_non_probands_w_variant
      cat("\t", length(non_proband_geno[non_proband_geno_het]), "\t", sep = "", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))
      
      # non_probands_w_variant
      cat(non_proband_names[non_proband_geno_het], sep = ",", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))
      cat("\t", sep = ",", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))

      # non_probands_w_variant_genotypes
      cat(non_proband_geno[non_proband_geno_het], sep = ",", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))

      non_proband_geno_hom = (non_proband_geno == "0/0")

      # N_non_probands_wo_variant
      cat("\t", length(non_proband_geno[non_proband_geno_hom]), "\t", sep = "", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))

      # non_probands_wo_variant
      cat(non_proband_names[non_proband_geno_hom], sep = ",", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))
      cat("\t", sep = ",", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))
  
      # non_probands_wo_variant_genotypes
      cat(non_proband_geno[non_proband_geno_hom], sep = ",", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))

      cat("\n", append = TRUE, file = paste0("intermediates/", args[1], "_valid_probands.txt"))

     }
  }
}
    


