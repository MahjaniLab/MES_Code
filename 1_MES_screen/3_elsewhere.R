# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Finds B* mirrored inheritence for A* cases

args = commandArgs(trailingOnly=TRUE)

# Gene A that we are comparing to other genes
geneA = args[1]

IN <- read.delim(file = "intermediates/1_valid_probands.passed.txt", header = TRUE, sep = "\t")
system(paste0("rm intermediates/", geneA, ".elsewhere.txt"))

# All families with this same gene
FAM_w_GENEA = subset(IN, Gene == geneA)
families = unique(FAM_w_GENEA$Family)

# For each family with A* inheritence in Gene A
for (fam in families){
  SP_FAMILY_all_GENES = subset(IN, (Family == fam))
  SP_FAMILY_w_GENEA = subset(IN, (Family == fam) & (Gene == geneA))

  # Parsing all variants in this family that exist in other genes
  for(var in 1:length(SP_FAMILY_all_GENES$Variant)){
    
    # Second variant is present in all other probands
    if(SP_FAMILY_all_GENES$N_probands_w_variant[var] == SP_FAMILY_all_GENES$N_probands[var]){

       # Now looking at the original A* variant
       for(var2 in 1:length(SP_FAMILY_w_GENEA$Variant)){     

         # No unaffected siblings with A* have B*
         var_B = unlist(strsplit(SP_FAMILY_all_GENES$Siblings_w_variant[var], ","))
         var_A = unlist(strsplit(SP_FAMILY_w_GENEA$Siblings_w_variant[var2], ","))

         if(!any(var_B %in% var_A)){

           # And other variant is not in parent where variant A came from
           if((SP_FAMILY_w_GENEA$Maternal_genotype[var2] != SP_FAMILY_all_GENES$Maternal_genotype[var]) && (SP_FAMILY_w_GENEA$Paternal_genotype[var2] != SP_FAMILY_all_GENES$Paternal_genotype[var])){

            cat(fam, "\t", SP_FAMILY_w_GENEA$Gene[var2], "\t", 
		SP_FAMILY_w_GENEA$Variant[var2], "\t", SP_FAMILY_all_GENES$Variant[var], "\t", 
		SP_FAMILY_all_GENES$Gene[var], "\n", 
		sep = "", file = paste0("intermediates/", geneA, ".elsewhere.txt"), append = TRUE)
           }
         }
        }
     }
  }
}
