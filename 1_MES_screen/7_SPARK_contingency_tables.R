# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Creates contingency tables for ASD/non-ASD individuals in SPARK based on A* and B* status

library(stringr)
library(rapport)
`%ni%` <- Negate(`%in%`)
args = commandArgs(trailingOnly=TRUE)

# Has ASD status info:
FAM <- read.delim(file = "SPARK.iWES_v1.mastertable.2022_02.tsv", header = TRUE)

unaffected = c(FAM$spid[FAM$asd == 1])
affected = c(FAM$spid[FAM$asd == 2])

IN <- read.delim(file = args[1], header = FALSE)
genes_split = unlist(str_split(args[1], "[.]"))
genes_split[1] = str_remove(genes_split[1], "intermediates/")
genes = genes_split[1:2]

AB_no_asd = unaffected[unaffected %ni% IN$V2]
AB_asd = affected[affected %ni% IN$V2]

matrix_cont = matrix(nrow = 4, ncol = 2)
matrix_cont[1,1:2] = c(length(AB_asd), length(AB_no_asd))

A_no_asd = vector()
A_asd = vector()

B_no_asd = vector()
B_asd = vector()

both_no_asd = vector()
both_asd = vector()

for(samp in affected){
   A <- unique(IN$V1[IN$V2 == samp])
   if(length(A) == 1){
     if(A[1] == genes[1]){ A_asd = c(A_asd, samp)}
     if(A[1] == genes[2]){ B_asd = c(B_asd, samp)}
   }
   if(length(A) == 2){
     both_asd = c(both_asd, samp)  
  }
}

for(samp in unaffected){
   A <- unique(IN$V1[IN$V2 == samp])
   if(length(A) == 1){
     if(A[1] == genes[1]){ A_no_asd = c(A_no_asd, samp)}
     if(A[1] == genes[2]){ B_no_asd = c(B_no_asd, samp)}
   }
   if(length(A) == 2){
     both_no_asd = c(both_no_asd, samp)
  }
}

matrix_cont[2,1:2] = c(length(unique(A_asd)), length(unique(A_no_asd)))
matrix_cont[3,1:2] = c(length(unique(B_asd)), length(unique(B_no_asd)))
matrix_cont[4,1:2] = c(length(unique(both_asd)), length(unique(both_no_asd)))

out_matrix_file = paste0("intermediates/", genes[1], "_", genes[2], ".table")
write.table(matrix_cont, file = out_matrix_file, sep = ",", row.names = FALSE, col.names = FALSE)
