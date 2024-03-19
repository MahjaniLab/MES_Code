# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Collapses cases where one family has more than one A* or B* in the same gene

args = commandArgs(trailingOnly=TRUE)
geneA = args[1]

ELSE <- read.delim(file = paste0("intermediates/", geneA, ".elsewhere.txt"), header = FALSE)
skip = vector()

system(paste0("rm intermediates/", geneA, ".elsewhere.collapse.1.txt"))
system(paste0("rm intermediates/", geneA, ".elsewhere.collapse.2.txt"))


for(line in 1:length(ELSE$V1)){
  if(!(line %in% skip)){
    sub = which((ELSE$V1 == ELSE$V1[line]) & (ELSE$V5 == ELSE$V5[line]) & (ELSE$V2 == ELSE$V2[line]) & (ELSE$V3 == ELSE$V3[line]))

    if(length(sub) == 1){
      cat(ELSE$V1[line], "\t", ELSE$V2[line], "\t", ELSE$V3[line], "\t", ELSE$V4[line], "\t", ELSE$V5[line], "\n", sep = "", file = paste0("intermediates/", geneA, ".elsewhere.collapse.1.txt"), append = TRUE)
    }
    if(length(sub) > 1){
      cat(ELSE$V1[line], "\t", ELSE$V2[line], "\t", ELSE$V3[line], "\t", paste(ELSE$V4[sub], collapse = ","), "\t", ELSE$V5[line], "\n", sep = "", file = paste0("intermediates/", geneA, ".elsewhere.collapse.1.txt"), append = TRUE)
      skip = c(skip, sub)
    }
  }
}

ELSE <- read.delim(file = paste0("intermediates/", geneA, ".elsewhere.collapse.1.txt"), header = FALSE)
skip = vector()

for(line in 1:length(ELSE$V1)){
  if(!(line %in% skip)){
    sub = which((ELSE$V1 == ELSE$V1[line]) & (ELSE$V5 == ELSE$V5[line]) & (ELSE$V2 == ELSE$V2[line]) & (ELSE$V4 == ELSE$V4[line]))

    if(length(sub) == 1){
      cat(ELSE$V1[line], "\t", ELSE$V2[line], "\t", ELSE$V3[line], "\t", ELSE$V4[line], "\t", ELSE$V5[line], "\n", sep = "", file = paste0("intermediates/", geneA, ".elsewhere.collapse.2.txt"), append = TRUE)
    }
    if(length(sub) > 1){
      cat(ELSE$V1[line], "\t", ELSE$V2[line], "\t", paste(ELSE$V3[sub], collapse = ","), "\t", ELSE$V4[line], "\t", ELSE$V5[line], "\n", sep = "", file = paste0("intermediates/", geneA, ".elsewhere.collapse.2.txt"), append = TRUE)
      skip = c(skip, sub)
    }
  }
}

system(paste0("rm intermediates/", geneA, ".elsewhere.collapse.1.txt"))
system(paste0("mv intermediates/", geneA, ".elsewhere.collapse.2.txt intermediates/", geneA, ".elsewhere.collapse.txt"))

