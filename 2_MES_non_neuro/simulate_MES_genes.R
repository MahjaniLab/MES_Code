# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Simulates inheritence of genes for null MES in all individuals

args = commandArgs(trailingOnly=TRUE)

# This is the simulation ID:
iteration = args[1]

# Mutation rate probabilities by gene
PROB = read.delim(file = "../infiles/gene_mut_stats_ultrarare.txt", header = FALSE)
# Infiles:
SPARK = read.delim(file = "../infiles/SPARK.ped", header = FALSE)
KGP = read.delim(file = "../infiles/1kgp.ped", header = FALSE, sep = " ")
KGP = subset(KGP, select = -V7)
KGP$V6 = "0"
PED_MERGE = rbind(SPARK,KGP)

# Step 1: Simulate any orphans:
ORPHANS = c(unique(SPARK$V2[(SPARK$V3 == "0") & (SPARK$V4 == "0")]), unique(KGP$V2[(KGP$V3 == "0") & (KGP$V4 == "0")]))
system(paste0("rm sampling_results/variant_assignment/", iteration, ".variants"))
for(indiv in ORPHANS){
  variants_pull = PROB$V1[rbinom(length(PROB$V1), 1, PROB$V2) == 1]
  cat(indiv, "\t", paste0(variants_pull, sep = ","), "\n", sep = "", file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), append = TRUE)
}

# Step 2: Simulate offspring from two orphans
VARIANTS = read.delim(file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), header = FALSE)
PARENTS_ORPHANS = c(unique(SPARK$V2[(SPARK$V3 %in% ORPHANS) & (SPARK$V4 %in% ORPHANS)]), unique(KGP$V2[(KGP$V3 %in% ORPHANS) & (KGP$V4 %in% ORPHANS)]))

for(kid in PARENTS_ORPHANS){
  father = unique(PED_MERGE$V3[PED_MERGE$V2 == kid])
  mother = unique(PED_MERGE$V4[PED_MERGE$V2 == kid])
  father_variants = unlist(strsplit(VARIANTS$V2[VARIANTS$V1 == father], split = ","))
  mother_variants = unlist(strsplit(VARIANTS$V2[VARIANTS$V1 == mother], split = ","))
  child_variants = unique(c(father_variants[rbinom(length(father_variants), 1, 0.5) == 1], mother_variants[rbinom(length(mother_variants), 1, 0.5) == 1]))
  cat(kid, "\t", paste0(child_variants, sep = ","), "\n", sep = "", file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), append = TRUE)
}

# Step 3: Simulate offspring from one orphan
VARIANTS = read.delim(file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), header = FALSE)
MOTHER_ORPHANS = c(unique(SPARK$V2[(SPARK$V3 %in% ORPHANS) & !(SPARK$V4 %in% ORPHANS)]), unique(KGP$V2[(KGP$V3 %in% ORPHANS) & !(KGP$V4 %in% ORPHANS)]))

for(kid in MOTHER_ORPHANS){
  father = unique(PED_MERGE$V3[PED_MERGE$V2 == kid])
  father_variants = unlist(strsplit(VARIANTS$V2[VARIANTS$V1 == father], split = ","))
  mother_variants = PROB$V1[rbinom(length(PROB$V1), 1, PROB$V2) == 1]
  child_variants = unique(c(father_variants[rbinom(length(father_variants), 1, 0.5) == 1], mother_variants[rbinom(length(mother_variants), 1, 0.5) == 1]))
  cat(kid, "\t", paste0(child_variants, sep = ","), "\n", sep = "", file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), append = TRUE)
}

FATHER_ORPHANS = c(unique(SPARK$V2[!(SPARK$V3 %in% ORPHANS) & (SPARK$V4 %in% ORPHANS)]), unique(KGP$V2[!(KGP$V3 %in% ORPHANS) & (KGP$V4 %in% ORPHANS)]))

for(kid in FATHER_ORPHANS){
  mother = unique(PED_MERGE$V4[PED_MERGE$V2 == kid])
  mother_variants = unlist(strsplit(VARIANTS$V2[VARIANTS$V1 == mother], split = ","))
  father_variants = PROB$V1[rbinom(length(PROB$V1), 1, PROB$V2) == 1]
  child_variants = unique(c(father_variants[rbinom(length(father_variants), 1, 0.5) == 1], mother_variants[rbinom(length(mother_variants), 1, 0.5) == 1]))
  cat(kid, "\t", paste0(child_variants, sep = ","), "\n", sep = "", file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), append = TRUE)
}

VARIANTS = read.delim(file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), header = FALSE)
ALL_SAMPLES = unique(c(PED_MERGE$V2))
NOT_ORPHANS = ALL_SAMPLES[!(ALL_SAMPLES %in% VARIANTS$V1)]

# Step 4: Fill in fringe cases of missing but ID'd parents
for(kid in  NOT_ORPHANS){
  mother = unique(PED_MERGE$V4[PED_MERGE$V2 == kid])
  if(!(mother %in% VARIANTS$V1)){
    mother_variants = PROB$V1[rbinom(length(PROB$V1), 1, PROB$V2) == 1]
    cat(mother, "\t", paste0(mother_variants, sep = ","), "\n", sep = "", file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), append = TRUE)
    df = data.frame(V1 = mother, V2 = paste0(mother_variants, sep = ","))
    VARIANTS = rbind(VARIANTS, df)
  }
  father = unique(PED_MERGE$V3[PED_MERGE$V2 == kid])
  if(!(father %in% VARIANTS$V1)){
    father_variants = PROB$V1[rbinom(length(PROB$V1), 1, PROB$V2) == 1]
    cat(father, "\t", paste0(father_variants, sep = ","), "\n", sep = "", file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), append = TRUE)
    df = data.frame(V1 = father, V2 = paste0(father_variants, sep = ","))
    VARIANTS = rbind(VARIANTS, df)
  }
}


# Step 5: Adress the final caes
VARIANTS = read.delim(file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), header = FALSE)
ALL_SAMPLES = unique(c(PED_MERGE$V2))
NOT_ORPHANS = ALL_SAMPLES[!(ALL_SAMPLES %in% VARIANTS$V1)]

for (kid in NOT_ORPHANS){
  mother = unique(PED_MERGE$V4[PED_MERGE$V2 == kid])
  father = unique(PED_MERGE$V3[PED_MERGE$V2 == kid])
  if((father != "0") && (mother != "0") && (father %in% VARIANTS$V1) && (mother %in% VARIANTS$V1)){
    father_variants = unlist(strsplit(VARIANTS$V2[VARIANTS$V1 == father], split = ","))
    mother_variants = unlist(strsplit(VARIANTS$V2[VARIANTS$V1 == mother], split = ","))
    child_variants = unique(c(father_variants[rbinom(length(father_variants), 1, 0.5) == 1], mother_variants[rbinom(length(mother_variants), 1, 0.5) == 1]))
    cat(kid, "\t", paste0(child_variants, sep = ","), "\n", sep = "", file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), append = TRUE)
  }
  #Just the father, simulate the mother
  if((father != "0") && (mother == "0") && (father %in% VARIANTS$V1)){
    father_variants = unlist(strsplit(VARIANTS$V2[VARIANTS$V1 == father], split = ","))
    mother_variants = PROB$V1[rbinom(length(PROB$V1), 1, PROB$V2) == 1]
    child_variants = unique(c(father_variants[rbinom(length(father_variants), 1, 0.5) == 1], mother_variants[rbinom(length(mother_variants), 1, 0.5) == 1]))
    cat(kid, "\t", paste0(child_variants, sep = ","), "\n", sep = "", file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), append = TRUE)
  }
  # Just the mother, simulate the father
  if((father == "0") && (mother != "0") && (mother %in% VARIANTS$V1)){
    mother_variants = unlist(strsplit(VARIANTS$V2[VARIANTS$V1 == mother], split = ","))
    father_variants = PROB$V1[rbinom(length(PROB$V1), 1, PROB$V2) == 1]
    child_variants = unique(c(father_variants[rbinom(length(father_variants), 1, 0.5) == 1], mother_variants[rbinom(length(mother_variants), 1, 0.5) == 1]))
    cat(kid, "\t", paste0(child_variants, sep = ","), "\n", sep = "", file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), append = TRUE)
  }
}

# Step 6: Do remaining samples
VARIANTS = read.delim(file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), header = FALSE)
DONE = VARIANTS$V1
SPARK_to_do = SPARK$V2[!(SPARK$V2 %in% DONE)]
KGP_to_do = KGP$V2[!(KGP$V2 %in% DONE)]
BioME_to_do = paste0("BIOME", c(1:16586))
AoU_to_do = paste0("AoU", c(1:156550))

TO_DO = c(SPARK_to_do,KGP_to_do,BioME_to_do,AoU_to_do)

for(sample in TO_DO){
  variants = PROB$V1[rbinom(length(PROB$V1), 1, PROB$V2) == 1]
  cat(sample, "\t", paste0(variants, sep = ","), "\n", sep = "", file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), append = TRUE)
}


## An additional script then corrects for MZ twins. Code:
# CORRECTION FOR MZ TWINS
#args = commandArgs(trailingOnly=TRUE)
#iteration = args[1]
#VARIANTS = read.delim(file = paste0("sampling_results/variant_assignment/", iteration, ".variants"), header = FALSE)
#MZ = read.delim(file = "../../double_offspring/infiles/vcfs/fam/MZ_pairs.txt", header = FALSE)
#for(twin in 1:length(MZ$V1)){
#  VARIANTS$V2[VARIANTS$V1 == MZ$V2[twin]] = VARIANTS$V2[VARIANTS$V1 == MZ$V1[twin]]
#}

#write.table(VARIANTS, paste0("sampling_results/variant_assignment/", iteration, ".mzfixed.variants"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# The remaining steps analyze the output as performed for real results.
