# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Calculates probability of co-occurance

args = commandArgs(trailingOnly=TRUE)
gene = args[1]

# Calculated from the SPARK parents. Gene and frequency of ultra-rare deleterious variant.
STAT <- read.delim(file = "parental_gene_mut_stats.txt", header = FALSE)

# Probability of success is the frequency of alleles in population
prob_succ = STAT$V2
# Ranked success based on the gene with the most deleterious variants (TTN).
gene_prob = STAT$V2 / max(STAT$V2)

IN <- read.delim(file = paste0("intermediates/", gene, ".elsewhere.collapse.txt"), header = FALSE)

# Calculates the number of families with B* inheritence given B* inheritence
tb = table(IN$V5)
n_A = as.integer(length(unique(IN$V1)))

system(paste0("rm intermediates/", gene, ".table"))
for (gene2 in rownames(tb)){
 if(tb[[gene2]] > 1){
    n_B = as.integer(tb[[gene2]])
    gene_IN = gene2
    sums = 0
    for(i in 1:length(STAT$V1)){
       sums = sums + dbinom(n_B, n_A, prob_succ[i])
    }
    cat(gene, "\t", n_A, "\t", gene2, "\t", n_B, "\t", sums, "\t", sums* gene_prob[STAT$V1 == gene2], 
	"\n", sep = "", file = paste0("intermediates/", gene, ".table"), append = TRUE)
  }
}


