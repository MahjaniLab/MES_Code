# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Performs the equilibrum test for SINGLE GENES, not pairs

# The input file is a merge of data from all cohorts. Header:
#Gene1,Gene2,N_AB_ASD,N_AB_non_ASD,N_A*_ASD,N_A*_non_ASD,N_B*_ASD,N_B*_non_ASD,N_A*B*_ASD,N_A*B*_non_ASD

library(stringr)
`%ni%` = Negate(`%in%`)

IN = read.csv(file = "merged_conting.csv", header = FALSE)
for(pair in 1:length(IN$V1)){
  total_samp = IN$V3[pair] + IN$V4[pair] + IN$V5[pair] + IN$V6[pair] + IN$V7[pair] + IN$V8[pair] + IN$V9[pair] + IN$V10[pair]
  ASD_total_samp = IN$V3[pair] + IN$V5[pair] + IN$V7[pair] + IN$V9[pair]
  No_ASD_total_samp = IN$V4[pair] + IN$V6[pair] + IN$V8[pair] + IN$V10[pair]

  freq_A = (IN$V5[pair] + IN$V6[pair] + IN$V9[pair] + IN$V10[pair]) / total_samp
  freq_B = (IN$V7[pair] + IN$V8[pair] + IN$V9[pair] + IN$V10[pair]) / total_samp

  # Probability of seeing at least X A* in ASD
  P_A_ASD = qbinom(0.95, ASD_total_samp, freq_A)
  P_B_ASD = qbinom(0.95, ASD_total_samp, freq_B)

  # Neither gene is overrepresented in ASD
  pass = (IN$V5[pair] < P_A_ASD) && (IN$V7[pair] < P_B_ASD) 

  # Outputs some statistics and the pass score for the pair
  cat(IN$V1[pair], qbinom(0.95, ASD_total_samp, freq_A), qbinom(0.05, ASD_total_samp, freq_A), IN$V5[pair], 
    IN$V2[pair], qbinom(0.95, ASD_total_samp, freq_B), qbinom(0.05, ASD_total_samp, freq_B), IN$V7[pair], pass, "\n", sep = ",")

}
