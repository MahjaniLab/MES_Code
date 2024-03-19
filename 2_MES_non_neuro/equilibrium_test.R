# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Performs equilibrium test

library(stringr)
`%ni%` = Negate(`%in%`)
args = commandArgs(trailingOnly=TRUE)

# The input file is a merge of data from all cohorts. Header:
#Gene1,Gene2,N_AB_ASD,N_AB_non_ASD,N_A*_ASD,N_A*_non_ASD,N_B*_ASD,N_B*_non_ASD,N_A*B*_ASD,N_A*B*_non_ASD

IN = read.csv(file = "../compare_deleterious_var_rates/merged_conting.csv", header = FALSE)
for(pair in 1:length(IN$V1)){
  total_samp = IN$V3[pair] + IN$V4[pair] + IN$V5[pair] + IN$V6[pair] + IN$V7[pair] + IN$V8[pair] + IN$V9[pair] + IN$V10[pair]
  ASD_total_samp = IN$V3[pair] + IN$V5[pair] + IN$V7[pair] + IN$V9[pair]
  No_ASD_total_samp = IN$V4[pair] + IN$V6[pair] + IN$V8[pair] + IN$V10[pair]

  freq_A = (IN$V5[pair] + IN$V6[pair] + IN$V9[pair] + IN$V10[pair]) / total_samp
  freq_B = (IN$V7[pair] + IN$V8[pair] + IN$V9[pair] + IN$V10[pair]) / total_samp

  p_geneA_ASD = pbinom((IN$V5[pair] + IN$V9[pair]), ASD_total_samp, freq_A)
  p_geneB_ASD = pbinom((IN$V7[pair] + IN$V9[pair]), ASD_total_samp, freq_B)

  p_geneA_non_ASD = pbinom((IN$V6[pair] + IN$V10[pair]), No_ASD_total_samp, freq_A)
  p_geneB_non_ASD = pbinom((IN$V8[pair] + IN$V10[pair]), No_ASD_total_samp, freq_B)

  # Probability of seeing at least X A*B* in ASD
  P_AB_ASD = pbinom(IN$V9[pair] -1, ASD_total_samp, freq_A*freq_B, lower.tail = FALSE)
  
  # Probability of seeing A*B* in ASD
  P_AB_no_ASD = pbinom(IN$V10[pair], No_ASD_total_samp, freq_A*freq_B, lower.tail = TRUE)

  P_all = pbinom((IN$V9[pair] + IN$V10[pair] -1), total_samp, freq_A*freq_B, lower.tail = FALSE)

  cat(IN$V1[pair], IN$V2[pair], P_AB_ASD*P_AB_no_ASD, P_all, P_all/(P_AB_ASD*P_AB_no_ASD), "\n", sep = ",")

}
