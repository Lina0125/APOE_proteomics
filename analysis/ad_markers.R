library(readr)
library(QuantPsyc)
library(dplyr)
library(writexl)

#----------------------Proteomics data load-------------------------------------
proteomics = read.csv("data/BF2SomaLogic/preprocessed/V5_with24.csv")

#----------------------------AD markers load------------------------------------
bf2 = read.csv("data/BF2SomaLogic/raw_data/alexa__20240623_175941.csv")
bf2 = subset(bf2,Study=="BF2" & Visit==0 & data_index==0)
tau_pet = bf2[,c("sid","tnic_cho_com_I_IV","tnic_cho_com_I_II",
                  "tnic_cho_com_V_VI","fnc_ber_com_composite", "ct_adsign_lr",
                  "mmse_score","mPACC_v1","mPACC_v2","mPACC_v3")]

#----------------------------Data processing------------------------------------
merged = merge(proteomics,tau_pet, by="sid",all.x=TRUE)
all_aptamers = grep('ANML$', names(merged), value = TRUE)
# Load APOE-associated proteins
apoe2protein = read.csv("data/BF2SomaLogic/Rresults/V5_e4vse3e3/Whole/apoe2protein.csv")
apoe2protein_e2 = read.csv("data/BF2SomaLogic/Rresults/V5_e2vse3e3/Whole/apoe2protein.csv")
apoe2protein = subset(apoe2protein, p_adjusted < 0.05)
apoe2protein_e2 = subset(apoe2protein_e2, p_adjusted < 0.05)
use_aptamers = unique(c(apoe2protein$apt_name, apoe2protein_e2$apt_name))

aptamers = c()
for (x in all_aptamers){
  parts = strsplit(x, "__")[[1]]
  result = parts[length(parts) - 1]
  if (result %in% use_aptamers){
    aptamers = c(aptamers, x)
  }
}

#--------------------------------Analysis---------------------------------------
Abneg_group = get_function_related(subset(merged, Abnormal_CSF_Ab42_Ab40_Ratio==0),
                                   aptamers, 
                                   c('age', 'gender_baseline_variable', 
                                     'plasma_ANML_mean'))
Abneg_group = map_symbol(Abneg_group)

Abpos_group = get_function_related(subset(merged, Abnormal_CSF_Ab42_Ab40_Ratio==1),
                                   aptamers, 
                                   c('age', 'gender_baseline_variable', 
                                     'plasma_ANML_mean'))
Abpos_group = map_symbol(Abpos_group)

#--------------------------Save and check results-------------------------------
write_xlsx(Abneg_group,"data/BF2SomaLogic/ad_markers_associations_Abneg.xlsx")
write_xlsx(Abpos_group,"data/BF2SomaLogic/ad_markers_associations_Abpos.xlsx")

resultsAbpos_fdr_any = Abpos_group %>% filter(if_any(contains("_fdr"), ~ .x < 0.05))
resultsAbneg_fdr_any = Abneg_group %>% filter(if_any(contains("_fdr"), ~ .x < 0.05))
