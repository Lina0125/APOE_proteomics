source('APOE_proteomics/analysis/setting.R')
#-----------------------------------Set analysis--------------------------------
dataset_identifier = 'ADNI_MS'
apoe_identifier = 'APOE4'
analysis = 'Aβ'

local(
  {
    loadData = main_driver(dataset_identifier, apoe_identifier, analysis)
    analysis_df <<- loadData$df
    apoe_col <<- loadData$apoe_col
    pathology_col <<- loadData$pathology_col
    protein_cols <<- loadData$protein_cols
    result_folder <<- loadData$result_folder
    cov_cols <<- loadData$cov_cols
    age_col <<- loadData$age_col
  } 
)

get_demographic(analysis_df, dataset_identifier, apoe_identifier, FALSE)

#-------------------------------Analysis process--------------------------------
apoe2protein = apoeOrAb2protein(analysis_df, protein_cols, apoe_col, cov_cols, dataset_identifier)
apoe2protein_adjab = apoeOrAb2protein(analysis_df, protein_cols, apoe_col, 
                                      c(cov_cols, pathology_col), dataset_identifier)
ab2protein = apoeOrAb2protein(analysis_df, protein_cols, pathology_col, cov_cols, dataset_identifier)
ab2protein_adjapoe = apoeOrAb2protein(analysis_df, protein_cols, pathology_col, 
                                      c(cov_cols, apoe_col), dataset_identifier)

write.csv(apoe2protein, file = paste0(result_folder,"apoe2protein.csv"), 
          row.names = FALSE)
write.csv(apoe2protein_adjab, 
          file = paste0(result_folder,"apoe2protein_adjab.csv"), 
          row.names = FALSE)
write.csv(ab2protein, file = paste0(result_folder,"ab2protein.csv"), 
          row.names = FALSE)
write.csv(ab2protein_adjapoe, 
          file = paste0(result_folder,"ab2protein_adjapoe.csv"), 
          row.names = FALSE)

#-------------------------------Mediation analysis------------------------------
sig_df = subset(apoe2protein, p_adjusted < 0.05)
nrow(sig_df)
protein_mediation = protein_mediate_analysis(data_df = analysis_df,
                                             cov_cols = cov_cols, 
                                             pathology_col = pathology_col, 
                                             apoe_col = apoe_col,
                                             proteins = sig_df$Protein_id,
                                             dataset_identifier = dataset_identifier)

ab_mediation = pathology_mediate_analysis(data_df = analysis_df, 
                                          cov_cols = cov_cols, 
                                          pathology_col = pathology_col, 
                                          apoe_col = apoe_col,
                                          proteins = sig_df$Protein_id,
                                          dataset_identifier = dataset_identifier)

# Save mediation
write.csv(protein_mediation, 
          file = paste0(result_folder,"protein_mediation.csv"), 
          row.names = FALSE)

write.csv(ab_mediation, 
          file = paste0(result_folder,"ab_mediation.csv"), 
          row.names = FALSE)

#----------------------APOE associated proteins in Aβ-negative------------------
##----------------------In the whole Aβ-negative--------------------------------
analysis_df_neg = subset(analysis_df, analysis_df[[pathology_col]] == 0)
get_demographic(analysis_df_neg, dataset_identifier, apoe_identifier)
neg_allage = apoeOrAb2protein(analysis_df_neg, protein_cols, apoe_col, cov_cols, dataset_identifier)
##-------------------------In the younger Aβ-negative---------------------------
neg_younger = subset(analysis_df_neg, analysis_df_neg[[age_col]] <= median(analysis_df_neg[[age_col]]))
get_demographic(neg_younger, dataset_identifier, apoe_identifier)
neg_younger_result = apoeOrAb2protein(neg_younger, protein_cols, apoe_col, cov_cols, dataset_identifier)
# In the older Aβ-
neg_older = subset(analysis_df_neg, analysis_df_neg[[age_col]] > median(analysis_df_neg[[age_col]]))
get_demographic(neg_older, dataset_identifier, apoe_identifier)
neg_older_result = apoeOrAb2protein(neg_older, protein_cols, apoe_col, cov_cols, dataset_identifier)

write.csv(neg_allage, file = paste0(result_folder,"neg_allage.csv"), row.names = FALSE)
write.csv(neg_younger_result, file = paste0(result_folder,"neg_younger.csv"), row.names = FALSE)
write.csv(neg_older_result, file = paste0(result_folder,"neg_older.csv"), row.names = FALSE)
##-----------------------APOE*age interaction in Aβ-negative--------------------
sig_df = subset(neg_allage, p_adjusted < 0.05)
analysis_df_neg$apoeNage = analysis_df_neg[[apoe_col]] * analysis_df_neg[[age_col]]
apoeNage_result = apoeOrAb2protein(analysis_df_neg, sig_df$Protein_id, 
                                   'apoeNage', c(cov_cols, apoe_col), dataset_identifier)

write.csv(apoeNage_result, file = paste0(result_folder,"neg_allage_interaction.csv"), row.names = FALSE)

#------------------------------AD in e3e3---------------------------------------
# pos_ine3e3 = apoeOrAb2protein(subset(initi, APOE==33), protein_cols, pathology_col, cov_cols)
# pos_ine3e3 = map_symbol(pos_ine3e3, dataset_identifier)
# write.csv(apoeNage_result, file = paste0(result_folder,"pos_ine3e3.csv"), row.names = FALSE)

#---------------------APOE associated proteins in Aβ-postive--------------------
analysis_df_pos = subset(analysis_df, analysis_df[[pathology_col]] == 1)
get_demographic(analysis_df_pos, dataset_identifier, apoe_identifier)
pos_allage = apoeOrAb2protein(analysis_df_pos, protein_cols, apoe_col, cov_cols, dataset_identifier)

write.csv(pos_allage, file = paste0(result_folder,"pos_allage.csv"), row.names = FALSE)

#----------------------------------APOE4 * AD-----------------------------------
apoe2protein = read.csv(paste0(result_folder,"apoe2protein.csv"))

sig_df = subset(apoe2protein, p_adjusted < 0.05)
analysis_df$apoe_ab = analysis_df[[apoe_col]] * analysis_df[[pathology_col]]
apoe_ab_interaction = apoeOrAb2protein(analysis_df, sig_df$Protein_id, 
                                       'apoe_ab', c(cov_cols, apoe_col, pathology_col),
                                       dataset_identifier)

write.csv(apoe_ab_interaction, file = paste0(result_folder,"apoe_ab_interaction.csv"), row.names = FALSE)

#--------------------------------APOE * protein---------------------------------
get_demographic(subset(analysis_df,
                       diagnosis_baseline_variable != 'AD'), dataset_identifier, apoe_identifier)

apoeNprotein_onAb = apoeNproteinOnAD(subset(analysis_df,
                                            diagnosis_baseline_variable != 'AD'),
                                     protein_cols, apoe_col, pathology_col, cov_cols, dataset_identifier)

apoeNprotein_onAb$standardized_beta = apoeNprotein_onAb$coef #Only for placeholder in result_summary.ipynb

write.csv(apoeNprotein_onAb, file = paste0(result_folder,"protein_regulation.csv"), row.names = FALSE)
