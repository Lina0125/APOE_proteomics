#========================Main sensitivity analysis==============================
source('APOE_proteomics/analysis/setting.R')
dataset_identifier = 'BF2SomaLogic'
walk(c('PC10', 'drug', 'wml', 'AbRatio', 'AbPET'), function(analysis){
  walk(c('APOE4', 'APOE2'), function(apoe_identifier){
    #--------------Set analysis-----------
    loadData = main_driver(dataset_identifier, apoe_identifier, analysis)
    analysis_df = loadData$df
    apoe_col = loadData$apoe_col
    pathology_col = loadData$pathology_col
    protein_cols = loadData$protein_cols
    result_folder = loadData$result_folder
    cov_cols = loadData$cov_cols
    age_col = loadData$age_col
    
    get_demographic(analysis_df, dataset_identifier, apoe_identifier)
    
    #------------Model 1, 2, 3-----------------------
    apoe2protein = apoeOrAb2protein(analysis_df, protein_cols, apoe_col, cov_cols, 
                                    dataset_identifier)
    apoe2protein_adjab = apoeOrAb2protein(analysis_df, protein_cols, apoe_col, 
                                          c(cov_cols, pathology_col), 
                                          dataset_identifier)
    ab2protein = apoeOrAb2protein(analysis_df, protein_cols, pathology_col, cov_cols, 
                                  dataset_identifier)
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
    
    #----------Mediation analysis---------------------
    apoe2protein = read.csv(paste0(result_folder,"apoe2protein.csv")) 
    
    sig_df = subset(apoe2protein, p_adjusted < 0.05)
    print(nrow(sig_df))
    
    if (analysis == 'AbRatio' || analysis == 'AbPET') {
      protein_mediation = protein_mediate_Abcontinuous(data_df = analysis_df,
                                                       cov_cols = cov_cols,
                                                       pathology_col = pathology_col,
                                                       apoe_col = apoe_col,
                                                       proteins = sig_df$Protein_id,
                                                       dataset_identifier = dataset_identifier)
      
      ab_mediation = Abcontinuous_mediate_analysis(data_df = analysis_df,
                                                   cov_cols = cov_cols,
                                                   pathology_col = pathology_col,
                                                   apoe_col = apoe_col,
                                                   proteins = sig_df$Protein_id,
                                                   dataset_identifier = dataset_identifier)
    } else {
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
    }
    
    #Save mediation
    write.csv(protein_mediation,
              file = paste0(result_folder,"protein_mediation.csv"),
              row.names = FALSE)
    
    write.csv(ab_mediation,
              file = paste0(result_folder,"ab_mediation.csv"),
              row.names = FALSE)
  })
})

#========================APOE * sex effect on early dysregulated proteins=======
source('APOE_proteomics/analysis/setting.R')
dataset_identifier = 'BF2SomaLogic'
apoe_identifier = 'APOE4'
analysis = 'Aβ'

loadData = main_driver(dataset_identifier, apoe_identifier, analysis)
analysis_df = loadData$df
apoe_col = loadData$apoe_col
pathology_col = loadData$pathology_col
protein_cols = loadData$protein_cols
result_folder = loadData$result_folder
cov_cols = loadData$cov_cols
age_col = loadData$age_col

get_demographic(analysis_df, dataset_identifier, apoe_identifier)

sex_list = list(
  male = 0,
  female = 1
)

##-----Load results for proteins associated with APOE at Abeta-negative stages--
neg_df = analysis_df[analysis_df[[pathology_col]] == 0, ]
neg = read_csv(glue('{result_folder}/neg_allage.csv'))
neg_sig = neg[neg$p_adjusted < 0.05, ]
protein_cols = neg_sig$Protein_id
##-------------APOE * sex effect on those proteins--------------
neg_df$apoe_sex = neg_df[[apoe_col]] * neg_df$gender_baseline_variable
apoe_sex_interaction = apoeOrAb2protein(neg_df, protein_cols, 
                                        'apoe_sex', c(cov_cols, apoe_col),
                                        dataset_identifier)

write.csv(apoe_sex_interaction, 
          file = glue("{result_folder}/apoe_sex_interaction.csv"), 
          row.names = FALSE)

##--------Sex stratification analysis-----------
table(neg_df$gender_baseline_variable)

sex_results = list()
walk(names(sex_list), function(sex){
  tmp_df = neg_df[neg_df$gender_baseline_variable == sex_list[[sex]],]
  print(get_demographic(tmp_df, dataset_identifier, apoe_identifier))
  
  re = apoeOrAb2protein(tmp_df, protein_cols, apoe_col, cov_cols,
                        dataset_identifier)
  re$Gender = sex
  
  sex_results[[sex]] <<- re
  write.csv(re, file = glue("{result_folder}/apoe2protein_neg_{sex}.csv"), 
            row.names = FALSE)
})

#==============================Medication and SPC25=============================
source('APOE_proteomics/analysis/setting.R')
dataset_identifier = 'BF2SomaLogic'
apoe_identifier = 'APOE4'
analysis = 'Aβ'

loadData = main_driver(dataset_identifier, apoe_identifier, analysis)
analysis_df = loadData$df
apoe_col = loadData$apoe_col
pathology_col = loadData$pathology_col
protein_cols = loadData$protein_cols
result_folder = loadData$result_folder
cov_cols = loadData$cov_cols
age_col = loadData$age_col

get_demographic(analysis_df, dataset_identifier, apoe_identifier)

spc25 = 'SPC25__seq.22782.80__ANML'
drug_dict = list(
  "drugs_platelet_inhibitors_baseline_variable"       = "Platelet inhibitors",
  "drugs_antidepressants_baseline_variable"           = "Antidepressants",
  "drugs_antiinflammatory_baseline_variable"          = "Antiinflammatory",
  "drugs_hypertension_cardioprotective_baseline_variable" = "Hypertension cardioprotective",
  "drugs_lipid_lowering_baseline_variable"            = "Lipid lowering",
  "drugs_cholinesterase_inhibitors_baseline_variable" = "Cholinesterase inhibitors"
  )

##-------------------------SPC25-APOE4 in Ab-negative ChEI non user-------------
non_user = analysis_df[(analysis_df$drugs_cholinesterase_inhibitors_baseline_variable == 0) & 
                         (analysis_df$Abnormal_CSF_Ab42_Ab40_Ratio == 0),]
res = apoeOrAb2protein(non_user, 
                 c(spc25), 
                 'apoe_class_num', 
                 c("age",
                   "gender_baseline_variable",
                   "plasma_ANML_mean"),
                 dataset_identifier)
#p = 3.818665e-174

##-----------------------------Mediation analysis-------------------------------
all_mediations = data.frame()
for (drug in names(drug_dict)) {
   cov_cols = c(drug, "age","gender_baseline_variable","plasma_ANML_mean")
   mediations = protein_mediate_analysis(data_df = analysis_df,
                                         cov_cols = cov_cols,
                                         pathology_col = pathology_col,
                                         apoe_col = apoe_col,
                                         proteins = c(spc25),
                                         dataset_identifier = dataset_identifier)
  
   mediations$drug = drug
   all_mediations = bind_rows(all_mediations, mediations)
 }

##------------------------------SPC25 ~ ChEI stratify by APOE4 carriers---------
chei_inapoe_neg = apoeOrAb2protein(analysis_df[analysis_df$apoe_class_num == 0,], c(spc25), 
                 'drugs_cholinesterase_inhibitors_baseline_variable',
                 c('Abnormal_CSF_Ab42_Ab40_Ratio', "age","gender_baseline_variable","plasma_ANML_mean"),
                 dataset_identifier) #p = 0.2989068

chei_inapoe_pos = apoeOrAb2protein(analysis_df[analysis_df$apoe_class_num == 1,], 
                                   c(spc25), 
                                   'drugs_cholinesterase_inhibitors_baseline_variable',
                                   c('Abnormal_CSF_Ab42_Ab40_Ratio',
                                     "age",
                                     "gender_baseline_variable",
                                     "plasma_ANML_mean"),
                                   dataset_identifier) #p = 0.001718218

##------------SPC25 change group by APOE4 * Ab * ChEI---------------------------
plot_df = analysis_df[,]
plot_df = plot_df %>% mutate(
  ab = case_when(
    is.na(Abnormal_CSF_Ab42_Ab40_Ratio) ~ NA_character_,
    Abnormal_CSF_Ab42_Ab40_Ratio == 0 ~ 'Aβ-',
    Abnormal_CSF_Ab42_Ab40_Ratio == 1 ~ 'Aβ+',
    TRUE ~ NA_character_
  ),
  ChEI = case_when(
    is.na(drugs_cholinesterase_inhibitors_baseline_variable) ~ NA_character_,
    drugs_cholinesterase_inhibitors_baseline_variable == 0 ~ 'ChEI-',
    drugs_cholinesterase_inhibitors_baseline_variable == 1 ~ 'ChEI+',
    TRUE ~ NA_character_
  ),
  apoe4 = case_when(
    apoe_class_num == 0 ~ 'ε4-',
    apoe_class_num == 1 ~ 'ε4+',
    TRUE ~ NA_character_
  )
)

model_apoe_drug_inter = lm(SPC25__seq.22782.80__ANML ~ apoe4 * ChEI * ab
                           + age + gender_baseline_variable + plasma_ANML_mean ,
                           data = plot_df)

p0 = proteinlevelby_auto(plot_df, spc25, 
                         c("age","gender_baseline_variable","plasma_ANML_mean", pathology_col, 'ChEI'), 
                         'apoe4', 'Adjustment for Aβ and ChEIs', 'Residual (SPC25)',
                         '', legend_flag=TRUE)

p1 = proteinlevelby_auto(plot_df, spc25, 
                         c("age","gender_baseline_variable","plasma_ANML_mean", 
                           apoe_col), 
                         'ab', 'Adjustment for APOE4', 'Residual (SPC25)',
                         '', legend_flag=TRUE)
p2 = proteinlevelby_auto(plot_df, spc25, 
                         c("age","gender_baseline_variable","plasma_ANML_mean", 
                           apoe_col,'drugs_cholinesterase_inhibitors_baseline_variable'), 
                         'ab', 'Adjustment for APOE4 and ChEIs', 'Residual (SPC25)',
                         '', legend_flag=TRUE)
p3 = proteinlevelby_auto(plot_df, spc25, 
                         c("age","gender_baseline_variable","plasma_ANML_mean", apoe_col, pathology_col), 
                         'ChEI', 'Adjustment for Aβ and APOE4', 'Residual (SPC25)',
                         '', legend_flag=TRUE) 

p0$plot+p1$plot+p2$plot+p3$plot+plot_layout(nrow = 1)
ggsave(glue('Figs/{dataset_identifier}/spc25_drug.svg'),
       p0$plot+p1$plot+p2$plot+p3$plot+plot_layout(nrow = 1),
       width = 20,
       height = 6,
       dpi = 300)
