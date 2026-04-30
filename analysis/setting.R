source('APOE_proteomics/analysis/src.R')

#-------------------------------Load data---------------------------------------
main_driver = function(dataset_identifier=c('GNPC',
                                            'BF2SomaLogic', 
                                            'ADNI', 
                                            'ADNI_MS',
                                            'BF2OLINK', 
                                            'UKBB', 
                                            'ROSMAP'), 
                       apoe_identifier=c('APOE4', 
                                         'APOE2'), 
                       analysis = c('Aβ', 
                                    'AD', 
                                    'PC10', 
                                    'AbRatio', 
                                    'drug', 
                                    'wml', 
                                    'AbPET')){
  if (dataset_identifier == "GNPC"){
    handle_gnpc(apoe_identifier)
  } else if (dataset_identifier == "BF2SomaLogic"){
    # Sensitivity analysis setting
    if(analysis %in% c('PC10')){
      handle_bf2somalogic_pc10(apoe_identifier)
      
    } else if(analysis == 'drug') {
      handle_bf2somalogic_drug(apoe_identifier)
      
    } else if(analysis == 'wml') {
      handle_bf2somalogic_wml(apoe_identifier)
      
    } else if(analysis %in% c('AbRatio')){
      handle_bf2somalogic_csfAbRatio(apoe_identifier)
      
    } else if(analysis == 'AbPET') {
      handle_bf2somalogic_AbPET(apoe_identifier)
      
    }
    else {
      handle_bf2somalogic_main(apoe_identifier, analysis)
    }
  } else if (dataset_identifier == "ADNI"){
    handle_adni(apoe_identifier, analysis)
    
  } else if (dataset_identifier == "ADNI_MS"){
    handle_adni_ms(apoe_identifier)
    
  } else if (dataset_identifier == "BF2OLINK"){
    handle_bf2olink(apoe_identifier, analysis)
    
  } else if (dataset_identifier == "UKBB"){
    handle_ukbb(apoe_identifier)
    
  } else if (dataset_identifier == 'ROSMAP'){
    handle_rosmap(apoe_identifier)
    
  } else {
    stop("Unknown dataset identifier")
  }
}

handle_gnpc = function(apoe_identifier){
  age_col = 'age_at_visit'
  cov_cols = c(age_col, 'sex', 'plasma_ANML_mean', 'contributor_code')
  pathology_col = 'diagnosis'
  protein_cols = grep('^seq', names(initial), value = TRUE)
  
  if (apoe_identifier == 'APOE4'){
    apoe_col = 'apoe_class_num'
    result_folder = '../results/e4vse3e3'
    analysis_df = subset(initial, 
                         APOE %in% c(33, 34, 43, 44) & contributor_code != 'C')
  } else {
    apoe_col = 'apoe_class_num_e2'
    result_folder = '../results/e2vse3e3/'
    analysis_df = subset(initial, 
                         APOE %in% c(22, 23, 32, 33) & contributor_code != 'C')
  }
  
  return(list(
    df = analysis_df,
    cov_cols = cov_cols,
    age_col = age_col,
    apoe_col = apoe_col,
    pathology_col = pathology_col,
    protein_cols = protein_cols,
    result_folder = result_folder
  ))
}

handle_bf2somalogic_pc10 = function(apoe_identifier){
  initial = read.csv('data/BF2SomaLogic/preprocessed/with24_addpc10.csv')
  
  get_demographic(initial, 'BF2SomaLogic', apoe_identifier, FALSE)
  
  cov_cols = c('age', 'gender_baseline_variable', 'plasma_ANML_mean', 
               'PC1', 'PC2', 'PC3', 'PC4', 'PC5')
  age_col = 'age'
  pathology_col = 'Abnormal_CSF_Ab42_Ab40_Ratio'
  
  if (apoe_identifier == 'APOE4'){
    apoe_col = 'apoe_class_num'
    analysis_df = subset(initial, apoe_genotype_baseline_variable %in% c(33, 34, 43, 44))
    result_folder = 'APOE_proteomics/results/BF2SomaLogic/e4vse3e3_pc10/'
  } else {
    apoe_col = 'apoe_class_num_e2'
    analysis_df = subset(initial, apoe_genotype_baseline_variable %in% c(22, 23, 32, 33))
    result_folder = 'APOE_proteomics/results/BF2SomaLogic/e2vse3e3_pc10/'
  }
  
  protein_cols = grep('ANML$', names(analysis_df), value = TRUE)
  
  return(list(
    df = analysis_df,
    cov_cols = cov_cols,
    age_col = age_col,
    apoe_col = apoe_col,
    pathology_col = pathology_col,
    protein_cols = protein_cols,
    result_folder = result_folder
  ))
}

handle_bf2somalogic_csfAbRatio = function(apoe_identifier){
  initial = read.csv('data/BF2SomaLogic/preprocessed/with24_7285prot.csv')
  additional_abeta_data = read_csv("data/BF2SomaLogic/preprocessed/additional_Abeta_data.csv")
  ab_ratiao_dict   = deframe(dplyr::select(additional_abeta_data, sid, csf_clinical_routine_Abeta42_40_ratio_x10))
  initial = initial %>%
    mutate(
      csf_clinical_routine_Abeta42_40_ratio_x10 = ab_ratiao_dict[sid],
    ) %>%
    ungroup()
  
  get_demographic(initial, "BF2SomaLogic", apoe_identifier, FALSE)
  
  age_col = 'age'
  pathology_col = 'csf_clinical_routine_Abeta42_40_ratio_x10'
  cov_cols = c('age', 'gender_baseline_variable', 'plasma_ANML_mean')
  
  if (apoe_identifier == 'APOE4'){
    apoe_col = 'apoe_class_num'
    analysis_df = subset(initial, apoe_genotype_baseline_variable %in% c(33, 34, 43, 44))
    result_folder = 'APOE_proteomics/results/BF2SomaLogic/e4vse3e3_AbRatio/'
  } else {
    apoe_col = 'apoe_class_num_e2'
    analysis_df = subset(initial, apoe_genotype_baseline_variable %in% c(22, 23, 32, 33))
    result_folder = 'APOE_proteomics/results/BF2SomaLogic/e2vse3e3_AbRatio/'
  }
  
  analysis_df = analysis_df[!is.na(analysis_df['csf_clinical_routine_Abeta42_40_ratio_x10']),]
  protein_cols = grep('ANML$', names(analysis_df), value = TRUE)
  
  return(list(
    df = analysis_df,
    cov_cols = cov_cols,
    age_col = age_col,
    apoe_col = apoe_col,
    pathology_col = pathology_col,
    protein_cols = protein_cols,
    result_folder = result_folder
  ))
}

handle_bf2somalogic_AbPET = function(apoe_identifier){
  initial = read.csv('data/BF2SomaLogic/preprocessed/with24_7285prot.csv')
  additional_abeta_data = read_csv("data/BF2SomaLogic/preprocessed/additional_Abeta_data.csv")
  abPET_dict   = deframe(dplyr::select(additional_abeta_data, sid, fnc_ber_com_composite))
  initial = initial %>%
    mutate(
      fnc_ber_com_composite = abPET_dict[sid],
    ) %>%
    ungroup()
  
  initial = initial[!is.na(initial['fnc_ber_com_composite']),]
  
  get_demographic(initial, "BF2SomaLogic", apoe_identifier, FALSE)
  
  age_col = 'age'
  pathology_col = 'fnc_ber_com_composite'
  cov_cols = c('age', 'gender_baseline_variable', 'plasma_ANML_mean')
  
  if (apoe_identifier == 'APOE4'){
    apoe_col = 'apoe_class_num'
    analysis_df = subset(initial, apoe_genotype_baseline_variable %in% c(33, 34, 43, 44))
    result_folder = 'APOE_proteomics/results/BF2SomaLogic/e4vse3e3_AbPET/'
  } else {
    apoe_col = 'apoe_class_num_e2'
    analysis_df = subset(initial, apoe_genotype_baseline_variable %in% c(22, 23, 32, 33))
    result_folder = 'APOE_proteomics/results/BF2SomaLogic/e2vse3e3_AbPET/'
  }
  
  protein_cols = grep('ANML$', names(analysis_df), value = TRUE)
  
  return(list(
    df = analysis_df,
    cov_cols = cov_cols,
    age_col = age_col,
    apoe_col = apoe_col,
    pathology_col = pathology_col,
    protein_cols = protein_cols,
    result_folder = result_folder
  ))
}

handle_bf2somalogic_wml = function(apoe_identifier){
  initial = read.csv('data/BF2SomaLogic/preprocessed/with24_7285prot.csv')
  
  get_demographic(initial, "BF2SomaLogic", apoe_identifier, FALSE)
  
  vascular = c('samseg_wmhs_WMH_total_mm3', 'icv_mm3')
  cov_cols = c('age', 'gender_baseline_variable', 'plasma_ANML_mean', vascular)
  age_col = 'age'
  pathology_col = 'Abnormal_CSF_Ab42_Ab40_Ratio'
  
  if (apoe_identifier == 'APOE4'){
    apoe_col = 'apoe_class_num'
    analysis_df = subset(initial, apoe_genotype_baseline_variable %in% c(33, 34, 43, 44))
    result_folder = 'APOE_proteomics/results/BF2SomaLogic/e4vse3e3_wml/'
  } else {
    apoe_col = 'apoe_class_num_e2'
    analysis_df = subset(initial, apoe_genotype_baseline_variable %in% c(22, 23, 32, 33))
    result_folder = 'APOE_proteomics/results/BF2SomaLogic/e2vse3e3_wml/'
  }
  
  protein_cols = grep('ANML$', names(analysis_df), value = TRUE)
  
  return(list(
    df = analysis_df,
    cov_cols = cov_cols,
    age_col = age_col,
    apoe_col = apoe_col,
    pathology_col = pathology_col,
    protein_cols = protein_cols,
    result_folder = result_folder
  ))
}

handle_bf2somalogic_drug = function(apoe_identifier){
  initial = read.csv('data/BF2SomaLogic/preprocessed/with24_7285prot.csv')
  
  get_demographic(initial, "BF2SomaLogic", apoe_identifier, FALSE)
  
  drugs = grep('^drug', names(initial), value = TRUE)
  cov_cols = c('age', 'gender_baseline_variable', 'plasma_ANML_mean', drugs)
  age_col = 'age'
  pathology_col = 'Abnormal_CSF_Ab42_Ab40_Ratio'
  
  if (apoe_identifier == 'APOE4'){
    apoe_col = 'apoe_class_num'
    analysis_df = subset(initial, apoe_genotype_baseline_variable %in% c(33, 34, 43, 44))
    result_folder = 'APOE_proteomics/results/BF2SomaLogic/e4vse3e3_drug/'
    } else {
    apoe_col = 'apoe_class_num_e2'
    analysis_df = subset(initial, apoe_genotype_baseline_variable %in% c(22, 23, 32, 33))
    result_folder = 'APOE_proteomics/results/BF2SomaLogic/e2vse3e3_drug/'
  }
  
  protein_cols = grep('ANML$', names(analysis_df), value = TRUE)
  
  return(list(
    df = analysis_df,
    cov_cols = cov_cols,
    age_col = age_col,
    apoe_col = apoe_col,
    pathology_col = pathology_col,
    protein_cols = protein_cols,
    result_folder = result_folder
  ))
}

handle_bf2somalogic_main = function(apoe_identifier, analysis){
  initial = read.csv('data/BF2SomaLogic/preprocessed/with24_7285prot.csv')
  
  get_demographic(initial, 'BF2SomaLogic', apoe_identifier, FALSE)
  
  cov_cols = c('age', 'gender_baseline_variable', 'plasma_ANML_mean')
  age_col = 'age'
  
  if (apoe_identifier == 'APOE4'){
    apoe_col = 'apoe_class_num'
    analysis_df = subset(initial, apoe_genotype_baseline_variable %in% c(33, 34, 43, 44))
    
    if (analysis == 'Aβ'){
      pathology_col = 'Abnormal_CSF_Ab42_Ab40_Ratio'
      result_folder = 'APOE_proteomics/results/BF2SomaLogic/e4vse3e3/'
    } else if (analysis == 'AD'){
      pathology_col = 'diagnosis_baseline_variable'
      analysis_df = subset(analysis_df, analysis_df[[pathology_col]] %in% c('CU', 'AD'))
      analysis_df[[pathology_col]] = ifelse(analysis_df[[pathology_col]] == 'CU', 0 ,1)
      result_folder = 'APOE_proteomics/results/BF2SomaLogic/e4vse3e3_AD/'
    }
  } else {
    apoe_col = 'apoe_class_num_e2'
    analysis_df = subset(initial, apoe_genotype_baseline_variable %in% c(22, 23, 32, 33))
    
    if (analysis == 'Aβ'){
      pathology_col = 'Abnormal_CSF_Ab42_Ab40_Ratio'
      result_folder = 'APOE_proteomics/results/BF2SomaLogic/e2vse3e3/'
    } else if (analysis == 'AD'){
      pathology_col = 'diagnosis_baseline_variable'
      analysis_df = subset(analysis_df, analysis_df[[pathology_col]] %in% c('CU', 'AD'))
      analysis_df[[pathology_col]] = ifelse(analysis_df[[pathology_col]] == 'CU', 0 ,1)
      result_folder = 'APOE_proteomics/results/BF2SomaLogic/e2vse3e3_AD/'
    }
  }
  protein_cols = grep('ANML$', names(analysis_df), value = TRUE)
  
  return(list(
    df = analysis_df,
    cov_cols = cov_cols,
    age_col = age_col,
    apoe_col = apoe_col,
    pathology_col = pathology_col,
    protein_cols = protein_cols,
    result_folder = result_folder
  ))
}

handle_adni = function(apoe_identifier, analysis){
  initial = read.csv('data/ADNI/preprocessed/adni_humanprotein_align.csv')
  initial$DX_bl = ifelse(initial$DX_bl %in% c('EMCI', 'LMCI'), 'MCI', initial$DX_bl)
  initial$AB_CSF_status = ifelse(initial$AB_CSF_status == 'AB_Neg', 0, ifelse(
    initial$AB_CSF_status == 'AB_Pos', 1, NA
  ))
  
  get_demographic(initial[initial$APOE != 24, ], 'ADNI', apoe_identifier, FALSE)
  
  cov_cols = c('AGE', 'PTGENDER', 'CSF_ANML_mean')
  age_col = 'AGE'
  
  if (apoe_identifier == 'APOE4'){
    apoe_col = 'apoe_class_num'
    analysis_df = subset(initial, APOE %in% c(33, 34, 43, 44))
    
    if (analysis == 'Aβ'){
      pathology_col = 'AB_CSF_status'
      result_folder = 'APOE_proteomics/results/ADNI/e4vse3e3/'
    } else if (analysis == 'AD'){
      pathology_col = 'DX_bl'
      analysis_df = subset(analysis_df, analysis_df[[pathology_col]] %in% c('CU', 'AD'))
      analysis_df[[pathology_col]] = ifelse(analysis_df[[pathology_col]] == 'CU', 0 ,1)
      result_folder = 'APOE_proteomics/results/ADNI/e4vse3e3_AD/'
    }
  } else {
    apoe_col = 'apoe_class_num_e2'
    analysis_df = subset(initial, APOE %in% c(22, 23, 32, 33))
    
    if (analysis == 'Aβ'){
      pathology_col = 'AB_CSF_status'
      result_folder = 'APOE_proteomics/results/ADNI/e2vse3e3/'
    } else if (analysis == 'AD'){
      pathology_col = 'DX_bl'
      analysis_df = subset(analysis_df, analysis_df[[pathology_col]] %in% c('CU', 'AD'))
      analysis_df[[pathology_col]] = ifelse(analysis_df[[pathology_col]] == 'CU', 0 ,1)
      result_folder = 'APOE_proteomics/results/ADNI/e2vse3e3_AD/'
    }
  }
  protein_cols = grep('ANML$', names(analysis_df), value = TRUE)
  
  return(list(
    df = analysis_df,
    cov_cols = cov_cols,
    age_col = age_col,
    apoe_col = apoe_col,
    pathology_col = pathology_col,
    protein_cols = protein_cols,
    result_folder = result_folder
  ))
}

handle_adni_ms = function(apoe_identifier){
  initial = read.csv('data/ADNI/preprocessed/adni_soma_ms_merged.csv')
  
  initial$DX_bl = ifelse(initial$DX_bl %in% c('EMCI', 'LMCI'), 'MCI', initial$DX_bl)
  initial$AB_CSF_status = ifelse(initial$AB_CSF_status == 'AB_Neg', 0, ifelse(
    initial$AB_CSF_status == 'AB_Pos', 1, NA
  ))
  
  get_demographic(initial[initial$APOE != 24, ], 'ADNI', apoe_identifier, FALSE)
  
  cov_cols = c('AGE', 'PTGENDER', 'CSF_TMT_mean')
  age_col = 'AGE'
  
  if (apoe_identifier == 'APOE4'){
    apoe_col = 'apoe_class_num'
    analysis_df = subset(initial, APOE %in% c(33, 34, 43, 44))
    pathology_col = 'AB_CSF_status'
    result_folder = 'APOE_proteomics/results/ADNI_MS/e4vse3e3/'
  } else {
    apoe_col = 'apoe_class_num_e2'
    analysis_df = subset(initial, APOE %in% c(22, 23, 32, 33))
    pathology_col = 'AB_CSF_status'
    result_folder = 'APOE_proteomics/results/ADNI_MS/e2vse3e3/'
  }
  
  protein_cols = grep('_MS$', names(analysis_df), value = TRUE)

  return(list(
    df = analysis_df,
    cov_cols = cov_cols,
    age_col = age_col,
    apoe_col = apoe_col,
    pathology_col = pathology_col,
    protein_cols = protein_cols,
    result_folder = result_folder
  ))
}

handle_bf2olink = function(apoe_identifier, analysis){
  initial = read.csv('data/BF2/preprocessed/with24_1391proteins.csv')
  
  get_demographic(initial, 'BF2OLINK', apoe_identifier, FALSE)
  
  cov_cols = c('age', 'gender_baseline_variable', 'CSF_NPX_mean')
  age_col = 'age'
  
  if (apoe_identifier == 'APOE4'){
    apoe_col = 'apoe_class_num'
    analysis_df = subset(initial, apoe_genotype_baseline_variable %in% c(33, 34, 43, 44))
    
    if (analysis == 'Aβ'){
      pathology_col = 'Abnormal_CSF_Ab42_Ab40_Ratio'
      result_folder = 'APOE_proteomics/results/BF2OLINK/e4vse3e3/'
    } else if (analysis == 'AD'){
      pathology_col = 'diagnosis_baseline_variable'
      analysis_df = subset(analysis_df, analysis_df[[pathology_col]] %in% c('CU', 'AD'))
      analysis_df$DX_bl_fac = ifelse(analysis_df[[pathology_col]] == 'CU', 0 ,1)
      pathology_col = 'DX_bl_fac'
      result_folder = 'APOE_proteomics/results/BF2OLINK/e4vse3e3_AD/'
    }
  } else {
    apoe_col = 'apoe_class_num_e2'
    analysis_df = subset(initial, apoe_genotype_baseline_variable %in% c(22, 23, 32, 33))
    
    if (analysis == 'Aβ'){
      pathology_col = 'Abnormal_CSF_Ab42_Ab40_Ratio'
      result_folder = 'APOE_proteomics/results/BF2OLINK/e2vse3e3/'
    } else if (analysis == 'AD'){
      pathology_col = 'diagnosis_baseline_variable'
      analysis_df = subset(analysis_df, analysis_df[[pathology_col]] %in% c('CU', 'AD'))
      analysis_df$DX_bl_fac = ifelse(analysis_df[[pathology_col]] == 'CU', 0 ,1)
      pathology_col = 'DX_bl_fac'
      result_folder = 'APOE_proteomics/results/BF2OLINK/e2vse3e3_AD/'
    }
  }
  protein_cols = grep('NPX$', names(analysis_df), value = TRUE)
  
  return(list(
    df = analysis_df,
    cov_cols = cov_cols,
    age_col = age_col,
    apoe_col = apoe_col,
    pathology_col = pathology_col,
    protein_cols = protein_cols,
    result_folder = result_folder
  ))
}

handle_ukbb = function(apoe_identifier){
  initial = read.csv('data/UKBB/preprocessed/ukbb_clean.csv')
  withdraw = read_table("data/UKBB/raw_data/withdraw105777_432_20251223.txt", col_names = FALSE)
  initial = initial[!initial$eid %in% withdraw$X1, ]
  
  get_demographic(initial, 'UKBB', apoe_identifier, FALSE)
  
  cov_cols = c('age', 'gender', 'Plasma_NPX_mean')
  protein_cols = grep('^CSF.*NPX$', names(initial), value = TRUE)
  
  if (apoe_identifier == 'APOE4') {
    apoe_col = 'apoe_class_num'
    analysis_df = subset(initial, APOE %in% c(33, 34, 43, 44))
    result_folder = 'APOE_proteomics/results/UKBB/e4vse3e3/'
  } else if (apoe_identifier == 'APOE2') {
    apoe_col = 'apoe_class_num_e2'
    analysis_df = subset(initial, APOE %in% c(22, 23, 32, 33))
    result_folder = 'APOE_proteomics/results/UKBB/e2vse3e3/'
    }
  
  return(list(
    df = analysis_df,
    cov_cols = cov_cols,
    age_col = NULL,
    apoe_col = apoe_col,
    pathology_col = NULL,
    protein_cols = protein_cols,
    result_folder = result_folder
  ))
}

handle_rosmap = function(apoe_identifier){
  initial = read.csv('data/ROSMAP/preprocessed/rosmap_apoe_mediators.csv')
  initial = initial %>% filter(
    apoe_genotype != 24,
    diagnosis_at_death %in% c('CU', 'AD')
  ) %>% mutate(
    AD_diagnosis = case_when(
      diagnosis_at_death == 'CU' ~ 0,
      diagnosis_at_death == 'AD' ~ 1
    )
  )
  
  get_demographic(initial, 'ROSMAP', apoe_identifier, FALSE)
  
  cov_cols = c('age_death_num', 'sex_fac', 'pmi')
  protein_cols = grep('__TPM$', names(initial), value = TRUE)
  
  if (apoe_identifier == 'APOE4') {
    apoe_col = 'apoe4_num'
    analysis_df = initial[initial$apoe_genotype %in% c(33, 34, 43, 44), ]
    result_folder = 'APOE_proteomics/results/ROSMAP/e4vse3e3/'
  } else if (apoe_identifier == 'APOE2') {
    apoe_col = 'apoe2_num'
    analysis_df = initial[initial$apoe_genotype %in% c(22, 23, 32, 33), ]
    result_folder = 'APOE_proteomics/results/ROSMAP/e2vse3e3/'
  }
  
  return(list(
    df = analysis_df,
    cov_cols = cov_cols,
    age_col = NULL,
    apoe_col = apoe_col,
    pathology_col = 'AD_diagnosis',
    protein_cols = protein_cols,
    result_folder = result_folder
  ))
}

#----------------------------------Pre-defined----------------------------------
apoe_list = list(
  APOE2 = 'e2vse3e3',
  APOE4 = 'e4vse3e3'
)