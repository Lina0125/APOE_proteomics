library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(ggsci)
library(mediation)
library(ggpp)
library(QuantPsyc)
library(car)
library(readxl)
library(stringr)
library(data.table)
library(table1)
library(interactions)
library(ggpubr)
library(RColorBrewer)
library(ggsci)
library(rstatix)
library(ggtext)
library(tidyr)
library(readr)
library(tibble)
library(glue)
library(openxlsx)
library(cowplot)

set.seed(123)

#---------------------------------Accessibility---------------------------------
get_demographic = function(df, dataset, apoe = 'APOE4', grouped = TRUE) {
  set_apoe_col = function(df, apoe, 
                          e4_col = "apoe_class_num", 
                          e2_col = "apoe_class_num_e2", 
                          e4_fallback_col = "apoe_class_num_e4") {
    if (apoe == 'APOE2') {
      df$apoe_col = factor(df[[e2_col]], levels = c(0, 1), labels = c('ε2-', 'ε2+'))
    } else {
      col = ifelse(e4_fallback_col %in% names(df), e4_fallback_col, e4_col)
      df$apoe_col = factor(df[[col]], levels = c(0, 1), labels = c('ε4-', 'ε4+'))
    }
    return(df)
  }
  
  # dataset-specific preprocessing
  if (dataset %in% c('BF2SomaLogic', 'BF2OLINK')) {
    df$apoe_genotype_baseline_variable = factor(df$apoe_genotype_baseline_variable)
    df$Abnormal_CSF_Ab42_Ab40_Ratio = factor(df$Abnormal_CSF_Ab42_Ab40_Ratio, 
                                             levels = c(0, 1), 
                                             labels = c('Aβ-', 'Aβ+'))
    df$gender_baseline_variable = factor(df$gender_baseline_variable, 
                                         levels = c(0, 1), 
                                         labels = c('Male', 'Female'))
    df$diagnosis_baseline_variable = factor(df$diagnosis_baseline_variable)
    df = set_apoe_col(df, apoe)
    
    formula = ~ Abnormal_CSF_Ab42_Ab40_Ratio + gender_baseline_variable + age + apoe_genotype_baseline_variable + diagnosis_baseline_variable
  } 
  else if (dataset == 'ADNI' || dataset == 'ADNI_MS') {
    df$PTGENDER = factor(df$PTGENDER, levels = c(1, 2), labels = c('Male', 'Female'))
    df$AB_CSF_status = factor(df$AB_CSF_status, levels = c(0, 1), labels = c('Aβ-', 'Aβ+'))
    df$APOE = factor(df$APOE)
    df = set_apoe_col(df, apoe)
    
    formula = ~ AB_CSF_status + PTGENDER + AGE + APOE + DX_bl
  } 
  else if (dataset == 'UKBB') {
    df$gender = factor(df$gender, levels = c(0, 1), labels = c('Male', 'Female'))
    df$APOE = factor(df$APOE)
    df = set_apoe_col(df, apoe)
    
    formula = ~ gender + age + APOE
  } else if (dataset == 'ROSMAP') {
    df$apoe_genotype = factor(df$apoe_genotype)
    df = set_apoe_col(df, apoe, e4_col = "apoe4_num", e2_col = "apoe2_num")
    
    formula = ~ sex_fac + age_death_num + apoe_genotype + pmi + diagnosis_at_death
  } else if (dataset %in% c('SEER CSF', 'SEER Plasma')) {
    df = set_apoe_col(df, apoe, e4_col = "apoe_class_num", e4_fallback_col = "apoe_class_num_e4")
    formula = ~ Diagnosis + as.factor(CSF.Ab.status) + as.factor(APOE) + as.factor(diagnosis_Gates_CU_CI)
  } 
  else {
    stop("Unsupported dataset")
  }
  
  if (grouped) {
    rhs = as.character(formula)[2]
    full_formula = as.formula(paste("~", rhs, "| apoe_col"))
    print(table1(full_formula, data = df, caption = dataset, overall = "Total"))
  } else {
    print(table1(formula, data = df, caption = dataset, overall = "Total"))
  }
}

map_symbol = function(input, dataset_identifier, bycol='Protein_id'){
  # Proteins are measured using multiple assays
  if(dataset_identifier == 'BF2OLINK'){
    symbol_map = read_csv("ref/olink_label_mapping.csv")
    input = merge(input, symbol_map, by= bycol, all.x = TRUE)
  } else if (dataset_identifier == 'ADNI_MS'){
    symbol_map = read_csv("../ADNI/data/preprocessed/ADNI_CSFTMTMS_labelMatch.csv")
    input = merge(input, symbol_map, by= bycol, all.x = TRUE)
  } else if (dataset_identifier == 'UKBB'){
    # Non duplicated gene symbols
    input$symbol = sapply(strsplit(input[[bycol]], "_"), `[`, 2)
    input$label = input$symbol
  } else if (dataset_identifier == 'ROSMAP'){
    # Non duplicated gene symbols
    input$symbol = sub("__TPM$", "", input[[bycol]])
    input$label = input$symbol
  } else {
    # Many proteins are measured using multiple assays
    symbol_map = read_excel("ref/symbol_mapping_across_platform.xlsx")
    input$apt_name = gsub(".*__(seq.\\d+.\\d+)__.*", "\\1", input[[bycol]])
    input = merge(input, symbol_map, by= "apt_name", all.x = TRUE)
  }
  
  return(input)
}

format_pval = function(p) {
  if (p < 0.001) {
    format(p, scientific = TRUE, digits = 2)
  } else {
    format(round(p, 3), nsmall = 3)
  }
}

#-----------------------------Analysis functions--------------------------------
qc_protein = function(qc_df, protein_col) {
  x = qc_df[[protein_col]]
  high_bound = mean(x, na.rm = TRUE) + 5 * sd(x, na.rm = TRUE)
  low_bound = mean(x, na.rm = TRUE) - 5 * sd(x, na.rm = TRUE)
  qc_df = qc_df[x <= high_bound & x >= low_bound, , drop = FALSE]
  return(qc_df)
}

model_to_rows = function(model, x_col, method){
  model_summary = summary(model)
  
  if (method == 'lm'){
    
    standardized = lm.beta(model)
    coef_table = model_summary$coefficients
    
    if (x_col %in% rownames(coef_table)) {
      X_coef = coef_table[x_col, "Estimate"]
      X_se = coef_table[x_col, "Std. Error"]
      X_t = coef_table[x_col, "t value"]
      X_p = coef_table[x_col, "Pr(>|t|)"]
      
      r_squared = model_summary$r.squared
      r_squared_adj = model_summary$adj.r.squared
      model_resid = model_summary$df[2]
      
      conf_int = confint(model, level = 0.95)
      X_95CI_lower = conf_int[x_col, 1]
      X_95CI_upper = conf_int[x_col, 2]
      
      result_row = data.frame(
        coef = X_coef,
        standardized_beta = standardized[x_col],
        se = X_se,
        t = X_t,
        p = X_p,
        R_squared = r_squared,
        R_squared_adj = r_squared_adj,
        resid = model_resid,
        CI95_lower = X_95CI_lower,
        CI95_upper = X_95CI_upper,
        AIC = AIC(model)
      )
    } else {
      result_row = data.frame(log = 'Unable to calculate beta')
    }
    
  } else if (method == 'logit'){
    coef_table = coef(model_summary)
    
    if (x_col %in% rownames(coef_table)) {
      conf_int = confint(model)
      
      coef_val = coef_table[x_col, "Estimate"]
      se_val = coef_table[x_col, "Std. Error"]
      z_val = coef_table[x_col, "z value"]
      p_val = coef_table[x_col, "Pr(>|z|)"]
      
      OR = exp(coef_val)
      CI_lower = exp(conf_int[x_col, 1])
      CI_upper = exp(conf_int[x_col, 2])
      
      result_row = data.frame(
        coef = coef_val,
        OddsRatio = OR,
        se = se_val,
        z = z_val,
        p = p_val,
        CI95_lower = CI_lower,
        CI95_upper = CI_upper,
        AIC = AIC(model)
      )
    } else {
      result_row = data.frame(log = 'Unable to calculate beta')
    }
  }
  
  return(result_row)
}

#----------------------------------Modelling------------------------------------
apoeOrAb2protein = function(dat, protein_cols, x_col, covs, dataset_identifier){
  result_list = list()
  for (protein_col in protein_cols) {
    cols = c(covs, x_col, protein_col)
    forqc = na.omit(dat[cols])
    qced = qc_protein(forqc, protein_col)
    
    model_str = paste(protein_col, '~', paste(c(x_col, covs), collapse = ' + '))
    print(model_str)
    
    formula = as.formula(model_str)
    model = lm(formula, data=qced)
    
    result_row = model_to_rows(model, x_col, 'lm')
    result_row$Protein_id = protein_col
      
    result_list[[protein_col]] = result_row
  }
  
  final_results = bind_rows(result_list)
  final_results$p_adjusted = p.adjust(final_results$p, method = 'fdr')

  final_results = map_symbol(final_results, dataset_identifier)
  
  return(final_results)
}

# APOE => continuous =>  binary
protein_mediate_analysis = function(data_df, cov_cols, pathology_col, apoe_col, proteins, dataset_identifier){
  mediation_results = list()
  for (protein_col in proteins) {
    cols = c(cov_cols, pathology_col, apoe_col, protein_col)
    forqc = na.omit(data_df[cols])
    setnames(forqc, protein_col, 'protein')
    qced = qc_protein(forqc, 'protein')
    
    print(paste("protein ~", apoe_col, "+", 
                paste(cov_cols, collapse = ' + ')))
    model_m = lm(as.formula(paste("protein ~", apoe_col, "+", 
                                  paste(cov_cols, collapse = ' + '))), data = qced)
    model_y = glm(as.formula(paste(pathology_col, "~", apoe_col, "+ protein +", 
                                   paste(cov_cols, collapse = ' + '))), 
                  data = qced, family = binomial())
    
    mediation_result = mediate(model.m = model_m, model.y = model_y, 
                               treat = apoe_col, mediator = 'protein', 
                               boot = TRUE, sims = 1000)
    
    result = data.frame(
      coef_x2m = coef(model_m)[[2]],
      coef_x2y = coef(model_y)[[2]],
      coef_m2y = coef(model_y)[[3]],
      acme = mediation_result$d.avg,
      acme_p = mediation_result$d.avg.p,
      ade = mediation_result$z.avg,
      ade_p = mediation_result$z.avg.p,
      total = mediation_result$tau.coef,
      total_p = mediation_result$tau.p,
      prop_m = mediation_result$n.avg,
      prop_m_p = mediation_result$n.avg.p,
      prop_m_ci_lower = mediation_result$n.avg.ci[1],
      prop_m_ci_upper = mediation_result$n.avg.ci[2]
    )
    mediation_results[[protein_col]] = result
  }
  
  final_results = do.call(rbind, mediation_results)
  final_results$p_adjusted = p.adjust(final_results$acme_p, method = 'fdr')
  final_results$ade_p_adjusted = p.adjust(final_results$ade_p, method = 'fdr')
  final_results$prop_m_p_adjusted = p.adjust(final_results$prop_m_p, method = 'fdr')
  final_results$Protein_id = rownames(final_results)
  
  final_results = map_symbol(final_results, dataset_identifier)
  
  return(final_results)
}

# APOE => binary => continuous
pathology_mediate_analysis = function(data_df, cov_cols, pathology_col, apoe_col, proteins, dataset_identifier){
  mediation_results = list()
  for (protein_col in proteins) {
    cols = c(cov_cols, pathology_col, apoe_col, protein_col)
    forqc = data_df[cols]
    forqc = na.omit(forqc)
    setnames(forqc, protein_col, 'protein')
    qced = qc_protein(forqc, 'protein')
    
    model_m = glm(as.formula(paste(pathology_col, "~", apoe_col, "+", paste(cov_cols, collapse = ' + '))), data = qced, family = binomial())
    model_y = lm(as.formula(paste("protein ~", apoe_col, "+", pathology_col, "+", paste(cov_cols, collapse = ' + '))), data = qced)
    mediation_result = mediate(model.m = model_m, model.y = model_y, treat = apoe_col, mediator = pathology_col, boot = TRUE, sims = 1000)
    
    result = data.frame(
      coef_x2m = coef(model_m)[[2]],
      coef_x2y = coef(model_y)[[2]],
      coef_m2y = coef(model_y)[[3]],
      acme = mediation_result$d.avg,
      acme_p = mediation_result$d.avg.p,
      ade = mediation_result$z.avg,
      ade_p = mediation_result$z.avg.p,
      total = mediation_result$tau.coef,
      total_p = mediation_result$tau.p,
      prop_m = mediation_result$n.avg,
      prop_m_p = mediation_result$n.avg.p,
      prop_m_ci_lower = mediation_result$n.avg.ci[1],
      prop_m_ci_upper = mediation_result$n.avg.ci[2]
    )
    mediation_results[[protein_col]] = result
  }
  
  final_results = do.call(rbind, mediation_results)
  final_results$p_adjusted = p.adjust(final_results$acme_p, method = 'fdr')
  final_results$ade_p_adjusted = p.adjust(final_results$ade_p, method = 'fdr')
  final_results$prop_m_p_adjusted = p.adjust(final_results$prop_m_p, method = 'fdr')
  final_results$Protein_id = rownames(final_results)
  
  final_results = map_symbol(final_results, dataset_identifier)
  
  return(final_results)
}

# APOE => continuous => continuous
protein_mediate_Abcontinuous = function(data_df, cov_cols, pathology_col, apoe_col, proteins, dataset_identifier){
  mediation_results = list()
  for (protein_col in proteins) {
    cols = c(cov_cols, pathology_col, apoe_col, protein_col)
    forqc = data_df[cols]
    forqc = na.omit(forqc)
    setnames(forqc, protein_col, 'protein')
    qced = qc_protein(forqc, 'protein')
    
    model_m = lm(as.formula(paste("protein ~", apoe_col, "+", paste(cov_cols, collapse = ' + '))), data = qced)
    model_y = lm(as.formula(paste(pathology_col, "~", apoe_col, "+ protein +", paste(cov_cols, collapse = ' + '))), data = qced)
    mediation_result = mediate(model.m = model_m, model.y = model_y, treat = apoe_col, mediator = 'protein', boot = TRUE, sims = 1000)
    
    result = data.frame(
      coef_x2m = coef(model_m)[[2]],
      coef_x2y = coef(model_y)[[2]],
      coef_m2y = coef(model_y)[[3]],
      acme = mediation_result$d.avg,
      acme_p = mediation_result$d.avg.p,
      ade = mediation_result$z.avg,
      ade_p = mediation_result$z.avg.p,
      total = mediation_result$tau.coef,
      total_p = mediation_result$tau.p,
      prop_m = mediation_result$n.avg,
      prop_m_p = mediation_result$n.avg.p,
      prop_m_ci_lower = mediation_result$n.avg.ci[1],
      prop_m_ci_upper = mediation_result$n.avg.ci[2]
    )
    mediation_results[[protein_col]] = result
  }
  
  final_results = do.call(rbind, mediation_results)
  final_results$p_adjusted = p.adjust(final_results$acme_p, method = 'fdr')
  final_results$ade_p_adjusted = p.adjust(final_results$ade_p, method = 'fdr')
  final_results$prop_m_p_adjusted = p.adjust(final_results$prop_m_p, method = 'fdr')
  final_results$Protein_id = rownames(final_results)
  
  final_results = map_symbol(final_results, dataset_identifier)
  
  return(final_results)
}

# APOE => continuous => continuous
Abcontinuous_mediate_analysis = function(data_df, cov_cols, pathology_col, apoe_col, proteins, dataset_identifier){
  mediation_results = list()
  for (protein_col in proteins) {
    cols = c(cov_cols, pathology_col, apoe_col, protein_col)
    forqc = data_df[cols]
    forqc = na.omit(forqc)
    setnames(forqc, protein_col, 'protein')
    qced = qc_protein(forqc, 'protein')
    
    model_m = lm(as.formula(paste(pathology_col, "~", apoe_col, "+", paste(cov_cols, collapse = ' + '))), data = qced)
    model_y = lm(as.formula(paste("protein ~", apoe_col, "+", pathology_col, "+", paste(cov_cols, collapse = ' + '))), data = qced)
    mediation_result = mediate(model.m = model_m, model.y = model_y, treat = apoe_col, mediator = pathology_col, boot = TRUE, sims = 1000)
    
    result = data.frame(
      coef_x2m = coef(model_m)[[2]],
      coef_x2y = coef(model_y)[[2]],
      coef_m2y = coef(model_y)[[3]],
      acme = mediation_result$d.avg,
      acme_p = mediation_result$d.avg.p,
      ade = mediation_result$z.avg,
      ade_p = mediation_result$z.avg.p,
      total = mediation_result$tau.coef,
      total_p = mediation_result$tau.p,
      prop_m = mediation_result$n.avg,
      prop_m_p = mediation_result$n.avg.p,
      prop_m_ci_lower = mediation_result$n.avg.ci[1],
      prop_m_ci_upper = mediation_result$n.avg.ci[2]
    )
    mediation_results[[protein_col]] = result
  }
  
  final_results = do.call(rbind, mediation_results)
  final_results$p_adjusted = p.adjust(final_results$acme_p, method = 'fdr')
  final_results$ade_p_adjusted = p.adjust(final_results$ade_p, method = 'fdr')
  final_results$prop_m_p_adjusted = p.adjust(final_results$prop_m_p, method = 'fdr')
  final_results$Protein_id = rownames(final_results)
  
  final_results = map_symbol(final_results, dataset_identifier)
  
  return(final_results)
}

apoeNproteinOnAD = function(dat, protein_cols, apoe_col, ab_col, covs, dataset_identifier){
  result_list = list()
  count = 0
  for (protein_col in protein_cols) {
    print(count)
    cols = c(covs, apoe_col, ab_col, protein_col)
    forqc = dat[cols]
    forqc = na.omit(forqc)
    setnames(forqc, protein_col, 'protein')
    setnames(forqc, apoe_col, 'apoe')
    qced = qc_protein(forqc, 'protein')
    
    qced$apoeNprotein = qced$apoe * qced$protein
    
    model = glm(
      as.formula(paste(ab_col, "~ apoe + protein + apoeNprotein +", 
                       paste(covs, collapse = ' + '))),
      data = qced,
      family = binomial()
    )
    
    result_row = model_to_rows(model, 'apoeNprotein', 'logit')
    result_row$Protein_id = protein_col
    
    result_list[[protein_col]] = result_row
    count = count + 1
  }
  final_results = bind_rows(result_list)
  final_results$p_adjusted = p.adjust(final_results$p, method = 'fdr')

  final_results = map_symbol(final_results, dataset_identifier)
  
  return(final_results)
}

get_function_related = function(dat, aptamers, cov_cols){
  results = list()
  for (aptamer in aptamers){
    result_row = c()
    for (adsign in c('tnic_cho_com_I_IV', 
                     'fnc_ber_com_composite', 
                     "ct_adsign_lr",
                     "mmse_score",
                     "mPACC_v1")){
      
      ondf = dat[c(aptamer, adsign, cov_cols)]
      ondf = na.omit(ondf)
      
      setnames(ondf, aptamer, 'protein')
      setnames(ondf, adsign, 'adsign')
      
      ondf = qc_protein(ondf, 'protein')
      
      formula = as.formula(paste('adsign', '~', paste(c('protein', cov_cols), collapse = ' + ')))
      
      model = lm(formula, data = ondf)
      
      result_row = c(result_row, lm.beta(model)['protein'], 
                     summary(model)$coefficients['protein', 'Pr(>|t|)'])
    }
    results[[aptamer]] = result_row
  }
  
  result_matrix = do.call(rbind, results)
  colnames(result_matrix) = c(
    'StdBeta_taupet','p_taupet', 'StdBeta_Abpet','p_Abpet','StdBeta_ADsig','p_ADsign',
    'StdBeta_mmse','p_mmse','StdBeta_mpacc','p_mpacc'
  )
  result_df = data.frame(Protein_id = rownames(result_matrix), result_matrix, row.names = NULL)
  
  result_df$p_fdr_taupet = p.adjust(result_df$p_taupet,method="fdr")
  result_df$p_fdr_Abpet = p.adjust(result_df$p_Abpet,method="fdr")
  result_df$p_fdr_ADsign = p.adjust(result_df$p_ADsign,method="fdr")
  result_df$p_fdr_mmse = p.adjust(result_df$p_mmse,method="fdr")
  result_df$p_fdr_mpacc = p.adjust(result_df$p_mpacc,method="fdr")
  
  return(result_df)
}
#----------------------------------Plots functions------------------------------
plot_betavsp = function(df, p_col='p_adjusted', beta_col='standardized_beta', label_col='label', 
                        title_, xlabel, ylabel="-log10(FDR)", sig_color, legend_flag=FALSE){
  
  df$neg_log10_p = -log10(df[[p_col]])
  df$neg_log10_p = ifelse(df$neg_log10_p > 300, 300, df$neg_log10_p)
  
  df$color = with(df, ifelse(neg_log10_p > -log10(0.05), "Sig", "Not Sig"))
  
  p = ggplot(df, aes(x = .data[[beta_col]], y = neg_log10_p, color=color, shape = color)) +
    geom_point(aes(color = color), size=3, alpha = 0.7) +
    geom_text_repel(aes(label = ifelse(neg_log10_p > -log10(0.05), .data[[label_col]], "")), #
                    size = 6, box.padding = 0,force = 0.5, force_pull = 2,direction = "y",
                    point.padding = 0, max.overlaps=5, nudge_x = 0.1, nudge_y = -0.1) +
    scale_color_manual(values = c("Sig" = sig_color, "Not Sig" = "grey"),
                       labels = c("Sig" = "P(FDR) < 0.05", "Not Sig" = ''), drop = TRUE) +
    scale_shape_manual(values = c("Sig" = 16, "Not Sig" = 16),
                       labels = c("Sig" = "P(FDR) < 0.05", "Not Sig" = ''), drop = TRUE) +
    #ylim(0, 70) +
    #coord_cartesian(xlim = c(NA, 1.52))  +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
    labs(x = xlabel, y = ylabel, title = title_)+
    theme_classic(base_size = 20) +
    theme(plot.title = element_text(size = 20),
          legend.position = ifelse(legend_flag == TRUE, 'right', 'none')
    )

  if (max(df$neg_log10_p) > 50) {
    p = p + scale_y_continuous(
      breaks = c(0, 50, 100, 150, 200, 250, 300),
      labels = c("0", "50", "100", "150", "200", "250", ">300")
    )
  }

  return(p)
}

plot_beta_vs_beta = function(df1, 
                             df2, 
                             title_, 
                             xlabel, 
                             ylabel, 
                             merge_by = 'Protein_id', 
                             method = 'spearman',
                             sig_cor='#DC0000A2'){
  # Proteins tested in both data set
  merged = merge(df1, df2, by = merge_by)
  merged$oncol = merged[[merge_by]] 
  
  sig_df1 = subset(df1, p_adjusted < 0.05)
  sig_df2 = subset(df2, p_adjusted < 0.05)
  
  merged$color = with(merged, ifelse(p_adjusted.x<0.05 & p_adjusted.y<0.05, 'sig', "nosig"))
  
  cor_test = cor.test(merged$standardized_beta.x, 
                      merged$standardized_beta.y, 
                      method = method)
  cor_value = round(cor_test$estimate, 3)
  
  print(cor_test)
  print(cor_test$p.value)
  p_label = if (cor_test$p.value == 0) {
    "p < 2.2e-16"
  } else {
    paste0("p = ", sprintf("%.2e", cor_test$p.value))
  }
  #merged = subset(merged, color != 'grey')
  p = ggplot(merged, aes(x = standardized_beta.x, y = standardized_beta.y, color = color)) +
    geom_point(size=4, shape = 16, stroke = 0) + 
    #scale_color_npg() +
    annotate("text", 
             x = -Inf, y = Inf,
             hjust = -0.05, vjust = 1.1,
             label = paste0("R = ", cor_value, ",\n", p_label),
             size = 6, color = "black") +
    scale_color_manual(values = c('sig' = sig_cor,"nosig" = "grey")) + #
    #geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.2) +
    #geom_smooth(method = "lm", se = FALSE, color = "Black", linewidth = 0.2)+
    geom_text_repel(aes(label = ifelse(color!='nosig', label.x, NA)), 
                    force = 0.5, force_pull = 2,direction = "y",
                    size = 6, box.padding = 0, point.padding = 0, 
                    max.overlaps=6, nudge_x = 0, nudge_y = 0,) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth=1)+
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey", linewidth=1)+
    theme_classic(base_size = 20) +
    theme(
      plot.title = element_text(size = 20),
      legend.position =  'none') + 
    labs(title = title_, x = xlabel, y = ylabel)
  
  return(p)
}

apoe4ad_interactionplot = function(dat, protein_col, apoe_col, ab_col, cov_cols,title_, ylabel, legend_flag=FALSE){
  
  cols = c(cov_cols, ab_col, apoe_col, protein_col)
  forqc = dplyr::select(dat, all_of(cols))
  qced = qc_protein(forqc, protein_col)
  qced = na.omit(qced)
  
  setnames(qced, protein_col, 'Protein')
  setnames(qced, apoe_col, 'APOE4')
  setnames(qced, ab_col, 'Ab')
  
  qced$Ab = ifelse(qced$Ab == 0, 'Ab-', 'Ab+')
  qced$APOE4 = ifelse(qced$APOE4 == 0, 'E4-', 'E4+')
  qced$Ab = ordered(qced$Ab, levels = c('Ab-','Ab+'))
  
  formula = as.formula(paste('Protein ~ APOE4 * Ab +', paste(cov_cols, collapse = ' + ')))
  full_model = lm(formula, data = qced)
  
  interact_plot(full_model, 
                pred = APOE4, 
                modx = Ab,
                vary.lty = TRUE,
                interval = TRUE,
                int.width = 0.95,
                facet.modx = TRUE,
                colors = c("#3686D3", "#DE2B25"),
                pred.labels = c('E4-', 'E4+'),
                legend.main = '',
                plot.points = FALSE,
                point.alpha = 0.4,
                partial.residuals = TRUE,
                modx.labels = c('Ab-', 'Ab+'),
                jitter = FALSE,
                x.label = "",
                y.label = ylabel)+
    labs(title = title_) +
    theme_classic(base_size = 20) +
    theme(plot.title = element_text(size = 20, face = 'bold'),
          legend.position = ifelse(legend_flag == TRUE, 'right', 'none')
    )
}

# Plot protein levels by APOE status, adjust for covariates
# t test for each group and p was adjusted for multiple comparisons
proteinlevelby_auto = function(dat, protein_col, cov_cols, group_col, title_, 
                               label_,legend_label, legend_flag=FALSE){
  cols = c(cov_cols, protein_col, group_col)
  forqc = dat[cols]
  forqc = na.omit(forqc)
  forqc$group = factor(forqc[[group_col]])
  forqc$protein = forqc[[protein_col]]
  qced = qc_protein(forqc, 'protein')
  
  model = lm(as.formula(paste('protein', "~", paste(cov_cols, collapse = ' + '))), data = qced)
  qced$residual = model$residuals
  
  colors_npg = rev(pal_npg("nrc")(length(unique(qced$group))))
  colors_npg[4] = '#CED094'
  
  stat.test = qced %>% 
    t_test(residual ~ group) %>%
    add_xy_position()
  
  if(nrow(stat.test) == 1){
    stat.test$p.adj = stat.test$p
    stat.test$p.adj.signif = stat.test$p.adj.signif = 
      ifelse(stat.test$p.adj <= 0.0001, "****",
             ifelse(stat.test$p.adj <= 0.001, "***",
                    ifelse(stat.test$p.adj <= 0.01, "**",
                           ifelse(stat.test$p.adj <= 0.05, "*", "ns"))))
  }
  
  print(stat.test)
  
  group_summary = qced %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean = mean(residual),
      sd = sd(residual),
      .groups = "drop"
    )
  
  p = ggplot(qced, aes(x = group, y = residual, color = group)) +
    geom_boxplot(outlier.shape=NA, alpha = 1) +
    stat_pvalue_manual(stat.test,
                       hide.ns = 'p.adj',
                       step.increase = 0.03,
                       label = "p.adj",
                       size = 4) +
    geom_jitter(aes(color = group), 
                width = 0.2, 
                alpha = 1, 
                size = 1.8) +
    scale_color_manual(values = colors_npg, 
                       guide = "none") +
    scale_fill_manual(values = colors_npg) +
    theme_classic(base_size = 20) +
    theme(plot.title = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          legend.position = ifelse(legend_flag == TRUE, 'right', 'none')
    ) +
    labs(title = title_, y = label_, x = '', fill= legend_label)
  return(list(plot = p, stats = stat.test, group_summary = group_summary))
}

plot_e4_vs_adine3e3 = function(df1, df2, title_, xlabel, ylabel){
  # Proteins tested in both data set
  merged = merge(df1, df2, by = "Protein_id")
  sig_df1 = subset(df1, p_adjusted< 0.05)
  sig_df2 = subset(df2, p_adjusted< 0.05)
  
  merged$significant = merged$Protein_id %in% sig_df1$Protein_id & merged$Protein_id %in% sig_df2$Protein_id
  
  merged = subset(merged, significant)
  
  merged$color = with(merged, ifelse((standardized_beta.x *standardized_beta.y) > 0, 
                                     "#DC0000A2", "#2166ACFF"))
  
  p = ggplot(merged, aes(x = standardized_beta.x, y = standardized_beta.y, 
                         color = color)) +
    geom_point(size=3, shape = 16, stroke = 0) + 
    scale_color_manual(values = c("#DC0000A2" = "#DC0000A2","#2166ACFF" = "#2166ACFF")) +
    geom_text_repel(aes(label = ifelse(color != 'grey', symbol.x, NA)), 
                    force = 0.5, force_pull = 2,direction = "y",
                    size = 6, box.padding = 0, point.padding = 0, 
                    max.overlaps=6, nudge_x = 0.05, nudge_y = 0.01) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
    theme_classic(base_size = 22) +
    theme(
      plot.title = element_text(size = 22),
      legend.position =  'none') + 
    labs(title = title_, x = xlabel, y = ylabel)
  
  return(p)
}

levelbyage_inneg  = function(df, cols, labels = NULL, 
                             group, ylabel, title_) {
  if (is.null(labels)) labels  = cols
  if (length(cols) != length(labels)) {
    stop("Unequal length of cols and labels")
  }
  
  df_long  = df %>%
    pivot_longer(
      cols = all_of(cols),
      names_to = "protein_name",
      values_to = "expression"
    ) %>%
    mutate(
      protein_name = factor(protein_name, levels = cols, labels = labels)
    )
  
  n_groups  = df_long %>% pull(.data[[group]]) %>% unique() %>% length()
  
  colors_npg  = rev(pal_npg("nrc")(n_groups))
  
  if (n_groups > 4) {
    colors_npg[4]  = "#CED094"
  }
  
  p  = ggplot(df_long,aes(x = age, y = expression, color = .data[[group]])) +
    geom_point(size=2, alpha=0.5) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 2) +
    labs(
      x = "Age", y = ylabel, color = "",
      title = title_) +
    scale_color_manual(values = colors_npg) +
    theme_classic(base_size = 20) +
    theme(plot.title = element_text(size = 20), legend.position = "right")
  
  # n_proteins  = nlevels(df_long$protein_name)
  # if (n_proteins > 1) {
  #   p  = p + facet_wrap(~ protein_name, scales = "free_y")
  # }
  
  return(p)
}

plot_pair_correlation = function(df, x, y,
                                 title_ = "",
                                 xlabel = x,
                                 ylabel = y,
                                 method = "spearman",
                                 legend_flag = FALSE) {
  
  cor_test = cor.test(df[[x]], df[[y]], method = method)
  cor_value = round(cor_test$estimate, 3)
  p_value = format_pval(cor_test$p.value)
  
  x_pos = min(df[[x]], na.rm = TRUE)
  y_pos = max(df[[y]], na.rm = TRUE)*0.99
  
  color = pal_npg("nrc")(10)[sample(1:10, 1)]
  
  p = ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(size = 3, alpha = 0.5, color = color) +
    geom_smooth(method = "lm", se = TRUE, linetype = "dashed",
                color = color) +
    annotate("text", 
             x = -Inf, y = Inf,
             hjust = -0.05, vjust = 1.1,
             label = paste0("R=", cor_value, "\nP=",p_value),
             hjust = 0, size = 6) +
    labs(title = title_, x = xlabel, y = ylabel) +
    theme_classic(base_size = 20) +
    theme(plot.title = element_text(size = 20, hjust=0.5),
          legend.position = ifelse(legend_flag == TRUE, 'right', 'none'))
  
  return(p)
}

plot_correlation_comprehensive_matrix = function(data_df, cols, labels = NULL) {
  library(GGally)
  
  df = data_df[, cols]
  
  if (!is.null(labels)) {
    colnames(df) = labels
  }
  
  p = GGally::ggpairs(
    df,
    upper = list(continuous = GGally::wrap("cor", method = "spearman")),
    lower = list(continuous = GGally::wrap("points", alpha = 0.6)),
    diag  = list(continuous = "densityDiag")
  )
  
  return(p)
}

display_venn = function(sets_list, labels=NULL, ...) {
  colors_npg = c("#56B4E9", "#E69F00", "#CD534CFF", "#999999","#009E73",
                 pal_npg("nrc")(10))[1:length(sets_list)]

  if (is.null(labels)) {
    labels = names(sets_list)
    if (is.null(labels)) {
      labels = paste0("Set", seq_along(sets_list))
    }
  }

  library(ggvenn)
  p = ggvenn(
    sets_list,
    fill_color = colors_npg,
    stroke_size = 0.5,
    fill_alpha = 0.5,
    set_name_size = -1,
    show_percentage = FALSE
  )

  return(p)
}
