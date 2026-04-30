source('APOE_proteomics/analysis/src.R')
#-------------------Plots: Volcano, beta vs beta---------------------------
##-------------------------------------Load data--------------------------------
apoe_identifier = 'APOE2'
dataset_identifier = 'BF2OLINK'
result_folder = glue('APOE_proteomics/results/{dataset_identifier}/{apoe_list[[apoe_identifier]]}')
phenotype = "AD" #Aβ

files = c(
  apoe2protein          = "apoe2protein.csv",
  apoe2protein_adjab    = "apoe2protein_adjab.csv",
  ab2protein            = "ab2protein.csv",
  ab2protein_adjapoe    = "ab2protein_adjapoe.csv",
  
  neg_allage            = "neg_allage.csv",
  pos_allage            = "pos_allage.csv",
  
  neg_younger_result    = "neg_younger.csv",
  neg_older_result      = "neg_older.csv",
  
  neg_interaction       = "neg_allage_interaction.csv",
  apoe_ab_interaction   = "apoe_ab_interaction.csv",
  
  # Only for GNPC
  ad_ine3e3            = "ad_ine3e3.csv",
  
  ab2protein            = "ad2protein.csv",
  ab2protein_adjapoe    = "ad2protein_adjapoe.csv",
  
  neg_allage            = "cu_allage.csv",
  pos_allage            = "ad_allage.csv",
  
  neg_younger_result    = "cu_younger.csv",
  neg_older_result      = "cu_older.csv",
  
  neg_interaction       = "cu_allage_interaction.csv",
  apoe_ab_interaction   = "apoe_ad_interaction.csv"
)

associations =
  files |>
  imap(function(fname, key) {
    f = file.path(result_folder, fname)
    if (!file.exists(f)) return(NULL)
    read.csv(f)
  }) |>
  compact()

##------------------------------Volcano-----------------------------------------
plot_lists = list(
  apoe2protein = list(xlab = glue('Standardized β ({apoe_identifier})'),
                      title = glue('{apoe_identifier}-associated proteins'),
                      color = '#DC0000A2'),
  ab2protein = list(xlab = glue('Standardized β (clinical AD)'),
                      title = glue('Without {apoe_identifier} adjustment'),
                      color = '#2166ACFF')
)

imap(plot_lists, function(feature, key){
  ggsave(
    glue("Figs/{dataset_identifier}/volcano_{apoe_identifier}_{key}.svg"),
    plot_betavsp(associations[[key]],
                 'p_adjusted', 
                 'standardized_beta',
                 'label',
                 feature$title, 
                 feature$xlab, 
                 "-log10(FDR)",
                 feature$color),
    width = 6,
    height = 6,
    dpi = 300,
    device = "svg"
  )
})
##---------------------------Beta vs. Beta--------------------------------------
###-------------------------Fig. 3f, APOE4 in CU vs. AD in ε3ε3-----------------
ggsave(
  glue('Figs/{dataset_identifier}/APOEinCU_vs_ADine3e3.svg'),
  plot_e4_vs_adine3e3(associations$neg_allage, 
                      associations$ad_ine3e3, 
                      glue('{apoe_identifier} in CU vs. AD in ε3ε3'),
                      glue("Standardized β ({apoe_identifier} in CU)"),
                      "Standardized β (AD in ε3ε3)"),
  width = 6,
  height = 6,
  dpi = 300
)
###-------------------Fig. 2d, 3d, Whole vc. CU---------------------------------
ggsave(
  glue('Figs/{dataset_identifier}/whole_vs_cu_{apoe_identifier}.svg'),
  plot_beta_vs_beta(associations$apoe2protein, 
                    associations$neg_allage, 
                    glue('{apoe_identifier}-associated proteins'), 
                    glue('Standardized β in whole cohort'), 
                    glue('Standardized β in CU ')),
  width = 6,
  height = 6,
  dpi = 300
)

dataset_identifier = 'ADNI_MS'
phenotype = 'Aβ'
###--Beta vs. Beta for APOE or Abeta with or withour adjustment for each other--
walk(c('APOE4', 'APOE2'), function(apoe_identifier){
  result_folder = glue('APOE_proteomics/results/{dataset_identifier}/{apoe_list[[apoe_identifier]]}')
  
  walk(c('APOE', phenotype), function(apoe_or_ab){
    if (apoe_or_ab == 'APOE'){
      apoe2protein = read_csv(glue("{result_folder}/apoe2protein.csv"))
      apoe2protein_adjab = read_csv(glue("{result_folder}/apoe2protein_adjab.csv"))
      
      p = plot_beta_vs_beta(apoe2protein, 
                            apoe2protein_adjab, 
                            glue('{apoe_identifier}-associated proteins'), 
                            glue('Standardized β (without {phenotype} adjustment)'), 
                            glue('Standardized β (with {phenotype} adjustment)'))
    } else {
      if (dataset_identifier == 'GNPC'){
        file = 'ad2protein'
      } else {
        file = 'ab2protein'
      }
      
      ab2protein = read_csv(glue("{result_folder}/{file}.csv"))
      ab2protein_adjapoe = read_csv(glue("{result_folder}/{file}_adjapoe.csv"))
      
      p = plot_beta_vs_beta(ab2protein, 
                            ab2protein_adjapoe, 
                            glue('{phenotype}-associated proteins'), 
                            glue('Standardized β (without {apoe_identifier} adjustment)'), 
                            glue('Standardized β (with {apoe_identifier} adjustment)'),
                            sig_cor = '#0072B2FF') 
    }
    
    ggsave(
      glue("Figs/{dataset_identifier}/betaVSbeta_for{apoe_or_ab}_{apoe_identifier}_{phenotype}.svg"),
      p,
      width = 6.5,
      height = 6.5,
      dpi = 300,
      device = "svg"
    )
  })
})
#----------------Fig. 2i, 3i: LDA-----------------------------------------------
# library(ggplot2)
# library(MASS)
# library(ggforce)
walk(apoe_list, function(apoe_folder){
  df = read_excel(glue('APOE_proteomics/results/GNPC/{apoe_folder}/lda.xlsx'))
  
  if(apoe_folder == 'e2vse3e3'){
    xlim_ = c(-13.5, 4)
    ylim_ = c(-8, 14)
  } else {
    xlim_ = c(-16, 4)
    ylim_ = c(-17, 10)
  }
  
  colors_npg = pal_npg("nrc")(length(unique(df$Group)))
  colors_npg[6] = 'lightgrey'
  
  background_group = "Ungrouped"
  
  df_sorted = df %>%
    mutate(ordering = ifelse(Group == background_group, 0, 1)) %>%
    arrange(ordering)
  
  print(table(df_sorted$Group))
  
  p = ggplot(df_sorted, aes(x = LD1, y = LD2, color = Group)) +
    geom_point(size = 4, alpha = 0.7) +
    stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.5) +
    geom_text(data = df_sorted %>% filter(Group != 'Ungrouped'),
              aes(label = Protein),
              vjust = -0.8, size = 4, check_overlap = TRUE) +
    labs(
      title = "LDA projection score of each group of proteins",
      x = "LD1",
      y = "LD2",
      color = ""
    ) +
    xlim(xlim_) +
    ylim(ylim_) +
    #scale_color_npg() +
    scale_color_manual(values = colors_npg) +
    scale_fill_manual(values = colors_npg) +
    theme_classic(base_size = 20) +
    theme(plot.title = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.position = 'top')
  print(p)
  ggsave(glue('Figs/GNPC/lda_{apoe_folder}.svg'),
         p,
         width = 8, 
         height = 6,
         dpi = 300
  )
})
#----------------Fig. 5e: S100A13*APOE4 effect on Abeta positivity--------------
dataset_identifier = 'BF2SomaLogic'
initial = read.csv('data/BF2Somalogic/preprocessed/with24_7285prot.csv')
nonAD = subset(initial, diagnosis_baseline_variable != 'AD' & apoe_genotype_baseline_variable %in% c(33, 34, 43, 44))
nonAD$apoe_genotype_baseline_variable = ifelse(nonAD$apoe_genotype_baseline_variable %in% c(22, 23, 32), '22/23', 
                                               nonAD$apoe_genotype_baseline_variable)
get_demographic(nonAD, dataset_identifier, 'APOE4')

apoeNs100a13_onAb = read_csv("data/BF2Somalogic/Rresults/V5_e4vse3e3/Whole/protein_regulation.csv")
apoeNs100a13_onAb = apoeNs100a13_onAb[apoeNs100a13_onAb$symbol == 'S100A13', ]

s100a13_inter = glm(Abnormal_CSF_Ab42_Ab40_Ratio ~ S100A13__seq.7223.60__ANML * apoe_class_num + age + gender_baseline_variable + plasma_ANML_mean, 
                     data = nonAD, family = binomial)
sim_slopes(s100a13_inter, pred = apoe_class_num, modx = S100A13__seq.7223.60__ANML, jnplot = TRUE)

fig5e = interact_plot(s100a13_inter, pred = apoe_class_num, modx = S100A13__seq.7223.60__ANML, 
                      plot.points = FALSE, 
                      legend.main = "S100A13",
                      main.title = 'APOE4 * protein effect on Aβ',
                      interval = TRUE,
                      x.label = "APOE4 carriers", 
                      colors = rev(pal_npg("nrc")(9)),
                      line.thickness = 1.5,
                      y.label = "Aβ Status") + 
  annotate("text", 
           x = -Inf, y = Inf,
           hjust = -0.05, vjust = 1.1,
           label = paste0("p(APOE4*S100A13) = ", format_pval(apoeNs100a13_onAb$p[1])),
           size = 6, color = "black") +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(size = 20),
        legend.position = 'top')
fig5e
ggsave(
  glue("Figs/{dataset_identifier}/fig5e.svg"),
  fig5e,
  width = 5.5, 
  height = 5.5,
  dpi = 300
)

# Source data
dat_plot = nonAD %>%
  dplyr::select(
    Abnormal_CSF_Ab42_Ab40_Ratio,
    apoe_class_num,
    S100A13__seq.7223.60__ANML,
    age,
    gender_baseline_variable,
    plasma_ANML_mean
  ) %>%
  dplyr::filter(stats::complete.cases(.))

modx_var = "S100A13__seq.7223.60__ANML"

m = mean(dat_plot[[modx_var]], na.rm = TRUE)
s = sd(dat_plot[[modx_var]], na.rm = TRUE)

gender_template = dat_plot$gender_baseline_variable

if (is.factor(gender_template)) {
  gender_val = factor(
    names(which.max(table(gender_template))),
    levels = levels(gender_template)
  )
} else if (is.numeric(gender_template) || is.integer(gender_template)) {
  gender_val = as.numeric(names(which.max(table(gender_template))))
} else {
  gender_val = names(which.max(table(gender_template)))
}

source_fig5e = expand.grid(
  apoe_class_num = c(0, 1),
  S100A13__seq.7223.60__ANML = c(m - s, m, m + s)
) %>%
  dplyr::mutate(
    S100A13_level = dplyr::case_when(
      abs(S100A13__seq.7223.60__ANML - (m - s)) < 1e-12 ~ "- 1 SD",
      abs(S100A13__seq.7223.60__ANML - m) < 1e-12       ~ "Mean",
      abs(S100A13__seq.7223.60__ANML - (m + s)) < 1e-12 ~ "+ 1 SD"
    ),
    age = mean(dat_plot$age, na.rm = TRUE),
    plasma_ANML_mean = mean(dat_plot$plasma_ANML_mean, na.rm = TRUE)
  )

source_fig5e$gender_baseline_variable = gender_val

source_fig5e = source_fig5e %>%
  dplyr::select(
    apoe_class_num,
    S100A13_level,
    S100A13__seq.7223.60__ANML,
    age,
    gender_baseline_variable,
    plasma_ANML_mean
  )

source_fig5e$predicted_probability = predict(
  s100a13_inter,
  newdata = source_fig5e,
  type = "response"
)

source_fig5e = source_fig5e %>%
  dplyr::mutate(
    APOE4_status = dplyr::if_else(apoe_class_num == 0, "Non-carriers", "Carriers")
  ) %>%
  dplyr::select(
    apoe_class_num,
    APOE4_status,
    S100A13_level,
    predicted_probability
  )

write.xlsx(source_fig5e,
  glue("SourceData/Fig5e_source_data.xlsx")
)

#----------Fig. 5b,d: Protein levels groups by APOE genotypes--------------------
dataset_identifier = 'BF2SomaLogic'
initial = read.csv('data/BF2Somalogic/preprocessed/with24_7285prot.csv')
get_demographic(subset(initial, Abnormal_CSF_Ab42_Ab40_Ratio == 0), 
                'BF2SomaLogic', 'APOE4', FALSE)

initial$apoe_genotype_baseline_variable = ifelse(initial$apoe_genotype_baseline_variable %in% c(22, 23, 32), '22/23', 
                                                 initial$apoe_genotype_baseline_variable)

apoe4_mediators = list(
  S100A13 = 'S100A13__seq.7223.60__ANML',
  SPC25 = 'SPC25__seq.22782.80__ANML',
  TBCA = 'TBCA__seq.12501.10__ANML'
)

apoe2_mediators = list(
  APOB = 'Apo_B__seq.2797.56__ANML',
  SNAP23 = 'SNP23__seq.20241.9__ANML',
  PCLAF = 'PAF__seq.19158.1__ANML'
)

plotlists = imap(list(APOE4=apoe4_mediators, APOE2=apoe2_mediators), function(lists, apoe){
  res_lists = imap(lists, function(col, name){
    res = proteinlevelby_auto(dat=initial[initial$Abnormal_CSF_Ab42_Ab40_Ratio==0,], 
                              protein_col=col, 
                              cov_cols=c('age', 'gender_baseline_variable', 'plasma_ANML_mean'), 
                              group_col='apoe_genotype_baseline_variable', 
                              title_=NULL, 
                              label_=name,
                              legend_label='', 
                              legend_flag=FALSE)
    
    res$stats$.y. = name
    res$group_summary$.y. = name
    res
  })
}) |> set_names(c('APOE4', 'APOE2'))

iwalk(plotlists, function(lists, name){
  p <- wrap_plots(map(lists, function(p) p$plot), nrow = 1) +
    plot_annotation(
      title = glue('{name}=>Aβ mediators in Aβ-'),
      theme = theme(
        plot.title = element_text(
          size = 19,
          hjust = 0.1,
          margin = margin(t = 5, r = 0, b = 5, l = 0)
        ),
        plot.margin = margin(t = 5, r = 1, b = 5, l = 1)
      )
    )
  
  ggsave(
    glue("Figs/{dataset_identifier}/{name}mediators_groupby_genetypes.svg"),
    plot = p,
    width = 12,
    height = 6,
    dpi = 300,
    device = "svg"
  )
})

flatted = flatten(plotlists)
stats = bind_rows(map(flatted, function(res) as.data.frame(res$stats)))
write.csv(stats[c(".y.","group1","group2","n1","n2","statistic","df","p","p.adj","p.adj.signif")], 
          file = 'Supplement/BF2_protein_groupbyAPOE.csv', row.names = FALSE)

levels = bind_rows(map(flatted, function(res) as.data.frame(res$group_summary)))
write.csv(levels, file = 'Supplement/BF2_protein_groupbyAPOE.csv', row.names = FALSE)

#------Extended Fig. 2: Sex-independent APOE effect at Aβ- stages---------------
iwalk(apoe_list, function(folder, apoe_identifier){
  apoe_sex_interaction = read_csv(glue("APOE_proteomics/results/BF2SomaLogic/{folder}/apoe_sex_interaction.csv"))
  ggsave(
    filename = glue("Figs/Other/BF2SomaLogic_{apoe_identifier}_sexInter.svg"),
    plot     = plot_betavsp(apoe_sex_interaction,
                            'p_adjusted', 
                            'standardized_beta',
                            'label',
                            glue('{apoe_identifier}*sex effects on early-altered proteins'), 
                            glue('Standardized β ({apoe_identifier}*Sex)'), 
                            "-log10(FDR)",
                            '#DC0000A2'),
    width    = 6,
    height   = 6,
    units    = "in",
    device   = "svg"
  )
  })

iwalk(apoe_list, function(folder, apoe_identifier){
  fuli_folder = glue("APOE_proteomics/results/BF2SomaLogic/{folder}")
  male = read_csv(glue("{fuli_folder}/apoe2protein_neg_male.csv"))
  female = read_csv(glue("{fuli_folder}/apoe2protein_neg_female.csv"))
  ggsave(
    filename = glue("Figs/BF2SomaLogic/{dataset_identifier}_{apoe_identifier}_betaMaleVsBetaFemale.svg"),
    plot     = plot_beta_vs_beta(male, 
                                 female,
                                 glue('Sex-stratified {apoe_identifier} effects in Aβ-'),
                                 'Standardized β in Aβ- men',
                                 "Standardized β in Aβ- women",
                                 'Protein_id'),
    width    = 6,
    height   = 6,
    units    = "in",
    device   = "svg"
  )
})


##-------------Export results----------------------
select_cols = function(df){
  df %>%
    select(apt_name, UniProt, EntrezGeneSymbol, label, coef, 
           standardized_beta, se, resid, p, p_adjusted)
}

wb = createWorkbook()
sex_results = map(names(apoe_list), function(apoe_identifier){
  result_folder = glue('data/BF2Somalogic/Rresults/V5_{apoe_list[[apoe_identifier]]}/Whole')
  
  map(names(sex_list), function(sex){
    tmp_df = read_csv(glue("{result_folder}/apoe2protein_neg_{sex}.csv"))
    
    sheet_inter = glue("{apoe_identifier}_In{toupper(sex)}")
    addWorksheet(wb, sheet_inter)
    writeData(wb, sheet_inter, select_cols(tmp_df))
  })
  
  inter_df = read_csv(glue("{result_folder}/apoe_sex_interaction.csv"))
  
  sheet_inter = glue("{apoe_identifier}xSex")
  addWorksheet(wb, sheet_inter)
  writeData(wb, sheet_inter, select_cols(inter_df))
})

saveWorkbook(
  wb,
  file = "SourceData/APOESexEffect.xlsx",
  overwrite = TRUE
)

#--------Extended Fig. 3: CSF APOE-proteomics signatures in TMT-MS (ADNI)-------
# See Plots: Volcano, beta vs beta section

#--Extended Fig. 4: Direct comparison of SomaLogic vs. TMT-MS in ADNI (CSF)-----
##----------------Extended Fig. 4a-----------
adni = read.csv('data/ADNI/preprocessed/adni_soma_ms_merged.csv')

anml_cols = colnames(adni)[grep('^(S100A13|TBCA|NEFL|LRRN1).*ANML$', colnames(adni))]
tmt_cols = colnames(adni)[grep('^(S100A13|TBCA|NEFL|LRRN1).*MS$', colnames(adni))]

residualize_into = function(df, y, x) {
  idx = complete.cases(df[[y]], df[[x]])
  out = rep(NA_real_, nrow(df))
  if (sum(idx) >= 2) {
    fit = lm(df[[y]][idx] ~ df[[x]][idx])
    out[idx] = resid(fit)
  }
  out
}

walk(anml_cols, function(col){
  adni[[col]] <<- residualize_into(adni, col, "CSF_ANML_mean")
})

walk(tmt_cols, function(col){
  adni[[col]] <<- residualize_into(adni, col, "CSF_TMT_mean")
})

resided = adni[c(anml_cols, tmt_cols)]

plots = map(c('S100A13', 'TBCA'), function(gene){
  anml = colnames(resided)[grep(glue('^{gene}.*ANML$'), colnames(resided))]
  tmt = colnames(resided)[grep(glue('^{gene}.*MS$'), colnames(resided))]
  
  plot_pair_correlation(resided, anml[[1]], tmt[[1]],
                        title_ = gene,
                        xlabel = 'SomaLogic',
                        ylabel = 'TMT-MS',
                        method = "spearman",
                        legend_flag = FALSE)
})

ggsave(
  glue("Figs/ADNI_MS/somavstmt.svg"),
  wrap_plots(plots) + plot_layout(nrow = 1, ncol = 2),
  width = 8,
  height = 5,
  dpi = 300,
  device = "svg"
)
##-------------Extended Fig. 4 b-d------------------
# See r2correlation.ipynb
corr_df = read_excel('~/Downloads/43587_2026_1123_MOESM13_ESM.xlsx')
corr_df = read_csv("Supplement/ADNI_CSFSoma_TMT_correlation_spearman.csv")
corr_df$p_adjusted = p.adjust(corr_df$p, method = 'fdr')
write.xlsx(corr_df, 'SourceData/ADNI_CSFSoma_TMT_correlation_spearman.xlsx')
corr_df$label = corr_df$symbol

threshold = 0.3

r = corr_df$r
data_in_range = r[abs(r)<= threshold]
count_in_range = length(data_in_range)
total_count = length(r)
proportion = count_in_range / total_count
proportion
# Create a simple density plot
p = ggplot(data = corr_df, aes(x = r)) +
  geom_histogram(fill = "#DC0000A2", alpha = 0.3, bins=20, color='black') +
  labs(x = "Spearman r", y = "Density", 
       title = "SomaLogic vs. TMT-MS (ADNI CSF)") +
  theme_classic(base_size = 18) +
  theme(plot.title = element_text(size = 18))+ 
  scale_x_continuous(breaks = seq(-0.3, 1, by = 0.2)) +
  #xlim(-0.31, 1)+
  geom_vline(xintercept = threshold, linetype = "dashed", color = "darkred", 
             linewidth = 1) +
  geom_vline(xintercept = -threshold, linetype = "dashed", color = "darkred", 
             linewidth = 1) + 
  annotate("text", size=6, x = 1, y = 300, 
           label = glue("[-{threshold}~{threshold}] {round(proportion*100, 1)}%"), 
           color = "black", hjust = 1.5)

p
ggsave("Figs/ADNI_MS/somsvstmt_corr_density.svg", 
       p,
       width = 6, 
       height = 6,
       dpi = 300)

ggsave("Figs/ADNI_MS/somavstmt_corr_volcano.svg", 
       plot_betavsp(corr_df,
                    'p_adjusted', 
                    'standardized_beta',
                    'label',
                    'SomaLogic vs. TMT-MS (ADNI CSF)',
                    'Spearman r',
                    "-log10(FDR)",
                    '#DC0000A2'),
       width = 6, 
       height = 6,
       dpi = 300)
#---Extended Data Fig. 4: APOE-proteomics changes in BioFINDER-2 CSF OLINK cohort--
# For Extended Data Fig. 4a,c, See Plots: Volcano, beta vs beta section

##-------Extended Data Fig. 5b,d----------------
initial = read.csv('data/BF2OLINK/preprocessed/V5_with24_1391proteins.csv')

plot_lists = list(
  APOE4 = 'SNAP25__Neurology_II_o2csf__NPX',
  APOE2 = 'LDLR__Cardiometabolic_o2csf__NPX'
)

plots = imap(plot_lists, function(col, apoe){
  if (apoe == 'APOE4'){
    cu_df = initial %>%
      filter(
        Abnormal_CSF_Ab42_Ab40_Ratio == 0 &
          apoe_genotype_baseline_variable %in% c(33, 34, 43, 44)
      )
    
    cu_df$apoe = ifelse(cu_df$apoe4 == 'E4-', 'ε3/ε3', 'ε4+')
  } else {
    cu_df = initial %>%
      filter(
        Abnormal_CSF_Ab42_Ab40_Ratio == 0 &
          apoe_genotype_baseline_variable %in% c(22, 23, 32, 33)
      )
    
    cu_df$apoe = ifelse(cu_df$apoe2 == 'E2-', 'ε3/ε3', 'ε2+')
  }
  
  levelbyage_inneg(cu_df, col, NULL, 'apoe', sub("__.*$", "", col),
                   NULL)
})

plots$APOE4 +
  plot_annotation(
    title = glue('Age-regulated APOE4-effect in Aβ-'),
    theme = theme(plot.title = element_text(size = 20))
  )

ggsave(
  glue("Figs/BF2OLINK/APOE4mediators_byage.svg"),
  plots$APOE4 +
    plot_annotation(
      title = glue('Age-regulated APOE4-effect in Aβ-'),
      theme = theme(plot.title = element_text(size = 20, hjust=0.3))
    ),
  width = 6,
  height = 6,
  dpi = 300,
  device = "svg"
)

ggsave(
  glue("Figs/BF2OLINK/APOE2mediators_byage.svg"),
  plots$APOE2 +
    plot_annotation(
      title = glue('Stable APOE2-effect across ages in Aβ-'),
      theme = theme(plot.title = element_text(size = 20, hjust=0.3))
    ),
  width = 6,
  height = 6,
  dpi = 300,
  device = "svg"
)

#----Extended Fig. 6: Direct comparison of plasma SomaLogic vs. CSF OLINK in BioFINDER-2-----
##---------------Extended Fig. 6 a,b-------------------
bf2soma = read.csv('data/BF2Somalogic/preprocessed/with24_7285prot.csv')
bf2olink = read.csv('data/BF2/preprocessed/V5_with24_1391proteins.csv')
bf2 = merge(bf2soma, bf2olink, by='sid')

soma_cols = colnames(bf2)[grepl('^(S100A13|TBCA|NFL|LRRN1|SIA8A).*ANML$', colnames(bf2))]
olink_cols = colnames(bf2)[grepl('^(S100A13|TBCA|NEFL|LRRN1|ST8SIA1).*NPX$', colnames(bf2))]

residualize_into = function(df, y, x) {
  idx = complete.cases(df[[y]], df[[x]])
  out = rep(NA_real_, nrow(df))
  if (sum(idx) >= 2) {
    fit = lm(df[[y]][idx] ~ df[[x]][idx])
    out[idx] = resid(fit)
  }
  out
}

walk(soma_cols, function(col){
  bf2[[col]] <<- residualize_into(bf2, col, "plasma_ANML_mean")
})

walk(olink_cols, function(col){
  bf2[[col]] <<- residualize_into(bf2, col, "CSF_NPX_mean")
})

resided = bf2[c(soma_cols, olink_cols)]

resided = resided %>% select(
  -SIA8A__seq.21663.149__ANML,
  -LRRN1_ECD__seq.11586.2__ANML
) %>%
  mutate(
    ST8SIA1__seq.21508.7__ANML= SIA8A__seq.21508.7__ANML,
    NEFL__seq.10082.251__ANML = NFL__seq.10082.251__ANML
    )

plots = map(c('S100A13', 'TBCA', 'NEFL'), function(gene){
  soma = colnames(resided)[grep(glue('^{gene}.*ANML$'), colnames(resided))]
  olink = colnames(resided)[grep(glue('^{gene}.*NPX$'), colnames(resided))]
  
  plot_pair_correlation(resided, soma[[1]], olink[[1]],
                        title_ = gene,
                        xlabel = 'SomaLogic (plasma)',
                        ylabel = 'OLINK (CSF)',
                        method = "spearman",
                        legend_flag = FALSE)
})

wrap_plots(plots) + plot_layout(nrow = 1, ncol = 3)
ggsave(
  glue("Figs/BF2SomaLogic/somavsolink_for_representatives.svg"),
  wrap_plots(plots) + plot_layout(nrow = 1, ncol = 3),
  width = 12,
  height = 6,
  dpi = 300,
  device = "svg"
)

ggsave(
  glue("Figs/BF2SomaLogic/ST8SIA1.svg"),
  plot_correlation_comprehensive_matrix(bf2[c(soma_cols, olink_cols)], 
                                        c('SIA8A__seq.21508.7__ANML',
                                         'SIA8A__seq.21663.149__ANML',
                                         'ST8SIA1__Inflammation_II_o2csf__NPX'), 
                                        labels = c('seq.21508.7', 
                                                   'seq.21663.149', 
                                                   'OLINK Inflammation')) 
  + ggplot2::ggtitle('ST8SIA1'),
  width = 4,
  height = 4,
  dpi = 300,
  device = "svg"
)

ggsave(
  glue("Figs/BF2SomaLogic/LRRN1.svg"),
  plot_correlation_comprehensive_matrix(bf2[c(soma_cols, olink_cols)], 
                                        c('LRRN1_CD__seq.11293.14__ANML',
                                          'LRRN1_ECD__seq.11586.2__ANML',
                                          'LRRN1__Inflammation_o2csf__NPX'), 
                                        labels = c('seq.11293.14', 
                                                   'seq.11586.2', 
                                                   'OLINK Inflammation')) 
  + ggplot2::ggtitle('LRRN1'),
  width = 4,
  height = 4,
  dpi = 300,
  device = "svg"
)

##---------------Extended Fig. 6 d,e-------------------
# See r2correlation.ipynb
soma_olink_cor = read.csv("Supplement/BF2_plasmaSoma_csfOLINK_correlation_spearman.csv")
soma_olink_cor$p_adjusted = p.adjust(soma_olink_cor$p, method = 'fdr')
write.xlsx(soma_olink_cor, "SourceData/BF2_plasmaSoma_csfOLINK_correlation_spearman.xlsx")
soma_olink_cor$label = soma_olink_cor$symbol

threshold = 0.3
r = soma_olink_cor$r
data_in_range = r[abs(r) <= 0.3]
count_in_range = length(data_in_range)
total_count = length(r)
proportion = count_in_range / total_count
proportion
# Create a simple density plot
p = ggplot(data = soma_olink_cor, aes(x = r)) +
  geom_histogram(fill = "#DC0000A2", alpha = 0.3, bins=20, color='black') +
  labs(x = "Spearman r", y = "Density", 
       title = "Plasma SomaLogic vs. CSF OLINK (BF2)") +
  theme_classic(base_size = 18) +
  theme(plot.title = element_text(size = 18))+ 
  scale_x_continuous(breaks = seq(-0.3, 1, by = 0.2)) +
  #xlim(-0.31, 1)+
  geom_vline(xintercept = threshold, linetype = "dashed", color = "darkred", 
             linewidth = 1) +
  geom_vline(xintercept = -threshold, linetype = "dashed", color = "darkred", 
             linewidth = 1) + 
  annotate("text", size=6, x = 1, y = 300, 
           label = glue("[-{threshold}~{threshold}] {round(proportion*100, 1)}%"), 
           color = "black", hjust = 1.5)
p

ggsave("Figs/other/somavsolink_corr_density.svg", 
       p,
       width = 6, 
       height = 6,
       dpi = 300)

ggsave("Figs/BF2SomaLogic/somavsolink_corr_volcano.svg", 
       plot_betavsp(soma_olink_cor,
                    'p_adjusted', 
                    'standardized_beta',
                    'label',
                    'Plasma SomaLogic vs. CSF OLINK (BF2)',
                    'Spearman r',
                    "-log10(FDR)",
                    '#DC0000A2'),
       width = 6, 
       height = 6,
       dpi = 300)

#-------Extended Data Fig. 7: APOE-proteomics changes in UKBB plasma OLINK cohort-----------------------
ukbb = read.csv('data/UKBB/preprocessed/ukbb_clean.csv')
withdraw = read_table("data/UKBB/raw_data/withdraw105777_432_20251223.txt", col_names = FALSE)
ukbb = ukbb[!ukbb$eid %in% withdraw$X1, ]

##-------------Extended Fig. 7a,b----------------
# See Plots: Volcano, beta vs beta section

##-------------Extended Fig. 7c------------------
cov_cols = c('age', 'gender', 'Plasma_NPX_mean')
cols = c('CSF_PLA2G7_NPX', 'CSF_BRK1_NPX', 'CSF_LDLR_NPX')
plots = map(cols, function(col){
  prot = sub("^CSF_(.*)_NPX$", "\\1", col)
  res = proteinlevelby_auto(ukbb, col, c(cov_cols), 'APOE', '', prot, TRUE)
  res$stats$.y. = prot
  res$group_summary$.y. = prot
  res
}) |> set_names(cols)

plots$CSF_PLA2G7_NPX$plot = plots$CSF_PLA2G7_NPX$plot + 
  labs(title = 'Protein level change group by APOE genotypes')

ggsave("Figs/UKBB/protein_level_byAPOE.svg", 
       wrap_plots(compact(map(plots, function(plot) plot$plot))),
       width = 12, 
       height = 6,
       dpi = 300)

group_summary = bind_rows(map(plots, function(plot) plot$group_summary))
stats = bind_rows(map(plots, function(plot) plot$stats))
write.xlsx(group_summary, 'SourceData/ukbb_group_sum.xlsx')
write.xlsx(stats, 'SourceData/ukbb_stats.xlsx')

##-------Extended Fig.7d---------------
ukbb = ukbb %>% mutate (
  apoe_plot = case_when(
    apoe == 'e2e2' ~ 'ε2/ε2',
    apoe == 'e2e3' ~ 'ε2/ε3',
    apoe == 'e2e4' ~ 'ε2/ε4',
    apoe == 'e3e3' ~ 'ε3/ε3',
    apoe == 'e3e4' ~ 'ε3/ε4',
    apoe == 'e4e4' ~ 'ε4/ε4',
    TRUE ~ NA_character_
  )
) 

ggsave("Figs/UKBB/protein_level_byage.svg", 
       levelbyage_inneg(ukbb, c('CSF_LDLR_NPX', 'CSF_PLA2G7_NPX'),
                        c('LDLR','PLA2G7'), 'apoe_plot', 
                        'OLINK NPX value',
                        'Changes in the levels of key proteins with age'),
       width = 12, 
       height = 6,
       dpi = 300)
#---------------------------Extended Fig. 8b-e: NEFL-------------------------------------
bf2soma_df = read.csv('data/BF2Somalogic/preprocessed/with24_7285prot.csv')
get_demographic(bf2soma_df, 'BF2SomaLogic', 'APOE4', FALSE)
bf2soma_df$APOE = ifelse(bf2soma_df$apoe_genotype_baseline_variable %in% c(22, 23), 
                         '22/23', bf2soma_df$apoe_genotype_baseline_variable)

bf2soma_df$Abnormal_CSF_Ab42_Ab40_Ratio = factor(bf2soma_df$Abnormal_CSF_Ab42_Ab40_Ratio,
                                                 levels = c(0,1), labels = c('Aβ-', 'Aβ+'))

SimoaNTK_df = read.csv('data/BF2Somalogic/preprocessed/drug_nfl_with24.csv')
SimoaNTK_df$PL_NFlight_pgmL_Simoa_UGOT_2022 = log2(SimoaNTK_df$PL_NFlight_pgmL_Simoa_UGOT_2022)
SimoaNTK_df$CSF_NFL_pgml_Imputed_NTK_2020 = log2(SimoaNTK_df$CSF_NFL_pgml_Imputed_NTK_2020)
SimoaNTK_df$APOE = ifelse(SimoaNTK_df$apoe_genotype_baseline_variable %in% c(22, 23), '22/23', 
                          SimoaNTK_df$apoe_genotype_baseline_variable)
SimoaNTK_df$Abnormal_CSF_Ab42_Ab40_Ratio = factor(SimoaNTK_df$Abnormal_CSF_Ab42_Ab40_Ratio,
                                                  levels = c(0,1), labels = c('Aβ-', 'Aβ+'))

adni_df = read.csv('data/ADNI/preprocessed/adni_humanprotein_V5align.csv')
get_demographic(adni_df, 'ADNI', 'APOE4', FALSE)
adni_df$APOE = ifelse(adni_df$APOE %in% c(22, 23), '22/23', adni_df$APOE)
adni_df$AB_CSF_status = ifelse(adni_df$AB_CSF_status == 'AB_Neg', 0, ifelse(
  adni_df$AB_CSF_status == 'AB_Pos', 1, NA
))
adni_df$AB_CSF_status = factor(adni_df$AB_CSF_status,
                               levels = c(0,1), labels = c('Aβ-', 'Aβ+'))

adni_tmt_df = read.csv('data/ADNI/preprocessed/adni_soma_ms_merged.csv')
get_demographic(adni_tmt_df, 'ADNI', 'APOE4', FALSE)
adni_tmt_df$APOE = ifelse(adni_tmt_df$APOE %in% c(22, 23), '22/23', adni_tmt_df$APOE)
adni_tmt_df$AB_CSF_status = ifelse(adni_tmt_df$AB_CSF_status == 'AB_Neg', 0, ifelse(
  adni_tmt_df$AB_CSF_status == 'AB_Pos', 1, NA
))
adni_tmt_df$AB_CSF_status = factor(adni_tmt_df$AB_CSF_status,
                                   levels = c(0,1), labels = c('Aβ-', 'Aβ+'))

bf2olink_df = read.csv('data/BF2OLINK/preprocessed/V5_with24_1391proteins.csv')
get_demographic(bf2olink_df, 'BF2OLINK', 'APOE4', FALSE)
bf2olink_df$APOE = ifelse(bf2olink_df$apoe_genotype_baseline_variable %in% c(22, 23), 
                          '22/23', bf2olink_df$apoe_genotype_baseline_variable)
bf2olink_df$Abnormal_CSF_Ab42_Ab40_Ratio = factor(bf2olink_df$Abnormal_CSF_Ab42_Ab40_Ratio,
                                                  levels = c(0,1), labels = c('Aβ-', 'Aβ+'))

nefl_lists = list(
  'SomaLogic (BF2)' = list(tissue='plasma', 
                           target_col='NFL__seq.10082.251__ANML',
                           df=bf2soma_df,
                           covs=c('age', 'gender_baseline_variable', 'plasma_ANML_mean'),
                           ab_col='Abnormal_CSF_Ab42_Ab40_Ratio',
                           apoe_cols=c('apoe_class_num', 'apoe_class_num_e2')),
  'Simoa (BF2)' = list(tissue='plasma',
                       target_col='PL_NFlight_pgmL_Simoa_UGOT_2022',
                       df=SimoaNTK_df,
                       covs=c('age', 'gender_baseline_variable'),
                       ab_col='Abnormal_CSF_Ab42_Ab40_Ratio',
                       apoe_cols=c('apoe_class_num', 'apoe_class_num_e2')),
  'SomaLogic (ADNI)' = list(tissue='CSF',
                            target_col='NEFL__seq.10082.251__ANML',
                            df=adni_df,
                            covs=c('AGE', 'PTGENDER', 'CSF_ANML_mean'),
                            ab_col='AB_CSF_status',
                            apoe_cols=c('apoe_class_num', 'apoe_class_num_e2')),
  'TMT-MS (ADNI)' = list(tissue='CSF',
                         target_col='NEFL_P07196_MS',
                         df=adni_tmt_df,
                         covs=c('AGE', 'PTGENDER', 'CSF_TMT_mean'),
                         ab_col='AB_CSF_status',
                         apoe_cols=c('apoe_class_num', 'apoe_class_num_e2')),
  'OLINK (BF2)' = list(tissue='CSF',
                       target_col='NEFL__Neurology_o2csf__NPX',
                       df=bf2olink_df,
                       covs=c('age', 'gender_baseline_variable', 'CSF_NPX_mean'),
                       ab_col='Abnormal_CSF_Ab42_Ab40_Ratio',
                       apoe_cols=c('apoe_class_num', 'apoe_class_num_e2')),
  'NTK (BF2)' = list(tissue='CSF',
                     target_col='CSF_NFL_pgml_Imputed_NTK_2020',
                     df=SimoaNTK_df,
                     covs=c('age', 'gender_baseline_variable'),
                     ab_col='Abnormal_CSF_Ab42_Ab40_Ratio',
                     apoe_cols=c('apoe_class_num', 'apoe_class_num_e2'))
)

tissue_lists = split(nefl_lists, map_chr(nefl_lists, "tissue"))

plots_by_tissue_group = imap(tissue_lists, function(platform_lists, tissue_name) {
  
  group_plots = map(c("APOE", "AB"), function(group_key) {
    
    panel_plots = imap(platform_lists, function(x, platform_name) {
      
      group_col = if (group_key == "APOE") "APOE" else x$ab_col
      additional_cov = if (group_key == "APOE") x$ab_col else x$apoe_cols
      
      res = proteinlevelby_auto(
        x$df,
        x$target_col,
        c(x$covs, additional_cov),
        group_col,
        platform_name,
        paste0("NEFL (", tissue_name, ")"),
        NULL
      )
      res$stats$.y. = glue('NEFL ({platform_name})')
      res$group_summary$.y. = glue('NEFL ({platform_name})')
      res
    })
    #Source data
    group_summary = bind_rows(imap(panel_plots, function(x, n) x$group_summary))
    stats = bind_rows(imap(panel_plots, function(x, n) x$stats))
    write.xlsx(group_summary, glue('SourceData/NEFL_sum_{tissue_name}_{group_key}.xlsx'))
    write.xlsx(stats, glue('SourceData/NEFL_stats_{tissue_name}_{group_key}.xlsx'))
    
    panel_plots
  })
  set_names(group_plots, c("APOE", "AB"))
})

# wrap_plots(imap(plots_by_tissue_group$plasma$APOE, function(x, n) x$plot))
# wrap_plots(imap(plots_by_tissue_group$plasma$AB, function(x, n) x$plot))
# wrap_plots(imap(plots_by_tissue_group$CSF$APOE, function(x, n) x$plot))
# wrap_plots(imap(plots_by_tissue_group$CSF$AB, function(x, n) x$plot))

iwalk(plots_by_tissue_group, function(apoeORab, tissue){
  iwalk(apoeORab, function(plots, apoe_or_ab){
    plots = wrap_plots(map(plots, function(x) x$plot))
    ggsave(
      glue("Figs/other/NEFL_{tissue}_{apoe_or_ab}.svg"),
      plots + plot_layout(nrow = 1) + theme(
        plot.margin = margin(1, 0.1, 1, 0.1)),
      width = 3*length(plots),
      height = 5,
      dpi = 300,
      device = "svg"
    )
  })
})

#--Supplementary Fig. 3_2: APOE-proteomics signatures across cohorts-------------
source('APOE_proteomics/analysis/betaVSbetaMatrix.R')