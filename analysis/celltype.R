install.packages("pkgbuild")
pkgbuild::has_build_tools(debug = TRUE)

install.packages("BiocManager")
BiocManager::install("EWCE")

library(EWCE)
library(ggdendro)
library(readxl)
library(readr)
library(dplyr)

#---------------------------Generate ctd using 81 cell types--------------------
expData = read.csv('celltype/ctd/wholebody_celltypes2ctd')
rownames(expData) = expData$Gene
expData = expData[, -1]

annotLevels = list(l1 = colnames(expData))

fNames_ALLCELLS = EWCE::generate_celltype_data(
  exp = expData,
  annotLevels = annotLevels,
  groupName = "Human",
  no_cores = 1,
  savePath = tempdir(),
  file_prefix = "ctd",
  as_sparse = TRUE,
  as_DelayedArray = FALSE,
  normSpec = FALSE,
  convert_orths = FALSE,
  input_species = "mouse",
  output_species = "human",
  non121_strategy = "drop_both_species",
  method = "homologene",
  force_new_file = TRUE,
  specificity_quantiles = TRUE,
  numberOfBins = 40,
  dendrograms = TRUE,
  return_ctd = FALSE,
  verbose = TRUE
)

wholebd = EWCE::load_rdata(fNames_ALLCELLS)

#EWCE::plot_ctd(wholebd, c('APOE','SPC25', 'TBCA', 'S100A13'), level = 1, metric = "specificity", show_plot = TRUE)

#------Cell type enrichment for early dysregulated proteins in APOE2 carriers----
library(readxl)
early_apoe = read_excel("celltype/data/early_apoe2.xlsx")
background = subset(early_apoe, bg_genesymbol != '')
up = subset(early_apoe, `early Regulation` == 'up' & early_genesymbol != '')
down = subset(early_apoe, `early Regulation` == 'down' & early_genesymbol != '')

bg_genes = background$bg_genesymbol
up_genes = as.vector(up$early_genesymbol)
down_genes = as.vector(down$early_genesymbol)

set.seed(123)
reps = 10000
annotLevel = 1

bg_genes = bg_genes[bg_genes != 'HES1']
up_results = EWCE::bootstrap_enrichment_test(sct_data = wholebd, 
                                                bg = bg_genes,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = up_genes,
                                                reps = reps,
                                                annotLevel = annotLevel)

down_results = EWCE::bootstrap_enrichment_test(sct_data = wholebd, 
                                             bg = bg_genes,
                                             sctSpecies = "human",
                                             genelistSpecies = "human",
                                             hits = down_genes,
                                             reps = reps,
                                             annotLevel = annotLevel)

View(down_results$results)
write.csv(up_results$results,'celltype/data/results/apoe2_early_up.csv', row.names = FALSE)
write.csv(down_results$results, 'celltype/data/results/apoe2_early_down.csv', row.names = FALSE)

#-------------------------------Bidirectional sets------------------------------
dat = read_csv("celltype/data/bidirectional.csv")

bg_genes = dat$bg
oppo = as.vector(dat$apoe2upstream_downinapoe2_apoe4downstream_upinad)
same = as.vector(dat$apoe2upstream_downinapoe2_apoe4downstream_downinad)
both_early = as.vector(dat$bothapoe2andapoe4early)
set.seed(123)
reps = 10000
annotLevel = 1

bg_genes = bg_genes[bg_genes != 'HES1']
oppo_results = EWCE::bootstrap_enrichment_test(sct_data = wholebd, 
                                             bg = bg_genes,
                                             sctSpecies = "human",
                                             genelistSpecies = "human",
                                             hits = oppo,
                                             reps = reps,
                                             annotLevel = annotLevel)

same_results = EWCE::bootstrap_enrichment_test(sct_data = wholebd, 
                                               bg = bg_genes,
                                               sctSpecies = "human",
                                               genelistSpecies = "human",
                                               hits = same,
                                               reps = reps,
                                               annotLevel = annotLevel)

both_early_results = EWCE::bootstrap_enrichment_test(sct_data = wholebd, 
                                               bg = bg_genes,
                                               sctSpecies = "human",
                                               genelistSpecies = "human",
                                               hits = both_early,
                                               reps = reps,
                                               annotLevel = annotLevel)

View(both_early_results$results)

#---------------------------BF2SomaLogic----------------------------------------
dat = read_csv("celltype/data/bf2soma.csv")
bg_genes = dat$bg
apoe4 = as.vector(dat$bf2soma_abneg_apoe4)
apoe2 = as.vector(dat$bf2soma_abneg_apoe2)
bothearly = as.vector(dat$bf2soma_both_early)
set.seed(123)
reps = 10000
annotLevel = 1

bg_genes = bg_genes[bg_genes != 'HES1']
apoe4_results = EWCE::bootstrap_enrichment_test(sct_data = wholebd, 
                                               bg = bg_genes,
                                               sctSpecies = "human",
                                               genelistSpecies = "human",
                                               hits = apoe4,
                                               reps = reps,
                                               annotLevel = annotLevel)

bg_genes = bg_genes[bg_genes != 'HES1']
apoe2_results = EWCE::bootstrap_enrichment_test(sct_data = wholebd, 
                                                bg = bg_genes,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = apoe2,
                                                reps = reps,
                                                annotLevel = annotLevel)

bothearly_results = EWCE::bootstrap_enrichment_test(sct_data = wholebd, 
                                                bg = bg_genes,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = bothearly,
                                                reps = reps,
                                                annotLevel = annotLevel)

View(apoe4_results$results)

#-------------------------Grouped proteins--------------------------------------
apoe2 = read_excel("ref/gnpc_apoe2associated_192apt_0828_previousFormat.xlsx")
apoe4 = read_excel("ref/gnpc_apoe4associated_357apt_0828_previousFormat.xlsx")

bg_genes = as.vector(unique(apoe2$bg_genesymbol))
set.seed(123)
reps = 10000
annotLevel = 1
bg_genes = bg_genes[bg_genes != 'HES1']

grouped_enrichment = function(df, bg_genes, genesymbols){
  grouped_results = data.frame()
  for (group in unique(df$group)) {
    for (reg in unique(df$Regulation)) {
      tmp = df[(df$group == group) & (df$Regulation == reg),]
      genes = as.vector(unique(tmp[[genesymbols]]))
      lb = paste0(group, '_', reg)
      
      if (length(genes) < 4) {
        next
      }
      
      results = EWCE::bootstrap_enrichment_test(sct_data = wholebd, 
                                                bg = bg_genes,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = genes,
                                                reps = reps,
                                                annotLevel = annotLevel)
      
      result = results$results
      result = result[result$p<0.05,]
      result$Cluster_Group_name = lb
        
      grouped_results = dplyr::bind_rows(grouped_results, result)
    }
  }
  
  return(grouped_results)
}

apoe2_results = grouped_enrichment(apoe2, bg_genes, 'APOE2_genesymbol')

apoe2_results_sorted = apoe2_results %>%
  group_by(Cluster_Group_name) %>%
  arrange(p, .by_group = TRUE)

apoe4_results = grouped_enrichment(apoe4, bg_genes, 'APOE4_genesymbol')

apoe4_results_sorted = apoe4_results %>%
  group_by(Cluster_Group_name) %>%
  arrange(p, .by_group = TRUE)

write_csv(apoe2_results_sorted,'celltype/data/results/apoe2_grouped.csv')

write_csv(apoe4_results_sorted,'celltype/data/results/apoe4_grouped.csv')
