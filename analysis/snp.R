library(TwoSampleMR)
library(purrr)
library(stringr)
library(biomaRt)
library(dplyr)

find_snps_in_genes = function(mediator_col, mediators, gwas, ensembl){
  
  gene_annotations = getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
    filters = "hgnc_symbol",
    values = mediators[[mediator_col]],
    mart = ensembl
  )
  
  gene_anno_clean = gene_annotations %>%
    filter(grepl("^([0-9]+|X|Y|MT)$", chromosome_name))
  
  map_dfr(unique(gene_anno_clean$hgnc_symbol), function(gene){
    
    tmp_pos = filter(gene_anno_clean, hgnc_symbol == gene)
    
    res = subset(gwas,
                 chr.outcome == tmp_pos$chr &
                   pval.outcome < 0.05 &
                   pos.outcome >= tmp_pos$start &
                   pos.outcome <= tmp_pos$end)
    
    if (nrow(res) > 0) {
      mutate(res, gene = gene)
    } else {
      print(paste0('No SNP found in ', gene, ' gene'))
      NULL
    }
  })
}
#----------------------------Load GWAS------------------------------------------
##--Load key mediators--
mediators = read.csv("ref/mediators.csv")

ensembl = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

##--Load Abeta GWAS--
abeta_outcome = read_outcome_data(
  filename = "data/Bf2Somalogic/preprocessed/GWAS/sums/GCST90129599_buildGRCh38.tsv",
  sep = ",",
  snp_col = "ID_38",
  beta_col = "beta",
  se_col = "SE",
  effect_allele_col = "effect_allele",
  other_allele_col = "",
  eaf_col = "effect_allele_frequency",
  pval_col = "p-value",
  samplesize_col = "N"
)

tmp = str_split_fixed(abeta_outcome$SNP, ":", 2)
abeta_outcome$chr.outcome = as.numeric(tmp[,1])
abeta_outcome$pos.outcome = as.numeric(tmp[,2])

##--Load AD GWAS--
ad_out = read_outcome_data(
  filename = "data/Bf2Somalogic/preprocessed/GWAS/sums/GCST002245_buildGRCh38_load.txt",
  sep = "\t",
  snp_col = "MarkerName",
  chr_col = "Chromosome",
  pos_col = "Position",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "Effect_allele",
  other_allele_col = "Non_Effect_allele",
  pval_col = "Pvalue",
)

#------------------------------Find SNPs on interested genes--------------------
abeta_apoe4mediators = find_snps_in_genes('abeta_apoe4mediators', mediators, abeta_outcome, ensembl)
abeta_apoe2mediators = find_snps_in_genes('abeta_apoe2mediators', mediators, abeta_outcome, ensembl)

ad_apoe4mediators = find_snps_in_genes('ad_apoe4mediators', mediators, ad_out, ensembl)
ad_apoe2mediators = find_snps_in_genes('ad_apoe2mediators', mediators, ad_out, ensembl)

#-----------------------------------Export--------------------------------------
e4_mediators = c('S100A13', 'TBCA', 'SPC25', 'CTF1', 'LRRN1') #'ARL2', 'BCDIN3D', 'HBQ1'
e2_mediators = c('APOB', 'PCLAF', 'SNAP23', 'WARS2')

table(abeta_apoe4mediators[abeta_apoe4mediators$gene %in% e4_mediators,]$gene)
table(abeta_apoe2mediators[abeta_apoe2mediators$gene %in% e2_mediators,]$gene)

table(ad_apoe4mediators[ad_apoe4mediators$gene %in% e4_mediators,]$gene)
table(ad_apoe2mediators[ad_apoe2mediators$gene %in% e2_mediators,]$gene)

export = function(ab_df, ad_df, mediators){
  ab = ab_df[ab_df$gene %in% mediators,]
  ab$GWAS = 'CSF Aβ42'
  ab$Study = 'GCST90129599'
  
  ad = ad_df[ad_df$gene %in% mediators,]
  ad$GWAS = 'AD'
  ad$Study = 'GCST002245'
  
  res = bind_rows(ab, ad)
  
  return(res)
}

e4 = export(abeta_apoe4mediators, ad_apoe4mediators, e4_mediators)
e2 = export(abeta_apoe2mediators, ad_apoe2mediators, e2_mediators)

library(openxlsx)

write.xlsx(
  list(
    'Key e4 mediators' = e4,
    'Key e2 mediators' = e2
  ),
  file = "Supplement/Supplement_Table_9_SNPs.xlsx"
)

