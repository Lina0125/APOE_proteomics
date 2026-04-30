source('APOE_proteomics/analysis/setting.R')
results_folder = 'APOE_proteomics/results'

get_last_token = function(x) str_extract(x, "[^-]+$")

plot_beta_vs_beta_matrix = function(df1, df2, title_, xlabel, ylabel, 
                                     merge_by = 'Protein_id', 
                                     label_col = 'symbol',
                                     method = 'pearson',
                                     sig_cor='#DC0000A2'){
  
  df1 = df1[order(df1[[merge_by]], df1$p_adjusted), ]
  df1 = df1[!duplicated(df1[[merge_by]]), ]
  
  df2 = df2[order(df2[[merge_by]], df2$p_adjusted), ]
  df2 = df2[!duplicated(df2[[merge_by]]), ]
  
  merged = merge(df1, df2, by = merge_by)
  
  merged$color = with(
    merged,
    ifelse(p_adjusted.x < 0.05 & p_adjusted.y < 0.05, "both",
           ifelse(p_adjusted.x < 0.05 & p_adjusted.y >= 0.05, "df1_only",
                  ifelse(p_adjusted.x >= 0.05 & p_adjusted.y < 0.05, "df2_only",
                         "nosig")))
  )
  
  ct=table(merged$color)
  
  cor_test = cor.test(
    merged$standardized_beta.x,
    merged$standardized_beta.y,
    method = method
  )
  cor_value = round(cor_test$estimate, 3)
  p_value = cor_test$p.value
  p_text = if (cor_test$p.value == 0) {
    "p < 2.2e-16"
  } else {
    paste0("p = ", sprintf("%.2e", p_value))
  }
  
  present=names(ct)[as.integer(ct) > 0]
  
  xr=range(merged$standardized_beta.x, finite = TRUE)
  yr=range(merged$standardized_beta.y, finite = TRUE)
  
  x_dot  = xr[2] - 0.04 * diff(xr)
  x_text=xr[2] - 0.10 * diff(xr)
  y0     = yr[1] + 0.01 * diff(yr)
  ystep  = 0.075 * diff(yr)
  
  count_df=data.frame(
    x_dot  = rep(x_dot,  length(present)),
    x_text = rep(x_text, length(present)),
    y      = y0 + ystep * (seq_along(present) - 1),
    n      = as.integer(ct[present]),
    color  = present
  )
  
  col_map=c(
    both     = sig_cor,
    df1_only = "#0072B2",
    df2_only = "#E69F00",
    nosig    = "grey"
  )
  
  p = ggplot(merged, aes(x = standardized_beta.x, y = standardized_beta.y, color = color)) +
    geom_point(size = 3) + 
    geom_smooth(method = "lm", se = TRUE, linetype = "dashed",
                color = 'darkgrey') +
    # R / p
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.05, vjust = 1.1,
      label = paste0("R = ", cor_value, ",\n", p_text),
      size = 5, color = "black"
    ) +
    
    # Count
    geom_point(
      data = count_df,
      aes(x = x_dot, y = y),
      inherit.aes = FALSE,
      size = 5,
      color = unname(col_map[count_df$color])
    ) +
    geom_text(
      data = count_df,
      aes(x = x_text, y = y, label = n),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 0.5,
      size = 4.5,
      color = "black"
    ) +
    
    scale_color_manual(
      values = c(
        both     = sig_cor,
        df1_only = "#0072B2",
        df2_only = "#E69F00",
        nosig    = "grey"
      )
    ) +
    
    geom_text_repel(
      aes(label = ifelse(color != 'nosig', merged[[label_col]], NA)),
      force = 0.5, force_pull = 2, direction = "y",
      size = 6, box.padding = 0, point.padding = 0,
      max.overlaps = 6
    ) +
    
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) +
    
    theme_classic(base_size = 20) +
    theme(
      plot.title = element_text(size = 20),
      legend.position = 'none'
    ) +
    labs(title = title_, x = xlabel, y = ylabel)
  
  return(p)
}

walk(names(apoe_list), function(apoe_identifier){
  
  apoe = apoe_list[[apoe_identifier]]
  
  walk(c('apoe2protein', 'apoe2protein_adjab'), function(file){
    if (file == 'apoe2protein') {
      data_list = list(
        'GNPC-plasma-SomaLogic' = glue('{results_folder}/GNPC/{apoe}'),
        'BF2-plasma-SomaLogic'  = glue('{results_folder}/BF2SomaLogic/{apoe}'),
        'ADNI-CSF-SomaLogic'    = glue('{results_folder}/ADNI/{apoe}'),
        'ADNI-CSF-TMT'          = glue('{results_folder}/ADNI_MS/{apoe}'),
        'BF2-CSF-OLINK'         = glue('{results_folder}/BF2OLINK/{apoe}'),
        'UKBB-plasma-OLINK'     = glue('{results_folder}/UKBB/{apoe}')
      )
      
      dfs = imap(data_list, \(folder, dataset_name){
        readr::read_csv(glue("{folder}/{file}.csv"), show_col_types = FALSE)
        
      })
      
    } else {
      data_list = list(
        'GNPC-plasma-SomaLogic' = glue('{results_folder}/GNPC/{apoe}'),
        'BF2-plasma-SomaLogic'  = glue('{results_folder}/BF2SomaLogic/{apoe}'),
        'ADNI-CSF-SomaLogic'    = glue('{results_folder}/ADNI/{apoe}'),
        'ADNI-CSF-TMT'          = glue('{results_folder}/ADNI_MS/{apoe}'),
        'BF2-CSF-OLINK'         = glue('{results_folder}/BF2OLINK/{apoe}')
      )
      
      dfs = imap(data_list, \(folder, dataset_name){
        readr::read_csv(glue("{folder}/{file}.csv"), show_col_types = FALSE)
      })
    }
    
    walk(c(TRUE, FALSE), function(common){
      if (common){
        common_symbols = reduce(
          map(dfs, ~ .x %>%
                distinct(symbol) %>%
                filter(!is.na(symbol)) %>%
                pull(symbol)),
          intersect
        )
        dfs = map(dfs, ~ .x %>% filter(symbol %in% common_symbols))
      }
      
      ds_names = names(dfs)
      n = length(ds_names)
      
      # Non-Diagonal panel
      diag_hist_cell = function(df){
        ggplot(df, aes(x = standardized_beta)) +
          geom_histogram(
            bins = 25,
            fill  = "#ffffff",
            color = "#000000",
            alpha = 0.5
          ) +
          geom_vline(
            xintercept = 0,
            linetype = "dashed",
            color = "grey40",
            linewidth = 0.4
          ) +
          theme_classic(base_size = 20) +
          theme(
            axis.title = element_blank(),
            axis.text  = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            plot.margin = margin(1, 1, 1, 1)
          )
      }
      
      # Diagonal panel
      offdiag_cell = function(x_name, y_name){
        a_last = get_last_token(x_name)
        b_last = get_last_token(y_name)
        
        merge_key = ifelse((a_last == "SomaLogic") && (b_last == "SomaLogic"),
                           "ref_name", "symbol")
        
        label_key = ifelse((a_last == "SomaLogic") && (b_last == "SomaLogic"),
                           "symbol.x", "symbol")
        
        p = plot_beta_vs_beta_matrix(
          dfs[[x_name]],
          dfs[[y_name]],
          '',
          x_name,
          y_name,
          merge_by = merge_key,
          label_col = label_key,
          sig_cor  = '#DC0000A2'
        )
        
        p + theme(
          plot.title = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          plot.margin = margin(0.1,0.1,0.1,0.1)
        )
      }
      
      # n×n matrix
      plot_list = vector("list", n*n)
      k = 1
      for(i in seq_len(n)){
        for(j in seq_len(n)){
          y_name = ds_names[i]
          x_name = ds_names[j]
          
          if(i == j){
            plot_list[[k]] = diag_hist_cell(dfs[[x_name]])
          } else {
            plot_list[[k]] = offdiag_cell(x_name = x_name, y_name = y_name)
          }
          k = k + 1
        }
      }
      
      title_ = glue("Pair correlation of {apoe_identifier}-effect on proteins across cohorts")
      
      if (file == 'apoe2protein_adjab') {
        title_ = glue("{title_} (with adj. Aβ/AD)")
      }
      
      big_matrix = patchwork::wrap_plots(plot_list, ncol = n) +
        patchwork::plot_annotation(
          title = title_,
          theme = theme(plot.title = element_text(size = 35, hjust = 0.5),
                        plot.margin = margin(t = -15, r = 0, b = 0, l = 0))
        )
      
      final_plot = ggdraw(big_matrix, xlim = c(-0.015, 1.01), ylim = c(-0.015, 1.01))+
        theme(
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA)
        )
      # Right rownames
      for(j in seq_along(ds_names)){
        final_plot = final_plot +
          draw_label(
            ds_names[j],
            x = (j - 0.5) / n,
            y = 0,
            angle = 0,
            #fontface = "bold",
            size = 28,
            hjust = 0.5,
            vjust = 1
          )
      }
      
      # Left rownames
      for(i in seq_along(ds_names)){
        final_plot = final_plot +
          draw_label(
            ds_names[i],
            x = 0,
            y = 1 - (i - 0.5) / n,
            angle = 90,
            #fontface = "bold",
            size = 28,
            hjust = 0.5,
            vjust = 0.5
          )
      }
      ggsave(
        ifelse(common, 
               glue("Figs/other/BetaVsBetaMatrix_{apoe_identifier}_{file}_onlyCommon.svg"),
               glue("Figs/other/BetaVsBetaMatrix_{apoe_identifier}_{file}.svg")),
        final_plot,
        width  = 5 * n,
        height = 5 * n,
        dpi = 300,
        device = "svg"
      )
      
      invisible(NULL)
    })
  })
})
