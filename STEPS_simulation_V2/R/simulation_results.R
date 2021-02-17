if(run.flag) {
  rm(list = ls())
  
  detachAllPackages <- function() {
    
    basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
    
    package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
    
    package.list <- setdiff(package.list,basic.packages)
    
    if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
    
  }
  
  detachAllPackages()
  
  
  if(!require('data.table'))  install.packages('data.table')
  if(!require('magrittr'))    install.packages('magrittr')
  if(!require('dplyr'))       install.packages('dplyr')
  if(!require('tidyr'))       install.packages('tidyr')
  if(!require('ggplot2'))     install.packages('ggplot2')
  
  
  
  all_output_files <- list.files('Output\\')
  all_output_files <- all_output_files[grepl('.csv', all_output_files)]
  class_res_files <- all_output_files[grepl("Classification", all_output_files) & grepl("Sce_i", all_output_files)]
  Cors_res_files <- all_output_files[grepl("Cors", all_output_files) & grepl("Sce_i", all_output_files)]
  Concordance_res_files <- all_output_files[grepl("Concordance", all_output_files) & grepl("Sce_i", all_output_files)]
  
  
  imp_fun <- function(files) {
    tmp_list <- rep(list(NA), length(files))
    for(i in 1:length(files)) {
      tmp_list[[i]] <- read.csv(paste0('Output\\', files[i]))
    }
    out_tab <- do.call('rbind', tmp_list)
    return(out_tab)
  }
  
  class_res <- imp_fun(class_res_files)
  Cors_res <- imp_fun(Cors_res_files)
  Concordance_res <- imp_fun(Concordance_res_files)
  
  colnames(Concordance_res) <- gsub("obserevd", "observed", colnames(Concordance_res))
  
  head(Concordance_res)
  
  
  ###### Concordance results
  Conc_res_plot <- Concordance_res %>% 
    dplyr::select(#-contains("_train"),
      -Iteration) %>% 
    pivot_longer(cols = c(concB.nldss.observed, concB.nldss.oracle,
                          # concB.lm1.observed, concB.lm1.oracle, 
                          # concB.lm2.observed, concB.lm2.oracle,
                          concW.nldss.observed, concW.nldss.oracle, 
                          concR.nldss.observed, concR.nldss.oracle
                          # , concW.lm1.observed, concW.lm2.observed
    ), names_to = "Method", values_to = "concordance_Between") %>% 
    separate(Method, c("Set", "Method", "oracleObserved")) %>%
    as.data.frame() %>% 
    group_by(Set, Method, oracleObserved, Obs, p, Noise_level, Noise_ratio, lambdaB, lambdaW, lambdaR, Rho, ntrees) %>%
    summarize("avg" = mean(concordance_Between),
              "SD" = sd(concordance_Between)) %>% 
    ungroup()
  Conc_res_plot <- Conc_res_plot %>% 
    mutate("Set2" = case_when(Set == "concB" ~ "Set B", 
                              Set == "concW" ~ "Set W", 
                              Set == "concR" ~ "Set R")) %>%
    filter(lambdaW != 10)
  Conc_res_plot$Set2 <- factor(Conc_res_plot$Set2, levels = c("Set B", "Set W", "Set R"), ordered = T)
  pBWconc <- Conc_res_plot %>%
    ggplot(aes(Noise_ratio, avg, color = paste( Method, oracleObserved) )) +
    geom_line(size = 1.1) +
    geom_point(size = 1.2) +
    geom_errorbar(aes(ymin = avg - SD, ymax = avg + SD ), width = 0.5,
                  position = position_dodge(0.05), size = 1.2) + 
    facet_grid (Set2 ~ lambdaW , scales = "free_y") + 
    ylab("Concordance percent") + xlab("Noise ratio") +
    theme_bw() +
    theme(text = element_text(face = "bold", size = 12), 
          legend.position = "bottom") + 
    scale_color_manual(name = "Method", values = c("#F8766D", "#00BFC4"), 
                       labels = c("STEPS", "Oracle STEPS")) + 
    guides(col=guide_legend(nrow=1,byrow=TRUE))
  
  
  
  
  ggsave("fig\\BWconcrdance.png", plot = pBWconc, width=10.5, height = 6, dpi=500)
  ggsave("fig\\tiff\\BWconcrdance.tiff", plot = pBWconc, width=10.5, height = 6, dpi=500, device = "tiff")
  
  
  
  ###### Correlations results
  
  cors_res_plot <- Cors_res%>% 
    dplyr::select( -Iteration) %>% 
    group_by(Method, index_type , Obs, p, 
             Noise_level, Noise_ratio,
             lambdaB, lambdaW, lambdaR, Rho, ntrees) %>%
    summarise_all(list("avg" = mean,
                       "SD" = sd))
  
  
  dat_STEPS <-  cors_res_plot %>%
    filter(index_type %in% c("", "mean_diff_from_median", "positive_mean_from_median"), 
           lambdaW != 10, Method %in% c("NLDSS", "NLDSS_oracle"))
  dat_Scores <-  cors_res_plot %>%
    filter(index_type %in% c("", "mean_diff_from_median", "positive_mean_from_median"), 
           lambdaW != 10, !(Method %in% c("NLDSS", "NLDSS_oracle")))
  pcors <- dat_STEPS %>% 
    ggplot(aes(Noise_ratio , avg, group = paste(Method, index_type), color = Method), size = 1.2) + 
    geom_line(aes(linetype = paste(index_type)), size = 1.1) +
    geom_point(size = 1.2) + 
    geom_line(data = dat_Scores, aes(color = Method), size = 1.1) +
    geom_point(data = dat_Scores, aes(color = Method), size = 1.2) + 
    facet_grid( ~ lambdaW) +
    theme_bw() +
    theme(text = element_text(face = "bold", size = 12), 
          legend.position = "bottom") + 
    ylab("Spearman correlation") + xlab("Noise ratio") + 
    scale_linetype_manual(name = "Index Aggregation", values = c("solid", "dashed"), 
                          labels = c("Mean diffrence from median", "Mean positive diffrence from median")) + 
    scale_color_manual(values = c("#F8766D", "#00BFC4", "#7CAE00", "#C77CFF"), 
                       labels = c("STEPS", "Oracle STEPS", "Score1", "Score2"))+ 
    guides(col=guide_legend(nrow=2,byrow=TRUE), linetype = guide_legend(override.aes = list(size = 0.5)))
  
  
  
  
  
  ggsave("fig\\correlation.png", plot = pcors, width=10.5, height = 6, dpi=500)
  ggsave("fig\\tiff\\correlation.tiff", plot = pBWconc, width=10.5, height = 6, dpi=500, device = "tiff")
  
  
  
  ##### Classification results
  head(class_res)
  class_res_plot <- class_res %>% 
    dplyr::select(-Iteration) %>%
    group_by(Obs, p,  Noise_level, Noise_ratio, 
             ntrees, Rho, lambdaB, lambdaW, lambdaR, Method, index_type) %>% 
    summarise_all(list("AUC_avg" = mean,
                       "SD" = sd)) %>% 
    ungroup()
  
  
  class_res_plot <- class_res %>% 
    # filter(index_type %in% c(""))
    dplyr::select(-Iteration) %>%
    group_by(Obs, p,  Noise_level, Noise_ratio, 
             ntrees, Rho, lambdaB, lambdaW, lambdaR, Method, index_type) %>% 
    summarise_all(list("AUC_avg" = mean,
                       "SD" = sd)) %>% 
    ungroup()
  
  
  
  
  dat_STEPS <-  class_res_plot %>%
    filter(index_type %in% c("", "mean_diff_from_median", "positive_mean_from_median"), 
           lambdaW != 10, Method %in% c("NLDSS", "NLDSS_oracle"))
  dat_Scores <-  class_res_plot %>%
    filter(index_type %in% c("", "mean_diff_from_median", "positive_mean_from_median"), 
           lambdaW != 10, !(Method %in% c("NLDSS", "NLDSS_oracle")))
  pclass <- dat_STEPS %>% 
    ggplot(aes(Noise_ratio, AUC_avg)) +
    geom_line(aes(linetype = index_type, color = Method), size = 1.1) +
    geom_point(aes(color = Method, group = paste(Method, index_type)), size = 1.2) +
    geom_line(data = dat_Scores, aes(Noise_ratio, AUC_avg, color = Method), size = 1.1) +
    geom_point(data = dat_Scores, aes(color = Method), size = 1.2) +
    facet_grid( ~ lambdaW  ) + 
    theme_bw() +
    theme(text = element_text(face = "bold", size = 12), 
          legend.position = "bottom") + 
    ylab("AUC") + xlab("Noise ratio") +
    scale_linetype_manual(name = "Index Aggregation", values = c("solid", "dashed"), 
                          labels = c("Mean diffrence from median", "Mean positive diffrence from median")) + 
    scale_color_manual(values = c("#F8766D", "#00BFC4", "#7CAE00", "#C77CFF"), 
                       labels = c("STEPS", "Oracle STEPS", "Score1", "Score2"))+ 
    guides(col=guide_legend(nrow=2,byrow=TRUE), linetype = guide_legend(override.aes = list(size = 0.5)))
  
  
  
  ggsave("fig\\classification.png", plot = pclass, width=10.5, height = 6, dpi=500)
  ggsave("fig\\tiff\\classification.tiff", plot = pBWconc, width=10.5, height = 6, dpi=500, device = "tiff")
  
  
}
