library(maftools)
library(dplyr)
library(stringr)
library(magrittr)
library(GenomicRanges)
library(ggpubr)


# Data Prep for the processing --------------------------------------------
rm(list = ls())
source('script/Utils.R')
source('script/functions/Map2Selector.R')
source('script/functions/Polisher.R')
source('script/functions/read_MAFs.R')
load("Rdata/NPC_cedat.Rdata")

my_theme = my_theme +
  theme(strip.background = element_blank())
box_fill = ggthemes::wsj_pal()(2)
# scales::show_col(box_fill)


# MAF of NPC data: filtering and polishing -------------------------------
load(file = 'Rdata/NPC_all_SNVs.Rdata') # npc_ftd, cf_gr from MafProcess.R

# iterating by Panels -----------------------------------------------------
# panel = 'Our_ori'
for(panel in c('Our', "Our_ori")){
  
 
  ## Capturing by panels --------------------
  
  ct_ftd = Map2Selector(gr = cf_gr, panel = panel) %>% as.data.frame()
  # save(ct_ftd, file = paste0("Rdata/ctDNA_NPC_PoN_",panel ," .Rdata"))

  ca_len = Map2Selector(panel = panel) %>% width %>% sum
  if(!dir.exists(paste0("fig/NPC/", panel)))
    dir.create(paste0("fig/NPC/", panel), recursive = TRUE)
  
  
  ## get bTMB_df ----------------
  
  bTMB_df = ct_ftd %>%
    mutate(Region = paste0("Exon_", Region)) %>%
    select(Tumor_Sample_Barcode, Region) %>%
    dplyr::rename(Sample = Tumor_Sample_Barcode) %>%
    group_by(Sample) %>%
    summarise(Num = length(Sample), bTMB = Num * 1e6 / ca_len) %>%
    left_join(x = ce_dat[, c("Patient", "Sample", "Group", "Response")],
              y = .,
              by = "Sample")
  bTMB_df[is.na(bTMB_df)] = 0
  
  if(panel == 'Our_ori'){
    ### a/b ratio -----------------------
    bTMB_ab = bTMB_df %>% 
      reshape2::dcast(formula = Patient + Response ~ Group, 
                      value.var = 'bTMB' ) %>% 
      mutate(abRatio = `Post-therapy` - `Pre-therapy`)
    
    bab = bTMB_ab %>% 
      ggplot(aes(x=Response, y=abRatio))+
      geom_boxplot(aes(fill=Response))+
      ggsci::scale_fill_aaas(alpha = 0.8, limits = c('Nonresponder', 'Responder'))+
      stat_compare_means(label.x = 1.25)+
      labs(x=NULL, y='Difference value of a minus b')+
      my_theme+
      theme(legend.position = 'none')
    ggsave(bab, filename = 'fig/NPC/Enricher-bTMB_abRatio.tiff',
           units = 'cm', height = 10, width = 10)
    save(bab, file = 'Rdata/fig/PO/Enricher-bTMB_abRatio.Rdata')
    
    ### ROC plot
    library(pROC)
    roc1 = roc(response = bTMB_ab$Response,
               predictor = bTMB_ab$abRatio %>% unlist)
    
    best = coords(roc1, 'best')[1,] %>% unlist
    roc_df = coords(roc1, transpose = FALSE)
    
    nrROC = ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
      geom_path(color = ggsci::pal_aaas(alpha = 0.7)(2)[2],
                size = 1.5) +
      geom_segment(aes(
        x = 0,
        y = 0,
        xend = 1,
        yend = 1
      ),
      colour = 'grey',
      linetype = 'dotdash') +
      labs(x = "1-Specificity",
           y = "Sensitivity") +
      annotate(geom = 'text', x = 0.6, y=0.15, hjust = 0,
               label = paste0("Sensitivity: ", round(best[['sensitivity']], 3), "\n",
                              "Specificity: ", round(best[['specificity']], 3), "\n",
                              "Threshold: ", round(best[['threshold']], 1), '\n',
                              "AUC: ", round(roc1$auc, 3)
               ))+
      my_theme
    
    ggsave(nrROC, filename = 'fig/NPC/Enricher-ROC-response.tiff', units = 'cm', height = 10, width = 10)
    save(nrROC, file = 'Rdata/fig/PO/NPC_Enricher-ROC-response.Rdata')
    }

  
  ## get KMR_df ------------------
  
  KMR_df = ct_ftd %>%
    mutate(Region = paste0("Exon_", Region)) %>%
    select(Tumor_Sample_Barcode, Region, VAF) %>%
    dplyr::rename(Sample = Tumor_Sample_Barcode) %>%
    # whether regions contains more than 1 muts
    group_by(Sample, Region) %>%
    summarise(mean = mean(VAF), len = length(VAF)) %>%
    ungroup() %>%
    dplyr::filter(!(mean < 1 & len == 1)) %>%
    group_by(Sample) %>%
    summarise(KMR = length(Region) * 1e6 / ca_len) %>%
    ungroup() %>%
    left_join(x = ce_dat[, c("Patient", "Sample", "Group", "Response")],
              y = .,
              by = "Sample") %>%
    as_tibble()
  
  KMR_df$Group = factor(
    KMR_df$Group,
    levels = c("Pre-therapy", "Post-therapy"),
    ordered = T
  )
  
  
  
  ## get VAF -----------------------
  
  VAF_df = ct_ftd %>%
    # ct_ftd %>%
    mutate(Region = paste0("Exon_", Region)) %>%
    select(Tumor_Sample_Barcode, VAF) %>%
    dplyr::rename(Sample = Tumor_Sample_Barcode) %>%
    left_join(x = ce_dat[, c("Patient", "Sample", "Group", "Response")],
              y = .,
              by = "Sample")
  VAF_df[is.na(VAF_df)] = 0
  pat_order = VAF_df %>% dplyr::filter(Group == "Pre-therapy") %>%
    group_by(Response, Patient) %>%
    summarise(mean = mean(VAF)) %>%
    arrange(Response, mean) %>% select(Patient)
  
  VAF_df$Patient = factor(VAF_df$Patient,
                          levels = pat_order$Patient,
                          ordered = T)
  VAF_df$Group = factor(
    VAF_df$Group,
    levels = c("Pre-therapy", "Post-therapy"),
    ordered = T
  )
  VAF_df$SupGroup = factor(
    paste0(VAF_df$Response, ":", VAF_df$Group),
    levels = c(
      "Responder:Pre-therapy",
      "Responder:Post-therapy",
      "Nonresponder:Pre-therapy",
      "Nonresponder:Post-therapy"
    ),
    ordered = T
  )
  
  ### and get the VAF polts
  if(1){
    source("D:/Projects/rProj/SunMoon/scripts/geom_split_violin.R")
    
    ## 
    pv = ggplot(VAF_df, aes(x = Patient, y = VAF, fill = SupGroup)) +
      geom_split_violin(trim = T, colour = NA) +
      geom_point(stat = 'summary',
                 fun = mean,
                 position = position_dodge(width = 0.5)) +
      ggsci::scale_fill_nejm() +
      stat_summary(
        fun.min = function(x) {
          quantile(x)[2]
        },
        fun.max = function(x) {
          quantile(x)[4]
        },
        geom = 'errorbar',
        color = 'black',
        width = 0.01,
        size = 0.5,
        position = position_dodge(width = 0.5)
      ) +
      my_theme +
      theme(legend.position = 'bottom') +
      guides(fill = guide_legend(nrow = 2))
    
    save(pv, file = paste0("Rdata/fig/PO/", panel, "_vaf_every.Rdata"))
    ggsave(
      plot = pv,
      filename = paste0("fig/NPC/", panel, "_vaf_every.jpeg"),
      width = 20,
      height = 15,
      units = 'cm'
    )
    
    ##
    pv2 = ggplot(data = VAF_df, aes(x = Group, y = VAF)) +
      geom_violin(aes(fill = Response), alpha = 0.9, color = NA) +
      geom_boxplot(fill = 'white', width = 0.1) +
      scale_fill_manual(limits = c("Responder", "Nonresponder"),
                        values = box_fill) +
      stat_compare_means(label.x = 1.3,
                         color = p_color,
                         label.y = 105) +
      facet_grid(. ~ Response) +
      my_theme +
      theme(legend.position = 'none',
            text = element_text(size = 15))
    save(pv2, file = paste0("Rdata/fig/PO/", panel, "_vaf_PRvsSD.Rdata"))
    ggsave(plot = pv2, filename = paste0("fig/NPC/", panel, "_vaf_PRvsSD.jpeg"))
    
    ##
    pv3 = ggplot(data = VAF_df, aes(x = Response, y = VAF)) +
      geom_boxplot(aes(fill = Response), alpha = 0.9) +
      scale_fill_manual(limits = c("Responder", "Nonresponder"),
                        values = box_fill) +
      facet_grid(. ~ Group) +
      stat_compare_means(
        method.args = list(alternative = "less"),
        color = p_color,
        label.x = 1.3
      ) + my_theme +
      theme(legend.position = 'none',
            text = element_text(size = 15))
    save(pv3, file = paste0("Rdata/fig/PO/", panel, "_vaf_AfterVsBefore.Rdata"))
    ggsave(plot = pv3, filename = paste0("fig/NPC/", panel, "_vaf_AfterVsBefore.jpeg"))
    
  }
  
  
  ## get mVAF_df -------------------
  
  mVAF_df = VAF_df %>% group_by(Group, Patient, Sample, Response) %>%
    summarise(mVAF = mean(VAF))
  
  
  ## iterating for box plots -------------
  for(var in c('bTMB', 'KMR', 'mVAF')){
    # var = 'bTMB'
    var_df = get(paste0(var, '_df')) %>% 
      dplyr::rename("VAR" = var) %>% 
      as.data.frame()
    

    p1 = ggplot(data = var_df, aes(x = Group, y = VAR)) +
      geom_boxplot(aes(fill = Response), alpha = 0.9) +
      facet_grid(. ~ Response) +
      scale_fill_manual(limits = c("Responder", "Nonresponder"),
                        values = box_fill) +
      stat_compare_means(
        paired = T,
        color = p_color,
        label.x = 1.3
      ) +
      labs(x=NULL, y=var) +
      my_theme +
      theme(legend.position = 'none',
            text = element_text(size = 15))
    save(p1, file = paste0("Rdata/fig/PO/", panel, "_", var, "_Group.Rdata"))
    ggsave(plot = p1, filename = paste0("fig/NPC/", panel, "/", var, "_Group.jpeg"))
    
    ##
    p2 = p1 + geom_point(
      aes(color = Response),
      alpha = 0.6,
      size = 4,
      shape = 16
    ) +
      geom_line(aes(group = Patient, color = Response)) +
      ggthemes::scale_color_tableau(limits = c("Responder", "Nonresponder"))
    save(p2, file = paste0("Rdata/fig/PO/", panel, "_", var, "_PRvsSD_dot.Rdata"))
    ggsave(plot = p2, filename = paste0("fig/NPC/", panel, "/", var, "_PRvsSD_dot.jpeg"))
    
    
    ## 
    p3 = ggplot(data = var_df, aes(x = Response, y = VAR)) +
      geom_boxplot(aes(fill = Response), alpha = 0.9) +
      stat_compare_means(label.x = 1.3, color = p_color) +
      labs(x = NULL, y=var) +
      scale_fill_manual(limits = c("Responder", "Nonresponder"),
                        values = box_fill) +
      facet_grid(. ~ Group) +
      my_theme +
      theme(legend.position = 'none',
            text = element_text(size = 15))
    
    save(p3, file = paste0("Rdata/fig/PO/", panel, "_", var, "_BeforevsAfter.Rdata"))
    ggsave(plot = p3, filename = paste0("fig/NPC/", panel, "/", var, "_BeforevsAfter.jpeg"))
    
    # paired-boxplot
    p4 = ggplot(data = var_df, aes(x = Group, y = VAR)) +
      geom_boxplot(
        fill = 'white',
        alpha = 0.1,
        size = 0.9,
        show.legend = FALSE,
        width = 0.2
      ) +
      geom_point(aes(fill = Response),
                 size = 5,
                 color = 'white',
                 shape = 21) +
      geom_line(aes(group = Patient, color = Response), lwd = 0.5) +
      ggthemes::scale_fill_wsj(limits = c("Responder", "Nonresponder")) +
      ggthemes::scale_color_wsj(limits = c("Responder", "Nonresponder")) +
      stat_compare_means(
        paired = TRUE,
        label.x = 1.3,
        label.y = max(var_df$VAR) * 1.05,
        size = 5,
        color = p_color
      ) +
      labs(x = NULL, y=var) +
      my_theme +
      theme(legend.position = 'bottom',
            text = element_text(size = 15))
    
    save(p4, file = paste0("Rdata/fig/PO/", panel, "_", var, "_pairBoxplot.Rdata"))
    ggsave(
      plot = p4,
      filename = paste0("fig/NPC/", panel, "/", var, "_pairBoxplot_bTMB.png"),
      width = 15,
      height = 15,
      units = 'cm'
    )

  }
  
  
  ## get ROC plots ----------
  
  ### get values ------------
  library(pROC)
  
  roc_df = NULL
  tss_df = NULL
  for (i in c("mVAF_df", "KMR_df", "bTMB_df")) {
    var = str_split(i, "_", simplify = T)[, 1]
    dat = get(i) %>% dplyr::select(Patient, Group, Response, all_of(var))
    dat$Res01 = ifelse(dat$Response == "Responder", 1, 0)
    
    for (grp in c("Pre-therapy", "Post-therapy", "Af-Be")) {
      if (grp != 'Af-Be') {
        dat_1 = dat %>% dplyr::filter(Group == grp)
        roc1 = roc(response = dat_1$Res01,
                   predictor = dat_1[, var] %>% unlist)
      } else {
        dat_b = dat %>% dplyr::filter(Group == "Pre-therapy") %>%
          arrange(Patient)
        dat_a = dat %>% dplyr::filter(Group == "Post-therapy") %>%
          arrange(Patient)
        dat_b$A_B = dat_a[, var] - dat_b[, var]
        roc1 = roc(response = dat_b$Res01,
                   predictor = dat_b$A_B %>% unlist)
      }
      
      best = coords(roc1, 'best')[1,] %>% unlist
      temp_df = coords(roc1, transpose = FALSE) %>%
        mutate(
          auc = round(roc1$auc, digits = 3),
          group = paste(var, "in", grp, "Treat: AUC=", auc),
          grp = grp,
          var = var
        )
      temp_ss = tibble(
        grp = grp,
        var = var,
        auc = round(roc1$auc, digits = 3),
        best_trd = best[['threshold']] %>% round(3),
        best_sen = best[['specificity']] %>% round(3),
        best_spe = best[['sensitivity']] %>% round(3)
      )
      
      roc_df = rbind(roc_df, temp_df)
      tss_df = rbind(tss_df, temp_ss)
    }
  }
  
  ### save the Threshold, Sensitivity and Specificity ----------
  tss_df %>%
    dplyr::filter(var != 'nSNVs') %>%
    arrange(desc(grp)) %>%
    select(
      Group = grp,
      Variables = var,
      AUC = auc,
      Threshold = best_trd,
      Sensitivity = best_sen,
      Specificity = best_spe
    ) %>% 
  openxlsx::write.xlsx(paste0("data/NPC_", panel, "_tss.xlsx"), overwrite = T)
  
  
  
  ### ROC plots -----------
  
  #### bTMB 
  bTMB_p = roc_df %>% 
    dplyr::filter(var == 'bTMB') %>% 
    mutate(group2 = str_replace_all(group, '(bTMB in )|( Treat)', "")) %>% 
    ggplot(aes(x = 1 - specificity, y = sensitivity)) +
    geom_path(aes(color = group2), size = 1.5) +
    ggsci::scale_color_aaas() +
    geom_segment(aes(
      x = 0,
      y = 0,
      xend = 1,
      yend = 1
    ),
    colour = 'grey',
    linetype = 'dotdash') +
    labs(x = "1-Specificity",
         y = "Sensitivity") +
    my_theme +
    # facet_grid(. ~ grp) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = c(0.7, 0.15),
      legend.title = element_blank(),
      legend.background = element_rect(
        fill = NULL,
        size = 0.5,
        linetype = "solid",
        colour = "black"
      )
    ) +
    guides(color = guide_legend(nrow = 3))
  save(bTMB_p, file = paste0("Rdata/fig/PO/", panel, "_bTMB_roc.Rdata"))
  ggsave(
    plot = bTMB_p,
    filename = paste0("fig/NPC/", panel, "_bTMB_roc.jpeg"),
    width = 11,
    height = 11,
    units = "cm"
  )
  #### series plots
  library(patchwork)
  plots = purrr::map(
    .x = c("Pre-therapy", "Post-therapy", "Af-Be"),
    .f = function(gg) {
      p_t = roc_df %>%
        dplyr::filter(var != "nSNVs") %>%
        dplyr::filter(grp == gg) %>%
        mutate(group2 = str_replace_all(group, paste0(" in ", gg, " Treat"), "")) %>%
        ggplot(aes(x = 1 - specificity, y = sensitivity)) +
        geom_path(aes(color = group2), size = 1.5) +
        ggsci::scale_color_aaas() +
        geom_segment(aes(
          x = 0,
          y = 0,
          xend = 1,
          yend = 1
        ),
        colour = 'grey',
        linetype = 'dotdash') +
        labs(x = "1-Specificity",
             y = "Sensitivity",
             title = paste0("ROC for ", gg)) +
        my_theme +
        # facet_grid(. ~ grp) +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.position = c(0.7, 0.15),
          legend.title = element_blank(),
          legend.background = element_rect(
            fill = NULL,
            size = 0.5,
            linetype = "solid",
            colour = "black"
          )
        ) +
        guides(color = guide_legend(nrow = 3))
      return(p_t)
    }
  )
  p3 = (plots[[1]]) | (plots[[2]]) | (plots[[3]])
  
  save(p3, file = paste0("Rdata/fig/PO/", panel, "_roc.Rdata"))
  ggsave(
    plot = p3,
    filename = paste0("fig/NPC/", panel, "_roc.jpeg"),
    width = 30,
    height = 12,
    units = "cm"
  )
  
  

}


