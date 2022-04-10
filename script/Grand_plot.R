library(ggplot2)
library(magrittr)
library(stringr)
library(patchwork)
library(ggridges)
library(RColorBrewer)
library(GenomicRanges)
library(dplyr)

## edited at 21-09-27, as all plots parameters had been chosen


# Mon Sep 27 09:43:52 2021 ------------------------------
#######################
#     Capture rate    #
#######################

# capture rate over all ---------------------------------------------------
source("script/functions/Map2Selector.R")
source('script/Utils.R')

# # test
# type = "TCGA"
# panel = "Our_ori"
for(type in c("TCGA")){
  for(panel in c("Our_ori", "Our", "Newman", "Burgener", "MSK")){
    se_len = Map2Selector(panel = panel) %>% width() %>% sum()
    
    load(paste0("Rdata/", panel, "_panel_", type, "_capRate.Rdata")) #resTemp, came from CaptureRatePlot.R
    load("Rdata/project_details.df.Rdata")
    
    resTemp2 = left_join(resTemp,
                         pro_df %>% select(Project, Organ = Primary.Site, World),
                         by = "Project") %>%
      dplyr::mutate(TMB = SNVs / se_len * 1e6, 
                    KMR = NumMutRegion / se_len * 1e6)
    
    
    resInfo = resTemp2 %>% group_by(Organ) %>%
      summarise(
        cap1_per = round(sum(NumMutRegion >= 1) / length(NumMutRegion) * 100, digits = 2),
        cap2_per = round(sum(NumMutRegion >= 2) / length(NumMutRegion) * 100, digits = 2),
        cap3_per = round(sum(NumMutRegion >= 3) / length(NumMutRegion) * 100, digits = 2),
        cap5_per = round(sum(NumMutRegion >= 5) / length(NumMutRegion) * 100, digits = 2),
        World = unique(World)
      ) %>%
      ungroup()
    
    ## subplot for capture percentage
    cMelt = resInfo %>%
      select(Organ, cap1_per, cap2_per, cap3_per, World) %>%
      mutate(rank = (0.5 * cap1_per + 0.3 * cap2_per + 0.2 * cap3_per)) %>% 
      arrange(desc(rank)) %>%
      select(-rank) %>% 
      mutate(Organ = factor(Organ, levels = rev(Organ), ordered = T)) %>%
      reshape2::melt(
        id.variable = c("Organ", "World"),
        variable.name = "CapTime",
        value.name = "Percentage"
      )
    
    # bar_color = c("#F0A500", "#1DB9C3")
    bar_color = ggsci::pal_d3()(2) %>% rev() #  %>% scales::show_col()

    p1_1 = ggplot(data = cMelt, aes(x = Organ, y = Percentage)) +
      annotate('rect', xmin = -Inf, xmax = Inf, ymin = 90, ymax = 100,
               fill = 'grey80') +
      geom_line(aes(group = CapTime, color = CapTime), linetype = 'dashed') +
      geom_point(aes(fill = World, shape = World)) +
      scale_shape_manual(
        limits = c("Y", "N"),
        labels = c("Yes", "No"),
        values = c(24, 21)
      ) +
      scale_fill_manual(
        limits = c("Y", "N"),
        labels = c("Yes", "No"),
        values = bar_color
      ) +
      scale_color_manual(
        limits = c("cap1_per", "cap2_per", "cap3_per"),
        values = c("#F85F73", "#1FAB89", "#3490DE"),
        labels = c("≥1 times", "≥2 times", "≥3 times")
      ) +
      scale_y_continuous(breaks = c(seq(0, 100, 20), 90)) +
      labs(x = NULL, y = "Captured Percentage(%)", fill = "World Top10",
           shape = 'World Top10',
           color = "Captured") +
      my_theme +
      theme(
        legend.background = element_rect(color = 'grey', fill = 'white'),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_text(hjust = 0.5)
      ) 
    
    save(p1_1, file = paste0("Rdata/fig/MC/", panel, "Organ_capRate_", type, ".Rdata"))
    ggsave(plot = p1_1,
           filename = paste0("fig/capRate/", panel, "Organ_capRate_", type, ".jpeg"))
    
    ### capture rate
    cc = cMelt %>% dplyr::filter(CapTime == 'cap1_per') %>%
      mutate(x_lab=rev(1:nrow(.))) 
    
    caRate = ggplot(cc, aes(x=x_lab, y=Percentage))+
      geom_area(fill="grey", alpha=0.5)+
      geom_line()+
      geom_point(aes(fill=World), shape=21, color = 'white', size= 3)+
      scale_fill_manual(
        limits = c("Y", "N"),
        labels = c("Yes", "No"),
        values = bar_color
      ) +
      labs(y="Capture Rate",x=NULL, fill="World Top10")+
      scale_x_continuous(breaks = cc$x_lab, 
                         labels = cc$Organ,
                         expand = c(0, 0.5))+
      coord_flip()+
      my_theme
    ggsave(caRate, filename = paste0("fig/capRate/", panel, "Organ_capRate_", type, "_only.jpeg"), dpi=600)
    save(caRate, file = paste0("Rdata/fig/MC/", panel, "Organ_capRate_", type, "_only.Rdata"))
    
    ### count total cancer numbers in capt times
    
    bMelt = cMelt %>% as_tibble() %>%
      group_by(CapTime, World) %>%
      summarise(
        "≥20%" = sum(Percentage >= 20),
        "≥40%" = sum(Percentage >= 40),
        "≥60%" = sum(Percentage >= 60),
        # "≥70%" = sum(Percentage >= 70 ),
        "≥80%" = sum(Percentage >= 80)
      ) %>%
      ungroup() %>%
      reshape2::melt(
        id.var = c("CapTime", "World"),
        value.name = 'Counts',
        variable.name = 'Percent'
      ) %>%
      dplyr::mutate(Percent = factor(Percent, levels = paste0("≥", seq(80, 20, -20), "%"), ordered = T))
    
    p1_2 = ggplot(data = bMelt, aes(x = CapTime, y = Counts, fill = World)) +
      geom_bar(stat = 'identity',
               width = 0.5,
               color = 'black') +
      # geom_text(aes(label = Counts, color=Percent), vjust=-0.5)+
      facet_grid(Percent ~ .) +
      scale_y_continuous(expand = c(0, 0)) +
      expand_limits(y = 28) +
      scale_x_discrete(limits = paste0("cap", c(1:3, 5), "_per"),
                       labels = paste0("≥", c(1:3, 5), " times")) +
      scale_fill_manual(limits = c("Y", "N"),
                            values = bar_color,
                            labels = c("Yes", "No")) +
      # ggsci::scale_color_npg()+
      labs(x = NULL, y = "Number of cancer types", fill = "World Top10") +
      my_theme +
      theme(legend.position = 'top')
    
    save(p1_2,
         file = paste0("Rdata/fig/MC/", panel, "Organ_capRate_bar_", type, ".Rdata"))
    ggsave(
      plot = p1_2,
      filename = paste0("fig/capRate/", panel, "Organ_capRate_bar_", type, ".jpeg")
    )
    
    # distribution for KMR --------------------------------------------
    ## boxplot
    resTemp2$Organ = factor(resTemp2$Organ,
                            levels = levels(cMelt$Organ),
                            ordered = T)
    #### KMR
    p2_1 = ggplot(data = resTemp2,
                  mapping = aes(y = KMR + 1, x = Organ)) +
      geom_boxplot(aes(fill = World), outlier.colour = 'grey') +
      scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50, 100),
                    labels = c(0, 2, 5, 10, 20, 50, 100)) +
      scale_fill_manual(
        limits = c("Y", "N"),
        labels = c("Yes", "No"),
        values = bar_color
      ) +
      labs(x = NULL, y = "KMR", fill = "World Top10") +
      my_theme +
      theme(
        legend.position = c(0.15, 0.8),
        legend.background = element_rect(colour = 'grey', fill = 'white'),
        axis.text.x = element_text(
          vjust = 0.5,
          hjust = 1,
          angle = 90
        ),
        legend.title = element_text(hjust = 0.5)
      )
    
    save(p2_1,
         file = paste0("Rdata/fig/MC/", panel, "Organ_boxplot_KMR_", type, ".Rdata"))
    ggsave(
      plot = p2_1,
      filename = paste0("fig/capRate/", panel, "Organ_boxplot_KMR_", type, ".jpeg")
    )
    
    ### mVAF
    p2_2 = ggplot(data = resTemp2,
                  mapping = aes(y = mVAF, x = Organ)) +
      geom_boxplot(aes(fill = World), outlier.colour = 'grey') +
      # scale_y_log10(breaks = c(1,2, 5, 10, 20, 50, 100),
      #               labels = c(0,2, 5, 10, 20, 50, 100))+
      scale_fill_manual(
        limits = c("Y", "N"),
        labels = c("Yes", "No"),
        values = bar_color
      ) +
      labs(x = NULL, y = "mVAF", fill = "World Top10") +
      my_theme +
      theme(
        legend.position = c(0.15, 0.85),
        legend.background = element_rect(colour = 'grey', fill = 'white'),
        axis.text.x = element_text(
          vjust = 0.5,
          hjust = 1,
          angle = 90
        ),
        legend.title = element_text(hjust = 0.5)
      )
    
    save(p2_2,
         file = paste0("Rdata/fig/MC/", panel, "Organ_boxplot_mVAF_", type, ".Rdata"))
    ggsave(
      plot = p2_2,
      filename = paste0("fig/capRate/", panel, "Organ_boxplot_mVAF_", type, ".jpeg")
    )
    
    ### TMB
    p2_3 = ggplot(data = resTemp2,
                  mapping = aes(y = TMB + 1, x = Organ)) +
      geom_boxplot(aes(fill = World), outlier.colour = 'grey') +
      scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50, 100),
                    labels = c(0, 2, 5, 10, 20, 50, 100)) +
      scale_fill_manual(
        limits = c("Y", "N"),
        labels = c("Yes", "No"),
        values = bar_color
      ) +
      labs(x = NULL, y = "TMB", fill = "World Top10") +
      my_theme +
      theme(
        legend.position = c(0.15, 0.8),
        legend.background = element_rect(colour = 'grey', fill = 'white'),
        axis.text.x = element_text(
          vjust = 0.5,
          hjust = 1,
          angle = 90
        ),
        legend.title = element_text(hjust = 0.5)
      )
    
    save(p2_3,
         file = paste0("Rdata/fig/MC/", panel, "Organ_boxplot_TMB_", type, ".Rdata"))
    ggsave(
      plot = p2_3,
      filename = paste0("fig/capRate/", panel, "Organ_boxplot_TMB_", type, ".jpeg")
    )
    
    
    # CDF plot for distribution -----------------------------------------------
    
    p3 = ggplot(data = resTemp2,
                mapping = aes(x = KMR + 1, group = Organ, color =
                                World)) +
      # stat_ecdf(geom = 'step')+
      stat_ecdf(geom = 'smooth') +
      scale_x_log10() +
      scale_color_manual(
        limits = c("Y", "N"),
        labels = c("Yes", "No"),
        values = bar_color
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      expand_limits(y = 1.05) +
      labs(x = "KMR",
           y = "Empirical cumulative distribution function",
           color = "World Top10") +
      # facet_grid(Project~.)+
      my_theme +
      theme(
        legend.position = c(0.8, 0.2),
        legend.background = element_rect(fill = 'white', color = 'grey'),
        axis.text.y.right = element_text()
      )
    save(p3, file = paste0("Rdata/fig/MC/", panel, "Organ_CDF_", type, ".Rdata"))
    ggsave(plot = p3,
           filename = paste0("fig/capRate/", panel, "Organ_CDF_", type, ".jpeg"))
    
    
    
    # mean barplot of 3 variables --------------------------------------
    mean_df = resTemp2 %>% group_by(Organ) %>%
      summarise(
        mean.KMR = mean(KMR),
        mean.mVAF = mean(mVAF),
        mean.TMB = mean(TMB),
        World = unique(World)
      )
    # use `fivenum` to count quartiles as trd
    trd_df = purrr::map_df(mean_df %>% select(-Organ, -World), fivenum) %>%
      ceiling() %>%
      dplyr::mutate(Group = paste0('Fivenum', 1:nrow(.)),
                    Trd = if_else(Group %in% c("Fivenum3", "Fivenum4"), "Yes", "No")) %>%
      reshape2::melt()
    save(trd_df, file = paste0("Rdata/", panel, "_ThreeVariables_OS_trd_df.Rdata"))
    ## trd plot
    p4_0 = trd_df %>%
      ggplot(aes(
        x = variable,
        y = value,
        fill = Trd,
        group = Group
      )) +
      geom_bar(color = 'black',
               alpha = 0.8,
               stat = 'identity',
               position = 'dodge') +
      scale_y_continuous(expand = c(0, 0)) +
      expand_limits(y = max(trd_df$value)*1.1) +
      scale_x_discrete(limits = paste0('mean.', c("KMR", 'TMB', 'mVAF')),
                       labels = paste0(c("KMR", 'TMB', 'mVAF'))) +
      labs(x = "Mean of variables for each TCGA datasets",
           y = "Celi values of Tukey's five number",
           group = 'Used as\nThreshold',
           fill = 'Used as\nThreshold') +
      geom_text(
        aes(label = value),
        size = 4,
        vjust = -0.5,
        position = position_dodge(0.9)
      ) +
      scale_fill_manual(limits = c("Yes", "No"),
                        values = c("#FF4848", "#316B83")) + #,
      my_theme +
      theme(
        legend.position = c(0.15, 0.8),
        legend.title = element_text(hjust = 0.5),
        legend.background = element_rect(color = 'grey', fill = 'white')
      )
    save(p4_0, file = paste0("Rdata/fig/MC/", panel, "quarters_trd_", type, ".Rdata"))
    ggsave(plot = p4_0,
           filename = paste0("fig/capRate/", panel, "quarters_trd_", type, ".jpeg"))
    
    
    ## mean KMR barplot
    KMR_trd = trd_df$value[trd_df$Trd == "Yes" &
                             trd_df$variable == 'mean.KMR']
    p4_1 = ggplot(mean_df, aes(x = reorder(Organ, mean.KMR), y = mean.KMR, fill = World)) +
      geom_bar(stat = 'identity') +
      geom_hline(yintercept = KMR_trd,
                 color = '#334756',
                 linetype = 'dashed') +
      annotate(
        'text',
        x = 2.5,
        y = KMR_trd * 1.05,
        label = KMR_trd,
        color = '#082032'
      ) +
      scale_fill_manual(
        limits = c("Y", "N"),
        labels = c("Yes", "No"),
        values = bar_color
      ) +
      labs(y = "Mean KMR", x = "Cancer types", fill = 'Wrold Top10') +
      my_theme +
      theme(
        legend.position = c(0.15, 0.8),
        legend.background = element_rect(colour = 'grey', fill = 'white'),
        axis.text.x = element_text(
          vjust = 0.5,
          hjust = 1,
          angle = 90
        ),
        legend.title = element_text(hjust = 0.5)
      )
    save(p4_1, file = paste0("Rdata/fig/MC/", panel, "meanKMR_", type, ".Rdata"))
    ggsave(plot = p4_1,
           filename = paste0("fig/capRate/", panel, "meanKMR_", type, ".jpeg"))
    
    ## mean mVAF: barplot
    mVAF_trd = trd_df$value[trd_df$Trd == "Yes" &
                              trd_df$variable == 'mean.mVAF']
    p4_2 = ggplot(mean_df, aes(x = reorder(Organ, mean.mVAF), y = mean.mVAF, fill = World)) +
      geom_bar(stat = 'identity') +
      geom_hline(yintercept = mVAF_trd,
                 color = '#334756',
                 linetype = 'dashed') +
      annotate(
        'text',
        x = 2.5,
        y = mVAF_trd * 1.05,
        label = mVAF_trd,
        color = '#082032'
      ) +
      scale_fill_manual(
        limits = c("Y", "N"),
        labels = c("Yes", "No"),
        values = bar_color
      ) +
      labs(y = "Mean mVAF", x = "Cancer types", fill = 'Wrold Top10') +
      my_theme +
      theme(
        legend.position = c(0.15, 0.8),
        legend.background = element_rect(colour = 'grey', fill = 'white'),
        axis.text.x = element_text(
          vjust = 0.5,
          hjust = 1,
          angle = 90
        ),
        legend.title = element_text(hjust = 0.5)
      )
    save(p4_2, file = paste0("Rdata/fig/MC/", panel, "meanmVAF_", type, ".Rdata"))
    ggsave(plot = p4_2,
           filename = paste0("fig/capRate/", panel, "meanmVAF_", type, ".jpeg"))
    
    ## mean TMB: barplot
    TMB_trd = trd_df$value[trd_df$Trd == "Yes" &
                             trd_df$variable == 'mean.TMB']
    p4_3 = ggplot(mean_df, aes(x = reorder(Organ, mean.TMB), y = mean.TMB, fill = World)) +
      geom_bar(stat = 'identity') +
      geom_hline(yintercept = TMB_trd,
                 color = '#334756',
                 linetype = 'dashed') +
      annotate(
        'text',
        x = 2.5,
        y = TMB_trd * 1.05,
        label = TMB_trd,
        color = '#082032'
      ) +
      scale_fill_manual(
        limits = c("Y", "N"),
        labels = c("Yes", "No"),
        values = bar_color
      ) +
      labs(y = "Mean TMB", x = "Cancer types", fill = 'Wrold Top10') +
      my_theme +
      theme(
        legend.position = c(0.15, 0.8),
        legend.background = element_rect(colour = 'grey', fill = 'white'),
        axis.text.x = element_text(
          vjust = 0.5,
          hjust = 1,
          angle = 90
        ),
        legend.title = element_text(hjust = 0.5)
      )
    save(p4_3, file = paste0("Rdata/fig/MC/", panel, "meanTMB_", type, ".Rdata"))
    ggsave(plot = p4_3,
           filename = paste0("fig/capRate/", panel, "meanTMB_", type, ".jpeg"))
    
    
  }
}


###############################################
#          for three original plots           #
###############################################
if(1){
  # rm(list = ls())
  panel = "Our_ori"
  load(paste0("Rdata/", panel, "_panel_TCGA_capRate.Rdata"), verbose = T)
  resPart = resTemp %>% dplyr::filter(Project %in% c("HNSC", "DLBC", "SARC"))
  # boxplot for original datasets
  # p5 = ggplot(data = resPart,
  #             mapping = aes(y = KMR+1, x = Project)) +
  #   geom_violin(aes(fill = Project)) +
  #   geom_boxplot(width = 0.1) +
  #   scale_y_log10(breaks = c(1, 3, 10 , 30),
  #                 labels = c(0, 3, 10 , 30))+
  #   labs(x = NULL, y = "KMR") +
  #   my_theme +
  #   theme(legend.position = 'none')
  # save(p5, file = "Rdata/fig/MC/ThreeDatasets_boxplot_.Rdata")
  # ggsave(plot=p5, filename = "fig/capRate/ThreeDatasets_boxplot_.jpeg")
  # 
  # 
  # # Mean VAF for original datasets
  # meanKMR_df = resPart %>% group_by(Project) %>%
  #   summarise(mean.KMR = mean(KMR))
  # p6 = ggplot(meanKMR_df, aes(x = Project, y = mean.KMR)) +
  #   geom_bar(stat = 'identity', width = 0.6, color='grey30', aes(fill = Project)) +
  #   labs(y = "Mean KMR", x = NULL) +
  #   scale_y_continuous(expand = c(0, 0)) +
  #   expand_limits(y = max(meanKMR_df$mean.KMR) + 1) +
  #   my_theme +
  #   theme(legend.position = 'none')
  # save(p6, file = "Rdata/fig/MC/ThreeDatasets_meanKMR_.Rdata")
  # ggsave(plot=p6, filename = "fig/capRate/ThreeDatasets_meanKMR_.jpeg")
  
  
  # capture rate of 1, 2, 3
  resInfo2 = resPart %>% group_by(Project) %>%
    summarise(
      cap1_per = round(sum(NumMutRegion >= 1) / length(NumMutRegion) * 100, digits = 2),
      cap2_per = round(sum(NumMutRegion >= 2) / length(NumMutRegion) * 100, digits = 2),
      cap3_per = round(sum(NumMutRegion >= 3) / length(NumMutRegion) * 100, digits = 2),
      # cap5_per = round(sum(NumMutRegion >= 5) / length(NumMutRegion) * 100, digits = 2),
    ) %>% reshape2::melt(id.var = "Project",
                         variable.name = "CapTimes",
                         value.name = "CapRate")
  
  p7 = ggplot(resInfo2, aes(x = Project, y = CapRate, fill = CapTimes)) +
    geom_bar(position = 'dodge',
             color = 'black',
             stat = 'identity',
             width = 0.6
             ) +
    geom_text(aes(label = round(CapRate, 1)), vjust = -0.5, size = 3,
              position = position_dodge(0.6)) +
    scale_fill_manual(
      limits = c("cap1_per", "cap2_per", "cap3_per"),
      values = c("#F85F73", "#1FAB89", "#3490DE"),
      labels = c("≥1 times", "≥2 times", "≥3 times")
    ) +
    labs(x=NULL, y = "Captured Percentage(%)", fill = "Captured") +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 20)) +
    expand_limits(y = 105) +
    my_theme
  save(p7, file = "Rdata/fig/MC/ThreeDatasets_bar_.Rdata")
  ggsave(plot=p7, filename = "fig/capRate/ThreeDatasets_bar_.jpeg")
  
}
