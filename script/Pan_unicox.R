library(dplyr)
library(stringr)
library(magrittr)
library(ggplot2)
library(survival)
library(survminer)

# main --------------------------------------------------------------------
# ## for test
# type = "TCGA"
# over = 'All'
# og = 'Lung'
# uper = 5
# status = "OS"
# var = "TMB"
# panel = "Our"
if (1) {
  rm(list = ls())
  a = Sys.time()
  
  source("script/functions/Map2Selector.R")
  source("script/Utils.R")
  source("script/functions/coxSurv.R")
  
  load("Rdata/project_details.df.Rdata")
  
  pro_df %<>% mutate(All = "Y")
  se_len = Map2Selector(panel = 'Our_ori') %>% width() %>% sum()
  resCox = NULL
  
  for (type in c("TCGA")) {

      for (panel in c("Our_ori", 'Our', "Newman", "Burgener", "MSK")) {
      # load data
      load(paste0("Rdata/", panel, "_ThreeVariables_OS_trd_df.Rdata"))
      
      load(paste0("Rdata/", panel, "_panel_", type, "_capRate.Rdata"))
      resTemp = resTemp %>%
        dplyr::mutate(TMB = SNVs / se_len * 1e6, 
                      KMR = NumMutRegion / se_len * 1e6)
      
      
      for (var in c("KMR", "TMB", "mVAF")) {
        # set variable to run
        ## get variable trd
        trd_vector = trd_df$value[trd_df$Trd == "Yes" &
                                    trd_df$variable == paste0('mean.', var)] %>%
          sort()
        
        
        for (over in c("All", "World", "Women")) {
          # study fileds indication
          group_df = pro_df[pro_df[, over] == "Y", ] %>%
            dplyr::filter(Source == type) %>%
            dplyr::select(Project = Project, Organ = Primary.Site)
          if (nrow(group_df) == 0)
            next
          dat = inner_join(resTemp, group_df, by = "Project")
          for (og in unique(dat$Organ)) {
            # set organ to search
            for (status in c("OS", "PFS")) {

              # set OS or PFS
              a_time = paste0(status, ".time")
              # set group var
              for (uper in trd_vector) {
                og_dat = dat %>% dplyr::filter(Organ == og)
                og_dat$Group = ifelse(og_dat[, var] > uper,
                                      "HRG",
                                      ifelse(og_dat[, var] == 0, "LRG", "MRG"))
                
                og_dat$Group = factor(
                  og_dat$Group,
                  levels = c("LRG", "MRG", "HRG"),
                  ordered = T
                )
                
                # make sure of multilevels
                if (max(og_dat[, var]) < uper)
                  next
                
                # doing cox
                try({
                  s = coxSurv(data = og_dat,
                              time = a_time,
                              status = status)
                  d = pairwise_survdiff(s[[1]]$formula,
                                        data = og_dat,
                                        p.adjust.method = "fdr")
                  t = d$p.value %>% as.data.frame() %>%
                    dplyr::mutate(Query = rownames(.)) %>%
                    reshape2::melt(
                      id.var = 'Query',
                      variable.name = 'Subject',
                      value.name = "Pvalue"
                    ) %>%
                    na.omit() %>%
                    transmute(
                      Panel = panel,
                      Over = over,
                      Source = type,
                      Variable = var,
                      Organ = og,
                      Status = status,
                      Compare = paste0(Query, " vs " , Subject),
                      Threshold = paste0("0~", uper),
                      Tlabel = if_else(uper == max(trd_vector),
                                           "UperQuertile", "Median"),
                      Pvalue_ori = Pvalue,
                      Pvalue = Pvalue %>% round(3),
                      Plabel = if_else(Pvalue == 0, "P<0.001", paste0("P=", Pvalue))
                    )
                  
                  # plot
                  file_path = paste0(
                    "fig/OS/curve/",
                    panel,
                    "/",
                    type,
                    "_",
                    status,
                    "/",
                    over,
                    "/",
                    var,
                    "/"
                  )
                  if (!dir.exists(file_path))
                    dir.create(file_path, recursive = T)
                  
                  sv = survfit(s[[1]]$formula, data = og_dat)
                  colors = RColorBrewer::brewer.pal(6, 'Set1')

                  p0 = ggsurvplot(sv, conf.int = T, alpha=0.8)
                  # edit
                  p_anno = paste(paste0(t$Compare, ": ", t$Plabel), collapse = "\n")
                  p_x = max(p0$data.survtable$time) * 0.8
                  
                  pa = p0$plot +
                    annotate(
                      "text",
                      color = '#082032',
                      size = 4,
                      family = 'serif',
                      x = p_x,
                      y = 0.85,
                      label = p_anno
                    ) +
                    scale_color_manual(limits = paste0("Group=", c("LRG", "MRG", "HRG")),
                                       values = colors[c(1, 2, 5)],
                                       labels = c("LRG", "MRG", "HRG"))+
                    scale_fill_manual(limits = paste0("Group=", c("LRG", "MRG", "HRG")),
                                      values = colors[c(1, 2, 5)],
                                      labels = c("LRG", "MRG", "HRG")) +
                    guides(fill='none') +
                    labs(color = NULL, 
                         title = paste0(status, " of ",  tolower(og), " cancer (", var, ": 0~", uper, ")"))+
                    theme(plot.title = element_text(hjust=0.5)) +
                    my_theme
                  
                  # save
                  save(
                    pa,
                    file = paste0(
                      "Rdata/fig/OS/", panel, "_", var, "_", type, "_", over, "_", og, "_0vs", uper, "_", status, ".Rdata"
                    )
                  )
                  ggsave(pa,
                         filename = paste0(file_path, og, "_0vs", uper, ".jpeg"))
                  
                  # concat
                  resCox = rbind(resCox, t)
                })
              }
            }
          }
        }
      }
    }
  }
  rownames(resCox) = 1:nrow(resCox)
  save(resCox, file = "Rdata/Pan_res_cox.Rdata")
  
  b = Sys.time() 
  print(b - a) # Time difference of 14.82709 mins
}


# tile plot -----------------------------------------------------------------

if (1) {
  # # # for test
  # type = 'TCGA'
  # over = "World"
  # lv = 'LRG vs HRG'
  
  # load data
  rm(list = ls())
  p_cut = 0.01
  source('script/Utils.R')
  load("Rdata/Pan_res_cox.Rdata", verbose = T)
  
  resCox = resCox %>% 
    mutate(Pgroup = if_else(Pvalue_ori < 0.01, 'P < 0.01', 'P ≥ 0.01'),
           Compare = str_replace_all(Compare, "High", "HRG"),
           Compare = str_replace_all(Compare, "Middle", "MRG"),
           Compare = str_replace_all(Compare, "Low", "LRG")
           )
  
  ## total comparison
  bar_panels = resCox %>% 
    dplyr::filter(Over == "World", Variable=="KMR") %>% 
    group_by(Tlabel, Status, Panel, Compare) %>% 
    summarise(Count = sum(Pvalue_ori < p_cut)) %>% 
    ungroup() %>%
    # rbind(tibble(
    #   Tlabel = "Median",
    #   Status = c("OS", "OS", "PFS", "PFS"),
    #   Panel = "Newman",
    #   Compare = rep(c("HRG vs MRG", "MRG vs LRG"), times = 2),
    #   Count = 0,
    # )) %>% 
    mutate(Panel = factor(Panel, levels = c("Our_ori", "Our", "MSK", "Newman", "Burgener"), ordered=T),
           Tlabel = paste0("0 ~ ", Tlabel)) %>% 
    ggplot(aes(x=Compare, y=Count, fill=Panel))+
    geom_bar(width=0.8, stat='identity', position = 'dodge')+
    facet_grid(Tlabel ~ Status)+
    geom_text(aes(label = Count), vjust=-0.5,size=3,  position = position_dodge(0.8))+
    labs(x=NULL, y="Number of significant cancer types", fill='Approaches')+
    scale_y_continuous(breaks = seq(0, 10, 2), expand = c(0, 0))+
    expand_limits(y=11)+
    ggsci::scale_fill_npg(limits = c("Our_ori","Our", "MSK", "Newman", "Burgener"),
                          labels = c("MaGICv1", "MaGICv2", "MSK-IMPACT", "Newman's", "Burgener's"))+
    my_theme +
    theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1),
          strip.background = element_blank())
  ggsave(bar_panels, filename = "fig/OS/Overview/World/zz_Panels.tiff", dpi=600)
  save(bar_panels, file = "Rdata/fig/OS/zz_Panel.Rdata")

  # color = c("#F38181", "#FCD2D1", "#90AACB")
  colors = RColorBrewer::brewer.pal(6, 'Set1')
  # scales::show_col(colors)
  color = colors[c(1, 2)]
  for (type in c("TCGA")) {
    #, "ICGC"
    for (panel in c("Our_ori", "Our", "Newman", "Burgener", "MSK")) {
      load(paste0("Rdata/", panel, "_ThreeVariables_OS_trd_df.Rdata"), verbose = T) # from grand_plot.R
      trd_df2 = trd_df %>% dplyr::filter(Trd == "Yes") %>%
        mutate(TrdLv = if_else(Group == "Fivenum3", 'a1', 'a2'))
      for (lv in c("HRG vs LRG", "MRG vs LRG", "HRG vs MRG")) {
        print(lv)
        
        for (over in c("All", "World")) {

          ov_cox = resCox %>%
            dplyr::filter(Source == type, Over == over, Compare == lv, Panel == panel) %>%
            dplyr::mutate(
              Threshold2 = if_else(
                Tlabel == 'Median',
                # as.integer(str_split(Threshold, "0~", simplify = T)[, 2])
                # %in% trd_df2$value[trd_df2$TrdLv == 'a1'],
                "0 ~ Median",
                "0 ~ Upper Quertile"
              ),
              # Pg = if_else(
              #   Pvalue_ori < 0.001,
              #   "UltSig",
              #   if_else(Pvalue_ori < 0.05, "NrmSig",
              #           "NS")
              # ),
              Variable = factor(
                Variable,
                levels = c("KMR", "TMB", "mVAF"),
                ordered = TRUE
              ),
              Plabel = if_else(Pvalue_ori < 0.01, "P<0.01", as.character(Pvalue))
            )
          
          if (nrow(ov_cox) == 0)
            next
          
          
          ## heatmap
          p1 = ov_cox %>%
            ggplot(aes(x = Variable, y = Organ)) +
            geom_tile(aes(fill = Pgroup), color = '#3F3351', alpha=0.7) +
            geom_text(aes(label = Plabel),
                      color = 'white',
                      size = 3) +
            scale_fill_manual(
              # limits = c('UltSig', "NrmSig", "NS"),
              # labels = c("P<0.001", "P<0.05", "NS"),
              values = color
            ) +
            labs(
              x = NULL,
              y = NULL,
              fill = NULL,
              title = lv
            ) +
            scale_x_discrete(expand = c(0, 0),
                             limits = c("KMR", 'TMB', 'mVAF'),
                             labels = c("KMR", "TMB", "mVAF")) +
            scale_y_discrete(expand = c(0, 0)) +
            facet_grid(Threshold2 ~ Status) +
            my_theme +
            theme(
              strip.background = element_blank(),
              panel.background = element_rect(fill = 'grey'),
              plot.title = element_text(hjust = 0.5)
            )
          if (!dir.exists(paste0('fig/OS/Overview/', over, "/")))
            dir.create(paste0('fig/OS/Overview/', over, "/"), recursive = T)
          save(p1,
               file = paste0("Rdata/fig/OS/",panel, "_", lv, "_", type, "_", over, ".Rdata"))
          
          ggsave(
            p1,
            filename = paste0(
              "fig/OS/Overview/",
              over,
              "/",
              panel, 
              "_",
              lv,
              "_",
              type,
              "_",
              over,
              ".jpeg"
            ),
            width = 18,
            height = 20,
            units = 'cm'
          )
          
          
          
          
          ## barplot
          ors = ifelse(over == 'World', 10, 24)
          bar_dat = ov_cox %>% group_by(Variable, Status, Threshold2, Pgroup) %>%
            summarise(Counts = n(), 
                      Percent = Counts / ors * 100) %>%
            ungroup() %>% 
            mutate(Pgroup = factor(Pgroup, ordered = T,
                                   levels = rev(c("P < 0.01", "P ≥ 0.01"))))
            # mutate(Pg = factor(Pg, levels = rev(c(
            #   'UltSig', "NrmSig", "NS"
            # )), ordered = T))
          
          p3 = ggplot(bar_dat, aes(x = Variable, y = Percent)) +
            scale_y_continuous(expand = c(0, 0),
                               breaks = seq(0, 100, 20)) +
            geom_bar(aes(fill = Pgroup), width=0.5, alpha=0.8,
                     stat = 'identity', position = 'stack') +
            labs(x = NULL,  y = "Cancer types (%)",
                 title = lv,
                 fill = NULL) +
            scale_fill_manual(
              limits = c("P < 0.01", "P ≥ 0.01"),
              # labels = c("P<0.001", "P<0.05", "NS"),
              values = color
            ) +
            scale_x_discrete(limits = c("KMR", 'TMB', 'mVAF'),
                             labels = c("KMR", "TMB", "mVAF")) +
            expand_limits(y = 100) +
            facet_grid(Threshold2 ~ Status) +
            my_theme +
            theme(strip.background = element_blank(),
                  plot.title = element_text(hjust = 0.5))
          
          save(p3,
               file = paste0(
                 "Rdata/fig/OS/",
                 panel, 
                 "_",
                 lv,
                 "_",
                 type,
                 "_",
                 over,
                 "_barplot.Rdata"
               ))
          
          ggsave(
            p3,
            filename = paste0(
              "fig/OS/Overview/",
              over,
              "/",
              panel, 
              "_",
              lv,
              "_",
              type,
              "_",
              over,
              "_barplot.jpeg"
            ),
            width = 18,
            height = 20,
            units = 'cm'
          )
        }
      }
    }
  }
  
  ## tile plot for TMB in Our_ori
  colors = RColorBrewer::brewer.pal(6, 'Set1')
  color = colors[c(1, 2)]
  ov_cox2 = resCox %>%
    dplyr::filter(Source == 'TCGA', Over == "World", Variable == 'TMB', Panel == 'Our_ori') %>%
    dplyr::mutate(
      Threshold2 = if_else(
        Tlabel == 'Median',
        # as.integer(str_split(Threshold, "0~", simplify = T)[, 2])
        # %in% trd_df2$value[trd_df2$TrdLv == 'a1'],
        "Threshold: HRG TMB>5",
        "Threshold: HRG TMB>11"
      ),
      # Pg = if_else(
      #   Pvalue_ori < 0.001,
      #   "UltSig",
      #   if_else(Pvalue_ori < 0.05, "NrmSig",
      #           "NS")
      # ),
      # Variable = factor(
      #   Variable,
      #   levels = c("KMR", "TMB", "mVAF"),
      #   ordered = TRUE
      # ),
      Plabel = if_else(Pvalue_ori < 0.01, "P<0.01", as.character(Pvalue))
    )
  
  ## heatmap
  p12 = ov_cox2 %>%
    ggplot(aes(x = Compare, y = Organ)) +
    geom_tile(aes(fill = Pgroup), color = '#3F3351', alpha=0.7) +
    geom_text(aes(label = Plabel),
              color = 'white',
              size = 3) +
    scale_fill_manual(
      limits = c('P < 0.01', "P ≥ 0.01"),
      values = color
    ) +
    labs(
      x = NULL,
      y = NULL,
      fill = NULL) +
    # scale_x_discrete(expand = c(0, 0),
    #                  limits = c("KMR", 'TMB', 'mVAF'),
    #                  labels = c("KMR", "TMB", "mVAF")) +
    scale_y_discrete(expand = c(0, 0)) +
    facet_grid(Threshold2 ~ Status) +
    my_theme_strip +
    theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1),
          panel.background = element_rect(fill = 'grey'),
          strip.background = element_blank())
  save(p12, file = paste0("Rdata/fig/OS/Our_ori_TMB_tile_world.Rdata"))
  
  ggsave(p12, filename = "fig/OS/Overview/World/Our_ori_TMB_tile_world.tiff", 
         width = 18, height = 20, units = 'cm')
  

  
  ## barplot
  
  bar_dat = ov_cox2 %>% 
    group_by(Compare, Status, Threshold2, Pgroup) %>%
    summarise(Counts = n(), 
              Percent = Counts / 10 * 100) %>%
    ungroup() %>%
    mutate(Pgroup = factor(Pgroup, ordered = T,
                           levels = rev(c("P < 0.01", "P ≥ 0.01"))))
  
  p32 = ggplot(bar_dat, aes(x = Compare, y = Percent)) +
    scale_y_continuous(expand = c(0, 0),
                       breaks = seq(0, 100, 20)) +
    geom_bar(aes(fill = Pgroup), width=0.5, alpha=0.8,
             stat = 'identity', position = 'stack') +
    labs(x = NULL, y = "Cancer types (%)",
         fill = NULL) +
    scale_fill_manual(
      limits = c('P < 0.01', "P ≥ 0.01"),
      values = color
    ) +
    # scale_x_discrete(limits = c("KMR", 'TMB', 'mVAF'),
    #                  labels = c("KMR", "TMB", "mVAF")) +
    expand_limits(y = 100) +
    facet_grid(Threshold2 ~ Status) +
    my_theme_strip +
    theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1),
          strip.background = element_blank())
  
  save(p32, file = "Rdata/fig/OS/Our_ori_TMB_bar_world.Rdata")
  
  ggsave(p32, filename = "fig/OS/Overview/World/Our_ori_TMB_bar_world.tiff",
         width = 18, height = 20, units = 'cm' )
  
}
