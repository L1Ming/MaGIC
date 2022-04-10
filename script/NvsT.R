
# load maf ----------------------------------------------------------------
load(file = 'Rdata/Prostate_all_SNVs.Rdata') # pst_ftd, p_gr2
load(file = 'Rdata/Normal_all_SNVs.Rdata') # nm_ftd, n_gr2
load(file = 'Rdata/exon_lib_.Rdata')

# load panels gr ----------------------------------------------------------

source("script/functions/Map2Selector.R")
selector = Map2Selector(panel = "Our_ori")
Catcher = Map2Selector(panel = "Our")
Newman = Map2Selector(panel = 'Newman')
Burgener = Map2Selector(panel = 'Burgener')
MSK = Map2Selector(panel = 'MSK')
load("panels/est_panels.Rdata", verbose = T) # from panels/established_4panels.R script
# Loading objects:
#   NCC_gr
#   plse_gr
#   F1Dx_gr
#   Gd360_gr



# count variables ---------------------------------------------------------

## define function
mapping = function(se_gr, a_gr, maf){
  overlap = findOverlaps(query = se_gr, subject = a_gr) %>% as.data.frame
  se_maf = maf[overlap$subjectHits,] %>% as.data.frame
  se_maf$Region = overlap$queryHits
  return(se_maf %>% select(Sample=Tumor_Sample_Barcode, VAF, Region))
}

P_Map = function(se_gr){
  mapping(se_gr, p_gr2, pst_ftd)
}

N_Map = function(se_gr){
  mapping(se_gr, n_gr2, nm_ftd)
}  

# calculate variables -----------------------------------------------------
# # test
# pnl = 'NCC_gr'
# data = N_maf
# RegionScale = 'NCC150'
# Group = 'Healthy'
# samples = N_sap

counting = function(data, panel_len, samples,
                    RegionScale, Group){
  
  samples = unique(samples) %>% sort
  
  a_df = data %>%   
    group_by(Sample, Region) %>% 
    summarise(mean = mean(VAF), len = length(VAF)) %>% 
    ungroup() %>% 
    dplyr::filter(!(mean < 1 & len == 1)) %>% 
    group_by(Sample) %>% 
    summarise(KMR = length(Region) / panel_len) %>%
    ungroup() 
    
  
  b_df = data %>% 
    group_by(Sample) %>% 
    summarise(mVAF = mean(VAF), 
              SNVs = length(Region),
              NMR = length(unique(Region)),
              bTMB = SNVs / panel_len
    )
  r_df = left_join(data.frame(Sample = samples),
                   a_df, by='Sample') %>% 
    left_join(b_df, by="Sample") %>% 
    mutate(RegionScale = RegionScale,
           Group = Group)
  r_df[is.na(r_df)] = 0
  
  return(r_df)
}



pheno_df = data.table::fread("D:/Data/NormalcfDNA/SraRunTable.txt") %>% 
  dplyr::filter(source_name != "Unhealthy centenarian_blood plasma") %>% 
  select(Sample = Run, Pheno = AGE) %>% 
  # mutate(Age = str_split(string = Group, pattern = " ", simplify = T)[, 1] %>% as.numeric())
  rbind(data.frame(Sample = paste0("P", 1:23),
                   Pheno = "Cancer"))
panels = c("exon_lib_", 
           "Newman", "Burgener", "NCC_gr", "plse_gr", "Gd360_gr", "F1Dx_gr", "selector",
           "Catcher")
pal_names = c("WES", "Newmans'", "Burgeners'", "NCC150",
              "Plasma SELECT", "Grand360", "F1CDx", 
              'MaGIC-Enricher', "MaGIC-Suppressor")

N_sap = pheno_df$Sample[pheno_df$Pheno != 'Cancer']
P_sap = pheno_df$Sample[pheno_df$Pheno == 'Cancer']

res_df = purrr::map2_df(.x=panels, .y = pal_names, .f = function(pnl, name){
  panel = get(pnl)
  panel_len = sum(width(panel)) / 1e6
  N_maf = N_Map(panel)
  P_maf = P_Map(panel)
  
  N_res = counting(N_maf, panel_len = panel_len, 
                   samples = N_sap, RegionScale = name, Group = "Healthy")
  P_res = counting(P_maf, panel_len = panel_len, 
                   samples = P_sap, RegionScale = name, Group = "Cancer")
  res = rbind(N_res, P_res)
})

save(res_df, file = 'Rdata/NvsProstat_ResDF.Rdata')

# Plot --------------------------------------------------------------------

library(ggplot2)
library(ggpubr)
options(scipen = 100)
 # from NMRs/countMutsInExon.R
load('Rdata/NvsProstat_ResDF.Rdata', verbose = T)
source('script/Utils.R')
## load group data


# add ROC plot and boxplot for Enricher-v1 and bTMB -----------------------

res_df1 = res_df %>% 
  dplyr::filter(RegionScale == "MaGIC-Enricher")

## boxplot
ntBox = ggplot(res_df1, aes(x=Group, y=bTMB, fill=Group))+
  geom_boxplot(outlier.colour = 'grey', width = 0.8)+
  ggsci::scale_fill_aaas(alpha = 0.8,
                         limits = c('Healthy', 'Cancer'))+
  ggpubr::stat_compare_means(label.x = 1.25)+
  my_theme+
  theme(legend.position = 'none')
ggsave(ntBox, filename = 'fig/NvsT/Enricher-Box.tiff', units = 'cm', height = 10, width = 11)
save(ntBox, file = 'Rdata/fig/PO/NvsT_Enricher-Box.Rdata')


## ROC plot
library(pROC)
roc1 = roc(response = res_df1$Group,
           predictor = res_df1$bTMB %>% unlist)

best = coords(roc1, 'best')[1,] %>% unlist
roc_df = coords(roc1, transpose = FALSE)

ntROC = ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
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

ggsave(ntROC, filename = 'fig/NvsT/Enricher-ROC.tiff', units = 'cm', height = 10, width = 10)
save(ntROC, file = 'Rdata/fig/PO/NvsT_Enricher-ROC.Rdata')

# other -------------------------------------------------------------------



res_df2 = res_df %>% 
  dplyr::filter(RegionScale %in% c("WES", 'MaGIC-Enricher', "MaGIC-Suppressor")) %>% 
  mutate(RegionScale = factor(RegionScale,
                              levels = c("WES", 'MaGIC-Enricher', "MaGIC-Suppressor"),
                              ordered = T),
         Group = factor(Group, levels = c("Healthy", "Cancer"), ordered = T)) %>% 
  left_join(pheno_df, by='Sample') # %>%
  # mutate(Pheno = factor(Pheno, levels = c("Cancer", paste0(c(25, 70, 100), ' years old')),
  #                       ordered = T))
  
rs = c('WES', "Enricher", "Suppressor")
names(rs) = c("WES", 'MaGIC-Enricher', "MaGIC-Suppressor")
## res DF
options(scipen = 0)
p1 = res_df2 %>%
  select(Sample, Group, bTMB, KMR, mVAF, RegionScale) %>%
  reshape2::melt(id.var = c("Sample", "Group", "RegionScale")) %>%
  mutate(variable = factor(variable, levels = c('mVAF', 'bTMB', 'KMR'))) %>% 
  ggplot(aes(x=Group, y=value, fill=Group))+
  geom_boxplot(alpha=0.8, outlier.colour = 'grey')+
  ggsci::scale_fill_aaas()+
  labs(x=NULL, y=NULL) +
  facet_grid(variable~RegionScale, scales = 'free', 
             labeller = labeller(RegionScale = rs))+
  scale_y_log10(expand = c(0, 0))+
  ggpubr::stat_compare_means(label.y.npc = 0.46, size = 3) +
  my_theme_strip +
  theme(legend.position = 'none', 
        text = element_text(size = 15))
save(p1, file = 'Rdata/fig/PO/NvsT.Rdata')
ggsave(p1, filename = 'fig/NvsT/WES12.tiff', width = 20, height = 20, units = 'cm')

# 
# ## WES, Selector and Catcher  
#   
# for(var in c("KMD", "bTMB", "mVAF")){
#   var_name = ifelse(var == 'KMD', "KMR per Mb", var) 
#   
#   p = res_df %>% 
#     # dplyr::filter(RegionScale %in% c("WES", "Catcher")) %>%
#     select(RegionScale, V=all_of(var), Group) %>% 
#     ggplot(mapping = aes(x=RegionScale, y=V, fill=Group))+
#     geom_boxplot(alpha=0.8, outlier.colour = 'grey')+
#     labs(x=NULL, fill=NULL, y=var_name) +
#     stat_compare_means(
#       # label.y.npc = 0.95,
#       # ref.group = 'WES',
#       label = 'p.signif',
#     ) +
#     ggsci::scale_fill_nejm() +
#     my_theme #+
#     # theme(legend.position = 'none',
#           # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
#   ggsave(p, filename = paste0('fig/NvsT/var of', var, '.jpeg'))
# }
# ## WES to other panels
# 
# for(var in c("KMD", "bTMB", "mVAF")){
#   var_name = ifelse(var == 'KMD', "KMR per Mb", var) 
#   
#   for(grp in c("Cancer", "Healthy")){
#     p = res_df %>% 
#       dplyr::filter(
#         # RegionScale %in% c("WES", 'Selector', 'Catcher'),
#                     Group == grp) %>%
#       select(RegionScale, V=all_of(var)) %>% 
#       mutate(Class = if_else(RegionScale == "WES", 'WES',
#                              if_else(RegionScale %in% c("Catcher", "Selector"), "Catcher",
#                                      "Other Panels"))) %>% 
#       ggplot(mapping = aes(x=RegionScale, y=V, fill=Class))+
#       geom_boxplot(alpha=0.8, outlier.colour = 'grey')+
#       labs(x=NULL, fill=NULL, y=var_name) +
#       stat_compare_means(
#         label.y.npc = 0.95,
#         ref.group = 'WES',
#         label = 'p.signif',
#       ) +
#       ggsci::scale_fill_nejm()+
#       my_theme+
#       theme(legend.position = 'none',
#             axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
#     
#     ggsave(p, filename = paste0('fig/NvsT/',grp, '_', var, '.jpeg'))
#   }
# }
# 
# 
# # ## WES among different age groups
# # for(var in c('SNVs', 'KMR', 'KMD', 'bTMB')){
# #   var_name = ifelse(var == 'KMD', "KMR per Mb", var)
# #   for(pnl in c("WES", "Selector", "Catcher")){
# # 
# #     p = res_df %>%
# #       dplyr::filter(RegionScale == pnl) %>%
# #       select(Pheno, Group, V = all_of(var)) %>%
# #       mutate(Group = factor(Group, levels = c("Healthy", "Cancer"))) %>%
# #       ggplot(aes(x = Pheno, y = V, fill = Group)) +
# #       geom_boxplot(alpha = 0.8) +
# #       geom_point(shape = 21, color='white', fill='green') +
# #       # scale_y_log10()+
# #       labs(x = NULL, fill = NULL, y = var) +
# #       stat_compare_means(
# #         # aes(group = Group),
# #         # hide.ns = T,
# #         label.y.npc = 0.95,
# #         label = 'p.signif',
# #         ref.group = 'Cancer',
# #         # method.args = list(alternative = 'greater')
# #       ) +
# #       ggthemes::scale_fill_wsj() +
# #       my_theme
# #     # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# #     ggsave(p, filename = paste0("fig/NvsT/panels/", var, " in ", pnl, ".jpeg"))
# #   }
# # }
# 
# 
# ## different variables ----
# 
# for(var in c('bTMB', 'KMD', 'mVAF')){
#   var_name = ifelse(var == 'KMD', "KMR per Mb", var)
# 
#   p = res_df %>%
#     dplyr::filter(!(RegionScale %in% c("Newmans'", "Burgeners'", 'Selector'))) %>%
#     select(RegionScale, Group, V = all_of(var)) %>%
#     mutate(Group = factor(Group, levels = c("Healthy", "Cancer"))) %>%
#     ggplot(aes(x=RegionScale, y=V, fill=Group))+
#     geom_boxplot(alpha=0.8)+
#     # scale_y_log10()+
#     labs(x=NULL, fill=NULL, y=var_name) +
#     stat_compare_means(
#       aes(group = Group),hide.ns = T, label.y.npc = 0.95,
#       label = 'p.signif',
#       method.args = list(alternative = 'greater')) +
#     ggthemes::scale_fill_wsj()+
#     my_theme
#   ggsave(p, filename = paste0("fig/NvsT/", var, ".jpeg"))
# 
# }
# 
# var = 'KMD'
# res_df %>%
#   dplyr::filter(!(RegionScale %in% c("Newmans'", "Burgeners'", 'Selector'))) %>%
#   select(RegionScale, Group, V = all_of(var)) %>%
#   compare_means(formula = V~Group, group.by = c("RegionScale"), alternative = 'greater') %>%
#   as.data.frame()
# 
# 
# ## SNVs
# p4 = res_df %>%
#   dplyr::filter(!(RegionScale %in% c("Newmans'", "Burgeners'", 'Selector'))) %>%
#   ggplot(aes(x=RegionScale, y=SNVs, fill=Group))+
#   geom_boxplot(alpha=0.8)+
#   scale_y_log10()+
#   labs(x=NULL, fill=NULL) +
#   stat_compare_means(label = 'p.signif') +
#   ggthemes::scale_fill_wsj()+
#   my_theme
# ggsave(p4, filename = 'fig/NvsT/SNVs.jpeg')
# 
# p5 = res_df %>%
#   dplyr::filter(!(RegionScale %in% c("Newmans'", "Burgeners'", 'Selector'))) %>%
#   ggplot(aes(x=RegionScale, y=KMR, fill=Group))+
#   geom_boxplot(alpha=0.8)+
#   scale_y_log10()+
#   labs(x=NULL, fill=NULL) +
#   stat_compare_means(label = 'p.signif') +
#   ggthemes::scale_fill_wsj()+
#   my_theme
# ggsave(p5, filename = 'fig/NvsT/KMR.jpeg')
# 
# 
# ##
# p6 = res_df %>%
#   dplyr::filter(RegionScale != "Initial Catcher") %>%
#   select(Sample, KMD, RegionScale, Group) %>%
#   # reshape2::dcast(formula =Sample ~ RegionScale, value.var = "KMD") %>%
#   ggplot(aes(x=RegionScale, y=KMD))+
#   geom_violin(aes(fill=RegionScale))+
#   geom_boxplot(width = 0.1, fill='white') +
#   ggsci::scale_fill_aaas(alpha = 0.6)+
#   facet_grid(Group~.)+
#   # labs(x=NULL, y=) +
#   geom_point(shape=21, color = 'white', fill="#064635", size = 2)+
#   geom_line(aes(group=Sample), color = '#064663', size = 0.3)+
#   stat_compare_means(label.y = 450) +
#   expand_limits(y=500) +
#   my_theme+
#   theme(legend.position = 'none')
# ggsave(p6, filename = "fig/NvsT/KMD_group.jpeg")
# 
# 
# p7 = res_df %>%
#   dplyr::filter(RegionScale != "Initial Catcher") %>%
#   select(Sample, bTMB, RegionScale, Group) %>%
#   ggplot(aes(x=RegionScale, y=bTMB))+
#   geom_violin(aes(fill=RegionScale))+
#   geom_boxplot(width = 0.1, fill='white') +
#   ggsci::scale_fill_aaas(alpha = 0.6)+
#   facet_grid(Group~.)+
#   geom_point(shape=21, color = 'white', fill="#064635", size = 2)+
#   geom_line(aes(group=Sample), color = '#064663', size = 0.3)+
#   stat_compare_means(label.y=800)+
#   expand_limits(y=950) +
#   my_theme+
#   theme(legend.position = 'none')
# ggsave(p7, filename = "fig/NvsT/bTMB_group.jpeg")
