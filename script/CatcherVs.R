
#####################################################################
# compare background noise rate among panels ----------------------
#####################################################################

source("script/functions/Map2Selector.R")
source('script/Utils.R')
our_ori = Map2Selector(panel = "Our_ori")
our = Map2Selector(panel = "Our")
newm = Map2Selector(panel = 'Newman')
burg = Map2Selector(panel = 'Burgener')
MSK = Map2Selector(panel = 'MSK')
load("panels/est_panels.Rdata") # from Established_4panels.R script
load("Rdata/exon_lib_.Rdata")
load("Rdata/Normal_all_SNVs.Rdata")
load("Rdata/over_norm.Rdata")

panel_df = map_df(
  .x = list("our_ori", "our", "newm", "burg", "MSK", "NCC_gr", "Gd360_gr", "F1Dx_gr", "plse_gr"),
  .f = function(gr_name) {
    # gr_name='our'
    gr = get(gr_name)
    over = findOverlaps(query = gr, subject = exon_lib_) %>% as.data.frame()
    norm_0 = nm_ftd[over_norm$queryHits,] %>%
      dplyr::mutate(
        exon_number = over_norm$subjectHits,
        width = width(exon_lib_[over_norm$subjectHits]),
        P.WBC = exon_lib_$P.WBC[over_norm$subjectHits],
        Symbol = exon_lib_$Symbol[over_norm$subjectHits],
      )
    
    ### count the average mutation number of each exon
    nmGrp_part = norm_0 %>% as_tibble() %>% 
      group_by(Tumor_Sample_Barcode, exon_number) %>%
      summarise(MutInExon = n(), 
                widthExon = unique(width), 
                P.WBC = unique(P.WBC)) %>%
      ungroup() %>% 
      group_by(exon_number) %>%
      summarise(
        MutsExon = mean(MutInExon),
        P.WBC = unique(P.WBC),
        MeanMuts = mean(MutInExon) / widthExon * 1000
      ) %>% # in Kb
      ungroup() %>% 
      distinct()
    ### add 0
    exNotIn = over$subjectHits[!(over$subjectHits %in% nmGrp_part$exon_number)]
    exNotIn_df = tibble(
      exon_number = exNotIn,
      P.WBC = exon_lib_$P.WBC[exNotIn],
      MutsExon = 0,
      MeanMuts = 0
    )
    
    ### add selector
    nmGrp1 = rbind(nmGrp_part, exNotIn_df) %>% 
      mutate(isSele = if_else(exon_number %in% over$subjectHits, "Yes", "No"))
    
    nmGrp1$panel = gr_name
    return(nmGrp1[nmGrp1$isSele == "Yes",])
  }
)

save(panel_df, file = 'Rdata/panel_df.Rdata')

## plot
if (1) {
  
  # --------------------------------------------------------------
  # load("Rdata/panel_df.Rdata")
  options(scipen = 100)
  h = seq(18, 48, length.out=7)
  p5 = panel_df %>% 
    dplyr::filter(MeanMuts < 100,
                  panel != "our_ori" 
    ) %>% 
    ggplot(aes(x = panel, y = MeanMuts+0.001, fill=panel)) +
    geom_jitter(width = 0.2, alpha=0.5, aes(fill=panel),shape=21)+
    geom_boxplot(fill='white', width=0.2) +
    ggpubr::stat_compare_means(comparisons = list(c('our', 'newm')),
                               mapping = aes(label = ..p.signif..),
                               method.args = list(alternative='less'),
                               label.y=h[1],
                               color=p_color)+
    ggpubr::stat_compare_means(comparisons = list(c('our', 'Gd360_gr')),
                               mapping = aes(label = ..p.signif..),
                               method.args = list(alternative='less'),
                               label.y=h[2],
                               color=p_color)+
    ggpubr::stat_compare_means(comparisons = list(c('our', 'plse_gr')),
                               mapping = aes(label = ..p.signif..),
                               method.args = list(alternative='less'),
                               label.y=h[3],
                               color=p_color)+
    ggpubr::stat_compare_means(comparisons = list(c('our', 'burg')),
                               mapping = aes(label = ..p.signif..),
                               method.args = list(alternative='less'),
                               label.y=h[4],
                               color=p_color)+
    ggpubr::stat_compare_means(comparisons = list(c('our', 'NCC_gr')),
                               mapping = aes(label = ..p.signif..),
                               method.args = list(alternative='less'),
                               label.y=h[5],
                               color=p_color)+
    ggpubr::stat_compare_means(comparisons = list(c('our', 'MSK')),
                               mapping = aes(label = ..p.signif..),
                               method.args = list(alternative='less'),
                               label.y=h[6],
                               color=p_color)+
    ggpubr::stat_compare_means(comparisons = list(c('our', 'F1Dx_gr')),
                               mapping = aes(label = ..p.signif..),
                               method.args = list(alternative='less'),
                               label.y=h[7],
                               color=p_color)+
    scale_fill_manual(limits = c("our", "newm", "Gd360_gr", "plse_gr",
                                 "burg", "NCC_gr", "MSK", "F1Dx_gr"),
                      values = c("#F05454", "#D1E8E4", "#D1E8E4", "#D1E8E4",
                                 "#D1E8E4", "#D1E8E4", "#D1E8E4", "#D1E8E4")) +
    scale_x_discrete(limits = c("our", "newm", "Gd360_gr", "plse_gr",
                                "burg", "NCC_gr", "MSK", "F1Dx_gr"),
                     labels = c("Suppressor", "Newmans'", "Grand360", "Plasma SELECT",
                                "Burgeners'", "NCC150", "MSK-IMPACT","F1CDx")) +
    labs(x = NULL, 
         y = expression(bar(m)~" of cfDNA in healthy donors")) +
    my_theme +
    coord_cartesian(ylim = c(0, 52)) +
    theme(legend.position = 'none',
          axis.text.x = element_text(hjust=1, vjust=1, angle=45, family = 'serif'))
  ggsave(plot = p5, filename = 'fig/Polish/exon_PvsNorm_sele_5.jpeg')
  save(p5, file = "Rdata/fig/BP/exon_PvsNorm_sele_5.Rdata")
  
  # non-zero m
  pb_dat = panel_df %>% 
    dplyr::filter(panel != 'our_ori') %>%
    group_by(panel) %>% 
    summarise(zero_len = sum(MeanMuts == 0),
              nonzero_len = sum(MeanMuts != 0), 
              pan_len = n(),
              zero_percent = zero_len / (zero_len + nonzero_len) * 100) %>% 
    ungroup %>% 
    mutate(group = if_else(panel == 'our', 'y', 'n'))
  pb_dat = pb_dat[match(c("our", "newm", "Gd360_gr", "plse_gr",
                          "burg", "NCC_gr", "MSK", "F1Dx_gr"),
                        pb_dat$panel), ] %>% 
    mutate(pLabel = c("Suppressor", "Newmans'", "Grand360", "Plasma SELECT",
                      "Burgeners'", "NCC150", "MSK-IMPACT","F1CDx"))
  
  p6 = ggplot(pb_dat, aes(y = reorder(pLabel, zero_percent), x=zero_percent, fill=group))+
    geom_bar(stat='identity', color='black', alpha=0.8, width = 0.8)+
    scale_fill_manual(
      limits = c('y' , 'n'),
      values = c("#FF4848", "#316B83")
    )+
    scale_x_continuous(expand = c(0, 0),
                       breaks = seq(0, 100, 20),
                       labels = paste0(seq(0, 100, 20), "%"))+
    expand_limits(x=105)+
    labs(y=NULL, x=expression("Percentage of exons'"~ bar(m)~ " are zero (%)"))+
    geom_text(aes(label = paste0(round(zero_percent), "%"), hjust=-0.3))+
    my_theme +
    theme(legend.position = 'none')
  
  ggsave(plot = p6, filename = 'fig/Polish/exon_PvsNorm_sele_6.jpeg')
  save(p6, file = "Rdata/fig/BP/exon_PvsNorm_sele_6.Rdata")
  
}



############################################################################
# selector vs background and signal ---------------------------------------
# Fri Oct 15 18:52:17 2021 ------------------------------
############################################################################

# load huada gene
hda_df = openxlsx::read.xlsx("data/Zhang_SupData3.xlsx",
                             startRow = 2) %>%
  dplyr::mutate(
    Width = Stop - Start,
    VAF = caseAF * 100,
    isCH = if_else(isCH == 'yes', "Yes", "No")
  ) %>%
  dplyr::filter(!is.na(Start), !is.na(Stop),
                Width == 1, Call != ".")

hda_gr = GRanges(
  seqnames = paste0("chr", hda_df$Chr),
  ranges = IRanges(start = hda_df$Stop, end = hda_df$Stop),
  isCH = hda_df$isCH
) %>% transHg38()

# hda_CH_gr = hda_gr[hda_gr$isCH == "Yes"]


## map huada to exonlib
over_hda = findOverlaps(query = hda_gr, subject = exon_lib_) %>% as.data.frame()

hda_ = hda_df[over_hda$queryHits,] %>%
  dplyr::mutate(
    exon_number = over_hda$subjectHits,
    P.WBC = exon_lib_$P.WBC[over_hda$subjectHits],
    width = width(exon_lib_[over_hda$subjectHits])
  )

### cfDNA percent vs WBC percent
hdGrp0 = hda_ %>%
  group_by(Sample, exon_number) %>%
  summarise(MutInExon = n(), widthExon = unique(width)) %>%
  group_by(exon_number) %>%
  summarise(
    MutsExon = mean(MutInExon),
    MeanMuts = mean(MutInExon) / widthExon * 1000)

panel_df2 = rbind(
  panel_df %>% as.data.frame() %>% 
    # dplyr::filter(panel != "our_ori") %>% 
    select(exon_number, MeanMuts, MutsExon, panel),
  hdGrp0 %>% mutate(panel = 'ctDNA')
) %>% mutate(Group = if_else(panel %in% c("our", "newm", "burg"),
                             'selector', panel))

# plot
if (1) {
  
  ## density --------------------------------------------------------------
  options(scipen = 100)
  colors = ggthemes::wsj_pal()(3)
  scales::show_col(colors)

  pd = panel_df2 %>% 
    dplyr::filter(panel %in% c("our", "our_ori", "ctDNA")) %>%
    dplyr::mutate(MeanMuts = if_else(MeanMuts == 0, 0.01, MeanMuts)) %>% 
    ggplot(aes(x = MeanMuts+0.01)) +
    geom_density(aes(fill = panel, alpha = panel), color = 'black') +
    scale_x_log10() +
    scale_fill_manual(
      limits = c('our', 'our_ori', 'ctDNA'),
      labels = c('Background in Suppressor', 
                 'Background in Enricher',
                 'Mutaion signals in ctDNA'),
      values = colors
    ) +
    scale_alpha_manual(
      limits = c('our', 'our_ori', 'ctDNA'),
      labels = c('Background in Suppressor',
                 'Background in Enricher',
                 'Mutaion signals in ctDNA'),
      values = c(0.9, 0.5, 0.5)
    ) +
    labs(x = expression(bar(m)), y='Density', fill=NULL, color=NULL) +
    my_theme +
    theme(
      legend.position = c(0.4, 0.8),
      legend.background = element_rect(color = 'grey', fill = 'white')
    ) +
    guides(alpha = 'none')
  pd
  ggsave(plot = pd, filename = "fig/Polish/bg_ct_se_density.jpeg")
  save(pd, file = "Rdata/fig/BP/bg_ct_se_density.Rdata")
}
