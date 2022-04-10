
###########################################################################
# NMR might be better than bTMB -------------------------------------------
###########################################################################

if(1){
  
  # data prep --------------
  load("Rdata/exon_lib_.Rdata")
  source('script/Utils.R')
  source('script/functions/transHg38.R')
  library(GenomicRanges)
  
  ## add classical CH genes (15 genes, in High-Inten ...)
  ## Paper title: High-intensity sequencing reveals the sources of plasma circulating cell-free DNA variants
  conCH = c(
    "DNMT3A", "TET2", "ASXL1", "PPM1D", "TP53",
    "JAK2", "RUNX1", "SF3B1", "SRSF2", "IDH1", 
    "IDH2", "U2AF1", "CBL", "ATM", "CHEK2")
  
  # load data from high-intensive cfDNA --------------------------------
  
  ## phenotypes
  hii_pheno = openxlsx::read.xlsx(
    "data/Razavi_Sup_Tab.xlsx",
    sheet = 1,
    startRow = 2
  ) %>% select(Tumor_Sample_Barcode = Patient.ID, Cancer = Tissue) %>%
    dplyr::filter(Cancer != 'Healthy') %>% 
    distinct()
  
  ## cfDNA
  hii_cf = openxlsx::read.xlsx(
    "data/Razavi_Sup_Tab.xlsx",
    sheet = 4,
    startRow = 2
  ) %>%
    dplyr::filter(str_detect(Tumor_Sample_Barcode, "^M")) %>%
    mutate(
      VAF = t_alt_count / t_depth * 100,
      isCH = if_else(Hugo_Symbol %in% conCH, "Yes", "No")) %>% 
    left_join(hii_pheno, by='Tumor_Sample_Barcode')
  
  ## wbc
  hii_wbc = openxlsx::read.xlsx(
    "data/Razavi_Sup_Tab.xlsx",
    sheet = 5,
    startRow = 2
  ) %>%
    dplyr::filter(str_detect(Tumor_Sample_Barcode, "^M")) %>% 
    mutate(VAF = n_alt_count / n_depth * 100) %>%
    left_join(hii_pheno, by = 'Tumor_Sample_Barcode')
  
  ## tissue
  hii_tis = openxlsx::read.xlsx(
    "data/Razavi_Sup_Tab.xlsx",
    sheet = 6,
    startRow = 2
  ) %>%
    dplyr::filter(str_detect(Tumor_Sample_Barcode, "^M")) %>% 
    mutate(VAF = t_alt_count / t_depth * 100) %>%
    left_join(hii_pheno, by='Tumor_Sample_Barcode')
  
  
  ## trans to gr objects ------------------
  hicf_gr = GRanges(
    seqnames = paste0("chr", hii_cf$Chromosome),
    ranges = IRanges(
      start = hii_cf$Start_Position,
      end = hii_cf$End_Position
    ),
    rowName = 1:nrow(hii_cf)
  ) %>% transHg38()
  
  hiw_gr = GRanges(
    seqnames = paste0("chr", hii_wbc$Chromosome),
    ranges = IRanges(
      start = hii_wbc$Start_Position,
      end = hii_wbc$End_Position
    )
  ) %>% transHg38()
  hit_gr = GRanges(
    seqnames = paste0('chr', hii_tis$Chromosome),
    ranges = IRanges(
      start = hii_tis$Start_Position,
      end = hii_tis$End_Position
    )
  ) %>% transHg38()
  
  
  # bTMB vs TMB in High-intensive -----------------------------------
  bTMB_df = hii_cf %>% group_by(Tumor_Sample_Barcode) %>% 
    summarise(Cancer=unique(Cancer),
              SNVs = n()) %>% 
    ungroup()
  
  TMB_df = hii_tis %>% group_by(Tumor_Sample_Barcode) %>% 
    summarise(Cancer = unique(Cancer),
              SNVs = n())
  
  vs_df = full_join(bTMB_df, TMB_df, 
                    by=c("Cancer", 'Tumor_Sample_Barcode'), 
                    suffix=c("_blood", "_tissue")) %>% 
    reshape2::melt(id.var=c("Cancer", 'Tumor_Sample_Barcode'), 
                   variable.name = c('Type'),
                   value.name = 'SNVs')
  # fill = c("#6D9886", "#C6D57E")
  fill = ggsci::pal_aaas(alpha = 0.6)(2) %>% rev()
  library(gghalves)
  TMB_vs = ggplot(data = vs_df, 
                  aes(x=Type, y=SNVs))+
    geom_half_violin(data=vs_df %>% dplyr::filter(Type=="SNVs_blood"), 
                     mapping = aes(x=Type, y=SNVs),fill=fill[1], 
                     # alpha=0.7, 
                     color=fill[1],
                     side = 'l')+
    geom_half_violin(data=vs_df %>% dplyr::filter(Type!="SNVs_blood"), 
                     mapping = aes(x=Type, y=SNVs), fill=fill[2], 
                     # alpha=0.7, 
                     color=fill[2],
                     side = 'r')+
    geom_half_boxplot(data=vs_df %>% dplyr::filter(Type=="SNVs_blood"), 
                      mapping = aes(x=Type, y=SNVs),width=0.1,
                      side = 'l')+
    geom_half_boxplot(data=vs_df %>% dplyr::filter(Type!="SNVs_blood"), 
                      mapping = aes(x=Type, y=SNVs), width=0.1,
                      side = 'r')+
    geom_line(aes(group=Tumor_Sample_Barcode, color=Cancer))+
    geom_point(aes(color=Cancer))+
    scale_y_log10()+
    ggpubr::stat_compare_means(label.x = 1.25, color=p_color, paired = T)+
    labs(x=NULL)+
    ggsci::scale_color_npg()+
    scale_x_discrete(limits = c('SNVs_blood', 'SNVs_tissue'),
                     labels = c("Blood","Tissue"))+
    theme(legend.position = c(0.85, 0.8),
          legend.background = element_rect(color='grey', fill='white'))+
    my_theme
  
  ggsave(TMB_vs, filename = "fig/Polish/bTMBvsTMB.split.jpeg")
  save(TMB_vs, file = 'Rdata/fig/DT/bTMBvsTMB.split.jpeg.Rdata')
  
  
  # signal_to_noise_ratio_discard
  # resource analysis -------------------------------------------------------
  over_cf = findOverlaps(hicf_gr, exon_lib_) %>% as.data.frame()
  over_wbc = findOverlaps(hiw_gr, exon_lib_) %>% as.data.frame()
  over_tis = findOverlaps(hit_gr, exon_lib_) %>% as.data.frame()
  ## cfDNA
  hii_cf_ = hii_cf[over_cf$queryHits,] %>%
    dplyr::mutate(exon_number = over_cf$subjectHits) %>% 
    group_by(Tumor_Sample_Barcode, exon_number) %>% 
    summarise(MutsInExon = n(),
              VAF = mean(VAF),
              Cancer = unique(Cancer)) %>% 
    ungroup() %>%
    dplyr::filter(MutsInExon == 1) %>% 
    distinct()
  ## wbc
  hii_wbc_ = hii_wbc[over_wbc$queryHits,] %>%
    dplyr::mutate(exon_number = over_wbc$subjectHits) %>% 
    group_by(Tumor_Sample_Barcode, exon_number) %>% 
    summarise(MutsInExon = n(),
              VAF = mean(VAF),
              Cancer = unique(Cancer)) %>% 
    ungroup() %>%
    # dplyr::filter(MutsInExon == 1) %>% 
    distinct()
  ## tissue
  hii_tis_ = hii_tis[over_tis$queryHits,] %>%
    dplyr::mutate(exon_number = over_tis$subjectHits) %>% 
    group_by(Tumor_Sample_Barcode, exon_number) %>% 
    summarise(MutsInExon = n(),
              VAF = mean(VAF),
              Cancer = unique(Cancer)) %>% 
    ungroup() %>%
    # dplyr::filter(MutsInExon == 1) %>% 
    distinct()
  
  
  # get matched exons in MutsInExon==1 ---------------------------------------
  tis_match = inner_join(
    x = hii_cf_,
    y = hii_tis_,
    by = c("Tumor_Sample_Barcode", "exon_number", "Cancer"),
    suffix = c(".blood", ".tissue")
  )
  
  wbc_match = inner_join(
    x = hii_cf_,
    y = hii_wbc_,
    by = c("Tumor_Sample_Barcode", "exon_number", "Cancer"),
    suffix = c(".blood", ".WBC")
  )
  
  unknow_source = anti_join(
    x = hii_cf_,
    y = tis_match,
    by = c("Tumor_Sample_Barcode", "exon_number", "Cancer")
  ) %>%
    anti_join(y = wbc_match,
              by = c("Tumor_Sample_Barcode", "exon_number", "Cancer"))
  
  purrr::map_chr(.x = list(tis_match, wbc_match, unknow_source), .f=nrow)
  # [1] "682"  "518"  "1330" # for cfDNA's MutsInExon == 1
  # [1] "51" "60" "98" # for cfDNA's MutsInExon > 1
  
  
  # to see the VAF distribution ---------------------------------------------
  
  VAF_df = rbind(
    tis_match %>% 
      select(Tumor_Sample_Barcode, exon_number, 
             MutsInExon = MutsInExon.blood, 
             VAF = VAF.blood, Cancer) %>% 
      mutate(Match = 'Tissue'),
    wbc_match %>%  
      select(Tumor_Sample_Barcode, exon_number,
             MutsInExon = MutsInExon.blood, 
             VAF = VAF.blood, Cancer) %>% 
      mutate(Match = 'WBC'),
    unknow_source %>% mutate(Match = 'Unknown')
  )
  
  ## boxplot 
  box = ggplot(VAF_df, aes(x=Match, y=VAF, fill=Match))+
    geom_violin(aes(fill=Match), color=NA, alpha=0.8)+
    geom_boxplot(width=0.1, fill='white')+
    scale_y_log10()+
    ggthemes::scale_fill_wsj(limits = c('Tissue', "WBC", "Unknown"),
                             labels = c("Tissue-matched", "WBC-matched", "Unknown-source"))+
    ggpubr::stat_compare_means(comparisons = list(c('Tissue', 'Unknown')),
                               mapping = aes(label = ..p.signif..),
                               label.y = log10(250),
                               color=p_color,
                               method.args = list(alternative='greater'))+
    ggpubr::stat_compare_means(comparisons = list(c('Tissue', 'WBC')), 
                               mapping = aes(label = ..p.signif..),
                               label.y=log10(100),
                               color=p_color,
                               method.args = list(alternative='greater'))+
    scale_x_discrete(limits = c('Tissue', "WBC", "Unknown"),
                     labels = c("Tissue-matched", "WBC-matched", "Unknown-source"))+
    labs(x="Mutations in ctDNA")+
    my_theme +
    theme(legend.position = 'none')
  
  ggsave(plot=box, filename = 'fig/Polish/boxplot_VAF_hii.jpeg')
  save(box, file = 'Rdata/fig/BP/boxplot_VAF_hii.Rdata')
  
  ## density plot
  den = ggplot(VAF_df, aes(VAF))+
    geom_density(alpha=0.54, mapping=aes(fill=Match), color=p_color)+
    scale_x_log10(expand = c(0, 0))+
    geom_vline(xintercept = 1, linetype='dashed', color='#6B4F4F')+
    scale_y_continuous(expand = c(0, 0)) +
    expand_limits(y=1.1)+
    labs(x="VAF of ctDNA", y='Density', fill='ctDNA of', color='ctDNA of')+
    ggthemes::scale_fill_wsj(limits = c('Tissue', "WBC", "Unknown"),
                             labels = c("Tissue-matched", "WBC-matched", "Unknown-source"))+
    ggthemes::scale_color_wsj(limits = c('Tissue', "WBC", "Unknown"),
                              labels = c("Tissue-matched", "WBC-matched", "Unknown-source"))+
    my_theme +
    theme(legend.position = c(0.65, 0.75),
          legend.background = element_rect(color='grey', fill = 'white'))
  ggsave(plot=den, filename = 'fig/Polish/density_VAF_hii.jpeg')
  save(den, file = 'Rdata/fig/BP/density_VAF_hii.Rdata')
  
  
  # number comparison: tissue vs blood -------------------------------------
  
  ## cfDNA
  hii_cf2 = hii_cf[over_cf$queryHits,] %>%
    mutate(exon_number = over_cf$subjectHits)
  hicGrp = hii_cf2 %>% group_by(Tumor_Sample_Barcode, exon_number) %>%
    summarise(MutsInExon = n(),
              Cancer = unique(Cancer)) %>% 
    ungroup() %>% 
    group_by(Cancer) %>% 
    summarise(one=sum(MutsInExon == 1),
              oneMore = sum(MutsInExon > 1),
              percent1 = one/n(), 
              percent2 = oneMore / n())
  
  ## tissue
  hii_tis2 = hii_tis[over_tis$queryHits,] %>%
    mutate(exon_number = over_tis$subjectHits)
  hitGrp = hii_tis2 %>% group_by(Tumor_Sample_Barcode, exon_number) %>%
    summarise(MutsInExon = n(),
              Cancer = unique(Cancer)) %>% 
    ungroup() %>% 
    group_by(Cancer) %>% 
    summarise(one=sum(MutsInExon == 1),
              oneMore = sum(MutsInExon > 1),
              percent1 = one/n(), 
              percent2 = oneMore / n())
  ## combine
  comb = rbind(hicGrp %>% select(Cancer, one, oneMore) %>% 
                 reshape2::melt(id.var="Cancer") %>% mutate(Source='ctDNA'),
               hitGrp %>% select(Cancer, one, oneMore) %>% 
                 reshape2::melt(id.var="Cancer") %>% mutate(Source='Tissue')) %>% 
    dplyr::rename(NumGroup = variable, SNVs = value)
  
  fill = c("#F6416C", "#68B0AB") %>%  rev
  cb_bar = ggplot(comb, aes(x=Cancer, y=SNVs, fill=NumGroup))+
    geom_bar(stat='identity', position='dodge', color=p_color)+
    scale_y_continuous(expand = c(0, 0))+
    facet_grid(Source~.)+#, 
    # scales = 'free',
    # labeller = labeller(NumGroup=c(ctDNA='One mutation in an exon', 
    #                                oneMore='>1 mutations in an exon')))+
    geom_text(aes(label=SNVs), position=position_dodge(0.9), vjust=-0.5, color='black')+
    scale_fill_manual(
      limits = c('one', 'oneMore'),
      values = fill,
      labels = c('1', ">1")) +
    expand_limits(y=1300)+
    # expand_limits(y=2000)+
    # scale_y_log10()+
    labs(x=NULL, y='Number of mutated regions', fill='Mutations\nin exons')+
    my_theme +
    theme(strip.background = element_rect(color = 'white', 
                                          fill = 'white'))
  ggsave(cb_bar, filename = 'fig/Polish/bloodvsTissue_numSNVs_bar.tiff', dpi = 600)
  save(cb_bar, file = 'Rdata/fig/BP/bloodvsTissue_numSNVs_bar.Rdata')
  
  ## percentage of 
  fpp = comb %>% group_by(Source, Cancer) %>% 
    summarise(Percent = 100 * SNVs[NumGroup == "oneMore"] / 
                (SNVs[NumGroup == "oneMore"] + SNVs[NumGroup == "one"])) %>% 
    ggplot(aes(x=Cancer, y=Percent, fill=Source)) +
    labs(x=NULL, y='Percentage (%)', fill=NULL) +
    scale_y_continuous(expand = c(0, 0))+
    expand_limits(y=10) +
    ggthemes::scale_fill_wsj()+
    geom_bar(stat='identity', alpha = 0.8,
             position = 'dodge', color=p_color) +
    geom_text(aes(label = round(Percent, 2)),
              vjust = -0.5,
              position = position_dodge(0.9))+
    my_theme +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(color='grey'))
  ggsave(fpp, filename = 'fig/Polish/OneMore_percent.jpeg')
  
  
  
  # stack bar of source ---------------------------------------------------
  ## before discard
  stack = VAF_df %>% 
    group_by(Cancer, Match) %>% 
    summarise(Exons = n())
  stack_bar = ggplot(stack, aes(x=Cancer, y=Exons, fill=Match))+
    geom_bar(stat = 'identity', position = 'fill', 
             width = 0.8,
             alpha = 0.8, color=p_color)+
    ggthemes::scale_fill_wsj(limits = c('Tissue', "WBC", "Unknown"),
                             labels = c("Tissue-matched", "WBC-matched", "Unknown-source"))+
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20))+
    labs(x=NULL, y="Percetage of ctDNA variances source (%)")+
    # expand_limits(y=1.03)+
    my_theme
  ggsave(stack_bar, filename = 'fig/Polish/PercentOfcfDNAOrgin.tiff', dpi=600)
  save(stack_bar, file = 'Rdata/fig/BP/PercentOfcfDNAOrgin.Rdata')
  ## discard VAF<1
  discard_df = VAF_df %>% 
    dplyr::filter(VAF > 1) %>% 
    group_by(Cancer, Match) %>% 
    summarise(Exons = n())
  discard_bar = ggplot(discard_df, aes(x=Cancer, y=Exons, fill=Match))+
    geom_bar(stat = 'identity', position = 'fill', 
             width = 0.8,
             alpha = 0.8, color=p_color)+
    ggthemes::scale_fill_wsj(limits = c('Tissue', "WBC", "Unknown"),
                             labels = c("Tissue-matched", "WBC-matched", "Unknown-source"))+
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20))+
    labs(x=NULL, y="Percetage of ctDNA variances source (%)")+
    # expand_limits(y=1.03)+
    my_theme
  ggsave(discard_bar, filename = 'fig/Polish/PercentOfcfDNAOrgin_discard.tiff', dpi=600)
  save(discard_bar, file = 'Rdata/fig/BP/PercentOfcfDNAOrgin_discard.Rdata')
  
  ## signal-to-noise ratio
  snr = rbind(discard_df %>% 
                dplyr::mutate(Group = 'Discard'), 
              stack %>% dplyr::mutate(Group = "Original")) %>% 
    reshape2::dcast(formula = Cancer + Group ~ Match, value.var = "Exons") %>% 
    mutate(lgFC = log2(Tissue / WBC),
           Group = factor(Group, levels = c("Original", "Discard"), ordered = T)) %>% 
    ggplot(aes(x=Cancer, y=lgFC))+
    geom_bar(aes(fill=Group), color = 'black',
             stat='identity', position = 'dodge')+
    ggsci::scale_fill_aaas(alpha = 0.7) +
    geom_hline(yintercept = 0, size = 1, color = p_color) +
    labs(x=NULL, y=expression(log[2]~"(signal-to-noise ratio)"), fill=NULL)+
    my_theme+
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(color = 'grey', fill='white'))
  save(snr, file = "Rdata/fig/BP/signal_to_noise_ratio_discard.Rdata")
  ggsave(snr, filename = "fig/Polish/signal_to_noise_ratio_discard.tiff", dpi=600)
  
  
  
  # cor plot ----------------------------------------------------------------
  
  ## NMR vs bTMB
  ### add info
  cf = hii_cf[over_cf$queryHits,] %>%
    mutate(exon_number = exon_lib_$exon_number[over_cf$subjectHits]) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(bTMB = n(),
              b_mVAF = mean(VAF),
              bNMR = length(unique(exon_number))) %>%
    ungroup()
  
  # cfc = hii_cf[over_cf$queryHits,] %>%
  #   mutate(exon_number = exon_lib_$exon_number[over_cf$subjectHits]) %>%
  #   dplyr::filter(isCH == "No") %>% 
  #   group_by(Tumor_Sample_Barcode) %>% 
  #   summarise(bcTMB = n(),
  #             bc_mVAF = mean(VAF),
  #             bcNMR = length(unique(exon_number))) %>%
  #   ungroup()
  
  cfv = hii_cf[over_cf$queryHits,] %>%
    mutate(exon_number = exon_lib_$exon_number[over_cf$subjectHits]) %>%
    dplyr::filter(VAF > 1) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(bvTMB = n(),
              bv_mVAF = mean(VAF),
              bvNMR = length(unique(exon_number))) %>%
    ungroup()
  
  
  ti = hii_tis[over_tis$queryHits,] %>%
    mutate(exon_number = exon_lib_$exon_number[over_tis$subjectHits]) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(tTMB = n(),
              t_mVAF = mean(VAF),
              tNMR = length(unique(exon_number))) %>%
    ungroup()
  
  d = cf %>%
    left_join(ti, by = 'Tumor_Sample_Barcode') %>%
    left_join(cfc, by = 'Tumor_Sample_Barcode') %>% 
    left_join(cfv, by = 'Tumor_Sample_Barcode') %>% 
    left_join(hii_pheno, by="Tumor_Sample_Barcode")
  d[is.na(d)] <- 0

  
  ##  cor heatmap
  a = d %>% 
    select(bTMB, mVAF = b_mVAF, NMR=bNMR, KMR=bvNMR, tNMR, tTMB) %>% 
    cor() %>% 
    pheatmap::pheatmap(display_numbers = T, border_color = 'white',
                       number_format = "%.3f", angle_col = 0,
                       treeheight_col = 10, treeheight_row = 10,
                       number_color = 'black', 
                       color = colorRampPalette(c("#FAEDC6", "#F05454"))(100))
  
  jpeg(filename = "fig/Polish/cor_pheatmap.jpeg", 
       width = 14, height = 11, units = 'cm', res = 600)
  print(a)
  dev.off()

  save(d, file = 'Rdata/fig/BP/cor_pheatmap.Rdata')
  load("Rdata/fig/BP/cor_pheatmap.Rdata")
  

  
  ## pearson bar plot
  fil_pp = ggsci::pal_nejm()(4) %>% rev
  scales::show_col(fil_pp)
  bar_pp = d %>%
    # bTMB, mVAF = b_mVAF, NMR=bNMR, KMR=bvNMR, tNMR, tTMB
    select(bTMB, mVAF=b_mVAF, NMR=bNMR, KMR=bvNMR, tTMB) %>% 
    cor %>% as.data.frame() %>% 
    select(tTMB) %>% mutate(Variable = rownames(.)) %>% 
    dplyr::filter(Variable != 'tTMB') %>%
    ggplot(aes(x=Variable, y=tTMB, fill=Variable))+
    geom_bar(stat = 'identity', width = 0.6,
             color=p_color)+
    scale_fill_manual(limits = c("bTMB", "mVAF", "NMR", "KMR"),
                      values = fil_pp)+
    scale_x_discrete(limits = c("bTMB", "mVAF", "NMR", "KMR"))+
    scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0))+
    geom_text(aes(label = round(tTMB, 3)), vjust=-0.3, size = 5)+
    expand_limits(y=1)+
    labs(x="Variables", y="Pearson correlation coefficient")+
    my_theme+
    theme(legend.position = 'none',
          text = element_text(size=15))
  
  ggsave(bar_pp, filename = "fig/Polish/NMR_TMB_linePoint.tiff", dpi=600)
  save(bar_pp, file = "Rdata/fig/BP/NMR_TMB_linePoint.Rdata")
}
