# Mon Oct 11 14:47:33 2021 ------------------------------

library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
library(GenomicRanges)
library(purrr)
source("script/functions/transHg38.R")
source('script/Utils.R')

# load data ---------------------------------------------------------------

load(file = "Rdata/Normal_all_SNVs.Rdata", verbose = T)

norm_gr = GRanges(
  seqnames = nm_ftd$Chromosome,
  ranges = IRanges(
    start = nm_ftd$Start_Position,
    end = nm_ftd$End_Position
  ),
  Donor = nm_ftd$Tumor_Sample_Barcode
)

# norm_gr$Donor %>% table # about 1M muts each one

load("Rdata/exon_lib_.Rdata")


# map to exonlib ----------------------------------------------------------

## map selector to it

source("script/functions/Map2Selector.R")
se_gr = Map2Selector(panel = 'Our_ori')
over_se = findOverlaps(query = se_gr, 
                       subject = exon_lib_, 
                       type='equal') %>% as.data.frame()


## map normal to exonlib
over_norm = findOverlaps(query = norm_gr, 
                         subject = exon_lib_) %>% as.data.frame()
norm_ = nm_ftd[over_norm$queryHits,] %>%
  dplyr::mutate(
    exon_number = over_norm$subjectHits,
    width = width(exon_lib_[over_norm$subjectHits]),
    P.WBC = exon_lib_$P.WBC[over_norm$subjectHits],
    Symbol = exon_lib_$Symbol[over_norm$subjectHits],
    isSele = if_else(exon_number %in% over_se$subjectHits, "Yes", "No")
  )
save(over_norm, file = "Rdata/over_norm.Rdata")

### count the average mutation number of each exon
nmGrp0 = norm_ %>% as_tibble() %>% 
  group_by(Tumor_Sample_Barcode, exon_number) %>%
  summarise(MutInExon = n(), widthExon = unique(width), 
            # isNPC = unique(isNPC),
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
index = exon_lib_[!(exon_lib_$exon_number %in% nmGrp0$exon_number)]
NotInDex = tibble(exon_number = index$exon_number,
                  MutsExon = 0,
                  P.WBC = index$P.WBC, 
                  MeanMuts = 0)
nmGrp = rbind(nmGrp0, NotInDex) %>% 
  mutate(isSele = if_else(exon_number %in% over_se$subjectHits, "Yes", "No")) %>% 
  distinct()


# plot
if (1) {
  # label selector  --------------------------------------------------------------
  barm = 1
  wbcc = 0.001
  
  nmGrp$Group = ifelse(nmGrp$MeanMuts < barm & 
                         nmGrp$P.WBC < wbcc, 'left', 'filter')
  nmGrp$MeanMuts_edit = ifelse(nmGrp$MeanMuts < 0.05, 0.05, nmGrp$MeanMuts)
  # colors = c("#F34D4D", "#8991BF")
  colors = ggthemes::wsj_pal()(2)
  # scales::show_col(colors)


  p1_1 = ggplot(data = nmGrp, aes(x = P.WBC, y = MeanMuts_edit)) +
    geom_point( color = 'grey'
      #shape = 21, fill = 'grey', color = 'white'
      ) +
    geom_point(data = nmGrp[nmGrp$isSele == "Yes", ],
               aes(x = P.WBC, y = MeanMuts_edit, fill=Group),
                # alpha=0.5,
               shape = 21, color = 'white', size = 2.5) +
    scale_fill_manual(limits = c("left", 'filter'),
                      values = colors) +
    scale_y_log10(breaks = c(0.05, 0.01, 0.1, 1, 10, 100, 1000),
                  labels = c(0, 0.01, 0.1, 1, 10, 100, 1000))+
    scale_x_continuous(breaks = seq(0, 1, 0.1))+
    labs(x = "Probability of mutations occur to an exon in WBC",
         y = expression(bar(m)~" of cfDNA in healthy donors")) +
    my_theme+
    theme(legend.position = 'none', panel.grid.minor.y = element_blank())
  ggsave(p1_1, filename = 'fig/Polish/exon_PvsNorm_sele.jpeg', 
         width = 14, height = 12, units = 'cm')
  save(p1_1, file = "Rdata/fig/BP/exon_PvsNorm_sele.Rdata")
  
  
  ## amplify the view
  p1_2  = p1_1 +
    coord_cartesian(xlim=c(0, 0.05))+
    scale_x_continuous(breaks = seq(0, 0.05, length.out=5))
  ggsave(p1_2, filename = 'fig/Polish/exon_Pvs Norm_sele_cartesian.jpeg',
         width = 14, height = 12, units = 'cm')
  save(p1_2, file = "Rdata/fig/BP/exon_PvsNorm_sele_cartesian.Rdata")
  
}


## get optimal Catcher ---------------------------------------------------------

index = nmGrp$exon_number[nmGrp$isSele == "Yes" & nmGrp$Group == 'left'] %>% 
  unique

strand(exon_lib_) = "*"
exon_lib_[index] %>% as_tibble() %>% 
  data.table::fwrite(file = "panels/Catcher.csv")

# exon_lib_[index] %>% length() # 811 exons
# exon_lib_[index] %>% reduce() %>% length() # 811 exons
# exon_lib_[index] %>% reduce() %>% width %>% sum # 364475 bp of 811 regions
