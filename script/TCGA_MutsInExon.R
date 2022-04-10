# count the MutsInExon
# Sun Oct 17 19:01:50 2021 ------------------------------

library(dplyr)
library(magrittr)
library(ggplot2)
library(stringr)
library(GenomicRanges)
load("Rdata/exon_lib_.Rdata")
source('script/functions/transHg38.R')
source("script/Utils.R")
# TCGA data ---------------------------------------------------------------

file = list.files("Rdata/Rmaf/", full.names = T)
projects = str_split(file, "[-_]", simplify = T)[, 2]

gp_df = NULL
for(i in seq_along(file)){
  prj = projects[i]
  file_path = file[i]
  maf_dat = load(file_path) %>% 
    get() %>% as.data.frame() %>%
    dplyr::select(contains(
      c('Tumor_Sample_Barcode', "chr", "start", "end", "Tumor_Seq_Allele2"),
      ignore.case = T
    )) %>%
    magrittr::set_colnames(c("Sample", "Chr", "Start", "End", "Allele")) %>%
    dplyr::mutate(PatientID = str_sub(Sample, 1, 12)) %>%
    dplyr::distinct()
  
  maf_gr = GRanges(
    seqnames = maf_dat$Chr,
    ranges = IRanges(start = maf_dat$Start, end = maf_dat$End)
  )
  
  over = findOverlaps(maf_gr, exon_lib_) %>% as.data.frame()
  maf = maf_dat[over$queryHits, ] %>% 
    mutate(exon_number = exon_lib_$exon_number[over$subjectHits],
           exon_width = width(exon_lib_[over$subjectHits]))
  
  mfGrp = maf %>% group_by(PatientID, exon_number) %>% 
    summarise(MutsInExon = n(),
              Project = prj,
              exon_width = unique(exon_width)) %>% 
    ungroup() %>%
    distinct()
  
  gp_df = rbind(gp_df, mfGrp)
}


load("Rdata/project_details.df.Rdata")
gp_df = left_join(gp_df,
                  pro_df %>% select(Project, Organ = Primary.Site),
                  by = "Project")

save(gp_df, file = "Rdata/MutsinExon.Rdata")

# load("Rdata/MutsinExon.Rdata")

a = gp_df %>% group_by(Organ) %>% 
  summarise(one=sum(MutsInExon == 1),
            oneMore = sum(MutsInExon > 1),
            percent1 = one / n() * 100, 
            percent2 = oneMore / n() * 100)

## plot percentage
b = a %>% select(contains('percent'), 'Organ') %>% 
  reshape2::melt(id.var='Organ') %>% 
  mutate(variable = factor(variable, levels = c("percent2", 'percent1'), ordered = T))

fill = c("#F6416C", "#68B0AB") %>% rev
scales::show_col(fill)
p_bar = ggplot() +
  geom_bar(
    data = b,
    mapping = aes(y = Organ, x = value, fill = variable),
    stat = 'identity',
    position = 'stack'
    # color = 'black'
  ) +
  geom_vline(xintercept = 95, linetype = 'dashed') +
  geom_text(
    data = a,
    mapping = aes(
      x = 70,
      y = Organ,
      label = round(percent1, 2)
    ),
    color = 'black',
    hjust = 0
  ) +
  scale_fill_manual(
    limits = c("percent1", 'percent2'),
    values = fill,
    labels = c("1", '>1')
  ) +
  scale_x_continuous(breaks = c(seq(0, 1, 0.25)) * 100)+
  labs(x = "Percentage (%)", y = NULL, fill = 'Mutations\n in Exon') +
  my_theme

ggsave(p_bar, filename = 'fig/Polish/MutsInExon_bar.jpeg')
save(p_bar, file = 'Rdata/fig/BP/MutsInExon_bar.Rdata')

## plot ridge
fill2 = c( "#A2CDCD",  "#FFE1AF", "#D57E7E", "#C6D57E")
scales::show_col(fill2)
p_ridge = ggplot(gp_df, aes(x = MutsInExon, y = Organ, fill = factor(stat(quantile)))) +
  ggridges::stat_density_ridges(
    bandwidth = 0.0466,
    geom = "density_ridges_gradient", 
    calc_ecdf = TRUE,
    color = "#444941",
    quantiles = 4, 
    quantile_lines = TRUE) +
  scale_fill_manual(name = "Quartiles", 
                    values = fill2,
                    limits = 1:4 %>% as.character(),
                    labels = 1:4 %>% as.character())+
  scale_x_log10()+
  labs(x="Mutations in an exon per patient", y='TCGA projects')+
  theme_bw()+
  theme(legend.position = c(0.9, 0.8), 
        legend.background = element_rect(fill='white'))+
  my_theme 
ggsave(p_ridge, filename = 'fig/Polish/MutsInExon_ridge.jpeg')
save(p_ridge, file = 'Rdata/fig/BP/MutsInExon_ridge.Rdata')

# ## concat them
# fp = (p_ridge | p_bar) + plot_layout(guides = "collect")
# ggsave(fp, filename = "Plot/fig3/b.tiff", dpi=600, height = 6, width = 8)