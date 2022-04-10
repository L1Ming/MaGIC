

load("Rdata/clin_pheno.Rdata")


# reads stat ------------------------------------------------------
fq_read = data.table::fread("data/stat/fq.reads.txt", header = F) %>% 
  select(1:4) %>% 
  magrittr::set_colnames(c("Sample", "ReadFR", "RawClean", "Counts")) %>% 
  mutate(Counts = str_split(Counts, pattern = ":", simplify = T)[, 2] %>% as.integer()) %>% 
  reshape2::dcast(Sample~RawClean+ReadFR, value.var = "Counts") %>% 
  select(Sample, `Raw_R1` = raw_R1, `Raw_R2` = raw_R2, 
         `Clean_R1` = clean_R1, `Clean_R2` = clean_R2)

# Reads of Mapped, Deduplicated, Size Selected -------------------------
abc = function(path){
  # path = '../stat/align.bam.reads.txt'
  sap_ls = map(LETTERS[1:5], paste0, 1:8) %>% unlist()
  read_stat = {
    readLines(path) %>%
      str_replace_all("^#", "") %>%
      str_split("\\s+\\(", simplify = T) %>%
      as_tibble()
  }[[1]] %>%
    str_split("\\+ 0 ", simplify = T) %>%
    as_tibble() %>%
    magrittr::set_colnames(c("Counts", 'Cols')) %>%
    mutate(
      Cols = if_else(Cols == "", "Sample", Cols),
      Sample = rep(sap_ls, each = 14),
      dd = 1:nrow(.)
    ) %>%
    dplyr::filter(!(dd %in% c(1:40 * 14))) %>%
    dplyr::filter(Cols != "Sample") %>%
    split(.$Sample) %>%
    purrr::map_df(
      .f = function(x) {
        x = x %>% as.data.frame()
        rownames(x) = x$Cols
        temp = x %>%
          select(1) %>%
          t %>% as.data.frame()
        rownames(temp) = x$Sample[1]
        return(temp)
      }
    ) %>%
    mutate(Sample = rownames(.)) %>%
    select(all_of(
      c(
        'Sample',
        'in total',
        'mapped',
        "paired in sequencing",
        "read1",
        "read2",
        "properly paired"
      )
    )) %>%
    magrittr::set_colnames(
      c(
        'Sample',
        "Total",
        "Mapped",
        "Paired",
        "Paired_R1",
        "Paired_R2",
        "Porperly_Paired"
      )
    )
  
  read_stat[, 2:ncol(read_stat)] = read_stat[, 2:ncol(read_stat)]  %>%
    apply(2, as.integer)
  read_stat = read_stat %>%
    mutate(`Porperly_Paired(%)` = Porperly_Paired / Paired * 100 %>% round(3))
  return(read_stat)
}

## run it
map_stat = abc('data/stat/align.bam.reads.txt')
dup_stat = abc('data/stat/rmDup.bam.reads.txt')
siz_stat = abc('data/stat/sizeSele.bam.reads.txt')


read_stat = left_join(ce_dat %>% select(Patient, Sample, Group),
                      fq_read,
                      by='Sample') %>%
  arrange(Sample) %>% 
  mutate(Group = if_else(Group == "Before", "Pre-therapy", "Post-therapy")) %>%
  select(Patient, Sample, Group, everything())


# depth among Enricher and Suppressor --------------------------

dep_e = data.table::fread('data/stat/Enricher_Depth_sizeSele.txt')
dep_s = data.table::fread('data/stat/Suppressor_Depth_sizeSele.txt')
library(GenomicRanges)
source('script/functions/Map2Selector.R')

e_gr = GRanges(seqnames = dep_e$V1, 
               ranges = IRanges(start =dep_e$V2, end = dep_e$V2),
               strand = '*',
               Depth = dep_e$V3,
               Sample = dep_e$V4)

s_gr = GRanges(seqnames = dep_s$V1, 
               ranges = IRanges(start =dep_s$V2, end = dep_s$V2),
               strand = '*',
               Depth = dep_s$V3,
               Sample = dep_s$V4)

ee = Map2Selector(e_gr, panel = 'Our_ori') %>% 
  as.data.frame() %>% select(Sample, Region, Depth) %>% 
  mutate(Region = stringr::str_replace(Region, 'Exon_', 'Region_')) %>% 
  reshape2::dcast(Region~Sample, value.var = 'Depth', fun.aggregate = mean, fill = 0)

ss = Map2Selector(s_gr, panel = 'Our') %>% 
  as.data.frame() %>% select(Sample, Region, Depth) %>% 
  mutate(Region = stringr::str_replace(Region, 'Exon_', 'Region_')) %>% 
  reshape2::dcast(Region~Sample, value.var = 'Depth', fun.aggregate = mean, fill = 0)


## size distribution ----------------------------------------------

bf = data.table::fread('data/sizeStat/size_beforeSele.txt') %>% 
  magrittr::set_colnames(c('Sample', 'Size')) %>% 
  mutate(Size = abs(Size)) %>% 
  group_by(Sample) %>% 
  summarise(
    Min = min(Size),
    Mean = mean(Size),
    Max = max(Size)) %>% 
  ungroup()

af = data.table::fread('data/sizeStat/size_AfterSele.txt') %>% 
  magrittr::set_colnames(c('Sample', 'Size')) %>% 
  mutate(Size = abs(Size)) %>% 
  group_by(Sample) %>% 
  summarise(
    Min = min(Size),
    Mean = mean(Size),
    Max = max(Size)) %>% 
  ungroup()



# save it ------------------------------------------------------------
openxlsx::write.xlsx(
  x = list(read_stat, map_stat, dup_stat, siz_stat, ee, ss, bf, af), 
  file = "data/stat/xlsx/reads.stat.xlsx", 
  sheetName = c("Total Reads", 
                "Alignment",
                "Deduplication",
                "Size selection",
                "Mean Depth of Enricher",
                "Mean Depth of Suppressor",
                "Size before select",
                "Size after select"
                ), 
  asTable = T, border = 'columns',
  zoom = 120, widths = 'auto',
  overwrite = T)


##
# 
## 