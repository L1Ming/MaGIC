# Tue Sep 28 18:13:00 2021 ------------------------------


# do parallel  ------------------------------------------------------------
if(1){
  ## define filepath
  chrs = c(paste0("D:/Data/1000G/Prep/ALL.chr", c(1:22), 
                  ".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.txt.gz"),
           "D:/Data/1000G/Prep/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.txt.gz",
           "D:/Data/1000G/Prep/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.txt.gz")
  names(chrs) = paste0("chr", c(1:22, "X", "Y"))
  
  # for test
  chrName='chrY'
  chrPath = chrs[1]
  library(parallel)
  library(foreach)
  library(doParallel)
  library(magrittr)
  
  cl = makeCluster(10)
  registerDoParallel(cl)
  PoN_DF = foreach(
    chrPath = chrs,
    .combine = plyr::rbind.fill,
    .packages = c(
      "plyr",
      "dplyr",
      "stringr",
      "rlang",
      "magrittr"
    ),
    .errorhandling = 'stop',
    .verbose = TRUE
  ) %dopar% {
    
    chrName = str_match(chrPath, "ALL\\.(\\w+)\\.")[, 2]
    dat = data.table::fread(chrPath, data.table = F)
    
    ## add cols to fix
    pon_df = dat %>% transmute(
      CHROM = paste0("chr", CHROM),
      POS = POS,
      REF = REF, 
      ALT = ALT,
      AC = AC %>% as.integer(),
      AF = AF %>% as.double(),
      AlleNum = 2505 - Double0,
      Percent = AlleNum / 2505
    )
    save(pon_df, file = paste0("Rdata/TempPoN/", chrName, "_transVCF.Rdata"))
    # return
    return(pon_df)
  }
  stopImplicitCluster()
  stopCluster(cl)
  save(PoN_DF, file = "Rdata/TempPoN/PoN_DF_allChrs.Rdata")
  # tuneR::play("C:/Users/DELL/Music/LegendNato.mp3")
}

# Thu Oct 14 08:24:45 2021 ------------------------------
## correct the chrY percent
## get the male number: 1233
chrY = data.table::fread("D:/Data/1000G/VCF/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz",
                  skip="#CHROM")
male_number = str_detect(colnames(chrY), pattern = "(^NA|^HG)") %>% sum # 1233
load("Rdata/TempPoN/PoN_DF_allChrs.Rdata")
PoN_DF$Percent = ifelse(PoN_DF$CHROM == "chrY", 
                        PoN_DF$Percent * 2505 / male_number,
                        PoN_DF$Percent)
save(PoN_DF, file = "Rdata/TempPoN/PoN_DF_allChrs.Rdata")

# visualization -----------------------------------------------------------


load("Rdata/TempPoN/PoN_DF_allChrs.Rdata")
# apply(PoN_DF[, c("Percent","AC", "AF")], 2, fivenum)
#           Percent   AC          AF
# [1,] -0.000399361    0 0.000000000
# [2,]  0.000000000    1 0.000199681
# [3,]  0.000399361    2 0.000399361
# [4,]  0.009185304   15 0.002995210
# [5,]  0.999600639 5008 1.000000000

# anyNA(PoN_DF[, c("Percent","AC", "AF")])
# PoN_DF[is.na(PoN_DF$AC), ]


PoN_DF2 = PoN_DF %>% na.omit() %>% 
  dplyr::filter(Percent <= 1)

# Wed Sep 29 16:41:12 2021 ------------------------------
## covert the PoN_DF to genomic ranges object
library(GenomicRanges)
pon_gr_all = GRanges(seqnames=PoN_DF2$CHROM,
                 ranges = IRanges(start=PoN_DF2$POS,
                                  end = PoN_DF2$POS),
                 REF = PoN_DF2$REF,
                 ALT = PoN_DF2$ALT,
                 Percent = PoN_DF2$Percent,
                 AC = PoN_DF2$AC,
                 AF = PoN_DF2$AF)

pon_gr = pon_gr_all[pon_gr_all$Percent >= 0.01 &
                      pon_gr_all$AC >= 1000]
save(pon_gr, file = "Rdata/pon_gr_hg19_1000G.Rdata")
source("script/functions/transHg38.R")
pon_gr = pon_gr %>% transHg38()
save(pon_gr, file="Rdata/pon_gr_hg38_1000G.Rdata")


# Get Exonlib.Rdata -------------------------------------------------------

# Load GTF files
G38 = rtracklayer::import("data/reference/gencode.v22.annotation.gtf")
# remove the version number suffix of Ensembl Gene Id
G38$gene_id = substr(G38$gene_id, start = 1, stop = 15)

ens2symbol = mcols(G38)[c("gene_id", "gene_name")] %>%
  magrittr::set_colnames(c("Gene_id", "Symbol")) %>%
  as.data.frame() %>%
  dplyr::distinct()

# filter and rearrange this file by exon
exons.join = G38[mcols(G38)$type == "exon"] %>%
  split(., mcols(.)$gene_id) %>%
  # remove the overlaps among exons from different transcripts
  IRanges::reduce()

exon_lib = exons.join %>% unlist()
exon_lib = exon_lib[seqnames(exon_lib) != "chrM"]
exon_lib$Gene_id = names(exon_lib)
names(exon_lib) = NULL
exon_df = exon_lib %>% as.data.frame() %>% distinct()
exon_df = left_join(exon_df, ens2symbol, by='Gene_id')
exon_lib = GRanges(seqnames = exon_df$seqnames, 
                   ranges = IRanges(start = exon_df$start, 
                                    end = exon_df$end),
                   Symbol = exon_df$Symbol)

# exon_lib = exons.join %>% unlist() %>% IRanges::reduce()

# here might be confusing, as the x$exon_number describe
#   the seq of exon in a gene(replaced by ExonofGene),
#   but here, we use the exon_number to describe the whole genome exon seq.
exon_lib$exon_number = 1:length(exon_lib)
save(exon_lib, file = "Rdata/exon_lib.Rdata")

# Thu Oct 14 08:20:59 2021 ------------------------------
###################################################
#   Get the probability of a exon is a CH-exon    #
###################################################

# Tue Oct 12 18:01:38 2021 ------------------------------

library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
library(GenomicRanges)
library(purrr)
library(readr)

# part0: Data split for limited RAM ------------------------------
if(1){
  
  file_split <- function(chr, chr_path, eachfile_lines_num = 100) {
    # title
    t = data.table::fread(chr_path, skip = "#CHROM", nrows = 0)
    
    c = file(chr_path, 'r')
    i <-
      1
    while (TRUE) {
      varname <- paste(chr, i, sep = "_")
      assign(varname, readLines(c, n = eachfile_lines_num))
      if (str_sub(get(varname)[1], 1, 2) != "##") {
        write.table(
          t,
          file = paste0("D:/Data/1000G/Temp/", varname, ".txt"),
          sep = '\t',
          quote = F,
          row.names = F
        )
      }
      write.table(
        get(varname),
        file = paste0("D:/Data/1000G/Temp/", varname, ".txt"),
        quote = F,
        row.names = F,
        col.names = F,
        append = T
      )
      
      if (length(get(varname)) < eachfile_lines_num) {
        break
      }
      else {
        i <- i + 1
      }
    }
    return(i)
  }
  
  
  chrs = paste0("chr", c(1:22, "X", "Y"))
  chr_df = data.frame(
    chrs = chrs,
    path = c(
      paste0(
        "D:/Data/1000G/VCF/ALL.chr",
        c(1:22),
        ".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
      ),
      "D:/Data/1000G/VCF/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz",
      "D:/Data/1000G/VCF/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz"
    )
  )
  for (chr in chrs) {
    chr_path = chr_df$path[chr_df$chrs == chr]
    file_split(chr, chr_path, eachfile_lines_num = 1000000)
  }
}

# # take the error
# chr = 'chr22'
# chr_path = chr_df$path[chr_df$chrs == chr]
# file_split(chr = chr, chr_path, eachfile_lines_num = 1000000)

# part1: read files and subset it by pon_ultrasig ------------------------------

if(1){
  
  memory.limit(300000)
  
  message("Start part1 at:\t", Sys.time())
  files = list.files("D:/Data/1000G/Temp/", full.names = T)
  load("Rdata/TempPoN/Pon_ultrasig.Rdata") # from GetExonLib.R
  
  chrs_all = {list.files("D:/Data/1000G/Temp/") %>% 
      str_split("[.]",simplify = T)}[, 1]
  
  chrs_done = {list.files("D:/Data/1000G/Count/subset/") %>% 
      str_split("[.]", simplify = T)}[, 1]
  chrs = setdiff(chrs_all, chrs_done)
  # # take the error
  # load("Rdata/TempPoN/Pon_ultrasig.Rdata")
  # chrs = c('chr22_1', "chr22_2")
  # chrs = 'chrY_1'
  for (chr in chrs) {
    try(expr = {
      chr_path = paste0("D:/Data/1000G/Temp/", chr, ".txt")
      a = data.table::fread(chr_path, skip = "#CHROM")
      chr_gr = GRanges(
        seqnames = paste0('chr', a$`#CHROM`),
        ranges = IRanges(start = a$POS, end = a$POS)
      )
      
      ## subset: in pon_ultrasig
      over_sig = findOverlaps(query = chr_gr,
                              subject = pon_ultrasig,
                              type = 'equal') %>% as.data.frame()
      a[over_sig$queryHits, ] %>%
        data.table::fwrite(file = paste0("D:/Data/1000G/Count/subset/", chr, ".txt.gz"))
      rm(a, chr_gr, over_sig)
    })
    gc()
  }
  rm(list = ls())
  gc()
  # protect the computer
  Sys.sleep(120)
}


## part 2, add exon info ------------------------------
if(1){
  # chr = "chr9_1"
  memory.limit(300000)
  
  message("Start part2 at:\t", Sys.time())
  
  load("Rdata/exon_lib.Rdata")
  source("script/functions/transHg38.R")
  chrs = {list.files("D:/Data/1000G/Temp/") %>% 
      str_split("[.]",simplify = T)}[, 1]
  
  for (chr in chrs) {
    try(expr = {
      a_sig = data.table::fread(
        file = paste0("D:/Data/1000G/Count/subset/", chr, ".txt.gz"),
        # file = "D:/Data/1000G/test/test.vcf.gz",
        skip = "#CHROM"
      )
      chr_sig = GRanges(
        seqnames = paste0('chr', a_sig$`#CHROM`),
        ranges = IRanges(start = a_sig$POS, end = a_sig$POS)
      ) %>%
        transHg38()
      over = findOverlaps(query = chr_sig, subject = exon_lib) %>% as.data.frame()
      a_sig[over$queryHits,] %>%
        dplyr::mutate(exon_number = over$subjectHits) %>%
        data.table::fwrite(file = paste0("D:/Data/1000G/Count/exon_info/", chr, "_add.txt.gz"))
      rm(a_sig, chr_sig, over)
    })
    gc()
  }
  
  rm(list=ls())
  gc()
  Sys.sleep(120)
  # }
  
  
  ## part 3, count number of exon ------------------------------
  # if(1){
  message("Start part3 at:\t", Sys.time())
  
  chrs = {list.files("D:/Data/1000G/Temp/") %>% 
      str_split("[.]",simplify = T)}[, 1]
  for (chr in chrs) {
    try(expr = {
      a_match = data.table::fread(
        file = paste0("D:/Data/1000G/Count/exon_info/", chr, "_add.txt.gz"),
        # file = "D:/Data/1000G/test/test.vcf.gz",
        skip = "#CHROM")
      
      b = a_match %>%
        dplyr::select(starts_with(c("HG", "NA")), exon_number) %>%
        dplyr::mutate(exon_number = exon_number %>% as.character())
      
      b2 = split(b, b$exon_number) %>%
        purrr::map_df(
          .f = function(mat) {
            row_x = apply(mat, 2, function(x) {
              sum(str_detect(x, "0|0", negate = T)) != 0
            })
            return(sum(row_x))
          }
        ) %>% t %>% as.data.frame() %>%
        dplyr::mutate(exon_number = rownames(.)) %>%
        dplyr::rename(MutNum = V1) %>%
        select(exon_number, MutNum)
      
      data.table::fwrite(b2, file = paste0("D:/Data/1000G/Count/exon_count/", chr, "_MutNum.csv"))
      rm(a_match, b, b2)
    })
    gc()
  }
  
  
  
  ## part4: concat Counts files ----------------------------------------------
  
  b = list.files("D:/Data/1000G/Count/exon_count/", full.names = T) %>% 
    purrr::map_df(.f = data.table::fread) %>% 
    # mutate(Percent = MutNum/2505) %>% 
    group_by(exon_number) %>% 
    summarise(MutNum = sum(MutNum),
              Percent = MutNum / 2505
              ) %>% 
    ungroup
  # fivenum(b$MutNum)
  # #[1]    0    0    1    1 2505
  data.table::fwrite(b, file = "D:/Data/1000G/Count/ExonMutCounts.csv")
}



