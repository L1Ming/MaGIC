library(GenomicRanges)

exonCount = data.table::fread( 
  "D:/Data/1000G/Count/ExonMutCounts.csv"
)# from GetPoN.R

load("Rdata/1000G_exonMutCounts.csv")
load(file = "Rdata/exon_lib.Rdata") # from GetPoN.R

# add exonMutCount info of WBC to exon_lib
exon_lib1 = exon_lib[exonCount$exon_number]
exonCount = exonCount[match(exon_lib1$exon_number, exonCount$exon_number), ]
exon_lib1$P.WBC = exonCount$Percent
exon_lib2 = exon_lib[setdiff(1:length(exon_lib), exonCount$exon_number)]
exon_lib2$P.WBC = 0
exon_lib_ = c(exon_lib1, exon_lib2) %>% sort() ;rm(exon_lib1, exon_lib2)
save(exon_lib_, file = 'Rdata/exon_lib_.Rdata')
