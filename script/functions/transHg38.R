transHg38 = function(gr, from = 'hg19') {
  library(GenomicRanges)
  library(dplyr)
  library(rtracklayer)
  
  chain <- import.chain(
    "data/reference/hg19ToHg38.over.chain"
    )
  ## only trans it when from is equal to hg19
  if (from == "hg19")
    gr = gr %>% liftOver(chain = chain) %>% unlist
  
  return(gr)
}

## test
# pon_gr = transHg38(pon_gr, from = 'hg19')
# save(pon_gr, file = "Rdata/pon_gr_hg38_1000G.Rdata")
