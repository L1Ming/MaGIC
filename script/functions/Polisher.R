Polisher = function(gr, hg='hg38'){
  library(GenomicRanges)
  
  hg = match.arg(hg, choices = c("hg38", "hg19"))
  
  ## load pon_gr
  if(hg=="hg38"){
    load("Rdata/pon_gr_hg38_1000G.Rdata")
  }else{
    load("Rdata/pon_gr_hg19_1000G.Rdata")
  }
  
  ## overlap and return filtered index
  overlap = findOverlaps(query = pon_gr, subject = gr, type='equal') %>% 
    as_tibble() %>% mutate(
      queryALT = pon_gr$ALT[queryHits],
      subjectALT = gr$ALT[subjectHits]
    ) %>% 
    dplyr::filter(queryALT == subjectALT)
  
  message("Overlaped ", length(overlap$subjectHits %>% unique),
          " of ", length(gr),
          " at ", round(100*length(overlap$subjectHits %>% unique)/length(gr), 3), "%")
  
  index = setdiff(1:length(gr), overlap$subjectHits)
  
  return(index)
}