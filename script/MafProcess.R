
# read maf ----------------------------------------------------------------
source('script/functions/read_MAFs.R')
pSaps = paste0("P", 1:23)
nSaps = paste0('SRR7170', 698:706)
load("Rdata/NPC_cedat.Rdata")

# ## for prostate cancer
pMaf = read_MAFs(samples = pSaps,
                 path_prefix = 'E:/PRJNA554329/ToDo/annotated/', multiThread = TRUE)
save(pMaf, file = "Rdata/pMaf.Rdata")

## for normal donors
nMaf = read_MAFs(samples = nSaps,
                 path_prefix = 'D:/Data/NormalcfDNA/annotated/', multiThread = TRUE)
save(nMaf, file = 'Rdata/nMaf.Rdata')

## for ctDNA in NPC
cMaf = read_MAFs(samples = ce_dat$Sample %>% sort, 
                 path_prefix = "D:/Data/NPC_sequencing/ctDNAseq/annotated/",
                 multiThread = T)
save(cMaf, file = 'Rdata/npcMaf.Rdata')

# load Maf and process it -------------------------------------------------

library(GenomicRanges)
load("Rdata/pMaf.Rdata")
load("Rdata/nMaf.Rdata")
load("Rdata/npcMaf.Rdata")
load("Rdata/exon_lib_.Rdata")
source("script/functions/Map2Selector.R")
source("script/functions/Polisher.R")


pMa = pMaf %>% 
  mutate(
    # Tumor_Sample_Barcode = Otherinfo8,
    Start_Position = Start_Position %>% as.numeric,
    End_Position = End_Position %>% as.numeric,
    # Filter=Otherinfo5,
    # AD=Otherinfo7 %>% as.numeric(),
    # DP=Otherinfo6 %>% as.numeric(),
    AD=AD %>% as.numeric(),
    DP=DP %>% as.numeric(),
    VAF=AD / DP * 100 %>% round(digits = 2)) %>%
  dplyr::filter(AD>1, Variant_Classification != 'Silent') %>% 
  as.data.frame()
pMa[is.na(pMa)] = '.'

nMa = nMaf %>% 
  mutate(
    # Tumor_Sample_Barcode = Otherinfo8,
    Start_Position = Start_Position %>% as.numeric,
    End_Position = End_Position %>% as.numeric,
    # Filter=Otherinfo5,
    # AD=Otherinfo7 %>% as.numeric(),
    # DP=Otherinfo6 %>% as.numeric(),
    AD=AD %>% as.numeric(),
    DP=DP %>% as.numeric(),
    VAF=AD / DP * 100 %>% round(digits = 2)) %>%
  dplyr::filter(AD>1, Variant_Classification != 'Silent') %>% 
  as.data.frame()
nMa[is.na(nMa)] = '.'


cMa = cMaf %>% 
  mutate(
    # Tumor_Sample_Barcode = Otherinfo8,
    Start_Position = Start_Position %>% as.numeric,
    End_Position = End_Position %>% as.numeric,
    # Filter=Otherinfo5,
    # AD=Otherinfo7 %>% as.numeric(),
    # DP=Otherinfo6 %>% as.numeric(),
    AD=AD %>% as.numeric(),
    DP=DP %>% as.numeric(),
    VAF=AD/DP*100 %>% round(digits = 2)) %>%
  dplyr::filter(AD>1, Variant_Classification != 'Silent') %>%  #, Filter == "PASS") %>% 
  as.data.frame()
cMa[is.na(cMa)] = "."

# get GRanges objects -----------------------------------------------------

## prostate cancers
p_gr = GRanges(seqnames = pMa$Chromosome,
               ranges = IRanges(start = pMa$Start_Position, 
                                end = pMa$End_Position),
               strand = "*",
               Tumor_Sample_Barcode = pMa$Tumor_Sample_Barcode,
               VAF = pMa$VAF,
               ALT = pMa$Reference_Allele)

## normal donors
n_gr = GRanges(seqnames = nMa$Chromosome,
               ranges = IRanges(start = nMa$Start_Position, 
                                end = nMa$End_Position),
               strand = "*",
               Tumor_Sample_Barcode = nMa$Tumor_Sample_Barcode,
               VAF = nMa$VAF,
               ALT = nMa$Reference_Allele)

## ctDNA in NPC
c_gr = GRanges(seqnames = cMa$Chromosome,
               ranges = IRanges(start = cMa$Start_Position, 
                                end = cMa$End_Position),
               strand = "*",
               Tumor_Sample_Barcode = cMa$Tumor_Sample_Barcode,
               VAF = cMa$VAF,
               ALT = cMa$Reference_Allele
)

# polish by PoN -----------------------------------------------------------
index_p = Polisher(p_gr)
index_n = Polisher(n_gr)
index_c = Polisher(c_gr)

pst_ftd =  pMa[index_p, ]
nm_ftd = nMa[index_n, ]
npc_ftd = cMa[index_c, ]

p_gr2 = p_gr[index_p]
n_gr2 = n_gr[index_n]
cf_gr = c_gr[index_c]
  
save(pst_ftd, p_gr2, file = 'Rdata/Prostate_all_SNVs.Rdata')
save(nm_ftd, n_gr2, file = 'Rdata/Normal_all_SNVs.Rdata')
save(npc_ftd, cf_gr, file = 'Rdata/NPC_all_SNVs.Rdata')
