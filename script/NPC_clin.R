# Title     : TODO
# Objective : TODO
# Created by: L1Ming
# Created on: 2021/5/14


# Encoding: GB2312

library(dplyr)
library(magrittr)
library(stringr)

# get sample info corresponding with patients

sam_pat = openxlsx::read.xlsx("D:/Data/NPC_sequencing/supfiles/样本及实验信息_表型补充0606.xlsx", sheet=2)
sam_pat = data.frame("Patient"=paste0("P_", sam_pat$ID顺序),
                     "PatientID"=sam_pat$住院ID,
                     "Sample"=sam_pat$接头,
                     "Group"=ifelse(str_split(sam_pat$序号,
                                              pattern = '-',
                                              simplify = TRUE)[,2] == 1,
                                    "Before",
                                    "After"))
str(sam_pat)



# get clinical data ------------------------------------------------------

clin_dat = openxlsx::read.xlsx("D:/Data/NPC_sequencing/supfiles/样本及实验信息_表型补充0606.xlsx", sheet=1, rows = 1:27)

clin_dat = clin_dat[, c("住院号", "性别", "年龄", "吸烟史", "饮酒史", "EBV", "T", "N", "M", "分期", "诱导", "同步", "化疗疗效")]

# rename it
colnames(clin_dat) = c("PatientID", "Sex", "Age", "Smoke", "Drink", "EBV", "T", "N", "M", "Stage", "Induce", "Synchronize", "Response")

# format it
clin_dat$Sex = ifelse(clin_dat$Sex == "男", 1, 0)
clin_dat$Smoke = ifelse(clin_dat$Smoke == "有", 1, 0)
clin_dat$Drink = ifelse(clin_dat$Drink == "有", 1, 0)
clin_dat$EBV = ifelse(clin_dat$EBV == "阴性", 0, clin_dat$EBV) %>% as.numeric()
clin_dat$Synchronize = clin_dat$Synchronize %>% zoo::na.fill("无")
clin_dat$Sync_orNot = ifelse(clin_dat$Synchronize != "无", 1, 0)
clin_dat$NaTuoZhu_only = ifelse(clin_dat$Synchronize == "尼妥珠", 1, 0)
clin_dat$NaTuoZhu_ShunBo = ifelse(str_detect(clin_dat$Synchronize, pattern = "顺铂"), 1, 0)
clin_dat$NaTuoZhu_XiMeiNa = ifelse(str_detect(clin_dat$Synchronize, pattern = "希美纳"), 1, 0)

# combind them
ce_dat = dplyr::left_join(x=sam_pat, y=clin_dat, by="PatientID")
save(ce_dat, file = "D:/Data/NPC_sequencing/wes/script/Rdata/clin_pheno.Rdata")
