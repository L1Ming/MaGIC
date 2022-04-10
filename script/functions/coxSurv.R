coxSurv <- function(data, time, status, Group="Group"){
  library(survival)
  library(survminer)
  # browser()
  fml = as.formula(paste0("Surv(", time, "," , status, ") ~ ", 
                          paste(Group, collapse = " + ")))
  
  fit = coxph(formula = fml, data = data)
  
  
  # get statical data
  sy = summary(fit)
  model.pval = sy$logtest['pvalue'] %>% round(3) # use log rank test
  HR = sy$conf.int[, 'exp(coef)'] %>% round(3)
  HR.confint.lower = sy$conf.int[,"lower .95"] %>% round(3)
  HR.confint.upper = sy$conf.int[,"upper .95"] %>% round(3)
  Pval = sy$coefficients[, "Pr(>|z|)"] %>% round(3)
  
  HR_df = data.frame(HR, HR.confint.lower, HR.confint.upper, Pval)
  
  # 
  # message(group)
  res = list(model=fit, model.pval=model.pval, HR_df=HR_df)
  return(res)
}
