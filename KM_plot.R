library(tidyverse)
library(survminer)
library(survival)
Train <- data.table::fread("data/CGGA.mRNAseq_693.RSEM-genes.20200506.txt") %>%
  column_to_rownames(., var = "Gene_Name") %>%
  t() %>%
  as.data.frame()
cli <- data.table::fread("data/CGGA.mRNAseq_693_clinical.20200506.txt") %>%
  select(.,CGGA_ID, OS, "Censor (alive=0; dead=1)") %>%
  column_to_rownames(., var = "CGGA_ID")
colnames(cli) <- c("OS.time", "OS")
identical(rownames(cli), colnames(Train))

Train <- log2(Train + 1)
cli$OS.time <- cli$OS.time/365

gene <- c("ABCC3", "ANXA1", "COL5A2", "EMP3", "G0S2", "GDF15", "TAGLN2", "VEGFA", "VEPH1")
plotp <- list()
for (name in gene) {
  rt = cbind(cli, Train[, name])
  rt$risk=ifelse(Train[, name]>median(Train[, name]),'high','low')
  diff=survdiff(Surv(OS.time, OS) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(OS.time, OS) ~ risk, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     legend.title="Risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(year)",
                     break.time.by = 1,
                     palette=c("darkblue", "orange"),
                     risk.table=TRUE,
                     risk.table.title=name,
                     risk.table.col = "strata",
                     risk.table.height=.25)
  plotp[[name]] <- surPlot
}
plotp[[9]]
dev.copy2pdf(file = "figures/day4/KM/KM_VEPH1.pdf")
