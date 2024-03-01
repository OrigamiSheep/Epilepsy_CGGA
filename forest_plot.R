options(scipen = 200)
outTab <- read.csv("data/unicox.csv", header = T, row.names = 1)
#读取输入文件
hrtable <- outTab
Sig=ifelse(hrtable$pvalue<0.001,"***",ifelse(hrtable$pvalue<0.01,"**",ifelse(hrtable$pvalue<0.05,"*","")))
tabletext <- cbind(c("Gene",hrtable$gene),
                   c("HR",format(round(as.numeric(hrtable$HR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(hrtable$HR.95L),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(hrtable$HR.95H),3),scientific = F,nsmall = 3)),
                   c("pvalue",formatC(as.numeric(hrtable$pvalue), format = "e", digits = 2)),
                   c("", formatC(Sig)))
tabletext

library(forestplot)
pdf("figures/day3/unicox_forest.pdf", width = 10, height = 12)
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(hrtable$HR)),#log2(HR)
           lower=c(NA,as.numeric(hrtable$HR.95L)), #log2(95%置信区间下限)
           upper=c(NA,as.numeric(hrtable$HR.95H)),#log2(95%置信区间上限)
           graph.pos=7,#图在表中的列位置
           graphwidth = unit(.25,"npc"),#图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),#box颜色
           boxsize=0.4,#box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,#不显示区间
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           xticks = c(0.5, 0.8, 1, 1.2, 1.4, 1.6),#横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab=expression("HR"),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),#第二行顶部加黑实线
                           "35" = gpar(lwd=2, col="black")),#最后一行底部加黑线，""中数字为nrow(tabletext) + 1
           txt_gp=fpTxtGp(label=gpar(cex=1.2),#各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           lineheight = unit(.75,"cm"),#固定行高
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
invisible(dev.off())
