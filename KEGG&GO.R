library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggh4x)
library(clusterProfiler)
load("data/WGCNA_gene.Rdata")

#Yes_Epilepsy
#KEGG
### 获得基因列表
gene <- Yes_gene
#基因名称转换，返回的是数据框
gene = bitr(gene, 
            fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")
head(gene)

EGG <- enrichKEGG(gene = gene$ENTREZID,
                  organism = 'hsa',
                  pvalueCutoff = 0.05)
kegg <- as.data.frame(EGG)

GO <-enrichGO(gene = gene$ENTREZID,
              'org.Hs.eg.db',
              ont = 'ALL',
              pvalueCutoff = 0.05)

go <- as.data.frame(GO)

df <- kegg %>% 
  mutate('Enrich factor' = Count/31*100) %>% 
  mutate('negative_log10_of_adjusted_p_value' = -log10(p.adjust)) %>%
  mutate('type' = "KEGG") %>%
  arrange(type,negative_log10_of_adjusted_p_value)


# order
df$term_name <- factor(df$Description,levels = unique(df$Description))
# plot save 8x15

ggplot(df,aes(y = negative_log10_of_adjusted_p_value,x = term_name)) +
  geom_col(aes(fill = type),width = 0.5,color = NA) +
  scale_fill_brewer(palette = 'Set1',name = 'Go type') +
  geom_line(aes(y = `Enrich factor`,group = type),
            color = 'grey75',
            size = 1) +
  geom_point(aes(y = `Enrich factor`),size = 5,
             color = '#FF9900') +
  scale_y_continuous(sec.axis = sec_axis(~.*0.01,name = 'Gene Ratio',
                                         labels = scales::label_percent())) +
  theme_bw() +
  #facet_wrap(~sample,scales = 'free') +
  xlab('') + ylab('-log10 adjusted Pvalue') +
  theme(axis.text = element_text(size = rel(1.2),color = 'black'),
        axis.title = element_text(size = rel(1.5),color = 'black'),
        strip.text = element_text(size = rel(1.2),color = 'black'),
        strip.background = element_rect(color = NA),
        strip.placement = 'outside') +
  coord_flip()
ggsave(file = "figures/day4/WGCNA_Yes_KEGG.pdf", width = 8, height = 6)



#GO
bp <- filter(go, ONTOLOGY == "BP") %>% mutate('Enrich factor' = Count/93*100) %>% 
  mutate('negative_log10_of_adjusted_p_value' = -log10(p.adjust)) %>%
  arrange(desc(negative_log10_of_adjusted_p_value))

cc <- filter(go, ONTOLOGY == "CC") %>% mutate('Enrich factor' = Count/101*100) %>% 
  mutate('negative_log10_of_adjusted_p_value' = -log10(p.adjust)) %>%
  arrange(desc(negative_log10_of_adjusted_p_value)) 

mf <- filter(go, ONTOLOGY == "MF") %>% mutate('Enrich factor' = Count/94*100) %>% 
  mutate('negative_log10_of_adjusted_p_value' = -log10(p.adjust)) %>%
  arrange(desc(negative_log10_of_adjusted_p_value))

df <- rbind(bp, cc, mf)
# order
df$term_name <- factor(df$Description,levels = unique(df$Description))

# plot save 8x15

ggplot(df,aes(y = negative_log10_of_adjusted_p_value,x = term_name)) +
  geom_col(aes(fill = ONTOLOGY),width = 0.5,color = NA) +
  scale_fill_brewer(palette = 'Set1',name = 'Go type') +
  geom_line(aes(y = `Enrich factor`,group = ONTOLOGY),
            color = 'grey75',
            size = 1) +
  geom_point(aes(y = `Enrich factor`),size = 5,
             color = '#FF9900') +
  scale_y_continuous(sec.axis = sec_axis(~.*0.01,name = 'Gene Ratio',
                                         labels = scales::label_percent())) +
  theme_bw() +
  # facet_wrap(~sample,scales = 'free') + #这里是分面
  xlab('') + ylab('-log10 adjusted Pvalue') +
  theme(axis.text = element_text(size = rel(1.2),color = 'black'),
        axis.title = element_text(size = rel(1.5),color = 'black'),
        strip.text = element_text(size = rel(1.2),color = 'black'),
        strip.background = element_rect(color = NA),
        strip.placement = 'outside') +
  coord_flip()
dev.copy2pdf(file = "figures/day4/WGCNA_Yes_GO.pdf")
