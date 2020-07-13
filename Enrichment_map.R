### LOAD LIBRARIES ###
library(tidyverse)
library(XGR)


### FIND SIGNIFICANT PATHWAYS (TRANSCRIPTOMIC) FOR EACH COMPARISON

load("total_genes.RData")

CTL_AD_GENES <- CTL_AD_FINAL[["GENES"]]
CTL_MCI_GENES <- CTL_MCI_FINAL[["GENES"]]
MCI_AD_GENES <- MCI_AD_FINAL[["GENES"]]

eTerm_CTL_AD_GENES <- vector("list")
eTerm_CTL_MCI_GENES <- vector("list")
eTerm_MCI_AD_GENES <- vector("list")

datalist <- c("MsigdbC2KEGG","MsigdbC5BP","MsigdbC5MF", "MsigdbC5CC")

for(i in datalist) {
  try({eTerm_CTL_AD_GENES <- xEnricherGenes(data=CTL_AD_GENES, ontology="MsigdbC5BP", test = "hypergeo", ontology.algorithm="none", true.path.rule=F, RData.location="http://galahad.well.ox.ac.uk/bigdata")
  eTerm_CTL_MCI_GENES[[i]] <- xEnricherGenes(data=CTL_MCI_GENES, ontology=i, test = "hypergeo", ontology.algorithm="none", true.path.rule=F, RData.location="http://galahad.well.ox.ac.uk/bigdata")
  eTerm_MCI_AD_GENES[[i]] <- xEnricherGenes(data=MCI_AD_GENES, ontology=i, test = "hypergeo", ontology.algorithm="none", true.path.rule=F, RData.location="http://galahad.well.ox.ac.uk/bigdata")
  print(i)}, silent = TRUE)}

## EXTRACT SIGNIFICANT PATHWAYS

pathways <- xEnrichViewer(
  eTerm_CTL_AD_GENES,
  sortBy = c("fdr"),
  top_num = "auto",
  decreasing = NULL,
  details = T)

## bar_ploting significant pathways


p <-ggplot(a, aes(x=reorder(PATHWAYS, logfdr), y=logfdr, label = P.VALUE)) + 
  geom_bar(stat="identity", fill = "steelblue") +
  geom_text(size = 2, position = position_stack(vjust=1.2))+
  xlab("Pathways") +
  ylab("FDR") +
  ggtitle("Molecular Functions: -log10(FDR)") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
p + theme(axis.text.y = element_text(hjust = 0.4, size = 6),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
) + coord_flip()






