rm(list = ls())
# 1. Install Packages
library(devtools)
install_local("C:/Users/runru/Desktop/seq/GSEA_enrichment/GSEA_enrichment/Xinpeng-Zhou-msigdf-master.zip")
# 2. 载入R包
library(ggplot2)
library(clusterProfiler)
library(tidyverse)
library(msigdf)
library(org.Hs.eg.db)
library(enrichplot)
# 3. 数据整理
all_gene <- read_csv("C:/Users/runru/Desktop/seq/GSEA_enrichment/GSEA_enrichment/RA-vs-OA-all_gene.csv", na = c("", "#NAME?", "Inf")) %>% dplyr::select(gene_id, log2FoldChange) %>% filter(!is.na(log2FoldChange))
# FC_4 <- read_csv("./data/RA-vs-OA-diff-pval-0.05-FC-4_gene.csv", na = c("","#NAME?"))
# FC_2 <- read_csv("./data/RA-vs-OA-diff-pval-0.05-FC-2_gene.csv", na = c("","#NAME?"))
id <- bitr(all_gene$gene_id, "SYMBOL", "ENTREZID", "org.Hs.eg.db") %>% mutate(gene_id = SYMBOL) %>% dplyr::select(gene_id, ENTREZID) %>% inner_join(all_gene) %>% dplyr::select(ENTREZID, log2FoldChange) %>% arrange(desc(log2FoldChange))
genelist <- id$log2FoldChange
names(genelist) <- id$ENTREZID
write.table(all_gene,"t1.csv",row.names=FALSE,col.names=TRUE,sep=",")
# 4. GSEA及作图
## Molecular Signatures Database (MSigDB)
## C1 positional gene sets
c1 <- msigdf.human %>% filter(collection == "c1") %>% dplyr::select(geneset, entrez) %>% as.data.frame()
gsea_c1 <- GSEA(genelist, TERM2GENE = c1, pvalueCutoff = 0.2)
gseaplot(gsea_c1, 6, title = gsea_c1@result[["ID"]][6])
gseaplot2(gsea_c1, 6, title = gsea_c1@result[["ID"]][6])
## C2 curated gene sets
c2 <- msigdf.human %>% filter(collection == "c2") %>% dplyr::select(geneset, entrez) %>% as.data.frame()
gsea_c2 <- GSEA(genelist, TERM2GENE = c2)
ridgeplot(gsea_c2, showCategory = 45, label_format = 30)
gseaplot(gsea_c2, 4, title = gsea_c2@result[["ID"]][4])
gseaplot2(gsea_c2, 4, title = gsea_c2@result[["ID"]][4])

## C5 GO gene set
c5 <- msigdf.human %>% filter(collection == "c5") %>% dplyr::select(geneset, entrez) %>% as.data.frame()
gsea_c5 <- GSEA(genelist, TERM2GENE = c5)
ridgeplot(gsea_c5, showCategory = 5, label_format = 30)
gseaplot(gsea_c5, 2, title = gsea_c5@result[["ID"]][2])
gseaplot2(gsea_c5, 2, title = gsea_c5@result[["ID"]][2])

## C7: immunologic signatures
c7 <- msigdf.human %>% filter(collection == "c7") %>% dplyr::select(geneset, entrez) %>% as.data.frame()
gsea_c7 <- GSEA(genelist, TERM2GENE = c7)
ridgeplot(gsea_c7, showCategory = 5, label_format = 30)
gseaplot(gsea_c7, 2, title = gsea_c7@result[["ID"]][2])
gseaplot2(gsea_c7, 2, title = gsea_c7@result[["ID"]][2])
