#### Plotting pathways analysis results using ClusterProfiler/ReactomePA
#TGFbeta1 response for pHLFs grown in HighDMEM vs lowDMEM vs Plasmax media

library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(ggupset)



## 1) KEGG pathway over-representation analysis (unranked list)
#Load gene lists
setwd("~/Documents/Chambers/Greg highvslow paper - Analysis2 Yale2022/Final-deseq2-analysis")
genes <- read.delim("lowDMEM-p0.05-fc1.5-symbols.txt", header=FALSE)
genes <- read.delim("highDMEM-p0.05-fc1.5-symbols.txt", header=FALSE)
genes <- read.delim("plasmax-p0.05-fc1.5-symbols.txt", header=FALSE)

#Convert gene symbols to entrez ID
anno <- AnnotationDbi::select(org.Hs.eg.db,
                              keys=genes$V1,
                              columns=c("ENTREZID"),
                              keytype="SYMBOL")
anno <- na.omit(anno)
de <- as.character(anno[,2])
head(de)

#Pathways analysis
x <- enrichKEGG(gene=de, pvalueCutoff=0.05, minGSSize = 10, organism='hsa')
head(x)

#Save analysis files
write.csv(x, file="lowDMEM-clusterprofiler-KEGG-unranked.csv")
write.csv(x, file="highDMEM-clusterprofiler-KEGG-unranked.csv")
write.csv(x, file="plasmax-clusterprofiler-KEGG-unranked.csv")

#Custom theme for publications
mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
          axis.text=element_text(face="bold")
    ))

#Visualise using dot plot
dotplot(x, showCategory=20) +
  ylab("KEGG term") +
  scale_color_gradient("Adjusted p-value", low="maroon1", high="maroon4") +
  mytheme +
  xlim(0,0.8)


#In case you want to do Reactome:
x <- enrichPathway(gene=de, pvalueCutoff=0.05, minGSSize = 10, readable=TRUE)
head(x)

#Visualise using dot plot
dotplot(x, showCategory=15) +
  ylab("") +
  scale_color_gradient("Adjusted p-value", low="maroon1", high="maroon4") +
  mytheme +
  xlim(0,0.08) +
  scale_y_discrete(label = function(y) {
    y %>% sub("Homo sapiens\r: ", "", .)
  })


## 2) KEGG pathway GSEA (ranked list)
setwd("~/Documents/Chambers/Greg highvslow paper - Analysis2 Yale2022/Final-deseq2-analysis")

genes <- read.table("Low DMEM ALL.tabular", header=TRUE)
genes <- read.table("High DMEM ALL.tabular", header=TRUE)
genes <- read.table("Plasmax ALL.tabular", header=TRUE)

#Extract log2foldchange as a numeric vector
geneList <- genes[,4]

#Add the entrezIDs as names in the vector
names(geneList) <- as.character(genes[,1])

#Order by log2Foldchange
geneList = sort(geneList, decreasing = TRUE)

#Run GSEA on kegg terms
y <- gseKEGG(geneList     = geneList,
             organism     = 'hsa',
             minGSSize    = 10,
             pvalueCutoff = 0.05,
             verbose      = FALSE)
head(y)

#Save analysis files
write.csv(y, file="lowDMEM-clusterprofiler-KEGG-GSEA.csv")
write.csv(y, file="highDMEM-clusterprofiler-KEGG-GSEA.csv")
write.csv(y, file="plasmax-clusterprofiler-KEGG-GSEA.csv")

#Custom theme for publications
mytheme2 = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
          axis.text=element_text(face="bold"), legend.position = "bottom", legend.box = "vertical"
    ))



#Plot using dotplot (DMEMlow)
dotplot(y, showCategory=20, x="NES", decreasing = TRUE, label_format = function(x) stringr::str_wrap(x, width=70)) + 
  ylab("") +
  scale_colour_gradient("Adjusted p-value", low="#66D2D6", high="#3c7d80", 
                        breaks=c(0.001,0.004,0.007)
  ) +
  scale_x_continuous(limit = c(-3, 3), breaks=seq(-3,3,1)) +
  geom_vline(xintercept = 0, linetype="dotted", color = "#444444", size=1) +
  mytheme2



#Plot using dotplot (DMEMhigh)
dotplot(y, showCategory=20, x="NES", decreasing = TRUE, label_format = function(x) stringr::str_wrap(x, width=60)) + 
  ylab("") +
  scale_colour_gradient("Adjusted p-value", low="#FBC740", high="#9e7a1e", 
                        breaks=c(0.001,0.004,0.007)
  ) +
  scale_x_continuous(limit = c(-3, 3), breaks=seq(-3,3,1)) +
  geom_vline(xintercept = 0, linetype="dotted", color = "#444444", size=1) +
  mytheme2



#Plot using dotplot (plasmax)
dotplot(y, showCategory=20, x="NES", decreasing = TRUE, label_format = function(x) stringr::str_wrap(x, width=150)) + 
  ylab("") +
  scale_colour_gradient("Adjusted p-value", low="#E56997", high="#7d3852", 
                        breaks=c(0.001,0.0055,0.01)
                        ) +
  scale_x_continuous(limit = c(-3, 3), breaks=seq(-3,3,1)) +
  geom_vline(xintercept = 0, linetype="dotted", color = "#444444", size=1) +
  mytheme2



#Visualisation of dotplots with Reactome:
#Run GSEA on Reactome terms
y <- gsePathway(geneList     = geneList,
                minGSSize    = 10,
                pvalueCutoff = 0.05,
                verbose      = FALSE)
head(y)

#Save analysis files
write.csv(y, file="lowDMEM-clusterprofiler-Reactome-GSEA.csv")
write.csv(y, file="highDMEM-clusterprofiler-Reactome-GSEA.csv")
write.csv(y, file="plasmax-clusterprofiler-Reactome-GSEA.csv")

#Plot using dotplot (DMEMlow) REACTOME
dotplot(y, showCategory=25, x="NES", decreasing = TRUE, label_format = function(x) stringr::str_wrap(x, width=150)) + 
  ylab("") +
  scale_colour_gradient("Adjusted p-value", low="#66D2D6", high="#3c7d80", 
                        breaks=c(0.0001,0.0006,0.0012)
  ) +
  scale_x_continuous(limit = c(-3, 3), breaks=seq(-3,3,1)) +
  geom_vline(xintercept = 0, linetype="dotted", color = "#444444", size=1) +
  mytheme2

#Plot using dotplot (DMEMhigh) REACTOME
dotplot(y, showCategory=25, x="NES", decreasing = TRUE, label_format = function(x) stringr::str_wrap(x, width=150)) + 
  ylab("") +
  scale_colour_gradient("Adjusted p-value", low="#FBC740", high="#9e7a1e", 
                        breaks=c(0.0005,0.002)
  ) +
  scale_x_continuous(limit = c(-3, 3), breaks=seq(-3,3,1)) +
  geom_vline(xintercept = 0, linetype="dotted", color = "#444444", size=1) +
  mytheme2

#Plot using dotplot (Plasmax) REACTOME
dotplot(y, showCategory=50, 
        x="NES", 
        #decreasing = TRUE, 
        label_format = function(x) stringr::str_wrap(x, width=150)) +
  ylab("") +
  scale_colour_gradient("Adjusted p-value", low="#E56997", high="#7d3852", 
                        breaks=c(0.00002,0.00009, 0.00016)
  ) +
  scale_x_continuous(limit = c(-3, 3), breaks=seq(-3,3,1)) +
  geom_vline(xintercept = 0, linetype="dotted", color = "#444444", size=1) +
  mytheme2
  

## 3) Looking at the intersection genes (those in all three medias)
setwd("~/Documents/Chambers/Greg highvslow paper - Analysis2 Yale2022/Final figures-Mar23-NOPlasmax")
genes <- read.delim("intersection-genes", header=FALSE)

#Convert gene symbols to entrez ID
anno <- AnnotationDbi::select(org.Hs.eg.db,
                              keys=genes$V1,
                              columns=c("ENTREZID"),
                              keytype="SYMBOL")
anno <- na.omit(anno)
de <- as.character(anno[,2])
head(de)

#Pathways analysis
#z <- enrichKEGG(gene=de, pvalueCutoff=0.05, minGSSize = 10, organism='hsa')
z <- enrichPathway(gene=de, pvalueCutoff=0.05, minGSSize = 10)

#Add gene symbol column
z.symbols <- setReadable(z, 'org.Hs.eg.db', 'ENTREZID')

#Save analysis file
write.csv(z, file="Reactome-unranked-intersection-noPlasmax.csv")
write.csv(z.symbols, file="Reactome-unranked-intersection-noPlasmax-symbols.csv")




#Visualise networks using cnetplot
## Top 5 pathways

#Remove the homo sapiens label from Reactome labels
z.symbols@result$Description <- gsub("Homo sapiens\r: ","", as.character(z.symbols@result$Description))

#cnet plot
cnetplot(z.symbols, showCategory = 5, foldChange = NULL, 
         node_label = "category",
         circular = TRUE, 
         colorEdge = TRUE,
         cex_label_category = 0.5,
         #cex_label_gene = 0.4
         cex_category = 1.5,
         shadowtext = "none"
         )

#p + theme(legend.position = "none")
   
#p + guides(size = guide_legend("Size"))

# Choose pathways 
#Choose the 10th, 13th and 14th pathway
selected_pathways <- z.symbols$Description[c(10, 13, 14, 8, 17)]
cnetplot(z.symbols, showCategory = selected_pathways, foldChange = NULL, 
         #circular = TRUE, 
         colorEdge = TRUE) 

#Other options for plots:
upsetplot(z.symbols)

edox2 <- pairwise_termsim(z.symbols)
treeplot(edox2)





