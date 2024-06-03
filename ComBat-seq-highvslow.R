###### ComBat-seq for batch effects (Zhang 2020)
#TGFbeta1 response for pHLFs grown in HighDMEM vs lowDMEM vs Plasmax media
#ComBat-Seq is an improved model based on the ComBat framework, which
#specifically targets RNA-Seq count data. It uses a negative binomial regression
#to model the count matrix, and estimate parameters representing the batch
#effects. Then it provides adjusted data by mapping the original data to an
#expected distribution if there were no batch effects. The adjusted data preserve
#the integer nature of count matrix.

BiocManager::install("sva")
library(sva)

setwd("~/Documents/Chambers/Greg highvslow paper - ComBat-seq Yale2022")

### Prepare the matrix
# Import the data
lowCounts <- read.delim("~/Documents/Chambers/Greg highvslow paper - ComBat-seq Yale2022/low glucose/Combined_counts_24hlow.tabular.txt", header=TRUE, comment.char="#", row.names=1)
highCounts <- read.delim("~/Documents/Chambers/Greg highvslow paper - ComBat-seq Yale2022/high glucose/Combined_counts_24hhigh.tabular.txt", header=TRUE, comment.char="#", row.names=1)
plasmaxCounts <- read.delim("~/Documents/Chambers/Greg highvslow paper - ComBat-seq Yale2022//Combined_counts_24hplasmax.tabular.txt", header=TRUE, comment.char="#", row.names=1)

# Convert to matrix
highCounts <- as.matrix(highCounts)
lowCounts <- as.matrix(lowCounts)
plasmaxCounts <- as.matrix(plasmaxCounts)

# Combine counts into one count matrix
count_matrix <- cbind(lowCounts, highCounts, plasmaxCounts)

### Define the variables for the batch adjustment:
batch <- c(rep(1,6), rep(2,8), rep(3,6))
group <- c(rep(1,3), rep(2,3), rep(1,4), rep(2,4), rep(1,3), rep(1,3))

### Compute the adjusted matrix
adjusted <- ComBat_seq(count_matrix, batch=batch, group=group)

#Save the file:
write.csv(adjusted, file="batch-adjusted-counts.csv")


#### DEseq2
library(DESeq2)
library(dplyr)
library(BiocManager)
library(org.Hs.eg.db)

#Run DEseq2 on the batch corrected counts using INTERACTION:

# Assign condition and media
(media <- factor(c(rep("low",6), rep("high",8), rep("plasmax",6))))
(condition <- factor(c(rep("MC",3), rep("TGF",3), rep("MC",4), rep("TGF",4), rep("MC",3), rep("TGF",3))))

# Create a coldata frame and instantiate the DESeqDataSet and set complex design
(coldata <- data.frame(row.names=colnames(adjusted), condition, media))
dds <- DESeqDataSetFromMatrix(countData=adjusted, colData=coldata, design=~ media + condition + media:condition)
dds

# Relevel so that high is the reference
dds$media <- relevel(dds$media, "high")

# Run the DESeq pipeline
dds <- DESeq(dds)
resultsNames(dds)







#Run DEseq2 on the non-batch corrected counts AND the batch corrected counts (to generate PCA plot):
# Combine non-corrected and corrected counts into one count matrix
combined_matrix <- cbind(count_matrix, adjusted)
#Save the file:
write.csv(combined_matrix, file="combined-counts-nonadj-adj.csv")
#Read back in with correct column names
combined_matrix_labelled <- read.csv("~/Documents/Chambers/Greg highvslow paper - ComBat-seq Yale2022//combined-counts-nonadj-adj.csv", header=TRUE, comment.char="#", row.names=1)

# Convert to matrix
combined_matrix <- as.matrix(combined_matrix_labelled)

#Set factors
(media <- factor(c(rep("low",6), rep("high",8), rep("plasmax",6), rep("adj.low",6), rep("adj.high",8), rep("adj.plasmax",6))))
(condition <- factor(c(rep("MC",3), rep("TGF",3), rep("MC",4), rep("TGF",4), rep("MC",3), rep("TGF",3),rep("MC",3), rep("TGF",3), rep("MC",4), rep("TGF",4), rep("MC",3), rep("TGF",3))))

# Create a coldata frame and instantiate the DESeqDataSet and set complex design
(coldata <- data.frame(row.names=colnames(combined_matrix), condition, media))
dds <- DESeqDataSetFromMatrix(countData=combined_matrix, colData=coldata, design=~ media + condition + media:condition)
dds

# Relevel so that high is the reference
dds$media <- relevel(dds$media, "high")

# Run the DESeq pipeline
dds <- DESeq(dds)
resultsNames(dds)





#### Get the differentially expressed gene lists:
###### Note: Haven't got the interaction working properly yet!

# Choose one:

##TGF effect
# Effect of TGF treatment in DMEMlow
res <- results(dds, list( c("condition_TGF_vs_MC","medialow.conditionTGF") ))
# Difference between DMEMlow and DMEMhigh TGF response
res <- results(dds, name="medialow.conditionTGF")

# Effect of TGF treatment in DMEMhigh
res <- results(dds, contrast=c("condition","TGF","MC"))

#Effect of TGF treatment in Plasmax
res <- results(dds, list( c("condition_TGF_vs_MC","mediaplasmax.conditionTGF") ))
# Difference between DMEMlow and Plasmax TGF response
res <- results(dds, name="mediaplasmax.conditionTGF")

##Baseline
# Difference between highDMEM and lowDMEM without treatment
res <- results(dds, contrast=c("media","high","low"))
# Difference between highDMEM and Plasmax without treatment
res <- results(dds, contrast=c("media","high","plasmax"))
#Difference between lowDMEM and Plasmax without treatment
res <- results(dds, contrast=c("media","low","plasmax"))



#Same analysis here for all scenarios:

# Merge with normalized count data
resdata<- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

## Annotation with gene symbols and filtering
#Create new data frame with annotation information
anno <- AnnotationDbi::select(org.Hs.eg.db,keys=resdata$Gene,
                              columns=c("SYMBOL"),
                              keytype="ENTREZID")

#Rename first column in resdata file from Gene to Ensembl
resdata.labelled <- resdata
colnames(resdata.labelled)
names(resdata.labelled)[names(resdata.labelled) == "Gene"] <- "ENTREZID"

#Bind the annotation information to the results data frame
results.annotated <- left_join(resdata.labelled, anno,by="ENTREZID")

#Change the order of the columns (move last column to start)
results.annotated <- results.annotated[,c(28,1:27)]

# Remove NAs
results.annotated <- na.omit(results.annotated)

#Filter for abs(Log2FoldChange)>=0. and padj<0.05
results.annotated.p0.05 <- filter(results.annotated, padj < 0.05)
results.annotated.p0.05.fc0.58 <- filter(results.annotated.p0.05, abs(log2FoldChange) >= 0.58)



## Write results to one of these files:

#lowDMEM
write.csv(results.annotated.p0.05.fc0.58, file="adjusted_lowDMEM_TGF_response.csv")
write.csv(results.annotated.p0.05.fc0.58, file="adjusted_lowvshighDMEM_TGF_interaction.csv")

#highDMEM
write.csv(results.annotated.p0.05.fc0.58, file="adjusted_highDMEM_TGF_response.csv")

#Plasmax
write.csv(results.annotated.p0.05.fc0.58, file="adjusted_Plasmax_TGF_response.csv")
write.csv(results.annotated.p0.05.fc0.58, file="adjusted_plasmaxvshighDMEM_TGF_interaction.csv")

#baseline
write.csv(results.annotated.p0.05.fc0.58, file="adjusted_baseline_highDMEMvslowDMEM.csv")
write.csv(results.annotated.p0.05.fc0.58, file="adjusted_baseline_highDMEMvsPlasmax.csv")
write.csv(results.annotated.p0.05.fc0.58, file="adjusted_baseline_lowDMEMvsPlasmax.csv")


#### Pathways analysis using ClusterProfiler
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggplot2)

#Load gene list
genes <- read.csv("adjusted_lowDMEM_TGF_response.csv", header=FALSE)
genes <- read.csv("adjusted_highDMEM_TGF_response.csv", header=FALSE)
genes <- read.csv("adjusted_Plasmax_TGF_response.csv", header=FALSE)

genes <- read.csv("adjusted_lowvshighDMEM_TGF_interaction.csv", header=FALSE)
genes <- read.csv("adjusted_plasmaxvshighDMEM_TGF_interaction.csv", header=FALSE)

de <- as.character(genes[,1])
head(de)

#Pathway enrichment using ReactomePA
x <- enrichPathway(gene=de, pvalueCutoff=0.05, readable=TRUE, minGSSize = 10, maxGSSize = 500)
head(x)

#Save analysis file
write.csv(x, file="pathways_results_lowDMEM.csv")
write.csv(x, file="pathways_results_highDMEM.csv")
write.csv(x, file="pathways_results_Plasmax.csv")

write.csv(x, file="pathways_results_lowvshighDMEM_TGF_interaction.csv")
write.csv(x, file="pathways_results_plasmaxvshighDMEM_TGF_interaction.csv")


#Custom theme for dot plot for publication
mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
          axis.text=element_text(face="bold")
    ))

#Visualise using dot plot
dotplot(x, showCategory=30) +
  ylab("Reactome term") +
  #scale_color_gradient("Adjusted p-value", low = "maroon1", high = "maroon4") +
  scale_color_gradient("Adjusted p-value", low = "steelblue1", high = "steelblue4") +
  xlim(0,0.1) +
  scale_y_discrete(label = function(y) {
    y %>% sub("Homo sapiens\r: ", "", .)
  }) +
  mytheme




#### Visualisation:

## Load libraries required
library(ggplot2)
library("gplots")
library(pheatmap)

#Custom theme for PCA plot for publication
mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", axis.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
          axis.text=element_text(face="bold")
    ))

#using rlog transformed data
rld <- rlog(dds)
pcaData <- plotPCA(rld, intgroup = c("media","condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#Plot PCA comparing before and after batch correction
ggplot(pcaData, aes(PC1, PC2, shape=condition, color=media)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  
  #customisable shapes and color labels
  #scale_color_manual(labels = c("Control", "AZD8055", "Rapamycin"),
                     #values = c("darkgoldenrod1", "steelblue1", "maroon3")) +
  labs(color="Treatment") +
  scale_shape_manual(labels = c("Media Control", "TGF-β1"),
                     values = c(16, 17)) +
  labs(shape="Media") +
  geom_point(size=2) +
  mytheme +
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2))






