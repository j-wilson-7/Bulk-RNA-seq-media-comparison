#### Plotting PCA from DEseq2 results (based on plotPCA function) and getting genes for each PC
#### to plot on a heatmap.
#TGFbeta1 response for pHLFs grown in HighDMEM vs lowDMEM vs Plasmax media

## Load libraries required
library(enrichplot)

library(enrichplot)
library(ReactomePA)
library(DOSE)
library(dplyr)
library(clusterProfiler)

library(ggplot2)
library("gplots")
library(pheatmap)
library(DESeq2)
library(BiocManager)
library(tibble)
library(org.Hs.eg.db)
library(msigdbr)


#### Carry out DEseq2 analysis using n of 3

setwd("~/Documents/Chambers/Greg highvslow paper - Analysis2 Yale2022/plotPCA")

### Prepare the matrix
# Import the data
lowCounts <- read.delim("~/Documents/Chambers/Greg highvslow paper - Analysis2 Yale2022/Counts/low glucose/Combined_counts_24hlow.tabular.txt", header=TRUE, comment.char="#", row.names=1)
highCounts <- read.delim("~/Documents/Chambers/Greg highvslow paper - Analysis2 Yale2022/Counts/high glucose/Combined_counts_24hhigh_n3.tabular.txt", header=TRUE, comment.char="#", row.names=1)
plasmaxCounts <- read.delim("~/Documents/Chambers/Greg highvslow paper - Analysis2 Yale2022/Counts/Combined_counts_24hplasmax.tabular.txt", header=TRUE, comment.char="#", row.names=1)

# Convert to matrix
highCounts <- as.matrix(highCounts)
lowCounts <- as.matrix(lowCounts)
plasmaxCounts <- as.matrix(plasmaxCounts)

# Combine counts into one count matrix
count_matrix <- cbind(lowCounts, highCounts, plasmaxCounts)

(media <- factor(c(rep("low",6), rep("high",6), rep("plasmax",6))))
(condition <- factor(c(rep("MC",3), rep("TGF",3), rep("MC",3), rep("TGF",3), rep("MC",3), rep("TGF",3))))

# Create a coldata frame and instantiate the DESeqDataSet and set complex design
(coldata <- data.frame(row.names=colnames(count_matrix), condition, media))
dds <- DESeqDataSetFromMatrix(countData=count_matrix, colData=coldata, design=~ media + condition + media:condition)
dds

# Relevel so that high is the reference
dds$media <- relevel(dds$media, "high")

# Run the DESeq pipeline
dds <- DESeq(dds)
resultsNames(dds)



#### PCA calculation

#using rlog transformed data
rld <- rlog(dds)
pcaData <- plotPCA(rld, intgroup = c("media","condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

object=rld
intgroup=c("media","condition")
ntop=500

  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assemble the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(object))
  

  
  #### Plot the PCA graph 
  
  #Custom theme for PCA plot for publication
  mytheme = list(
    theme_classic()+
      theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", axis.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
            axis.text=element_text(face="bold")
      ))
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="media", shape="condition", group ="media")) +  
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed() +
    labs(color="Treatment") +
    scale_shape_manual(labels = c("Media Control", "TGF-β1"),
                       values = c(16, 17)) +
    labs(shape="Media") +
    geom_point(size=2) +
    mytheme +
    guides(colour = guide_legend(order = 1), 
           shape = guide_legend(order = 2)) +
    scale_color_manual(labels = c("DMEMhigh", "DMEMlow" 
                                  #,"Plasmax"
                                  ),
    values = c("darkgoldenrod1", "steelblue1"
               #,"maroon3"
               )) +
    stat_ellipse(linetype=3)
    

 
library(ggforce)


  #### Get the top contributing PCA genes
  
  #Get the PC data
  loadings <- as.data.frame(pca$rotation)
  
  ## Annotation with gene symbols
  
  #Change row names into column
  loadings <- tibble::rownames_to_column(loadings, "ENTREZID")
  #loadings$ENTREZID <- rownames(loadings)
  
  #Create new data frame with annotation information
  anno <- AnnotationDbi::select(org.Hs.eg.db,keys=loadings$ENTREZID,
                                columns=c("SYMBOL"),
                                keytype="ENTREZID")
  
  #Bind the annotation information to the results data frame
  loadings.anno <- left_join(loadings, anno,by="ENTREZID")
  
  # Remove NAs
  loadings.anno <- na.omit(loadings.anno)
  
  #Change the order of the columns (move last column to start)
  loadings.anno <- loadings.anno[,c(20,1:19)]

  #Save as csv file
  write.csv(loadings.anno, file="PCAresults-high-low-plasmax-2.csv")  
  
  
## GSEA on PC1 and PC2
  
  #Take PC1 and PC2 and the gene labels
  absPC <- loadings.anno[,c("ENTREZID","PC1", "PC2")]
  #rownames(absPC) <- loadings.anno[,1]
  
  #Get the absolute values
  #absPC$PC1 <- abs(absPC$PC1)
  #absPC$PC2 <- abs(absPC$PC2)
  
  #Get into right format for ClusterProfiler - named vector
  PC1 <- absPC[,2]
  PC2 <- absPC[,3]
  
  names(PC1) <- as.character(absPC[,1])
  names(PC2) <- as.character(absPC[,1])
  
  #Order eigenvalues from high to low
  PC1.ordered <- sort(PC1, decreasing = TRUE)
  PC2.ordered <- sort(PC2, decreasing = TRUE)
  
  #Prepare the MSigDB database for GSEA
  msig <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, entrez_gene)
  head(msig)
  
  #GSEA using ClusterProfiler with MSigDB
  #Set scoreType = pos if using absolute values, otherwise remove
  GSEAresultsPC1 <- GSEA(PC1.ordered, 
                         TERM2GENE = msig, 
                         #scoreType = "pos"
                         )
  head(GSEAresultsPC1)
  
  GSEAresultsPC2 <- GSEA(PC2.ordered, 
                         TERM2GENE = msig, 
                         #scoreType = "pos"
                         )
  head(GSEAresultsPC2)
  
  #Save the results
  write.csv(GSEAresultsPC1, file="PCA-GSEA-PC1.csv")
  write.csv(GSEAresultsPC2, file="PCA-GSEA-PC2.csv")
  
 
  ## KEGG GSEA
  KEGGresultsPC1 <- gseKEGG(geneList = PC1.ordered,
                 organism     = 'hsa',
                 minGSSize    = 1,
                 pvalueCutoff = 0.05,
                 #scoreType = "pos"
                 )

  KEGGresultsPC2 <- gseKEGG(geneList = PC2.ordered,
                            organism     = 'hsa',
                            minGSSize    = 1,
                            pvalueCutoff = 0.05,
                            #scoreType = "pos"
                            )
  
  
  ## Reactome GSEA
  ReactomeresultsPC1 <- gsePathway(PC1.ordered,
                                   organism = "human",
                                   minGSSize = 10)
  
  ReactomeresultsPC2 <- gsePathway(PC2.ordered,
                                   organism = "human",
                                   minGSSize = 10)
  
  write.csv(ReactomeresultsPC1, file="PCA-GSEA-Reactome-PC1.csv")
  write.csv(ReactomeresultsPC2, file="PCA-GSEA-Reactome-PC2.csv")
  
  #Remove the homo sapiens label from Reactome labels
  ReactomeresultsPC1@result$Description <- gsub("Homo sapiens\r: ","", as.character(ReactomeresultsPC1@result$Description))
  ReactomeresultsPC2@result$Description <- gsub("Homo sapiens\r: ","", as.character(ReactomeresultsPC2@result$Description))
  
  #Get the absolute normalised enrichment score for PC2 because there is one negative one
  ReactomeresultsPC2@result$NES <- abs(ReactomeresultsPC2@result$NES)
  
  #Prepare publication theme
  mytheme = list(
    theme_classic()+
      theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
            axis.text=element_text(face="bold")
      ))
  
  #Plot using barplot (PC1)
  ggplot(ReactomeresultsPC1[1:9]) +
    mytheme +
    aes(x = reorder(Description, NES), y = NES) +
    coord_flip() +
    geom_col(position = "dodge", fill="burlywood1") +
    xlab("") +
    ylab("Absolute NES") +
    coord_flip() +
    ylim(0, 3)
  
 #Plot using barplot (PC2) with upregulated and downregulated marked
  ggplot(ReactomeresultsPC2) +
    mytheme +
    aes(x = reorder(Description, NES), y = NES, fill = NES > 2.14) +
    coord_flip() +
    geom_col(position = "dodge") +
    scale_fill_manual(labels = c("Upregulated", "Downregulated"), values = c("#f6b293","#9ed4f0")) +
    xlab("") +
    ylab("Absolute NES") +
    labs(fill = "") +
    coord_flip() +
    ylim(0, 3)
  
  #Plot using barplot (PC2) without all bars same colour
  ggplot(ReactomeresultsPC2[1:4]) +
    mytheme +
    aes(x = reorder(Description, NES), y = NES) +
    coord_flip() +
    geom_col(position = "dodge", fill="burlywood1") +
    xlab("") +
    ylab("Absolute NES") +
    coord_flip() +
    ylim(0, 3)
  
  
### Plotting a 3D PCA 

library(plot3D)
  
# assemble the data for the plot
d3 <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3], group=group, intgroup.df, name=colnames(object))

#plot
scatter3D(d3$PC1, d3$PC2, d3$PC3,
          clab = c(""),
          pch = 19,
          cex = 0.5,
          bty = 'f',
          colkey = FALSE
          )
  
  
  
  
  