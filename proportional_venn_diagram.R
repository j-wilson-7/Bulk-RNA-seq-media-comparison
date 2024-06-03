#### Producing proportional and non-proportional venn diagrams
#TGFbeta1 response for pHLFs grown in HighDMEM vs lowDMEM vs Plasmax media

### Prepare the data 

#Get the list of filtered up and downregulated genes for the venn diagram
setwd("~/Documents/Chambers/Greg highvslow paper - Analysis2 Yale2022/Final-deseq2-analysis")
low.genes <- read.table("~/Documents/Chambers/Greg highvslow paper - Analysis2 Yale2022/Final-deseq2-analysis/Low DMEM ALL.tabular", header=TRUE, comment.char="#", row.names=1)
high.genes <- read.table("~/Documents/Chambers/Greg highvslow paper - Analysis2 Yale2022/Final-deseq2-analysis/High DMEM ALL.tabular", header=TRUE, comment.char="#", row.names=1)
plasmax.genes <- read.table("~/Documents/Chambers/Greg highvslow paper - Analysis2 Yale2022/Final-deseq2-analysis/Plasmax ALL.tabular", header=TRUE, comment.char="#", row.names=1)

# lowDMEM
low.genes <- na.omit(low.genes)
low.genes.p0.05 <- filter(low.genes, X7.P.adj < 0.05)
low.genes.p0.05.fc0.58 <- filter(low.genes.p0.05, abs(X3.log2.FC.) >= 0.58)
write.csv(low.genes.p0.05.fc0.58, file="lowDMEM.fc1.5.padj0.05.csv")

# highDMEM
high.genes <- na.omit(high.genes)
high.genes.p0.05 <- filter(high.genes, X7.P.adj < 0.05)
high.genes.p0.05.fc0.58 <- filter(high.genes.p0.05, abs(X3.log2.FC.) >= 0.58)
write.csv(high.genes.p0.05.fc0.58, file="highDMEM.fc1.5.padj0.05.csv")

# plasmax
plasmax.genes <- na.omit(plasmax.genes)
plasmax.genes.p0.05 <- filter(plasmax.genes, X7.P.adj < 0.05)
plasmax.genes.p0.05.fc0.58 <- filter(plasmax.genes.p0.05, abs(X3.log2.FC.) >= 0.58)
write.csv(plasmax.genes.p0.05.fc0.58, file="plasmax.fc1.5.padj0.05.csv")




## Visualise using VennDiagram (non-proportional)
library(VennDiagram) 

# Triple venn
grid.newpage()                  # Move to new plotting page
draw.triple.venn(area1 = 1505,    # lowDMEM
                 area2 = 4265,    # highDMEM
                 area3 = 2654,    # Plasmax
                 n12 = 1279,      #intersection including n123 also
                 n23 = 1771,
                 n13 = 792,
                 n123 = 726,
                 euler.d = TRUE,
                 scaled = TRUE,
                 fill = c("lightblue","green","yellow"))   #change colours here

# Double venn (high vs low DMEM)
grid.newpage()                    # Create new plotting page
draw.pairwise.venn(area1 = 2784,    # Draw pairwise venn diagram
                             area2 = 5544,
                             cross.area = 1279,
                             fill = c("lightblue", "yellow"),
                             category.names = NULL)




## Visualise using Eulerr (proportional)
library(eulerr)
fit1 <- euler(c("A" = 1505, "B" = 4265, "C" = 2654,
                "A&B" = 1279, "A&C" = 792, "B&C" = 1771,
                "A&B&C" = 726))

# A=DMEMlow B=DMEMhigh C=plasmax
fit1 <- euler(c("A" = 160, "B" = 1941, "C" = 817,
                "A&B" = 553, "A&C" = 66, "B&C" = 1045,
                "A&B&C" = 726))

col <- c("#66D2D6", "#FBC740", "#E56997")

plot(fit1,
     quantities = TRUE,
     #legend = list(side = "right", font=1, labels = c("DMEMlow", "DMEMhigh", "Plasmax"), alpha=1, cex=1.2),
     labels = c("DMEMlow", "DMEMhigh","Plasmax"),
     fills = list(fill = col, alpha=0.6),
     )



## Proportional venn high and low DMEM only
fit1 <- euler(c("A" = 226, "B" = 2986,
                "A&B" = 1279))

col <- c("#66D2D6", "#FBC740")

plot(fit1,
     quantities = TRUE,
     #legend = list(side = "right", font=1, labels = c("DMEMlow", "DMEMhigh"), alpha=1, cex=1.2),
     labels = c("DMEMlow", "DMEMhigh"),
     fills = list(fill = col, alpha=0.6),
)

