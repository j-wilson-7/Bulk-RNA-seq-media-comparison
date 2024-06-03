## Enhanced Volcano Plots
#TGFbeta1 response for pHLFs grown in HighDMEM vs lowDMEM vs Plasmax media

#Download the packages required
#BiocManager::install('EnhancedVolcano')
#if (!require("pacman")) install.packages("pacman")

#Load packages into R session
library(EnhancedVolcano)
pacman::p_load_gh("trinker/textshape")

## Load DEseq2 results with annotated gene symbols
setwd("~/Documents/Chambers/Greg highvslow paper - Analysis2 Yale2022/Final-deseq2-analysis")
genes <- read.table("Low DMEM ALL.tabular", header=TRUE)
genes <- read.table("High DMEM ALL.tabular", header=TRUE)
genes <- read.table("Plasmax ALL.tabular", header=TRUE)

#Make the gene symbols into row columns
res <- column_to_rownames(genes, loc=2)

#Set publication theme
mytheme2 = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
          axis.text=element_text(face="bold"), legend.position = "none"
    ))



### lowDMEM
#Set colour of plot
color = c('darkgrey','darkgrey','darkgrey','#66D2D6') #dmemlow
#Plot volcano plot logFC2 and p-adjvalue0.05
p1 <- EnhancedVolcano(
  #Load data
  res,
  x = 'X3.log2.FC.',
  y = 'X7.P.adj',
  
  #Labelling
  lab = rownames(res),
  selectLab = '',
  cutoffLineType = 'dashed',

  title = '1,505 differentially expressed genes', #lowdmem
  subtitle = 'Fold change cutoff 1.5; adjusted p-value cutoff 0.05',
  #legendlabels = c('Not sig.', 'FC 1.5 & p-adj 0.05'),
  caption = '',
  #caption = bquote(~Log[2]~ "Fold change cutoff, 1.5; adjusted p-value cutoff, 0.05"),
  
  #Cutoffs
  pCutoff = 5*10e-03,
  FCcutoff = 0.58,
  
  #Size colour shape of points
  col = color,
  colAlpha=0.7)

#Custom axis tick marks
p1 +
  ggplot2::coord_cartesian(xlim=c(-10, 10)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-10, 10, 2)) +
  mytheme2 +
  theme(plot.subtitle = element_text(color='black', hjust=0.5, face='italic'), plot.title = element_text(color='black', hjust=0.5))







### highDMEM
color = c('darkgrey','darkgrey','darkgrey','#FBC740') #dmemhigh

#Plot volcano plot logFC2 and p-adjvalue0.05
p1 <- EnhancedVolcano(
  #Load data
  res,
  x = 'X3.log2.FC.',
  y = 'X7.P.adj',
  
  #Labelling
  lab = rownames(res),
  selectLab = '',
  cutoffLineType = 'dashed',
  
  title = '4,265 differentially expressed genes',   #highdmem
  subtitle = 'Fold change cutoff 1.5; adjusted p-value cutoff 0.05',
  caption = '',
  
  #Cutoffs
  pCutoff = 5*10e-03,
  FCcutoff = 0.58,
  
  #Size colour shape of points
  col = color,
  colAlpha=0.7)


#Custom axis tick marks
p1 +
  ggplot2::coord_cartesian(xlim=c(-10, 10)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-10, 10, 2)) +
  mytheme2 +
  theme(plot.subtitle = element_text(color='black', hjust=0.5, face='italic'), plot.title = element_text(color='black', hjust=0.5))








### Plasmax
color = c('darkgrey','darkgrey','darkgrey','#E56997') #plasmax

#Plot volcano plot logFC2 and p-adjvalue0.05
p1 <- EnhancedVolcano(
  #Load data
  res,
  x = 'X3.log2.FC.',
  y = 'X7.P.adj',
  
  #Labelling
  lab = rownames(res),
  selectLab = '',
  cutoffLineType = 'dashed',

  title = '2,654 differentially expressed genes',   #plasmax
  subtitle = 'Fold change cutoff 1.5; adjusted p-value cutoff 0.05',
  caption = '',
  
  #Cutoffs
  pCutoff = 5*10e-03,
  FCcutoff = 0.58,
  
  #Size colour shape of points
  col = color,
  colAlpha=0.7)
  

#Custom axis tick marks
p1 +
  ggplot2::coord_cartesian(xlim=c(-10, 10)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-10, 10, 2)) +
  mytheme2 +
  theme(plot.subtitle = element_text(color='black', hjust=0.5, face='italic'), plot.title = element_text(color='black', hjust=0.5))
  

  
  
  
  
