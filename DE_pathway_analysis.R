##########################################
## libraries
##########################################

library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(limma)
library(colorRamp2)
library(fgsea)
library(enrichR)
library(msigdbr)
library(ggrepel)
library(ComplexHeatmap)

######################################
## functions
######################################

volcanoPlot <- function(feature, coef, pval, padj, pos.cutoff, neg.cutoff, x.lab, padj.label, cutoff){
  
  data <- data.frame(feature = feature,
                     coef = coef,
                     pval = pval,
                     FDR = padj)
  
  if( padj.label == FALSE){
    
    data$diffexpressed <- "NO"
    data$diffexpressed[data$coef > 0 & data$pval < cutoff] <- paste(paste("Pval < ", cutoff, sep=""), "Coef > 0", sep=", ")
    data$diffexpressed[data$coef < 0 & data$pval < cutoff] <- paste(paste("Pval < ", cutoff, sep=""), "Coef < 0", sep=", ")
    
    mycolors <- c( "#EF8A62","#67A9CF", "#999999")
    names(mycolors) <- c(paste(paste("Pval < ", cutoff, sep=""), "Coef > 0", sep=", "),
                         paste(paste("Pval < ", cutoff, sep=""), "Coef < 0", sep=", "),
                         "NO")
    
    data$delabel <- NA
    data <- data[order(data$pval, decreasing = FALSE), ]
    id_pos <- data[data$coef > 0 , "feature"][1:pos.cutoff]
    id_neg <- data[data$coef < 0 , "feature"][1:neg.cutoff]
    id <- c(id_pos, id_neg)
    
    for(j in 1:length(id)){
      k <- which(data$feature == id[j])
      data$delabel[k] <- data[k, ]$feature
    }
    
  }else{
    
    data$diffexpressed <- "NO"
    data$diffexpressed[data$coef > 0 & data$FDR < cutoff] <- paste(paste("FDR < ", cutoff, sep=""), "Coef > 0", sep=", ")
    data$diffexpressed[data$coef < 0 & data$FDR < cutoff] <- paste(paste("FDR < ", cutoff, sep=""), "Coef < 0", sep=", ")
    
    mycolors <- c( "#EE6677","#4477AA", "#999999")
    names(mycolors) <- c(paste(paste("FDR < ", cutoff, sep=""), "Coef > 0", sep=", "),
                         paste(paste("FDR < ", cutoff, sep=""), "Coef < 0", sep=", "),
                         "NO")
    
    data$delabel <- NA
    data <- data[order(data$FDR, decreasing = FALSE), ]
    id_pos <- data[data$coef > 0 , "feature"][1:pos.cutoff]
    id_neg <- data[data$coef < 0 , "feature"][1:neg.cutoff]
    id <- c(id_pos, id_neg)
    
    for(j in 1:length(id)){
      k <- which(data$feature == id[j])
      data$delabel[k] <- data[k, ]$feature
    }
    
  }
  
  ggplot(data=data, aes(x=coef, y=-log10(pval), col= diffexpressed)) +
    geom_point(size = 3) + theme_minimal() +
    ylab("-log10 P value") +
    xlab(x.lab) +
    scale_colour_manual(values = mycolors) +
    theme(
      axis.text.x=element_text(size=10, face = "bold"),
      axis.title=element_text(size=10, face = "bold"),
      axis.text.y=element_text(size=10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position="bottom",
      legend.text = element_text(size = 9, face="bold"),
      legend.title = element_blank()) +
    geom_text_repel(aes(label= delabel),
                    size = 2.2,
                    color = "black",
                    min.segment.length = 0,
                    na.rm = TRUE, direction = "both",
                    seed = 2356,
                    fontface= "bold",
                    max.overlaps = 50)
}

####################################
## load CCA data
####################################

se <- readRDS(file.path("C:/lbr-pcs-immunotherapy/source-data/", "SE_CCA.rds"))

# Extract Clinical, rna, and annotation data 
clin <- data.frame(colData(se)) 
expr_tpm <- assays(se)[["expr_tpm"]] 
annot <- data.frame(rowData(se))  

clin$stage <- with(clin, ifelse(stage %in% c("IIIA", "IIIB", "IIIC"), "III",
                                ifelse(stage %in% c("IVA", "IVB"), "IV", stage)))

clin$sex <- ifelse(clin$sex == "F", "Female", "Male")

clin$psc_ibd[clin$psc_ibd == 'y'] <- 'PSC/IBD'
clin$psc_ibd[clin$psc_ibd == 'n'] <- 'non-PSC/IBD'

clin$TMB <- ifelse(clin$tmb > 10, "TMB > 10", "TMB < 10")
clin$msi <- ifelse(clin$msi == "y", "MMRd", "MMRp")

#####################################
## data pre-processing
#####################################

annot_proteincoding <- annot[annot$gene_type == "protein_coding",] 
expr_tpm<- expr_tpm[rownames(expr_tpm) %in% rownames(annot_proteincoding),]

r <- as.numeric(apply(expr_tpm, 1, function(i) sum(round(i, 6) == round(log2(1), 6))))
remove <- which(r > dim(expr_tpm)[2] * 0.5) 
expr_tpm <- expr_tpm[-remove, ] 

#####################################
## differential expression analysis
#####################################

design <- model.matrix(~ clin$psc_ibd)
fit <- lmFit(expr_tpm, design)
fit <- eBayes(fit)

top.table <- topTable(fit, sort.by = "P", n = Inf)

#write.csv(top.table, file = file.path("C:/lbr-pcs-immunotherapy/results/", "DE_analysis.csv"), row.names = T)

## Create volcano plot

jpeg(file = file.path("C:/lbr-pcs-immunotherapy/results/", "DE_FDR_0.05.jpeg"), 
     width = 800, height = 700, res = 150)

p <- volcanoPlot(feature = rownames(top.table),
                 coef = top.table$logFC,
                 pval = top.table$P.Value,
                 padj = top.table$adj.P.Val,
                 pos.cutoff = 25,
                 neg.cutoff = 20,
                 x.lab = expression(log[2] ~ 'fold change'),
                 padj.label = TRUE,
                 cutoff = 0.05)

p

dev.off()

## Create heatmpa uisng top 100 significant genes
df <- top.table[1:100, ]
expr <- expr_tpm[rownames(expr_tpm) %in% rownames(df), ]
expr <- t(scale(t(expr)))

col = list( "PSC/IBD" = c( "non-PSC/IBD" = "#67A9CF" , "PSC/IBD" = "#EF8A62" ) ,
            Location = c("dCCA" = "#cc4c02", "GBC" = "#88419d", "iCCA" = "#238443", 
                         "pCCA" = "#525252"),
            Sex = c( "Female" = "#fee090" , "Male" = "#4477AA" ),
            TMB = c("TMB < 10" = "#969696", "TMB > 10" = "#b2182b"), 
            Stage = c("II" = "#35978f" , "III" = "#bf812d" , "IV" = "#8073ac"),
            MMR = c("MMRp" = "#4d004b", "MMRd" = "#fa9fb5")) 

# Create the heatmap annotation
ha <- HeatmapAnnotation(
  "PSC/IBD" = clin$psc_ibd,
  Location = clin$location,
  Sex = clin$sex,
  TMB = clin$TMB,
  Stage = clin$stage,
  MMR = clin$msi,
  col = col,
  annotation_name_gp = gpar(fontsize = 8),
  simple_anno_size = unit(0.5, "cm")
)

# Combine the heatmap and the annotation
jpeg(file=file.path("C:/lbr-pcs-immunotherapy/results/", "heatmap_top_DE_genes.jpeg"),
     width = 900, height = 1200, res = 150)

ht_data <- Heatmap(expr, name = "expr",
                   top_annotation = ha,
                   show_row_names = FALSE,
                   #row_names_gp = gpar(fontsize = 6),
                   column_names_gp = gpar(fontsize = 8),
                   show_column_names = TRUE,
                   colorRamp2(c(-3, 0, 3), c("#045a8d", "white", "#800026")))

draw(ht_data,
     column_title_gp = gpar(fontsize = 8, fontface = "bold"),
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side="bottom")

dev.off()

#######################################
## MSigDB Hallmark: GSEA Analysis
#######################################
# rank statistics
ranks <- top.table$logFC
#ranks <- sign(top.table$logFC) * (-log10(top.table$P.Value))
names(ranks) <- rownames(top.table)
ranks <- sort(ranks, decreasing = T)

# load downloaded Hallmark pathway
pathwaysH <- gmtPathways(file.path("C:/lbr-pcs-immunotherapy/source-data/", "h.all.v2023.2.Hs.symbols.gmt"))

fgseaRes <- fgsea(pathwaysH, 
                  ranks, 
                  minSize=15, 
                  maxSize = 500)
fgseaRes <- fgseaRes[!is.na(fgseaRes$padj), ]

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) 

fgseaResTidy$leadingEdge <- sapply(fgseaResTidy$leadingEdge, paste, collapse = ";")
fgseaResTidy <- fgseaResTidy[order(fgseaResTidy$padj), ]
fgseaResTidy$pathway <- substr(fgseaResTidy$pathway, 10, nchar(fgseaResTidy$pathway))


#write.csv(fgseaResTidy, file = file.path("C:/lbr-pcs-immunotherapy/results/", "pathway_association_Hallmark.csv"), row.names = FALSE)

## create dotplot

jpeg(file=file.path("C:/lbr-pcs-immunotherapy/results/", "dotplot_hallmark_FDR_0.05.jpeg"),
     width = 1000, height = 700, res = 150)

p <- ggplot(fgseaResTidy[fgseaResTidy$padj < 0.05, ], 
            aes(y=reorder(pathway, size),x= NES)) + 
  geom_point(aes(color=padj, size=size)) + 
  scale_color_gradientn(colours = c("#084594","#b10026")) + 
  labs(x='normalized enrichment score', y=NULL ) + 
  theme(
    axis.text.x=element_text(size=10),
    axis.title=element_text(size=10, face = "bold"),
    axis.text.y=element_text(size=8, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.text = element_text(size = 8, face="bold"))

p

dev.off()

#######################################
## MSigDB KEGG: GSEA Analysis
#######################################

pathwaysKEGG <- gmtPathways(file.path("C:/lbr-pcs-immunotherapy/source-data/", "c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt"))

fgseaRes_kegg <- fgsea(pathwaysKEGG, 
                       ranks, 
                       minSize=15, 
                       maxSize = 500)

fgseaResTidy_kegg <- fgseaRes_kegg %>%
  as_tibble() %>%
  arrange(desc(NES)) 

fgseaResTidy_kegg$leadingEdge <- sapply(fgseaResTidy_kegg$leadingEdge, paste, collapse = ";")
fgseaResTidy_kegg <- fgseaResTidy_kegg[order(fgseaResTidy_kegg$padj), ]
fgseaResTidy_kegg $pathway <- substr(fgseaResTidy_kegg$pathway, 14, nchar(fgseaResTidy_kegg$pathway))


#write.csv(fgseaResTidy_kegg, file = file.path("C:/lbr-pcs-immunotherapy/results/", "pathway_association_KEGG_MEDICUS.csv"), row.names = FALSE)

#######################################
## MSigDB GO: GSEA Analysis
#######################################

pathwaysGO <- gmtPathways(file.path("C:/lbr-pcs-immunotherapy/source-data/", "c5.go.v2023.2.Hs.symbols.gmt"))
fgseaRes_go <- fgsea(pathwaysGO, 
                     ranks, 
                     minSize=15, 
                     maxSize = 500)

fgseaRes_go <- fgseaRes_go[!is.na(fgseaRes_go$padj), ]
fgseaResTidy_go <- fgseaRes_go %>%
  as_tibble() %>%
  arrange(desc(NES)) 

fgseaResTidy_go$leadingEdge <- sapply(fgseaResTidy_go$leadingEdge, paste, collapse = ";")
fgseaResTidy_go <- fgseaResTidy_go[order(fgseaResTidy_go$padj), ]
fgseaResTidy_go $pathway <- substr(fgseaResTidy_go$pathway, 6, nchar(fgseaResTidy_go$pathway))

#write.csv(fgseaResTidy_go, file = file.path("C:/lbr-pcs-immunotherapy/results/", "pathway_association_GO.csv"), row.names = FALSE)







