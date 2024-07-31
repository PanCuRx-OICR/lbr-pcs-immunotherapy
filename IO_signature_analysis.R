## libraries

library(PredictioR)
library(GSVA)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)
library(MultiAssayExperiment)

############################################
## load data
############################################
dat <- readRDS(file.path("C:/lbr-pcs-immunotherapy/source-data", "SE_CCA.rds"))
expr <- as.data.frame(assay(dat))
clin <- as.data.frame(colData(dat))
annot <- rowData(dat)

clin$stage <- with(clin, ifelse(stage %in% c("IIIA", "IIIB", "IIIC"), "III",
                                ifelse(stage %in% c("IVA", "IVB"), "IV", stage)))

clin$sex <- ifelse(clin$sex == "F", "Female", "Male")

clin$psc_ibd[clin$psc_ibd == 'y'] <- 'PSC/IBD'
clin$psc_ibd[clin$psc_ibd == 'n'] <- 'non-PSC/IBD'

clin$TMB <- ifelse(clin$tmb > 10, "TMB > 10", "TMB < 10")
clin$msi <- ifelse(clin$msi == "y", "MMRd", "MMRp")

###########################################################################
## load signature data from repo: https://github.com/bhklab/SignatureSets
###########################################################################
## load signature's information
signature <- read.csv( file.path("C:/lbr-pcs-immunotherapy/source-data", "signature_information.csv"))
signature$method <- as.character( signature$method )

## load curated signature's data
GeneSig_list <- list.files(file.path("C:/lbr-pcs-immunotherapy/source-data", "signature"))
GeneSig_list <- GeneSig_list[order(GeneSig_list)]
signature <- signature[order(signature$signature), ]

######################################################
## compute signature score
######################################################

AllGeneSig <- lapply(1:length(GeneSig_list), function(i){ 
  
  print(GeneSig_list[i])
  
  load(file.path("C:/lbr-pcs-immunotherapy/source-data", "signature", GeneSig_list[i]))
  
  if(signature$method[i] == "GSVA"){
    
    geneSig <- geneSigGSVA(dat.icb = expr,
                           sig = sig,
                           sig.name = signature$signature[i],
                           missing.perc = 0.5,
                           const.int = 0.001,
                           n.cutoff = 15,
                           sig.perc = 0.8,
                           study = "CCA")
    
    if(sum(!is.na(geneSig)) > 0){
      geneSig <- geneSig[1,]
    }     
    
  }
  
  if(signature$method[i] == "Weighted Mean"){
    
    geneSig <- geneSigMean(dat.icb = expr,
                           sig = sig,
                           sig.name = signature$signature[i],
                           missing.perc = 0.5,
                           const.int = 0.001,
                           n.cutoff = 15,
                           sig.perc = 0.8,
                           study = "CCA")
    
  }
  
  
  if(signature$method[i] == "ssGSEA"){
    
    geneSig <- geneSigssGSEA(dat.icb = expr,
                             sig = sig,
                             sig.name = signature$signature[i],
                             missing.perc = 0.5,
                             const.int = 0.001,
                             n.cutoff = 15,
                             sig.perc = 0.8,
                             study = "CCA")
    
    if(sum(!is.na(geneSig)) > 0){
      geneSig <- geneSig[1,]
    }     
    
  }
  
  if(signature$method[i] == "Specific Algorithm" & signature$signature[i] == "COX-IS_Bonavita"){
    
    geneSig <- geneSigCOX_IS(dat.icb = expr,
                             sig = sig,
                             sig.name = signature$signature[i],
                             missing.perc = 0.5,
                             const.int = 0.001,
                             n.cutoff = 15,
                             sig.perc = 0.8,
                             study = "CCA")
    
  }
  
  if(signature$method[i] == "Specific Algorithm" & signature$signature[i] == "IPRES_Hugo"){
    
    geneSig <- geneSigIPRES(dat.icb = expr,
                            sig = sig,
                            sig.name = signature$signature[i],
                            missing.perc = 0.5,
                            const.int = 0.001,
                            n.cutoff = 15,
                            sig.perc = 0.8,
                            study = "CCA")
    
  }
  
  
  if(signature$method[i] == "Specific Algorithm" & signature$signature[i] == "PredictIO_Bareche"){
    
    geneSig <- geneSigPredictIO(dat.icb = expr,
                                sig = sig,
                                sig.name = signature$signature[i],
                                missing.perc = 0.5,
                                const.int = 0.001,
                                n.cutoff = 15,
                                sig.perc = 0.8,
                                study = "CCA")
    
  }
  
  if(signature$method[i] == "Specific Algorithm" & signature$signature[i] == "IPS_Charoentong"){
    
    geneSig <- geneSigIPS(dat.icb = expr,
                          sig = sig,
                          sig.name = signature$signature[i],
                          missing.perc = 0.5,
                          const.int = 0.001,
                          n.cutoff = 15,
                          sig.perc = 0.8,
                          study = "CCA")
    
  }
  
  if(signature$method[i] == "Specific Algorithm" & signature$signature[i] == "IPSOV_Shen"){
    
    geneSig <- geneSigIPSOV(dat.icb = expr,
                            sig = sig,
                            sig.name = signature$signature[i],
                            missing.perc = 0.5,
                            const.int = 0.001,
                            n.cutoff = 15,
                            sig.perc = 0.8,
                            study = "CCA")
    
  }
  
  if(signature$method[i] == "Specific Algorithm" & signature$signature[i] == "PassON_Du"){
    
    geneSig <- geneSigPassON(dat.icb = expr,
                             sig = sig,
                             sig.name = signature$signature[i],
                             missing.perc = 0.5,
                             const.int = 0.001,
                             n.cutoff = 15,
                             sig.perc = 0.8,
                             study = "CCA")
    
  }
  
  geneSig 
  
})

AllGeneSig <- do.call(rbind, AllGeneSig)
rownames(AllGeneSig) <- signature$signature
remove <- which(is.na(rowSums(AllGeneSig)))
if(length(remove) > 0){
  AllGeneSig <- AllGeneSig[-remove, ]
}

#save(AllGeneSig, file=file.path("C:/lbr-pcs-immunotherapy/results", "sig.RData"))

##################################################
## signatures and PSC/IBD: Wilcox test
##################################################
# Wilcox test
wilcox.res <- lapply(1:nrow(AllGeneSig), function(i){
  
  fit <- wilcox.test(AllGeneSig[i, ] ~ clin$psc_ibd, exact = FALSE)
  expr_PSC_IBD <- AllGeneSig[i, ][names(AllGeneSig[i, ]) %in% clin[clin$psc_ibd == "PSC/IBD", "id"]]
  expr_non_PSC_IBD <- AllGeneSig[i, ][names(AllGeneSig[i, ]) %in% clin[clin$psc_ibd == "non-PSC/IBD", "id"]]
  
  data.frame(signature = rownames(AllGeneSig)[i],
             pval = round(fit$p.value, 3),
             median_PSC_IBD = median(expr_PSC_IBD, na.rm = TRUE),
             median_non_PSC_IBD = median(expr_non_PSC_IBD, na.rm = TRUE))
  
})

wilcox.res.IOsignature <- do.call(rbind, wilcox.res)
#write.csv(wilcox.res.IOsignature, file = file.path("C:/lbr-pcs-immunotherapy/results", "IOsignature_wilcox.csv"), 
#          row.names = T)

## create violin plot 

for(i in 1:nrow(AllGeneSig)){
  
  df <- data.frame( y = AllGeneSig[i, ],
                    group = clin$psc_ibd)
  
  dir <- file.path("C:/lbr-pcs-immunotherapy/results", "signature_PSC_IBD")
  jpeg(paste(dir, paste(rownames(AllGeneSig)[i], ".jpg", sep=""), sep="/"), 
       width = 4, height = 4, units = 'in', res = 300)
  
  p <- ggplot(df, aes(x=group, y=y, fill = group)) + 
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.07, fill = "white")+
    xlab("")+
    ylab("")+
    ggtitle(strsplit(rownames(AllGeneSig)[i], "_")[[1]][1]) +
    scale_fill_manual(values=c("#67A9CF", "#EF8A62"))+
    theme(axis.text.x=element_text(size=10,  face="bold"),
          axis.title=element_text(size=10,face="bold"),
          axis.text.y=element_text(size=10),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill=FALSE)
  
  print(p)
  
  dev.off()
  
}

#################################################################
## signatures and TMB: Wilcox test
#################################################################
# Wilcox test

wilcox.res <- lapply(1:nrow(AllGeneSig), function(i){
  
  fit <- wilcox.test(AllGeneSig[i, ] ~ clin$TMB, exact = FALSE)
  
  expr_tmb_less <- AllGeneSig[i, ][names(AllGeneSig[i, ]) %in% clin[clin$TMB == "TMB < 10", "id"]]
  expr_tmb_great <- AllGeneSig[i, ][names(AllGeneSig[i, ]) %in% clin[clin$TMB == "TMB > 10", "id"]]
  
  data.frame(signature = rownames(AllGeneSig)[i],
             pval = round(fit$p.value, 3),
             median_TMB_less = median(expr_tmb_less, na.rm = TRUE),
             median_TMB_great = median(expr_tmb_great, na.rm = TRUE))
  
})

wilcox.res <- do.call(rbind, wilcox.res)
#write.csv(wilcox.res, file = file.path("C:/lbr-pcs-immunotherapy/results", "IOsignature_TMB.csv"), 
#          row.names = T)

## create violin plot 
for(i in 1:nrow(AllGeneSig)){
  
  df <- data.frame( y = AllGeneSig[i, ],
                    group = clin$TMB)
  
  dir <- file.path("C:/lbr-pcs-immunotherapy/results", "signature_TMB")
  jpeg(paste(dir, paste(rownames(AllGeneSig)[i], ".jpg", sep=""), sep="/"), 
       width = 4, height = 4, units = 'in', res = 300)
  
  p <- ggplot(df, aes(x=group, y=y, fill = group)) + 
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.07, fill = "white")+
    xlab("")+
    ylab("")+
    ggtitle(strsplit(rownames(AllGeneSig)[i], "_")[[1]][1]) +
    scale_fill_manual(values=c("#969696", "#b2182b"))+ 
    theme(axis.text.x=element_text(size=10,  face="bold"),
          axis.title=element_text(size=10,face="bold"),
          axis.text.y=element_text(size=10),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill=FALSE)
  
  print(p)
  
  dev.off()
  
}

#################################################################
## Selected IO biomarkers (gene level) and PSC/IBD: Wilcox test
#################################################################

genes <- c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "HAVCR2", "LAG3", "TIGIT")
sub_expr <- expr[rownames(expr) %in% genes, ]

# Wilcox test
wilcox.res <- lapply(1:nrow(sub_expr), function(i){
  
  fit <- wilcox.test(as.numeric(sub_expr[i, ]) ~ clin$psc_ibd, exact = FALSE)
  expr_PSC_IBD <- as.numeric(sub_expr[i, ][names(sub_expr[i, ]) %in% clin[clin$psc_ibd == "PSC/IBD", "id"]])
  expr_non_PSC_IBD <- as.numeric(sub_expr[i, ][names(sub_expr[i, ]) %in% clin[clin$psc_ibd == "non-PSC/IBD", "id"]])
  
  data.frame(signature = rownames(sub_expr)[i],
             pval = round(fit$p.value, 3),
             median_PSC_IBD = median(expr_PSC_IBD, na.rm = TRUE),
             median_non_PSC_IBD = median(expr_non_PSC_IBD, na.rm = TRUE))
  
})

wilcox.res <- do.call(rbind, wilcox.res)
# write.csv(wilcox.res, file = file.path("C:/lbr-pcs-immunotherapy/results", "genes_PSC_IBD.csv"), row.names = T)

## create violin plot 
for(i in 1:nrow(sub_expr)){
  
  df <- data.frame( y = as.numeric(sub_expr[i, ]),
                    group = clin$psc_ibd)
  
  dir <- file.path("C:/lbr-pcs-immunotherapy/results", "genes_PSC_IBD")
  jpeg(paste(dir, paste(rownames(sub_expr)[i], ".jpg", sep=""), sep="/"), 
       width = 4, height = 4, units = 'in', res = 300)
  
  p <- ggplot(df, aes(x=group, y=y, fill = group)) + 
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.07, fill = "white")+
    xlab("")+
    ylab("")+
    ggtitle(rownames(sub_expr)[i]) +
    scale_fill_manual(values=c("#67A9CF", "#EF8A62"))+
    theme(axis.text.x=element_text(size=10,  face="bold"),
          axis.title=element_text(size=10,face="bold"),
          axis.text.y=element_text(size=10),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill=FALSE)
  
  print(p)
  
  dev.off()
  
}


#############################################
## Heatmap for IO bimarkers 
#############################################

## Part 1: IO genes and selected signatures
genes <- c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "HAVCR2", "LAG3", "TIGIT" )
sig_name <- c("Chemokine_Messina", "IFN_Ayers", "IRG_Ayers",
              "Inflammatory_Thompson", "TIS_Damotte", "CYT_Rooney")

sub_expr <- expr[rownames(expr) %in% genes, ]
sub_AllGeneSig <- AllGeneSig[rownames(AllGeneSig) %in% sig_name, ]

df <- rbind(sub_expr, sub_AllGeneSig)
df <- t(scale(t(df)))

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
jpeg(file=file.path("C:/lbr-pcs-immunotherapy/results", "heatmap_selected_signature_genes.jpeg"),
     width = 900, height = 800, res = 150)

ht_data <- Heatmap(df, name = "expr",
                   top_annotation = ha,
                   show_row_names = TRUE, 
                   show_column_names = TRUE,
                   row_names_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 8),
                   colorRamp2(c(min(df), 0, max(df)), 
                              c("#045a8d", "white", "#800026")))

draw(ht_data,
     column_title_gp = gpar(fontsize = 8, fontface = "bold"),
     #merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side="bottom")

dev.off()

## Part 2: sensitive-IO signatures (i.e., gene sets) and genes. Include sensitive signatures from original published paper

genes <- c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "HAVCR2", "LAG3", "TIGIT" )
sig_name <- signature[signature$association == "sensitive", "signature"]

sub_expr <- expr[rownames(expr) %in% genes, ]
sub_AllGeneSig <- AllGeneSig[rownames(AllGeneSig) %in% sig_name, ]

df <- rbind(sub_expr, sub_AllGeneSig)
df <- t(scale(t(df)))

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
jpeg(file=file.path("C:/lbr-pcs-immunotherapy/results", "heatmap_sensitive_IOSignature_genes.jpeg"),
     width = 1000, height = 1200, res = 150)

ht_data <- Heatmap(df, name = "expr",
                   top_annotation = ha,
                   show_row_names = TRUE,
                   row_names_gp = gpar(fontsize = 6),
                   column_names_gp = gpar(fontsize = 8),
                   show_column_names = TRUE,
                   colorRamp2(c(-3, 0, 3), c("#045a8d", "white", "#800026")))

draw(ht_data,
     column_title_gp = gpar(fontsize = 8, fontface = "bold"),
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side="bottom")


dev.off()

## Part 3: IO signatures (i.e., gene sets) and genes (with Wilcox P-value < 0.1)

genes <- c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "HAVCR2", "LAG3", "TIGIT" )
sig_name <- wilcox.res.IOsignature[wilcox.res.IOsignature$pval < 0.1, "signature"]

sub_expr <- expr[rownames(expr) %in% genes, ]
sub_AllGeneSig <- AllGeneSig[rownames(AllGeneSig) %in% sig_name, ]

df <- rbind(sub_expr, sub_AllGeneSig)
df <- t(scale(t(df)))

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
jpeg(file=file.path("C:/lbr-pcs-immunotherapy/results", "heatmap_IOSignature_genes_Wilcox_0.1.jpeg"),
     width = 1000, height = 800, res = 150)

ht_data <- Heatmap(df, name = "expr",
                   top_annotation = ha,
                   show_row_names = TRUE,
                   row_names_gp = gpar(fontsize = 6),
                   column_names_gp = gpar(fontsize = 8),
                   show_column_names = TRUE,
                   colorRamp2(c(-3, 0, 3), c("#045a8d", "white", "#800026")))

draw(ht_data,
     column_title_gp = gpar(fontsize = 8, fontface = "bold"),
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side="bottom")


dev.off()


## Part 4: IO signatures (i.e., gene sets) and genes (with Wilcox P-value < 0.05)

genes <- c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "HAVCR2", "LAG3", "TIGIT" )
sig_name <- wilcox.res.IOsignature[wilcox.res.IOsignature$pval < 0.05, "signature"]

sub_expr <- expr[rownames(expr) %in% genes, ]
sub_AllGeneSig <- AllGeneSig[rownames(AllGeneSig) %in% sig_name, ]

df <- rbind(sub_expr, sub_AllGeneSig)
df <- t(scale(t(df)))

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
jpeg(file=file.path("C:/lbr-pcs-immunotherapy/results", "heatmap_IOSignature_genes_Wilcox_0.05.jpeg"),
     width = 1000, height = 600, res = 150)

ht_data <- Heatmap(df, name = "expr",
                   top_annotation = ha,
                   show_row_names = TRUE,
                   row_names_gp = gpar(fontsize = 6),
                   column_names_gp = gpar(fontsize = 8),
                   show_column_names = TRUE,
                   colorRamp2(c(-3, 0, 3), c("#045a8d", "white", "#800026")))

draw(ht_data,
     column_title_gp = gpar(fontsize = 8, fontface = "bold"),
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side="bottom")


dev.off()


## Part 4: All IO signatures (i.e., gene sets) and genes (no cut-off based on Wilcox test and/or original association from publication)

genes <- c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "HAVCR2", "LAG3", "TIGIT" )

sub_expr <- expr[rownames(expr) %in% genes, ]
sub_AllGeneSig <- AllGeneSig[rownames(AllGeneSig) %in% signature$signature, ]

df <- rbind(sub_expr, sub_AllGeneSig)
df <- t(scale(t(df)))

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
jpeg(file=file.path("C:/lbr-pcs-immunotherapy/results", "heatmap_IOSignature_genes.jpeg"),
     width = 1000, height = 1200, res = 150)

ht_data <- Heatmap(df, name = "expr",
                   top_annotation = ha,
                   show_row_names = TRUE,
                   row_names_gp = gpar(fontsize = 6),
                   column_names_gp = gpar(fontsize = 8),
                   show_column_names = TRUE,
                   colorRamp2(c(-3, 0, 3), c("#045a8d", "white", "#800026")))

draw(ht_data,
     column_title_gp = gpar(fontsize = 8, fontface = "bold"),
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side="bottom")


dev.off()



