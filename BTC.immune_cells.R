#! /usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggdendro)
library(ggpubr)
library(cowplot)
library(FactoMineR)
library(sf)

'%ni%' <- function(x,y)!('%in%'(x,y))
font_size = 18
NO_AUTO = "non-PSC/IBD"
YES_AUTO = "PSC/IBD"


basedir <- "~/Documents/scripts/bitbucket/"
setwd(basedir)

source(paste0(basedir, "/BTC.functions.R"))
cell_metatypes <- fread('autoimmune-disorders/source-data/cell_metatypes.txt')


#### IMPORT ####

immune_matrix_raw <- fread('autoimmune-disorders/source-data/estimation_matrix.csv')

immune_cohort <- fread('autoimmune-disorders/source-data/samples.tsv')
immune_cohort$TMB_bin[immune_cohort$TMB <= 10] <- "TMB < 10"
immune_cohort$TMB_bin[immune_cohort$TMB > 10] <- "TMB > 10"
immune_cohort$final_status[immune_cohort$final_status == "No"] <- NO_AUTO
immune_cohort$final_status[immune_cohort$final_status == "Yes"] <- YES_AUTO

clinical_raw <- fread("autoimmune-disorders/source-data/clinical_data.txt")
clinical_raw$location[clinical_raw$location == 5] <- NA
clinical_raw$location[clinical_raw$location == 1] <- 2
clinical_raw$MSI[clinical_raw$MSI == "F"] <- "MMRp"
clinical_raw$MSI[clinical_raw$MSI == "T"] <- "MMRd"

immune_cohort <-  left_join(immune_cohort, clinical_raw, by=c("Donor"="ID"))

immune_matrix_split <- separate(data = immune_matrix_raw, col = cell_type, into = c("cell_type","software"), sep = "_")
immune_matrix <- left_join(immune_matrix_split, cell_metatypes, by=c("cell_type"="cell_type"))

immune_cohort$Do_nor <- gsub("BTC","BTC_",immune_cohort$Donor)

### descriptions ####

immune_matrix <- immune_matrix %>% 
  filter(software %in% c( "CIBERSORT-ABS"))

immune_melt_desc <- melt(immune_matrix)

ggplot(immune_melt_desc, aes(x = value, y = cell_type, fill = value)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")

immune_matrix$sd <- apply(immune_matrix[,-c(1,2, 36)],1, sd, na.rm = TRUE)
immune_matrix$var <- apply(immune_matrix[,-c(1,2, 36)],1, var, na.rm = TRUE)
immune_matrix$zero <- apply(immune_matrix[,-c(1,2, 36)] == 0, 1, sum)
immune_matrix$median <- apply(immune_matrix[,-c(1,2, 36)] , 1, median)

immune_matrix_t <- as.data.frame(t(immune_matrix[,-c(1, 2, 36:41)]))
names(immune_matrix_t) <- immune_matrix$cell_type
immune_matrix_t$samples <- row.names(immune_matrix_t)
immune_matrix_t <- left_join(immune_matrix_t, immune_cohort, by=c("samples"="Do_nor"))

t.test(immune_matrix_t$`Macrophage M2`[immune_matrix_t$final_status == YES_AUTO],
       immune_matrix_t$`Macrophage M2`[immune_matrix_t$final_status == NO_AUTO])

t.test(immune_matrix_t$`Macrophage M1`[immune_matrix_t$final_status == YES_AUTO],
       immune_matrix_t$`Macrophage M1`[immune_matrix_t$final_status == NO_AUTO])

#### CIBERSORT HEATMAP ####

TIL_META <- c("CD4+","CD8+","M", "NK") #, 

immune_matrix_filtered <- immune_matrix %>% 
  filter(software %in% c( "CIBERSORT-ABS")  & cell_metatype %in% TIL_META & !is.na(cell_metatype))
immune_melt <- melt(immune_matrix_filtered )

immune_melt_IO <- left_join(immune_cohort, immune_melt,  by=c("Do_nor"="variable"))
immune_melt_IO <- immune_melt_IO[!is.na(immune_melt_IO$value),]


immune_melt_IO$log_score <- log(immune_melt_IO$value)
immune_melt_IO$log_score[immune_melt_IO$log_score == -Inf] <- log(0.0001)


## cluster
CIBERSORT_dcast <- dcast(immune_melt_IO[,c("cell_type", "Donor", "log_score")], formula=Donor ~cell_type, fun.aggregate = mean)
immune_matrix_CIBERSORT_mt <- as.matrix(CIBERSORT_dcast[,-1])
rownames(immune_matrix_CIBERSORT_mt) <- CIBERSORT_dcast$Donor


CIBERSORT_dcast_cohort <- left_join(CIBERSORT_dcast, immune_cohort, by=c('Donor'='Donor'))

allOfit_pca <- PCA(immune_matrix_CIBERSORT_mt, graph = F)

allOfit_pca_labs <- cbind.data.frame("donor"=CIBERSORT_dcast_cohort$Donor,
                                     "autoimmune"=CIBERSORT_dcast_cohort$final_status,
                                     allOfit_pca$ind$coord
)



allOfit_pca_contribs<- cbind.data.frame(
  "Var"=row.names(allOfit_pca$var$coord),
  allOfit_pca$var$coord)

allOfit_pca_contribs$Varf <- 
  factor(allOfit_pca_contribs$Var,
         levels=allOfit_pca_contribs$Var[order(-allOfit_pca_contribs$Dim.1)])

allOfit_pca_contribs_sub <- allOfit_pca_contribs[allOfit_pca_contribs$Varf %in% c("Macrophage M1","T cell CD4+ naive", "NK cell activated"),]

png(paste0("autoimmune-disorders/results/cibersort.pca.png"), width = 1000, height = 1200, bg = "transparent")

  ggplot() +
    geom_point(data=allOfit_pca_labs,aes(x=Dim.1,y=Dim.2,color=autoimmune),size=4) + 
    theme_classic(base_size=26) +
    labs(color="PSC/IBD",
         x=paste("PC1 (",round(allOfit_pca$eig[1,2],digits=2),"%)",sep=""),
         y=paste("PC2 (",round(allOfit_pca$eig[2,2],digits=2),"%)",sep="")) +
    coord_sf() + 
    scale_color_manual( values=c("#67A9CF","#EF8A62" ), na.translate = F) +
    
    geom_segment(data=allOfit_pca_contribs_sub, aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2), arrow = arrow(length = unit(0.5, "cm"))) + 
    geom_text(data=allOfit_pca_contribs_sub, aes(x = Dim.1, y = Dim.2, label=Varf), hjust = -0.05, nudge_x = 0.05, size=5)
                              
        
   dev.off()
    
  
#### PLOT HEAT ####
  
   cibersort_median <- median(immune_melt_IO$log_score)
   
   immune_melt_IO$donor_clustered <- factor(x = immune_melt_IO$Donor,
                                            levels = allOfit_pca_labs$donor[order(-allOfit_pca_labs$Dim.1)], 
                                            ordered = TRUE)
   
 #  immune_melt_IO$cell_type_clustered <- factor(x = immune_melt_IO$cell_type,
 #                                               levels = names(CIBERSORT_dcast[,-1])[order.dendrogram(CIBERSORT_dendro_t)], 
 #                                               ordered = TRUE)
   

   this_theme <- theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       
                       panel.background = element_rect(fill='transparent'), #transparent panel bg
                       plot.background = element_rect(fill='transparent', color='transparent'), #transparent plot bg
                       legend.background = element_rect(color = NA,fill=NA), #transparent legend bg
                       legend.box.background = element_rect(color = NA,fill=NA),
                       plot.title = element_blank(),
                       axis.ticks.x =element_blank(),
                       axis.title.x = element_blank(),
                       
   )
   
 
   
  heatmap <- 
     ggplot(immune_melt_IO , aes(x=donor_clustered, y=cell_type, fill=log_score)) + 
     geom_tile(  color="white",) +
     facet_grid(cell_metatype~.,space="free",scales="free", switch="y")+
     labs(x="Sample",y="",fill="CIBERSORT\nscore (log)") + 
     theme_bw(base_size = font_size)+
     scale_fill_gradient2(midpoint = cibersort_median, low = "cornflowerblue", mid = "white", high = "red" ,na.value="white") +
     this_theme + 
     theme_classic(base_size = 22)+  theme(
       panel.grid.minor.y = element_blank(),
       legend.spacing = unit(0, 'cm'),
       panel.background = element_rect(fill='transparent'), #transparent panel bg
       plot.background = element_rect(fill='transparent', color='transparent'), #transparent plot bg
       legend.background = element_rect(color = NA,fill=NA), #transparent legend bg
       legend.box.background = element_rect(color = NA,fill=NA), #transparent legend panel
       plot.margin = unit(c(0,1,0,1), "lines"),
       
       #  legend.title = element_blank(),
       axis.title.y=element_blank(),
       axis.title.x=element_blank(),
       axis.text.x=element_blank(),
       axis.ticks.x=element_blank(),
       panel.grid.major = element_blank(), #remove major gridlines
       panel.grid.minor = element_blank(), #remove minor gridlines
       strip.text.y.left = element_text(angle = 0)
     )  
   
   
   TMB_PLOT <- 
     ggplot(immune_melt_IO, aes(x=donor_clustered, y="TMB", fill=TMB_bin)) + 
     geom_tile() +
     theme_bw(base_size = font_size)+
     guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
     # guides(fill="none")+
     labs(y="",fill='') + 
     scale_fill_manual( values = c("#969696", "#B2182B")) +
     this_theme + 
     theme(
       legend.title = element_blank(),
       strip.text = element_blank(),
       panel.grid.major.x = element_blank(), 
       plot.margin = unit(c(0,0,0,0), "lines"),
       axis.text.x =element_blank()
     )
   
   
   mmrd_plot <- 
     ggplot(immune_melt_IO, aes(x=donor_clustered,y="MMRd",fill=MSI)) + 
     geom_tile() +
     labs(fill='') + 
     scale_fill_manual( values=c("#FA9FB5","#4D004B" ), na.translate = F) +
     
     theme_bw(base_size = font_size)+
     guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
     this_theme+  theme(  legend.title = element_blank(), axis.text.x =element_blank(),axis.title.y=element_blank() )   
   
   location_plot <- 
     ggplot(immune_melt_IO, aes(x=donor_clustered,y="Location",fill=as.factor(location))) + 
     geom_tile() +
     theme_bw(base_size = font_size)+
     labs(fill='') + 
     
     scale_fill_manual(labels=c("eCCA", "pCCA", "iCCA", "Gall"), values=c("#CC4C02","#525252","#238443","#88419D"), na.translate = F) +
     this_theme+  theme(  legend.title = element_blank(),axis.text.x =element_blank(), axis.title.y=element_blank())    +guides(fill=guide_legend(nrow=1,byrow=TRUE))
   
   
   pscibd_plot <-  
     ggplot(immune_melt_IO, aes(x=donor_clustered,y="PSC/IBD",fill=as.factor(final_status))) + 
     geom_tile() +
     labs(fill='') + 
     
     theme_bw(base_size = font_size)+
     scale_fill_manual( values=c("#67A9CF","#EF8A62" ), na.translate = F) +
     
     this_theme+  theme(  legend.title = element_blank(),axis.text.x =element_blank(), axis.title.y=element_blank())    +guides(fill=guide_legend(nrow=1,byrow=TRUE))
   
   
   specimen_plot <-  
     ggplot(immune_melt_IO, aes(x=donor_clustered,y="specimen",fill=as.factor(specimen))) + 
     geom_tile() +
     labs(fill='') + 
     
     theme_bw(base_size = font_size)+
     scale_fill_manual( values=c("#67A9CF","#EF8A62" ), na.translate = F) +
     
     this_theme+  theme(  legend.title = element_blank(),axis.text.x =element_blank(), axis.title.y=element_blank())    +guides(fill=guide_legend(nrow=1,byrow=TRUE))
   
   png(paste0("autoimmune-disorders/results/cibersort.heatmap.png"), width = 1200, height = 550) #, bg = "transparent"
   gap = -0.005
   skinny = 0.015
   
   plot_grid(
     pscibd_plot, NULL, 
     location_plot, NULL, 
     
     mmrd_plot, NULL,  
     TMB_PLOT, NULL, 
     heatmap, 
     
     
     rel_heights = c( 
       skinny, gap, 
       skinny, gap, 
       
       skinny, gap,
       skinny, gap,
       0.1) ,
     ncol = 1,  align = 'hv', axis = "rltb")
   
   dev.off()
   
#### bimodality ####
   
library(mclust)

pc1_psc <- allOfit_pca_labs$Dim.1[allOfit_pca_labs$autoimmune == YES_AUTO]
   
 v <- rep(1, 8)
 for (i in 1:8) {
   fit_mclust <- Mclust(pc1_psc,G=i)
   v[i] <- max(fit_mclust$BIC,na.rm = T)
 }
 
 mclust_BIC <- cbind.data.frame(v,"k"=seq(1,length(v),1))
 ggplot(mclust_BIC,aes(x=k,y=v))+
   geom_point() + 
   geom_line() +
   geom_vline(xintercept = mclust_BIC$k[mclust_BIC$v == min(mclust_BIC$v)],
              linetype="longdash",color="gray",alpha=0.75)+
   theme_bw(base_size=15) +labs(y="Bayesian Information Criterion (BIC)",x="Number of groups G")

 
 
 fit_mclust <- Mclust(pc1_psc, G=2)
 pc1_psc_df <- as.data.frame(pc1_psc)
 pc1_psc_df$mclusters <- fit_mclust$classification  
 
 ggplot(data=pc1_psc_df,aes(x=pc1_psc, fill=as.factor(mclusters ))) +
   geom_histogram()+ 
   theme_bw() 
   
