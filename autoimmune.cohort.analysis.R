#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(ggplot2)
library(jsonlite)
library(tidyr)
library(cowplot)

'%ni%' <- function(x,y)!('%in%'(x,y))
fontSize = 28

#basedir <- paste(Sys.getenv(c("BASE_DIR")), sep='/')
basedir <- "~/Documents/scripts/bitbucket/"
setwd(basedir)

source(paste0(basedir, "/BTC.functions.R"))

stage_plotting <- fread("autoimmune-disorders/source-data/stage_numeric.txt")
cosmic3.meta <- fread('Documents/scripts/github/lbr-pcs-immunotherapy/source-data/cosmic3.meta.txt', header = T)
actionable_list <- fread("autoimmune-disorders/source-data/vogel_lamarca_kendre.actionable.txt", header=T)
actionable_list <- rbind.data.frame(actionable_list, cbind.data.frame("gene"="TP53", "actionable_type"="point"))


NO_AUTO = "non-PSC/IBD"
YES_AUTO = "PSC/IBD"
PSC_colors = c("#67A9CF","#EF8A62" )

#### LOAD COHORT ####

  samples_raw <- fread("autoimmune-disorders/source-data/samples.tsv")
  samples_raw$TMB_bin[samples_raw$TMB <= 10] <- "TMB < 10"
  samples_raw$TMB_bin[samples_raw$TMB > 10] <- "TMB > 10"
  HYPERMUTATORS <- c("BTC0025", "BTC0045", "BTC0013")
  
  clinical_raw <- fread("autoimmune-disorders/source-data/clinical_data.txt")
  clinical_raw$location[clinical_raw$location == 5] <- NA
  clinical_raw$location[clinical_raw$location == 1] <- 2
  # 0=Male, 1=Female; 
  # Other autoimmune disease (0=PSC only, 1=UC/PSC, 2=Crohn's/PSC, 3=other)
  # Primary Location (0=Extrahepatic/distal, 1=Hilar, 2=Perihilar/Klatskin, 3=Intrahepatic, 4=Gallbladder, 5=Unknown)
  
  treatments_raw <- fread("autoimmune-disorders/source-data/treatments.txt")
  names(treatments_raw)[5] <- "chemo_pre_sequencing"
  summary_raw <- fread("temp_data/BTC.summary.csv")
  
  # translate some binaries
  summary_raw$ploidy_bin[summary_raw$ploidy < 2.25] <- "diploid"
  summary_raw$ploidy_bin[summary_raw$ploidy >= 2.25] <- "polyploid"
  summary_raw$mmr_bin <- "MMRp"
  summary_raw$mmr_bin[summary_raw$mmr_score >= 4] <- "MMRD"

  ## `stage_plotting` codes numeric values for the roman-numeral staging
  clinical <- left_join(clinical_raw, stage_plotting, by=c("stage"="stage"))

  samples <-  left_join(samples_raw, clinical, by=c("Donor"="ID"))
  samples <-  left_join(samples, treatments_raw, by=c("Tumor"="BTC_ID"))
  samples <- left_join(samples, summary_raw, by=c("Tumor"="tumour", "normal"="normal"))
  
  samples$final_status[samples$final_status == "No"] <- NO_AUTO
  samples$final_status[samples$final_status == "Yes"] <- YES_AUTO
  
 # names(MDAnderson_raw)[c(1,3)] <- c("donor","actionable_type")

  MDA.co <- fread('autoimmune-disorders/source-data/MDAnderson.cohort.txt')
  MDA.co$donor <- paste0('MDA',MDA.co$patient_id)
  MDA.co$final_status <- 'PSC/IBD'
  MDA.co$cohort <- "MDA"
  MDA.co$NGS <- "TAR"
  #Other autoimmune disease (0=PSC only, 1=UC/PSC, 2=Crohn's/PSC, 3=other)
  #Primary Location (0=Extrahepatic/distal, 1=Hilar, 2=Perihilar/Klatskin, 3=Intrahepatic, 4=Gallbladder, 5=Unknown)
  MDA.co$location[MDA.co$location == 5] <- NA
  MDA.co$location[MDA.co$location == 1] <- 2
  
  MDA.var.raw <- fread('autoimmune-disorders/source-data/MDAnderson.variants.txt')
  MDA.var.raw$donor <- paste0('MDA',MDA.var.raw$patient_number)
  
  mutation_type.actionable_type <- fread("legresley/source-data/mutation_type.actionable_type.txt")
  #fusion must be confirmed with RNA
  mutation_type.actionable_type <- mutation_type.actionable_type[mutation_type.actionable_type$actionable_type != "fusion"] 
  
  fusions_raw <- fread('temp_data/BTC.star.fusions.txt', header=F)
  names(fusions_raw)[1] <- 'donor'
  fusions <- separate(fusions_raw, V3, c("gene","gene2"), sep = "--")
  
  all_somatic_variants_raw <- fread("temp_data/BTC.somatic.variants.csv")
  
  COSMIC3_raw <-  fread('autoimmune-disorders/source-data/BTC.signatures.COSMIC3.csv')
  COSMIC3_raw$sampleName <- gsub(x = COSMIC3_raw$sampleName, pattern = "_SBS", replacement = "" )
  row.names(COSMIC3_raw) <- COSMIC3_raw$sampleName
  

  
#### PSC/IBD ACTIONABLE ####
  
  
  psc_actionable <- samples %>% filter(final_status == YES_AUTO) %>% select(Tumor, normal, location) %>% 
    left_join(all_somatic_variants_raw, by=c("Tumor"="tumour", "normal"="normal"))

  PMH.var <- get_actionable_mutations(psc_actionable, actionable_list, mutation_type.actionable_type )
  PMH.var$mutation_type[!is.na(PMH.var$aa_context)] <- get_first_aa(PMH.var)
  
  PMH.var <- PMH.var %>% select(donor, gene, actionable_type, mutation_type)
  
  # add MSI
  actionable_MSI <- summary_raw[summary_raw$mmr_score >= 4, c("donor","mmr_bin","mmr_bin","mmr_bin")] 
  names(actionable_MSI) <- c("donor","gene" ,"actionable_type", "mutation_type")
  actionable_MSI$gene <- 'MSI'
  actionable_MSI$actionable_type <- 'HIGH'
  actionable_MSI$mutation_type <- 'MSI-H'
  
  PMH.var <- rbind.data.frame(PMH.var, actionable_MSI)
  
  #need to add samples with no hits back in
  PMH.var <- samples %>% filter(final_status == YES_AUTO) %>% select(Donor, location, NGS) %>% left_join(PMH.var, by=c("Donor"="donor"))
  PMH.var$cohort <- 'PMH'

  MDA.var <- MDA.var.raw %>% filter(gene %in% c(actionable_list$gene, "TP53")) %>% select(donor, gene, actionable_type, variant)
  MDA.var <- MDA.co %>% select(donor,  location, NGS) %>% left_join(MDA.var, by=c("donor"="donor"))
  MDA.var$cohort <- 'MDA'
  
  names(PMH.var) <- names(MDA.var) 
  actionable_variants <- rbind.data.frame(PMH.var, MDA.var)
  
  length(unique(actionable_variants$donor[actionable_variants$gene == "TP53"]))
  length(unique(actionable_variants$donor))
  
  actionable_variants <- actionable_variants[actionable_variants$gene != "TP53",]
  
  actionable_variants <- make_oncoplot_factor_order(actionable_variants)
  
  actionable_variants_frequency <- get_gene_frequencies(actionable_variants)
  actionable_variants_frequency$frequency[is.na(actionable_variants_frequency$gene)] <- NA
  actionable_samples <- unique(actionable_variants[,c('donor_factor','location', 'cohort', 'NGS')])
  
  
  my_theme <- theme_classic(base_size = fontSize)+  theme(
    text=element_text(family="Arial"),
    legend.title = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.spacing = unit(0, 'cm'),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color='transparent'), #transparent plot bg
    legend.background = element_rect(color = NA,fill=NA), #transparent legend bg
    legend.box.background = element_rect(color = NA,fill=NA) #transparent legend panel
  )
  
  
  TOP_fz_plot <- 
    ggplot(actionable_variants_frequency, aes(x=frequency*100, y=gene_factor)) + 
    geom_bar(stat = "identity")+
    scale_y_discrete(limits=rev, breaks=na.omit(unique(actionable_variants$gene_factor)))+ 
    scale_fill_manual(values=PSC_colors) +
    #guides(fill="none")+  
    labs(x="Freq. (%)") + 
    my_theme+
    theme(
      legend.position = 'bottom',
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      plot.margin = unit(c(0,1,0,-1), "lines"),
      axis.title.y =element_blank(),
      axis.text.y =element_blank(),
      axis.ticks.y =element_blank(),
    ) 
  
  my_theme <- my_theme +  theme(
    
    axis.title.x =element_blank(),
    axis.ticks.x =element_blank(),

    plot.margin = unit(c(0,1,0,1), "lines"),
    
    
    axis.text.x=element_blank(),
    legend.direction="horizontal",
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
  )
  
  
  location_plot <- 
    ggplot(actionable_samples, aes(x=donor_factor,y="Location",fill=as.factor(location))) + 
    geom_tile() +
    guides(fill="none")+
   scale_fill_manual(labels=c("eCCA", "pCCA", "iCCA", "Gall"), values=c("#CC4C02","#525252","#238443","#88419D"), na.translate = F) +
    my_theme+  theme(   axis.title.y=element_blank())    +guides(fill=guide_legend(nrow=1,byrow=TRUE))
  
  my_theme <- my_theme +  theme(strip.text = element_blank())
  
  cohort_plot <- 
    ggplot(actionable_samples, aes(x=donor_factor,y="Cohort",fill=as.factor(cohort))) + 
    geom_tile() +
    guides(fill="none")+
    my_theme+  theme(   axis.title.y=element_blank())    +guides(fill=guide_legend(nrow=1,byrow=TRUE))
  
  
  TOP_VAR_PLOT <- 
    ggplot(actionable_variants , aes(x=donor_factor, y=gene_factor, fill=actionable_type)) + 
      geom_tile(aes(fill=actionable_type), na.rm = TRUE) +
      labs(x="Sample",y="",fill="Mutation Type") + 
      my_theme + 
      scale_fill_discrete( na.value = "white", na.translate = F) +
      theme( legend.position = 'bottom',
        axis.text.y = element_text(face = "italic") )+
   #    guides(fill="none")+
      scale_y_discrete(limits=rev, breaks=na.omit(unique(actionable_variants$gene_factor)) )
  
  
  gap = -0.08
  skinny = 0.12
  width_gap = -0.4
  
  png("autoimmune-disorders/results/actionable_onco.png", width = 2000, height = 800) #, bg = "transparent"

  plot_grid(location_plot, NULL, NULL, NULL, NULL, NULL,
            cohort_plot, NULL, NULL, NULL, NULL, NULL,
            TOP_VAR_PLOT, NULL, TOP_fz_plot, # NULL, NULL, NULL,

            rel_heights = c(skinny, gap,skinny, gap, 0.7  ), 
            rel_widths = c(1,width_gap,0.8 ),
            ncol = 3,  align = 'hv', axis = "rlbt")
  
  dev.off()
  
  # actionable pcs/ibd tallies and stats
  
  psc_tally <- actionable_samples %>% group_by(NGS) %>% tally()
  sum(psc_tally$n)
  
  gene_tally <-  actionable_variants %>% group_by(gene) %>% tally()
  variant_tally <-  actionable_variants %>% group_by(gene, actionable_type, variant) %>% tally()
  variant_tally %>% filter(gene == 'KRAS') %>% arrange(variant, -n)
  variant_tally %>% filter(gene == 'ERBB2') %>% arrange(variant, -n)
  variant_tally %>% filter(gene == 'BRAF') %>% arrange(variant, -n)
  variant_tally %>% filter(gene == 'FGFR2') %>% arrange(variant, -n)
  
#  write.table(
#    top_variants_by_sample %>% filter(gene %in% c("TP53","KRAS") & final_status == YES_AUTO) %>% select(Tumor, donor, final_status, gene, copy_number, mutation_class, mutation_type),
#    file = "autoimmune-disorders/results/TP53andKRAS.txt",
#    row.names = FALSE, quote = FALSE, sep = "\t"
#  )
  
#### DRIVERS ####
  
  samples %>% group_by(final_status) %>% tally()
  cohort_size = length(samples$donor)
  VARIANT_FREQUENCY_CUTOFF = 0.05
  
  TOP_GENE_COUNT = round(VARIANT_FREQUENCY_CUTOFF * cohort_size)
  
  driver_zhang_raw <- fread("legresley/source-data/zhang.list.txt", header=T)
  driver_zhang_raw <- driver_zhang_raw[driver_zhang_raw$oncokb == 'YES',]
  driver_zhang <- driver_zhang_raw$gene
  
  all_somatic_variants <- samples %>% select(Tumor, normal, final_status) %>% 
    left_join(all_somatic_variants_raw, by=c("Tumor"="tumour", "normal"="normal"))
  
  variants_tally <- all_somatic_variants %>% group_by(gene)  %>% 
    filter(gene %in% driver_zhang) %>% tally()

  top_onco_variants <- variants_tally %>% filter(n > TOP_GENE_COUNT & !is.na(gene))

  top_variants_by_sample <- all_somatic_variants %>% 
    filter(gene %in% top_onco_variants$gene )
  
  no_hit_samples <- samples$Donor[samples$Donor %ni% unique(top_variants_by_sample$donor)]

  no_hit_rows <- cbind.data.frame(unique(all_somatic_variants[all_somatic_variants$donor %in% no_hit_samples, c(1:5)]),  t(rep(NA, length(names(all_somatic_variants)) -5 )))
  names(no_hit_rows) <- names(all_somatic_variants)
  top_variants_by_sample <- rbind.data.frame(top_variants_by_sample, no_hit_rows)
  
  top_variants_by_sample <- make_oncoplot_factor_order(top_variants_by_sample)
  #set x-axis of plot by variant counts
  samples$donor_factor    <- factor(samples$donor,  levels = levels(top_variants_by_sample$donor_factor))
  
  #make sure all the samples in the cohort are still in data.frame
  cohort_size == length(levels(samples$donor_factor))

  top_variants_by_sample_by_pathway <- left_join(top_variants_by_sample, driver_zhang_raw, by=c("gene"="gene"))
  top_variants_by_sample_by_pathway$pathway[top_variants_by_sample_by_pathway$pathway == ''] <- 'Other'
  top_variants_by_sample_by_pathway$pathway[is.na(top_variants_by_sample_by_pathway$pathway)] <- 'Other'
  
  top_variants_by_sample_by_pathway$pathway[top_variants_by_sample_by_pathway$pathway == "Epigenetic Modifiers"] <- "Epigenetic\nModifiers"
  top_variants_by_sample_by_pathway$pathway[top_variants_by_sample_by_pathway$pathway == "Kinase-RAS/RAF"] <- "Kinase\nRAS/RAF"
  
  top_variants_by_sample_by_pathway$pathway_factor    <- 
    factor(top_variants_by_sample_by_pathway$pathway,  levels = c("p53", "Cell cycle", "Kinase\nRAS/RAF", 
                                                                       "Epigenetic\nModifiers", "PI3K","TGF-B", "WNT","Notch","Myc", "Other"))
  
  TOP_variant_frequency <- get_gene_frequencies_final_status_pathway(top_variants_by_sample_by_pathway)
  TOP_variant_frequency$frequency[is.na(TOP_variant_frequency$gene)] <- NA


#### SIGNATURES ####

COSMIC3 <- COSMIC3_raw[,-c("Residuals", "AbsResiduals","RSS", "N.Mutations", "sampleName")]
COSMIC3 <- COSMIC3 / COSMIC3_raw$N.Mutations

# removed signatures w.o any calls
COSMIC3_called_cols <- COSMIC3[,colSums(COSMIC3) > 0]
COSMIC3_called <- COSMIC3[,..COSMIC3_called_cols]
COSMIC3_called$tumour <- COSMIC3_raw$sampleName
COSMIC3_joined <- samples %>% select(Tumor, Donor, final_status, location, TMB_bin) %>% left_join(COSMIC3_called, by=c("Tumor"="tumour"))

COSMIC3_melt <- melt(COSMIC3_joined, id.vars = c("Tumor",  "Donor", "final_status","location", "TMB_bin"))
COSMIC3_melt <- left_join(COSMIC3_melt, cosmic3.meta, by=c("variable"="signature_name"))
COSMIC3_melt <- COSMIC3_melt[COSMIC3_melt$Donor %ni% HYPERMUTATORS,]

COSMIC3_sample_meta <- COSMIC3_melt %>% group_by(Donor, final_status , meta, TMB_bin ) %>% summarise(meta_sum = sum(value))
COSMIC3_sample_meta_order <- COSMIC3_sample_meta %>% group_by(meta ) %>% summarise(meta_mean = mean(meta_sum))
COSMIC3_sample_meta_order <- COSMIC3_sample_meta_order[order(-COSMIC3_sample_meta_order$meta_mean),]

COSMIC3_melt$meta_factor    <- factor(COSMIC3_melt$meta,  levels = COSMIC3_sample_meta_order$meta)

COSMIC3_melt <- COSMIC3_melt[COSMIC3_melt$meta_factor != "unknown",]
COSMIC3_melt$variable <- gsub("Signature.","SBS", COSMIC3_melt$variable)

COSMIC3_meta <- COSMIC3_melt %>% group_by(final_status, meta_factor, Donor, TMB_bin) %>% summarise(meta_sum = sum(value))
COSMIC3_melt <- left_join(COSMIC3_melt, COSMIC3_meta)

COSMIC3_melt_sig <- COSMIC3_melt %>% filter(meta_factor %in% c('aflatoxin', 'APOBEC', 'SBS17' ))

png("autoimmune-disorders/results/signatures_heatmap.png", width = 1200, height = 750) #, bg = "transparent"

ggplot(COSMIC3_melt, aes(x=Donor, y=variable, fill=log(value))) + 
  facet_grid(meta_factor    ~ final_status, scales="free", space="free", switch="y") +
  labs(fill="log %")+
  geom_tile(color="white")+ 
  scale_fill_gradient2(midpoint = -4, low = "white", mid = "cornflowerblue", high = "red" ,na.value="white")+ 
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

dev.off()
  

png("autoimmune-disorders/results/signatures.barplot.png", width = 1200, height = 900, bg = "transparent")

  ggplot(COSMIC3_meta, aes(y=meta_factor,fill=final_status, x=meta_sum)) + 
    #facet_grid(meta_factor    ~ ., scales="free", space="free") +
  #  geom_violin() + 
    geom_boxplot() +
    labs(fill='')+
    scale_fill_manual(values=PSC_colors)+ 
    theme_classic(base_size = 26)+  theme(
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
    #  axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      strip.text.y = element_blank(),
      legend.position = 'top',
    )  

  dev.off()

#### DO PLOT ####
  
  my_theme <- theme_classic(base_size = fontSize)+  theme(
    text=element_text(family="Arial"),
    legend.title = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.spacing = unit(0, 'cm'),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color='transparent'), #transparent plot bg
    legend.background = element_rect(color = NA,fill=NA), #transparent legend bg
    legend.box.background = element_rect(color = NA,fill=NA) #transparent legend panel
  )

  
  
TOP_fz_plot <- 
  ggplot(TOP_variant_frequency, aes(x=frequency*100, y=gene_factor, fill=final_status)) + 
  geom_bar(stat = "identity")+
  facet_grid(pathway_factor ~ ., scales="free", space="free", switch="y") +
  
  scale_y_discrete(limits=rev, breaks=na.omit(unique(top_variants_by_sample$gene_factor)))+ 
  scale_fill_manual(values=PSC_colors) +
#  guides(fill="none")+  
  labs(x="Freq. (%)") + 
  my_theme+
  theme(
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        plot.margin = unit(c(0,1,0,-1), "lines"),
        axis.title.y =element_blank(),
        axis.text.y =element_blank(),
        axis.ticks.y =element_blank(),
        strip.text.y = element_blank()
  ) 
  
  my_theme <- my_theme +  theme(
    
    axis.title.x =element_blank(),
    axis.ticks.x =element_blank()
  )
  
  

TMB_PLOT <- 
  ggplot(samples, aes(x=donor_factor, y=TMB, fill=TMB_bin)) + 
  geom_bar(stat = "identity")+
 # guides(fill="none")+
  labs(y="TMB",fill='') + 
  facet_grid(. ~ final_status, scales="free_x", space="free_x") +
  scale_fill_manual( values=c("#969696","#B2182B")) +
  my_theme + 
  theme(
    strip.text = element_blank(),
    panel.grid.major.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "lines"),
    axis.text.x =element_blank()
  )+ scale_y_log10() 

my_theme <- my_theme +  theme(
  plot.margin = unit(c(0,1,0,1), "lines"),
  
  
  axis.text.x=element_blank(),
  legend.direction="horizontal",
  panel.grid.major = element_blank(), #remove major gridlines
  panel.grid.minor = element_blank(), #remove minor gridlines
)




location_plot <- 
  ggplot(samples, aes(x=donor_factor,y="Location",fill=as.factor(location))) + 
  facet_grid(. ~ final_status, scales="free_x", space="free_x") +
  geom_tile() +
  scale_fill_manual(labels=c("eCCA", "pCCA", "iCCA", "Gall"), values=c("#CC4C02","#525252","#238443","#88419D"), na.translate = F) +
  my_theme+  theme(   axis.title.y=element_blank())    +guides(fill=guide_legend(nrow=1,byrow=TRUE))

my_theme <- my_theme +  theme(strip.text.x = element_blank())

TOP_VAR_PLOT <- 
  ggplot(top_variants_by_sample_by_pathway , aes(x=donor_factor, y=gene_factor, fill=mutation_type)) + 
    geom_tile(aes(fill=mutation_type), na.rm = TRUE) +
    labs(x="Sample",y="",fill="Mutation Type") + 
    facet_grid(pathway_factor ~ final_status, scales="free", space="free", switch="y") +
    my_theme + 
    scale_fill_discrete( na.value = "white", na.translate = F) +
    theme(  axis.text.y = element_text(face = "italic"), strip.text.y.left = element_text(angle = 0) )+
   # guides(fill="none")+
    scale_y_discrete(limits=rev, breaks=na.omit(unique(top_variants_by_sample$gene_factor)) )

variant_class_legend =  get_plot_component(TOP_VAR_PLOT, pattern = "guide-box-right" )
TOP_VAR_PLOT = TOP_VAR_PLOT + guides(fill="none")

my_theme <- my_theme +  theme( axis.title.y=element_blank())


sex_plot <- 
  ggplot(samples, aes(x=donor_factor,y="Sex",fill=as.factor(sex))) + 
  facet_grid(. ~ final_status, scales="free_x", space="free_x") +
  geom_tile() +
  scale_fill_manual(labels=c("male","female"), values=c("#3690C0","#FEE090")) +
  my_theme 

stage_plot <- 
  ggplot(samples, aes(x=donor_factor,y="Stage",fill=stage_bin)) + 
  facet_grid(. ~ final_status, scales="free_x", space="free_x") +
  geom_tile() +
  scale_fill_manual(breaks =  c("I", "II", "III", "IV"), values=c("grey","#35978F","#BF812D", "#8073AC")) +
  my_theme+  theme(    legend.spacing.y = unit(0.5, 'cm') )   +
  guides(fill = guide_legend(byrow = TRUE))

ploidy_plot <-  
  ggplot(samples, aes(x=donor_factor,y="Ploidy",fill=ploidy_bin)) + 
  facet_grid(. ~ final_status, scales="free_x", space="free_x") +
  geom_tile() +
  
  scale_fill_manual( values=c("#ffb14e","#ea5f94"), na.translate = F) +
  
  my_theme 


mmrd_plot <- 
  ggplot(samples, aes(x=donor_factor,y="MMRd",fill=mmr_bin)) + 
  facet_grid(. ~ final_status, scales="free_x", space="free_x") +
  geom_tile() +
  scale_fill_manual( values=c("#FA9FB5", "#4D004B" ), na.translate = F) +
  my_theme


sigs_plot <- 
ggplot(COSMIC3_melt_sig, aes(x=Donor, y=meta_factor, fill=log(meta_sum))) + 
  facet_grid(.    ~ final_status, scales="free", space="free", switch="y") +
  labs(fill="log %")+
  geom_tile(color="white")+ 
  guides(fill = guide_colourbar(
    theme = theme(legend.key.width = unit(10, "lines")), 
    ticks = TRUE, ticks.colour = "black")
  ) +
  scale_fill_gradient2(midpoint = -4, 
                       breaks = c(log(0.01),  log(0.1)), 
                       labels = c('1%', '10%'), 
                       low = "white", 
                       mid = "cornflowerblue",
                       high = "red" ,
                       na.value="white") +
  my_theme

gap = -0.095
skinny = 0.14
width_gap = -0.55

png("autoimmune-disorders/results/oncoplot.png", width = 1600, height = 1100* 1.4) #, bg = "transparent"

plot_grid(
          location_plot, NULL, NULL, NULL, NULL, NULL,
          sex_plot, NULL, NULL, NULL, NULL, NULL,
          stage_plot, NULL, NULL, NULL, NULL, NULL,
          ploidy_plot , NULL, NULL,NULL, NULL, NULL,
          mmrd_plot, NULL, NULL, NULL, NULL, NULL,
          
          TMB_PLOT, NULL, NULL, NULL, NULL, NULL,
          sigs_plot, NULL, NULL, NULL, NULL, NULL,
          TOP_VAR_PLOT, NULL, TOP_fz_plot, # NULL, NULL, NULL,
          variant_class_legend, NULL, NULL,
          
          rel_heights = c(skinny,gap,
                          skinny, gap,
                          skinny, gap,
                          skinny,gap,  
                          skinny,gap,
                          
                          0.2, gap, 
                          0.2, gap,
                          0.7* 1.4,  0.2  ), 
          rel_widths = c(1,width_gap,0.8 ),
          ncol = 3,  align = 'hv', axis = "rlbt")

dev.off()


  cosmic3_meta_lm <- lm(meta_sum ~  meta_factor * final_status , data = COSMIC3_meta)
  summary(cosmic3_meta_lm)
  
  cosmic3_lm <- lm(value ~  variable * final_status , data = COSMIC3_melt)
  summary(cosmic3_lm)
  
  
  cosmic3_lm <- lm(value ~  variable * TMB_bin , data = COSMIC3_melt)
  summary(cosmic3_lm)
  
  cosmic3_meta_lm <- lm(meta_sum ~  meta_factor *  TMB_bin , data = COSMIC3_meta)
  summary(cosmic3_meta_lm)
  
  
  #### STATS ####
  
  ## what are the top mutated genes in PSC/IBD?
  
  TOP_variant_frequency %>% filter(final_status == YES_AUTO) %>% arrange(desc(frequency))
  
  
  ## Which drivers are enriched in PCS/IBD ?
  
  
  
  ## ZEROES
  for(this_gene_name in c("BAP1", "NRAS", "ARID2", "TGFBR1")){
    cat("\n",this_gene_name,"\n")
    this_gene <- top_variants_by_sample_by_pathway %>% filter(gene == this_gene_name) %>% group_by(final_status) %>% tally()
    
    this_gene <- rbind.data.frame(this_gene, cbind.data.frame("final_status"="PSC/IBD","n"=0))
    
    this_gene$no_mutation <- totals$n - this_gene$n
    print(fisher.test(this_gene[,-1]))
  }
  
  # LOW FREQUENCIES
  totals <- samples %>% group_by(final_status) %>% tally()
  
  for(this_gene_name in c("MYC","RBM10","FBXW7", "KMT2C")){
    cat("\n",this_gene_name)
    this_gene <- top_variants_by_sample_by_pathway %>% filter(gene == this_gene_name) %>% group_by(final_status) %>% tally()
    
    this_gene$no_mutation <- totals$n - this_gene$n
    print(chisq.test(this_gene[,-1]))
    Xsq <-chisq.test(this_gene[,-1])
    Xsq$expected
    Xsq$observed
    print(fisher.test(this_gene[,-1]))
  }
  
  ## Do PCS/IBD cases have higher TMB?
  
  mean(samples$TMB[samples$final_status == YES_AUTO]) / mean(samples$TMB[samples$final_status == NO_AUTO])
  
  quantile(samples$TMB[samples$final_status == YES_AUTO], probs = c(0,0.25,0.5,0.75,1))
  quantile(samples$TMB[samples$final_status == NO_AUTO], probs = c(0,0.25,0.5,0.75,1))
  
  wilcox.test(samples$TMB[samples$final_status == YES_AUTO],
              samples$TMB[samples$final_status == NO_AUTO])
  
  wilcox.test(samples$TMB[samples$final_status == YES_AUTO & samples$mmr_bin == "MMRp" & samples$Donor != "BTC0025"],
              samples$TMB[samples$final_status == NO_AUTO & samples$mmr_bin == "MMRp" & samples$Donor != "BTC0025"])
  
  
  samples$TMB_FDA[samples$TMB >= 10] <- "True"
  samples$TMB_FDA[samples$TMB < 10] <- "False"
  
  FDA_TMB_tally <- samples %>% group_by(TMB_FDA, final_status) %>% tally()
  FDA_TMB_TABLE <- dcast(FDA_TMB_tally, formula=TMB_FDA ~ final_status, fill = 0)
  row.names(FDA_TMB_TABLE) <- FDA_TMB_TABLE$TMB_FDA
  FDA_TMB_TABLE <- FDA_TMB_TABLE[,-1]
  
  fisher.test(FDA_TMB_TABLE)
  