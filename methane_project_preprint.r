library(tidyverse)
library(ggrepel)
library(phyloseq)
library(vegan)
library(viridis) #Viridis color library
library(microbiome) # data analysis and visualisation
library(plotly)
library("IRdisplay")
library(ggpmisc) #to use stat_poly_eq
library(ggpubr)
library("gridExtra")
library("factoextra")
library(rioja) # plotting poackages for tabular bubbleplots
library(reshape2)
library(RColorBrewer) # nice color options
library(randomcoloR) #generate random colors
library(ggthemes)
# set size of the plot
options(repr.plot.width=14, repr.plot.height=12)

save.image()

theme_glab <- function(base_size = 30,
                    base_family = "",
                    base_line_size = base_size / 180,
                    base_rect_size = base_size / 180) {
   
    font <- "Helvetica" #assign font family up front
   
    theme_bw(base_size = base_size,
                base_family = base_family,
                base_line_size = base_line_size) %+replace%
    theme(
        legend.background =  element_blank(),
        legend.title =       element_text(color = rgb(100, 100, 100, maxColorValue = 255),
                                          size = rel(0.65),
                                         hjust = 0),
        legend.text =        element_text(color = rgb(100, 100, 100, maxColorValue = 255),
                                          size = rel(0.65)),
        legend.key.size =    unit(0.8, "lines"),
     
      plot.title = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        hjust = 0),
       
      axis.title = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        size = rel(0.65)),
      axis.text = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        size = rel(0.65)),
       
      plot.caption = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        size = rel(0.7),
        hjust = 1),
       
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.border = element_rect(fill = NA, colour = rgb(100, 100, 100, maxColorValue = 255)),

     
      complete = TRUE
    )
}

#Uploading Phyloseq Object (Basili et al., 2024 https://doi.org/10.1371/journal.pone.0308756)

cr_meth<-readRDS("cr17_18_normalized_physeq.rds")
cr_meth

# Methane cycling prokaryotes subset
meth_prok <- subset_taxa(cr_meth,

                            Order == "Methanobacteriales" | #Archaea
                            Order == "Methanococcales" | #Archaea
                            Order == "Methanomicrobiales" | #Archaea
                            Order == "Methanosarcinales" | #Archaea
                            Order == "Methanocellales" | #Archaea
                            Order == "Methanopyrales" | #Archaea
                            Order == "Methanomassiliicoccales" | #Archaea
                            Class == "Methanomicrobia" | #Archaea
                            Phylum == "Methanobacteriota" | #Archaea
                            Genus == "Candidatus Methanophagales" | #Archaea
                            Genus == "Candidatus Methanoxibalbensis" | #Archaea
                            Genus == "Candidatus Methanospirare" | #Archaea
                            Order == "Methanophagales" | #Archaea
                            Class == "Methanosarcinia" | #Archaea
                            Class == "ANME-1" | #Archaea
                            Class == "Methanococci" | #Archaea
                            Class == "Methanopyri" | #Archaea
                            Order == "Methanofastidiosales" | #Archaea
                            Class == "Methanocellia" | #Archaea
                            Class == "Methanonatronarchaeia" | #Archaea
                            Genus == "Candidatus Methanoliparia" | #Archaea
                            Genus == "Candidatus Methanodesulfokores" | #Archaea
                            Class == "Methanomethylia" | #Archaea
                            Class == "Methanomethylicia" | #Archaea
                            Order == "Methanomethyliales" | #Archaea
                            Order == "Methanonatronoarchaeales" | #Archaea
                            Order == "Methanonatronaechaeales" | #Archaea
                            Order == "Candidatus Methanomethyliales" | #Bacteria
                            Order == "Methylococcales" | #Bacteria
                            Order == "Methylacidiphilales" | #Bacteria
                            Phylum == "Methylomirabilota" | #Bacteria
                            Family == "Methylomirabilaceae" | #Bacteria
                            Family == "Methylophagaceae" | #Bacteria
                            Family == "Methyloligellaceae" | #Bacteria
                            Family == "Methylophilaceae" | #Bacteria
                            Genus == "Methyloterrigena" | #Bacteria
                            Genus == "Methylarcula" | #Bacteria
                            Genus == "Methylomusa" | #Bacteria
                            Genus == "Methylobrevis" | #Bacteria
                            Family == "Methylopilaceae" | #Bacteria
                            Genus == "Candidatus Methyloumidiphilus" | #Bacteria
                            Genus == "Methyloterricola" | #Bacteria
                            Genus == "Candidatus Methylospira" | #Bacteria
                            Genus == "Methylomagnum" | #Bacteria
                            Genus == "Methylotetracoccus" | #Bacteria
                            Genus == "Methyloparacoccus" | #Bacteria
                            Genus == "Methylobacterium" | #Bacteria
                            Genus == "Methylosinus" | #Bacteria
                            Genus == "Methylorubrum" | #Bacteria
                            Genus == "Methylocapsa" | #Bacteria
                            Genus == "Methylocella" | #Bacteria
                            Genus == "Methylocystis" | #Bacteria
                            Genus == "Methyloferula" | #Bacteria
                            Genus == "Methylorosula" | #Bacteria
                            Genus == "Methylovirgula" | #Bacteria
                            Genus == "Marinosulfonomonas"  #Bacteria
)
meth_prok
sum(readcount(meth_prok))
sum(readcount(meth_prok))/sum(readcount(cr_meth))

#cr_env_data<-read.csv("../../matteo_run/bac_arc_marging2/dataset/DCO_BMS_clean_dataset_Selci_CR17_CR18.csv",header=TRUE,sep=",",row.names=3)
#sweet Baby Jesus, make it works without errors...
#dataset updated at 17 Jan 2023
cr_env_data<-read.csv("DCO_BMS__envdata CR17_CR18_preprint.csv",header=TRUE,sep=",",row.names=3)
head(cr_env_data)

cr_meth_prok<-merge_phyloseq(meth_prok,sample_data(cr_env_data))
cr_meth_prok

cr_meth_prok2<-cr_meth_prok

cr_meth_prok_ra = transform_sample_counts(cr_meth_prok2, function(x){x / sum(x)}) 
cr_meth_prok_ra

is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))

otu_table(cr_meth_prok_ra)[is.nan(otu_table(cr_meth_prok_ra))] <- 0

#removing sites with 0% abundance
cr_meth_prok_ra2 = subset_samples(cr_meth_prok_ra, sample_names(cr_meth_prok_ra) != "HAS")
cr_meth_prok_ra2 = subset_samples(cr_meth_prok_ra2, sample_names(cr_meth_prok_ra2) != "QHS1")
cr_meth_prok_ra2 = subset_samples(cr_meth_prok_ra2, sample_names(cr_meth_prok_ra2) != "QHS2")

## Agglomerate at a specific taxonomic level at the Genus level
prok_ra_genus = tax_glom(cr_meth_prok_ra2, "Genus", NArm = F)
prok_ra_family = tax_glom(cr_meth_prok_ra2, "Family", NArm = F)
prok_ra_order = tax_glom(cr_meth_prok_ra2, "Order", NArm = F)
prok_ra_class = tax_glom(cr_meth_prok_ra2, "Class", NArm = F)
prok_ra_phyla = tax_glom(cr_meth_prok_ra2, "Phylum", NArm = F)
prok_ra_kingdom = tax_glom(cr_meth_prok_ra2, "Kingdom", NArm = F)

plot_bar(prok_ra_phyla, fill="Phylum", x="ID2", 
         title = "Taxonomic level:Phylum") + scale_y_continuous(labels=c(0,25,50,75,100)) +
theme_glab() + labs(x="") + 
theme(plot.title = element_text(size=8),legend.position = "bottom",axis.text.x = element_text(angle = 90, 
vjust = 0.4, hjust=1.2),axis.text=element_text(size=8),
legend.key.size = unit(4, 'mm'),legend.text = element_text(size=8))

#############################################################################################################

plot_bar(prok_ra_class, fill="Class", x="ID2", 
         title = "Taxonomic level:Class") + scale_y_continuous(labels=c(0,25,50,75,100)) +
theme_glab() + labs(x="") +
theme(plot.title = element_text(size=8),legend.position = "bottom",axis.text.x = element_text(angle = 90, 
vjust = 0.4, hjust=1.2),axis.text=element_text(size=8),
legend.key.size = unit(4, 'mm'),legend.text = element_text(size=8))

#############################################################################################################

plot_bar(prok_ra_order, fill="Order", x="ID2", 
         title = "Taxonomic level:Order") + scale_y_continuous(labels=c(0,25,50,75,100)) +
theme_glab() + labs(x="") +
theme(plot.title = element_text(size=8),legend.position = "bottom",axis.text.x = element_text(angle = 90, 
vjust = 0.4, hjust=1.2),axis.text=element_text(size=8),
legend.key.size = unit(4, 'mm'),legend.text = element_text(size=8))
#ggsave("plot/barplot/mathane_groups_orders.svg",width=16,height=12)
#############################################################################################################

plot_bar(prok_ra_family, fill="Family", x="ID2", 
         title = "Taxonomic level:Family") + scale_y_continuous(labels=c(0,25,50,75,100)) +
theme_glab() + labs(x="") +
theme(plot.title = element_text(size=8),legend.position = "bottom",axis.text.x = element_text(angle = 90, 
vjust = 0.4, hjust=1.2),axis.text=element_text(size=8),
legend.key.size = unit(4, 'mm'),legend.text = element_text(size=8))


samples_order_geounit2<-list("CWF","CWS","CWF2","CWS2","LWF","LWS","LWF2","LWS2","PSS","PSS2","RSS","SIS",
"CIF","CIS","CZF","CZS","ERF","ERS","ERS2","HAF","HAS","MCF","MCS","QHS1","QHF2",
"QHS2","SM","EPF","EPS","ESF","SR","BCF","BCF2","BSF","BSS","BWF","CHF","CHS",
"CHS2","CLF","CLS","CVF","CVS","CVS2","LBS","LBS2","LHF","LHS","LPF","LPS","RRF"
,"RRS","SCF","SCS","SCF2","SCS2","GES","GES2","LEF","LES","LEF2","LES2","PXF",
"PXS","RCF","RCS","XFF","YRF","YRF2","BQ","BRF1","BRS1","BRF2","BRS2","CYF","CYS",
"ETF","ETS","FAF","FAS","HN","MTF","PFF","PFS","PLS","QNF","QNS","RV","SLF","SLS",
"STS","TCF","TCS","VCS")

order_methane<-plot_bar(prok_ra_order, fill="Order", x="ID2", 
         title = "Taxonomic level:Order") +
theme_glab() + labs(y="Relative Abundance (%)") +
scale_fill_manual(na.value="black",values=c(colorRampPalette(RColorBrewer::brewer.pal(11,"PRGn"))(15))) +
scale_y_continuous(labels=c(0,25,50,75,100)) +
theme(plot.title = element_text(size=18),legend.position = "bottom",axis.text.x = element_text(angle = 90, 
vjust = 0.4, hjust=1.2),axis.text=element_text(size=15),
legend.key.size = unit(8, 'mm'),legend.text = element_text(size=18)) 

order_methane$data$ID2 <- factor(order_methane$data$ID2, levels = samples_order_geounit2)

print(order_methane)
#ggsave("plot/barplot/mathane_geounit2_orders.svg",width=16,height=12)

tax_table(prok_ra_family)[18,5]<-"Order Methanomassiliicoccales; Family Unknown"
tax_table(prok_ra_family)[20,5]<-"Order Methanobacteriales; Family Unknown"
tax_table(prok_ra_family)[21,5]<-"Order Methanofastidiosales; Family Unknown"
tax_table(prok_ra_family)[32,5]<-"Order ANME-1; Family Unknown"
tax_table(prok_ra_family)[34,5]<-"Class Methanomicrobia; Family Unknown"
tax_table(prok_ra_family)[37,5]<-"Order Methanosarcinales; Family Unknown"
tax_table(prok_ra_family)[38,5]<-"Order Methanocellales; Family Unknown"
tax_table(prok_ra_family)[39,5]<-"Order Methanomicrobiales; Family Unknown"
tax_table(prok_ra_family)

family_methane<-plot_bar(prok_ra_family, fill="Family", x="ID2", 
         title = "Taxonomic level:Family") +
theme_glab() + labs(y="Relative Abundance (%)") +
scale_fill_manual(na.value="black",values=c(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(39))) +
scale_y_continuous(labels=c(0,25,50,75,100)) +
theme(plot.title = element_text(size=18),legend.position = "bottom",axis.text.x = element_text(angle = 90, 
vjust = 0.4, hjust=1.2),axis.text=element_text(size=15),legend.text = element_text(size=12)) +
theme(legend.position = "bottom")

family_methane$data$ID2 <- factor(family_methane$data$ID2, levels = samples_order_geounit2)

print(family_methane)
#ggsave("plot/barplot/mathane_geounit2_orders.svg",width=16,height=12)

# Saving a .csv with combined taxonomy and ASVs abundances for the selected taxonomic level
summary_object <-subset_samples(prok_ra_family, sample_type== "fluid")
summary_otu <- as.matrix(t(otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
summary_comb <- cbind(summary_otu, summary_tax)
#write.csv(summary_comb, "../summary_ra_abundance_241113/fluid_ch4_family_summary.csv")

cr_meth_prok_ra3 = transform_sample_counts(cr_meth, function(x){x / sum(x)}) 

# Methane cycling prokaryotes
meth_prok2 <- subset_taxa(cr_meth_prok_ra3, 
                            Order == "Methanobacteriales" | #Archaea
                            Order == "Methanococcales" | #Archaea
                            Order == "Methanomicrobiales" | #Archaea
                            Order == "Methanosarcinales" | #Archaea
                            Order == "Methanocellales" | #Archaea
                            Order == "Methanopyrales" | #Archaea
                            Order == "Methanomassiliicoccales" | #Archaea
                            Class == "Methanomicrobia" | #Archaea
                            Phylum == "Methanobacteriota" | #Archaea
                            Genus == "Candidatus Methanophagales" | #Archaea
                            Genus == "Candidatus Methanoxibalbensis" | #Archaea
                            Genus == "Candidatus Methanospirare" | #Archaea
                            Order == "Methanophagales" | #Archaea
                            Class == "Methanosarcinia" | #Archaea
                            Class == "ANME-1" | #Archaea
                            Class == "Methanococci" | #Archaea
                            Class == "Methanopyri" | #Archaea
                            Order == "Methanofastidiosales" | #Archaea
                            Class == "Methanocellia" | #Archaea
                            Class == "Methanonatronarchaeia" | #Archaea
                            Genus == "Candidatus Methanoliparia" | #Archaea
                            Genus == "Candidatus Methanodesulfokores" | #Archaea
                            Class == "Methanomethylia" | #Archaea
                            Class == "Methanomethylicia" | #Archaea
                            Order == "Methanomethyliales" | #Archaea
                            Order == "Methanonatronoarchaeales" | #Archaea
                            Order == "Methanonatronaechaeales" | #Archaea
                            Order == "Candidatus Methanomethyliales" | #Bacteria
                            Order == "Methylococcales" | #Bacteria
                            Order == "Methylacidiphilales" | #Bacteria
                            Phylum == "Methylomirabilota" | #Bacteria
                            Family == "Methylomirabilaceae" | #Bacteria
                            Family == "Methylophagaceae" | #Bacteria
                            Family == "Methyloligellaceae" | #Bacteria
                            Family == "Methylophilaceae" | #Bacteria
                            Genus == "Methyloterrigena" | #Bacteria
                            Genus == "Methylarcula" | #Bacteria
                            Genus == "Methylomusa" | #Bacteria
                            Genus == "Methylobrevis" | #Bacteria
                            Family == "Methylopilaceae" | #Bacteria
                            Genus == "Candidatus Methyloumidiphilus" | #Bacteria
                            Genus == "Methyloterricola" | #Bacteria
                            Genus == "Candidatus Methylospira" | #Bacteria
                            Genus == "Methylomagnum" | #Bacteria
                            Genus == "Methylotetracoccus" | #Bacteria
                            Genus == "Methyloparacoccus" | #Bacteria
                            Genus == "Methylobacterium" | #Bacteria
                            Genus == "Methylosinus" | #Bacteria
                            Genus == "Methylorubrum" | #Bacteria
                            Genus == "Methylocapsa" | #Bacteria
                            Genus == "Methylocella" | #Bacteria
                            Genus == "Methylocystis" | #Bacteria
                            Genus == "Methyloferula" | #Bacteria
                            Genus == "Methylorosula" | #Bacteria
                            Genus == "Methylovirgula" | #Bacteria
                            Genus == "Marinosulfonomonas"  #Bacteria
)
meth_prok2

cr_meth_prok_ra3<-merge_phyloseq(meth_prok2,sample_data(cr_env_data))
cr_meth_prok_ra3

#% of methane related ASV
#ASV methane cycling: 2918
#ASV before filtering: 54354 
(2918/54354)*100

save.image()

meth_ra_order2 = tax_glom(cr_meth_prok_ra3, "Order", NArm = F)
meth_ra_family2 = tax_glom(cr_meth_prok_ra3, "Family", NArm = F)

#removing sites with 0% abundance
meth_ra_order3 = subset_samples(meth_ra_order2, ! sample_names(meth_ra_order2) %in% c("HAS","QHS1","QHS2"))
meth_ra_family3 = subset_samples(meth_ra_family2, ! sample_names(meth_ra_family2) %in% c("HAS","QHS1","QHS2"))

# Saving a .csv with combined taxonomy and ASVs abundances for the selected taxonomic level
summary_object <-meth_ra_order3
summary_otu <- as.matrix(t(otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
summary_comb <- cbind(summary_otu, summary_tax)
#write.csv(summary_comb, "../summary_ra_abundance/real_ra_meth_order_summary.csv")

sum(otu_table(meth_ra_order3)[17,2:15])

write.csv(sample_data(meth_ra_order3), "../sample_data_meth_ra_order3_16S.csv")

meth_ra_order3

order_methane3<-plot_bar(meth_ra_order3, fill="Order", x="ID2", 
         title = "Taxonomic level:Order") +
scale_fill_manual(na.value="black",values=c(colorRampPalette(RColorBrewer::brewer.pal(11,"PRGn"))(16))) +
theme_glab() + labs(x="") + labs(y="Relative Abundance (%)") +
scale_y_continuous(labels=c(0,10,20,30,40)) +
theme(plot.title = element_text(size=8),legend.position = "bottom",axis.text.x = element_text(angle = 90, 
vjust = 0.4, hjust=1.2),axis.text=element_text(size=14),
legend.key.size = unit(4, 'mm'),legend.text = element_text(size=14))

order_methane3$data$ID2 <- factor(order_methane3$data$ID2, levels = samples_order_geounit2)

print(order_methane3)

#ggsave("plot/barplot/mathane_geounit2_orders.svg",width=16,height=12)

meth_ra_family3

tax_table(meth_ra_family3)

tax_table(meth_ra_family3)[18,5]<-"Order Methanomassiliicoccales; Family Unknown"
tax_table(meth_ra_family3)[22,5]<-"Order Methanobacteriales; Family Unknown"
tax_table(meth_ra_family3)[23,5]<-"Order Methanofastidiosales; Family Unknown"
tax_table(meth_ra_family3)[31,5]<-"Order ANME-1; Family Unknown"
tax_table(meth_ra_family3)[32,5]<-"Class Methanomicrobia; Family Unknown"
tax_table(meth_ra_family3)[37,5]<-"Order Methanosarcinales; Family Unknown"
tax_table(meth_ra_family3)[38,5]<-"Order Methanocellales; Family Unknown"
tax_table(meth_ra_family3)[39,5]<-"Order Methanomicrobiales; Family Unknown"
tax_table(meth_ra_family3)

family_methane3<-plot_bar(meth_ra_family3, fill="Family", x="ID2", 
         title = "Taxonomic level:Family") +
scale_fill_manual(na.value="black",values=c(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(39))) +
theme_glab() + labs(x="") +
scale_y_continuous(labels=c(0,10,20,30,40)) + labs(y="Relative Abundance (%)") +
theme(plot.title = element_text(size=8),legend.position = "bottom",axis.text.x = element_text(angle = 90, 
vjust = 0.4, hjust=1.2),axis.text=element_text(size=14),
legend.key.size = unit(4, 'mm'),legend.text = element_text(size=14))

family_methane3$data$ID2 <- factor(family_methane3$data$ID2, levels = samples_order_geounit2)

print(family_methane3)

#ggsave("plot/barplot/mathane_geounit2_orders.svg",width=16,height=12)

ggarrange(order_methane3,order_methane, nrow=2, ncol=1, widths = c(2,2),
  heights = c(1,2), common.legend = TRUE)

#ggsave("../figures/barplot_order.svg", width=14, height=14)

ggarrange(family_methane3,family_methane, nrow=2, ncol=1, widths = c(2,2),
  heights = c(1,2), common.legend = TRUE)

#ggsave("../figures/barplot_order.svg", width=14, height=14)

save.image()

set.seed(1000) # to make the analysis reproducible
options(repr.plot.width=16, repr.plot.height=12)

cr_meth_prok2

cr_meth_prok3 <- subset_samples(cr_meth_prok2, 
                 sample_names(cr_meth_prok2) != "QHS1") #contain NAs

cr_meth_prok3 <- subset_samples(cr_meth_prok3, 
                 sample_names(cr_meth_prok3) != "QHS2") #contain NAs

cr_meth_prok3 <- subset_samples(cr_meth_prok3, 
                 sample_names(cr_meth_prok3) != "HAS") #contain only 0

cr_meth_prok3 <- subset_samples(cr_meth_prok3, 
                 sample_names(cr_meth_prok3) != "XFF") #contain only 0

cr_meth_prok3 <- subset_samples(cr_meth_prok3, 
                 sample_names(cr_meth_prok3) != "RRS") #contain only 0

cr_meth_prok3

sample_data(cr_meth_prok3)$ID2

#Uploading methane aerobic and anaerobic groups
aerobic_anaerobic_list<-read.csv("../../matteo_run/bac_arc_marging2/csv/summary/cr_meth_aerobic_anaerobic.csv",header=TRUE,sep=",")

ch4_clump_dt1<-read.csv("../../matteo_run/bac_arc_marging2/dataset/ch4_clumped_compilation.csv",header=TRUE,sep=",")
#ch4_clump_dt2<-read.csv("../../matteo_run/bac_arc_marging2/dataset/ch4_clumped_compilation2.csv",header=TRUE,sep=",")

cr_methane_cycle<-read.csv("../../../cr17_cr18_metagenomes/mifaser/analysis_in_progress/jupyter/cr_methane_cycle.csv",header=TRUE,sep=",")
cr_methane_cycle<-cr_methane_cycle[,-1]
head(cr_methane_cycle)

#subsetting only for mcr and mmo
cr_methane_cyc<-cr_methane_cycle[c(2,4,6,7),c(3,6:75)]
row.names(cr_methane_cyc)<-cr_methane_cyc[,1]
cr_methane_cyc<-t(cr_methane_cyc[,-1])
cr_methane_cyc<-cbind(cr_methane_cyc,rownames(cr_methane_cyc))
colnames(cr_methane_cyc)<-c("mcra","mdh","mmo","pmo-amo","ID1")
rownames(cr_methane_cyc)<- NULL
cr_methane_cyc <- data.frame(cr_methane_cyc) %>%
  select(ID1, everything())
head(cr_methane_cyc)

i<-c(2:5) 

cr_methane_cyc[ , i] <- apply(cr_methane_cyc[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))
head(cr_methane_cyc)

ch4_clump_dt2.1<-list(ch4_clump_dt1,aerobic_anaerobic_list) %>%
                    reduce(left_join, by="ID1")
ch4_clump_dt2.1<-list(ch4_clump_dt2.1,cr_methane_cyc) %>%
                    reduce(left_join, by="ID1")
ch4_clump_dt2.1$mcr<-as.numeric(ch4_clump_dt2.1$mcr)
ch4_clump_dt2.1$mmo<-as.numeric(ch4_clump_dt2.1$mmo)

ch4_clump_dt2.1

subset(ch4_clump_dt2.1, ID1 == "BQF")

subset(ch4_clump_dt2.2, ID1 == "BQF")

ch4_clump_dt2.2<-list(ch4_clump_dt2.1,cr_env_data[,c(3,21:136)]) %>%
                    reduce(left_join, by="ID1")
ch4_clump_dt2.2

ch4_clump_dt2.1

save.image()

#WHITICAR PLOT + METHANOGE CYCLE EC Numbers for mcr and mmo 

ggplot(ch4_clump_dt2.1, aes(x=delta_d,y=delta_13c)) + 
geom_point(data = subset(ch4_clump_dt2.1, reference != c("Deeply sourced seeps")),alpha=0.5,
    aes(shape=reference,fill=samplename2),size=6,stroke=.3) + 
geom_point(data = subset(ch4_clump_dt2.1, reference %in% c("Deeply sourced seeps")),size=8,
           shape=21,fill="red2",stroke=.3) + 
scale_shape_manual(values = c(22:25)) +
scale_fill_viridis(discrete=TRUE) +
#geom_text_repel(aes(label= ID1), size=3) + 
scale_x_continuous(limits = c(-455,-45),breaks=seq(-450, -50, by = 50)) + 
scale_y_reverse(limits = c(-5, -115),breaks=seq(0, -120, by = -20)) +
#xlim(-450,-50) + 
#ylim(-5, -115) +
#guides(fill = guide_legend(override.aes = list(shape = 21) ),
#            shape = guide_legend(override.aes = list(fill = "black"))) +
theme_glab() +
labs(title="", x=expression(δD[CH[4]]*"(‰)"), y=expression(δ^{13}*"C"[CH[4]]*"(‰)"))
#ggsave("plot/geochem/ch4D_ch4_d13_mmo.svg",width=14,height=8)

ggplot(ch4_clump_dt2.2, aes(x=delta_13c,y=delta_13ch3d)) + 
geom_point(data = subset(ch4_clump_dt2.2, reference != c("Deeply sourced seeps")),
    aes(shape=reference,fill=samplename2), alpha=.4, size=10,stroke=.3) + 
geom_point(data = subset(ch4_clump_dt2.2, reference %in% c("Deeply sourced seeps")),size=10,
           shape=21,fill="red2",stroke=.3) + 
scale_shape_manual(values = c(22:25)) +
scale_fill_viridis(discrete=TRUE) +
#geom_text(data = subset(ch4_clump_dt2.2, reference %in% c("Deeply sourced seeps")),aes(label= ID1), size=4, hjust=1.5) + 
#scale_x_reverse() + 
#scale_y_reverse() +
xlim(-80,0) + 
#ylim(-1.5, 8) +
guides(fill = guide_legend(override.aes = list(shape = 21)),
            shape = guide_legend(override.aes = list(fill = "black"))) +
theme_glab() +
theme(legend.position = "right") +
labs(title="", x=expression(δ^{13}*"C"[CH[4]]*"(‰)"), y=expression("\u0394"^{13}*CH[3]*D*"(‰)"))

#ggsave("plot/geochem/figures_paper/clumped_isotopes_13c.svg",width=10,height=10)

ggplot(ch4_clump_dt2.2, aes(x=delta_d,y=delta_13ch3d)) + 
geom_point(data = subset(ch4_clump_dt2.2, reference != c("Deeply sourced seeps")),
    aes(shape=reference,fill=samplename2), alpha=.4, size=10,stroke=.3) + 
geom_point(data = subset(ch4_clump_dt2.2, reference %in% c("Deeply sourced seeps")),size=10,
           shape=21,fill="red2",stroke=.3) + 
scale_shape_manual(values = c(22:25)) +
scale_fill_viridis(discrete=TRUE) +
#geom_text(data = subset(ch4_clump_dt2.2, reference %in% c("Deeply sourced seeps")),aes(label= ID1), size=4, hjust=1.5) + 
#scale_x_reverse() + 
#scale_y_reverse() +
#xlim(-80,0) + 
#ylim(-1.5, 8) +
guides(fill = guide_legend(override.aes = list(shape = 21)),
            shape = guide_legend(override.aes = list(fill = "black"))) +
theme_glab() +
theme(legend.position = "right") +
labs(title="", x=expression(δD[CH[4]]*"(‰)"), y=expression("\u0394"^{13}*CH[3]*D*"(‰)"))

#ggsave("plot/geochem/figures_paper/clumped_isotopes_D.svg",width=10,height=10)

ch4_clump_dt2.2F<-subset(ch4_clump_dt2.2, sample_type %in% c("F"))
ch4_clump_dt2.2S<-subset(ch4_clump_dt2.2, sample_type %in% c("S"))

ggplot(ch4_clump_dt2.2F, aes(x=delta_13c,y=delta_13ch3d)) + 
geom_point(data = subset(ch4_clump_dt2.2, reference != c("Deeply sourced seeps")),
    aes(shape=reference),fill="darkgrey", alpha=.6, size=8,stroke=.3) + 
geom_point(data = subset(ch4_clump_dt2.2F, reference %in% c("Deeply sourced seeps")),aes(size=aerobic*100),fill="#2c7bb6",
           shape=21,stroke=.3, alpha=.8) + 
scale_size(range = c(5,25)) +
scale_shape_manual(values = c(22:25)) +
scale_fill_viridis(discrete=T) +
#geom_text(data = subset(ch4_clump_dt2.2F, reference %in% c("Deeply sourced seeps")),aes(label= ID1), size=4, hjust=1.5) + 
#scale_x_reverse() + 
#scale_y_reverse() +
xlim(-80,0) + 
#ylim(-1.5, 8) +
guides(fill = guide_legend(override.aes = list(shape = 21)),
            shape = guide_legend(override.aes = list(fill = "black"))) +
theme_glab() +
theme(legend.position = "none") +
labs(title="", x=expression(δ^{13}*"C"[CH[4]]*"(‰)"), y=expression("\u0394"^{13}*CH[3]*D*"(‰)"))

#ggsave("../../matteo_run/bac_arc_marging2/plot/geochem/figures_paper/clumped_aerobes_fluid.svg",width=10,height=10)

ggplot(ch4_clump_dt2.2F, aes(x=delta_13c,y=delta_13ch3d)) + 
geom_point(data = subset(ch4_clump_dt2.2, reference != c("Deeply sourced seeps")),
    aes(shape=reference),fill="darkgrey", alpha=.6, size=8,stroke=.3) + 
geom_point(data = subset(ch4_clump_dt2.2F, reference %in% c("Deeply sourced seeps")),aes(size=anaerobic*100),fill="#fdae61",
           shape=21,stroke=.3, alpha=.8) + 
scale_size(range = c(5,25)) +
scale_shape_manual(values = c(22:25)) +
scale_fill_viridis(discrete=T) +
#geom_text(data = subset(ch4_clump_dt2.2, reference %in% c("This study")),aes(label= ID1), size=4, hjust=1.5) + 
#scale_x_reverse() + 
#scale_y_reverse() +
xlim(-80,0) + 
#ylim(-1.5, 8) +
guides(fill = guide_legend(override.aes = list(shape = 21)),
            shape = guide_legend(override.aes = list(fill = "black"))) +
theme_glab() +
theme(legend.position = "none") +
labs(title="", x=expression(δ^{13}*"C"[CH[4]]*"(‰)"), y=expression("\u0394"^{13}*CH[3]*D*"(‰)"))

#ggsave("../../matteo_run/bac_arc_marging2/plot/geochem/figures_paper/clumped_anaerobes_fluid.svg",width=10,height=10)

summary(ch4_clump_dt2.2F$mmo)

summary(log((ch4_clump_dt2.2F$mmo)))

log10(30332)

subset(ch4_clump_dt2.2F, ID1 == "BQF")$

ggplot(ch4_clump_dt2.2F, aes(x=delta_13c,y=delta_13ch3d)) + 
geom_point(data = subset(ch4_clump_dt2.2, reference != c("Deeply sourced seeps")),
    aes(shape=reference),fill="darkgrey", alpha=.6, size=8,stroke=.3) + 
geom_point(data = subset(ch4_clump_dt2.2F, reference %in% c("Deeply sourced seeps")),aes(size=log10(as.numeric(mmo))),fill="#2c7bb6",
           shape=21,stroke=.3, alpha=.8) + 
scale_size(range = c(1,25)) +
scale_shape_manual(values = c(22:25)) +
scale_fill_viridis(discrete=T) +
#geom_text(data = subset(ch4_clump_dt2.2F, reference %in% c("Deeply sourced seeps")),aes(label= ID1), size=4, hjust=1.5) + 
#scale_x_reverse() + 
#scale_y_reverse() +
xlim(-80,0) + 
#ylim(-1.5, 8) +
guides(fill = guide_legend(override.aes = list(shape = 21)),
            shape = guide_legend(override.aes = list(fill = "black"))) +
theme_glab() +
theme(legend.position = "none") +
labs(title="mmo", x=expression(δ^{13}*"C"[CH[4]]*"(‰)"), y=expression("\u0394"^{13}*CH[3]*D*"(‰)"))

#ggsave("../../matteo_run/bac_arc_marging2/plot/geochem/figures_paper/clumped_mmo_fluid.svg",width=10,height=10)

ggplot(ch4_clump_dt2.2F, aes(x=delta_13c,y=delta_13ch3d)) + 
geom_point(data = subset(ch4_clump_dt2.2, reference != c("Deeply sourced seeps")),
    aes(shape=reference),fill="darkgrey", alpha=.6, size=8,stroke=.3) + 
geom_point(data = subset(ch4_clump_dt2.2F, reference %in% c("Deeply sourced seeps")),aes(size=log10(as.numeric(mdh))),fill="#2c7bb6",
           shape=21,stroke=.3, alpha=.8) + 
scale_size(range = c(1,25)) +
scale_shape_manual(values = c(22:25)) +
scale_fill_viridis(discrete=T) +
#geom_text(data = subset(ch4_clump_dt2.2F, reference %in% c("Deeply sourced seeps")),aes(label= ID1), size=4, hjust=1.5) + 
#scale_x_reverse() + 
#scale_y_reverse() +
xlim(-80,0) + 
#ylim(-1.5, 8) +
guides(fill = guide_legend(override.aes = list(shape = 21)),
            shape = guide_legend(override.aes = list(fill = "black"))) +
theme_glab() +
theme(legend.position = "none") +
labs(title="mdh", x=expression(δ^{13}*"C"[CH[4]]*"(‰)"), y=expression("\u0394"^{13}*CH[3]*D*"(‰)"))

#ggsave("../../matteo_run/bac_arc_marging2/plot/geochem/figures_paper/clumped_mdh_fluid.svg",width=10,height=10)

ggplot(ch4_clump_dt2.2F, aes(x=delta_13c,y=delta_13ch3d)) + 
geom_point(data = subset(ch4_clump_dt2.2, reference != c("Deeply sourced seeps")),
    aes(shape=reference),fill="darkgrey", alpha=.6, size=8,stroke=.3) + 
geom_point(data = subset(ch4_clump_dt2.2F, reference %in% c("Deeply sourced seeps")),aes(size=log10(as.numeric(pmo.amo))),fill="#2c7bb6",
           shape=21,stroke=.3, alpha=.8) + 
scale_size(range = c(1,25)) +
scale_shape_manual(values = c(22:25)) +
scale_fill_viridis(discrete=T) +
#geom_text(data = subset(ch4_clump_dt2.2F, reference %in% c("Deeply sourced seeps")),aes(label= ID1), size=4, hjust=1.5) + 
#scale_x_reverse() + 
#scale_y_reverse() +
xlim(-80,0) + 
#ylim(-1.5, 8) +
guides(fill = guide_legend(override.aes = list(shape = 21)),
            shape = guide_legend(override.aes = list(fill = "black"))) +
theme_glab() +
theme(legend.position = "none") +
labs(title="pmo.amo", x=expression(δ^{13}*"C"[CH[4]]*"(‰)"), y=expression("\u0394"^{13}*CH[3]*D*"(‰)"))

#ggsave("../../matteo_run/bac_arc_marging2/plot/geochem/figures_paper/clumped_pmoamo_fluid_legend.svg",width=10,height=10)


ggplot(ch4_clump_dt2.2F, aes(x=delta_13c,y=delta_13ch3d)) + 
geom_point(data = subset(ch4_clump_dt2.2, reference != c("Deeply sourced seeps")),
    aes(shape=reference),fill="darkgrey", alpha=.6, size=8,stroke=.3) + 
geom_point(data = subset(ch4_clump_dt2.2F, reference %in% c("Deeply sourced seeps")),aes(size=log10(as.numeric(mcra))),fill="#fdae61",
           shape=21,stroke=.3, alpha=.8) + 
scale_size(range = c(1,25)) +
scale_shape_manual(values = c(22:25)) +
scale_fill_viridis(discrete=T) +
#geom_text(data = subset(ch4_clump_dt2.2, reference %in% c("This study")),aes(label= ID1), size=4, hjust=1.5) + 
#scale_x_reverse() + 
#scale_y_reverse() +
xlim(-80,0) + 
#ylim(-1.5, 8) +
guides(fill = guide_legend(override.aes = list(shape = 21)),
            shape = guide_legend(override.aes = list(fill = "black"))) +
theme_glab() +
theme(legend.position = "none") +
labs(title="mcra", x=expression(δ^{13}*"C"[CH[4]]*"(‰)"), y=expression("\u0394"^{13}*CH[3]*D*"(‰)"))

#ggsave("../../matteo_run/bac_arc_marging2/plot/geochem/figures_paper/clumped_mcra_fluid.svg",width=10,height=10)

save.image()

ch4_clump_dt2.2S

ggplot(ch4_clump_dt2.2S, aes(x=delta_13c,y=delta_13ch3d)) + 
geom_point(data = subset(ch4_clump_dt2.2, reference != c("Deeply sourced seeps")),
    aes(shape=reference),fill="darkgrey", alpha=.6, size=8,stroke=.3) + 
geom_point(data = subset(ch4_clump_dt2.2S, reference %in% c("Deeply sourced seeps")),aes(size=log10(mmo)),fill="#2c7bb6",
           shape=21,stroke=.3, alpha=.8) + 
scale_size(range = c(1,25)) +
scale_shape_manual(values = c(22:25)) +
scale_fill_viridis(discrete=T) +
#geom_text(data = subset(ch4_clump_dt2.2S, reference %in% c("Deeply sourced seeps")),aes(label= ID1), size=4, hjust=1.5) + 
#scale_x_reverse() + 
#scale_y_reverse() +
xlim(-80,0) + 
#ylim(-1.5, 8) +
guides(fill = guide_legend(override.aes = list(shape = 21)),
            shape = guide_legend(override.aes = list(fill = "black"))) +
theme_glab() +
theme(legend.position = "none") +
labs(title="mmo", x=expression(δ^{13}*"C"[CH[4]]*"(‰)"), y=expression("\u0394"^{13}*CH[3]*D*"(‰)"))

#ggsave("../../matteo_run/bac_arc_marging2/plot/geochem/figures_paper/clumped_mmo_sediment.svg",width=10,height=10)

ggplot(ch4_clump_dt2.2S, aes(x=delta_13c,y=delta_13ch3d)) + 
geom_point(data = subset(ch4_clump_dt2.2, reference != c("Deeply sourced seeps")),
    aes(shape=reference),fill="darkgrey", alpha=.6, size=8,stroke=.3) + 
geom_point(data = subset(ch4_clump_dt2.2S, reference %in% c("Deeply sourced seeps")),aes(size=log10(as.numeric(mdh))),fill="#2c7bb6",
           shape=21,stroke=.3, alpha=.8) + 
scale_size(range = c(1,25)) +
scale_shape_manual(values = c(22:25)) +
scale_fill_viridis(discrete=T) +
#geom_text(data = subset(ch4_clump_dt2.2S, reference %in% c("Deeply sourced seeps")),aes(label= ID1), size=4, hjust=1.5) + 
#scale_x_reverse() + 
#scale_y_reverse() +
xlim(-80,0) + 
#ylim(-1.5, 8) +
guides(fill = guide_legend(override.aes = list(shape = 21)),
            shape = guide_legend(override.aes = list(fill = "black"))) +
theme_glab() +
theme(legend.position = "none") +
labs(title="mdh", x=expression(δ^{13}*"C"[CH[4]]*"(‰)"), y=expression("\u0394"^{13}*CH[3]*D*"(‰)"))

#ggsave("../../matteo_run/bac_arc_marging2/plot/geochem/figures_paper/clumped_mdh_sediment.svg",width=10,height=10)

ggplot(ch4_clump_dt2.2S, aes(x=delta_13c,y=delta_13ch3d)) + 
geom_point(data = subset(ch4_clump_dt2.2, reference != c("Deeply sourced seeps")),
    aes(shape=reference),fill="darkgrey", alpha=.6, size=8,stroke=.3) + 
geom_point(data = subset(ch4_clump_dt2.2S, reference %in% c("Deeply sourced seeps")),aes(size=log10(as.numeric(pmo.amo))),fill="#2c7bb6",
           shape=21,stroke=.3, alpha=.8) + 
scale_size(range = c(1,25)) +
scale_shape_manual(values = c(22:25)) +
scale_fill_viridis(discrete=T) +
#geom_text(data = subset(ch4_clump_dt2.2S, reference %in% c("Deeply sourced seeps")),aes(label= ID1), size=4, hjust=1.5) + 
#scale_x_reverse() + 
#scale_y_reverse() +
xlim(-80,0) + 
#ylim(-1.5, 8) +
guides(fill = guide_legend(override.aes = list(shape = 21)),
            shape = guide_legend(override.aes = list(fill = "black"))) +
theme_glab() +
theme(legend.position = "none") +
labs(title="pmo.amo", x=expression(δ^{13}*"C"[CH[4]]*"(‰)"), y=expression("\u0394"^{13}*CH[3]*D*"(‰)"))

#ggsave("../../matteo_run/bac_arc_marging2/plot/geochem/figures_paper/clumped_pmoamo_sediment.svg",width=10,height=10)


ggplot(ch4_clump_dt2.2S, aes(x=delta_13c,y=delta_13ch3d)) + 
geom_point(data = subset(ch4_clump_dt2.2, reference != c("Deeply sourced seeps")),
    aes(shape=reference),fill="darkgrey", alpha=.6, size=8,stroke=.3) + 
geom_point(data = subset(ch4_clump_dt2.2S, reference %in% c("Deeply sourced seeps")),aes(size=log10(as.numeric(mcra))),fill="#fdae61",
           shape=21,stroke=.3, alpha=.8) + 
scale_size(range = c(1,25)) +
scale_shape_manual(values = c(22:25)) +
scale_fill_viridis(discrete=T) +
#geom_text(data = subset(ch4_clump_dt2.2, reference %in% c("This study")),aes(label= ID1), size=4, hjust=1.5) + 
#scale_x_reverse() + 
#scale_y_reverse() +
xlim(-80,0) + 
#ylim(-1.5, 8) +
guides(fill = guide_legend(override.aes = list(shape = 21)),
            shape = guide_legend(override.aes = list(fill = "black"))) +
theme_glab() +
theme(legend.position = "none") +
labs(title="mcra", x=expression(δ^{13}*"C"[CH[4]]*"(‰)"), y=expression("\u0394"^{13}*CH[3]*D*"(‰)"))

#ggsave("../../matteo_run/bac_arc_marging2/plot/geochem/figures_paper/clumped_mcra_sediment.svg",width=10,height=10)

ggplot(data=methane_df, aes(x=trench,y=rc_ra)) + 
    geom_point(size=10,aes(fill=geologic_unit2),alpha=.6,shape=21,stroke=.2) +
#    geom_point(aes(x=province,y=as.numeric(co2_d13),size=as.numeric(co2)),
#               fill="gray",alpha=.5,shape=21,stroke=.2) +
#scale_fill_viridis(discrete=TRUE) +
#scale_size_continuous(limits = c(0,14000), range = c(10,25), breaks = c(1,100,1000,10000)) +
    scale_fill_manual(values = c("#f03b20","#3182bd", "#238443", "#bd0026"),labels=c("OFB","QD","SR","VR")) +
#    scale_x_discrete(name="",labels=c("QD","SR","OFB", "VR")) +
    labs(x="Distance from trench (Km)", y="Rc/Ra") +

    theme_glab() +
    theme(axis.text.x = element_text(colour = "grey3", size=20,angle=30, vjust =0.5, hjust=.5)) +
  theme(aspect.ratio = 1/2) 
ggsave("../figures/geochemistry/Supplementary_Figure1.svg", width=14,height=10)

ggplot(data=methane_df %>% arrange(-ch4)) + 
    geom_point(aes(x=factor(geologic_unit2,level=c('Quaternary deposits','Sedimentary rocks',
                                             'Ocean floor basalts','Volcanic rocks')),
                   y=ch4_d13,size=(ch4*16.04),fill=geologic_unit2),alpha=.6,shape=21,stroke=.2) +
#    geom_point(aes(x=province,y=as.numeric(co2_d13),size=as.numeric(co2)),
#               fill="gray",alpha=.5,shape=21,stroke=.2) +
    scale_size_continuous(limits = c(0,14000), range = c(10,25), breaks = c(1,100,1000,10000)) +
    scale_fill_manual(values = c("#f03b20","#3182bd", "#238443", "#bd0026"),labels=c("QD","SR","OFB", "VR")) +
    scale_x_discrete(name="",labels=c("QD","SR","OFB", "VR")) +

    scale_y_reverse() +
    labs(size=expression(CH[4]*-(mg/Kg)), x="", y=expression(δ^{13}*"C")) +
    theme_glab() +
    theme(axis.text.x = element_text(colour = "grey3", size=20,angle=30, vjust =0.5, hjust=.5)) +
  theme(aspect.ratio = 1/2) + ylim(-15,-85)


ggsave("../figures/methane_geologic_unit.svg", height=10, width=14)

summary(methane_df$ch4)
summary(methane_df$co2)
summary(methane_df$h2)

ggtern::ggtern(data=methane_df,aes(x=co2/5,y=he*1000, z=ch4*10)) +
geom_point(aes(fill=geologic_unit2),shape=21, size=7, stroke=0.3, alpha=1)+
scale_fill_viridis(discrete=T) +
#scale_fill_manual(values=c("#5ec962","#fde725", "#21918c", "#21918c"))+
#ggtern::tern_limit(T=1.15,L=1.15,R=1.15)+
#guides(fill = guide_legend(override.aes=list(shape=21)))+
labs( x       = expression(CO[2]/5),
      y       = "He*1000",
      z       = expression(CH[4]*"*10"),
    fill="")+
theme_glab()

# For Rendering the Lines, use Segment geometry
lines <- data.frame(x = c(0.5, 0, 0.5), 
                    y = c(0.5, 0.5, 0), 
                    z = c(0, 0.5, 0.5), 
                    xend = c(1.5, 1, 1)/3, 
                    yend = c(1.5, 1, 1)/3, 
                    zend = c(1.5, 1, 1)/3)

# Giggenbach plot of major cations in geothermal fluids

ggtern::ggtern(data=methane_df,aes(x=ca*40.1,y=mg*24.3, z=(na*23)+(k*39.1))) + #
geom_point(aes(fill=province2),shape=21, size=7, stroke=0.3, alpha=1)+
scale_fill_viridis(discrete=T) +
#scale_fill_manual(values=c("#5ec962","#fde725", "#21918c", "#21918c"))+
#ggtern::tern_limit(T=1.15,L=1.15,R=1.15)+
#guides(fill = guide_legend(override.aes=list(shape=21)))+
labs( x       = "Ca",
      y       = "Mg",
      z       = "Na+K",
    fill="") +
geom_segment(data = lines, aes(x=x, y=y, z=z, xend=xend, yend=yend, zend=zend), 
             color = 'black', size = .4) +
theme_glab()
#svg("../figures/ternary_cation.svg", width=10, height=10)


# Giggenbach plot of major water types in geothermal environemnts

ggtern::ggtern(data=subset(methane_df, sample_type == "fluid"),aes(x=so4*96,y=cl*35.45, z=dic*61)) +
geom_point(aes(fill=province2),shape=21, size=7, stroke=0.3, alpha=1)+
scale_fill_viridis(discrete=T) +
#geom_text(aes(label= ID1), size=4) + 
#scale_fill_manual(values=c("#5ec962","#fde725", "#21918c", "#21918c"))+
#ggtern::tern_limit(T=1.15,L=1.15,R=1.15)+
#guides(fill = guide_legend(override.aes=list(shape=21)))+
labs( x       = expression(SO[4]),
      y       = expression(Cl),
      z       = "DIC",
    fill="") +
geom_segment(data = lines, aes(x=x, y=y, z=z, xend=xend, yend=yend, zend=zend), 
             color = 'black', size = .4) +
theme_glab()
#dev.off()

ggtern::ggtern(data=methane_df,aes(x=he3,y=he4_g, z=ch4)) +
geom_point(aes(fill=geologic_unit2),shape=21, size=7, alpha=1, stroke=0.3)+
scale_fill_viridis(discrete=T) +
#scale_shape_manual(values = c(21,22)) +
#scale_fill_manual(values=c("#5ec962","#fde725", "#21918c", "#21918c"))+
#ggtern::tern_limit(T=1.15,L=1.15,R=1.15)+
#guides(fill = guide_legend(override.aes=list(shape=c(21))))+
labs( x       = "3He",
      y       = "4He",
      z       = "CH4",
    fill="")+
theme_glab() 
#best

save.image()

rast_out<-read.csv("dataset/rast_contigs_cr17_18.csv",header=TRUE, sep=",")
row.names(rast_out)<-rast_out[,1]
head(rast_out)

rast_out.m<-as.matrix(rast_out[2:10])
rast_out.m_otu<- otu_table(rast_out.m, taxa_are_rows=F)

rast_out.m_taxa<-data.frame(cbind(colnames(rast_out.m)))
colnames(rast_out.m_taxa)<-c("gene")
row.names(rast_out.m_taxa)<-rast_out.m_taxa[,1]
rast_out.m_taxa<-as.matrix(rast_out.m_taxa)
rast_out.m_taxa<-tax_table(rast_out.m_taxa)

rast_physeq<-phyloseq(sample_data(rast_out),
                     otu_table(rast_out.m_otu),
                     rast_out.m_taxa)
rast_physeq

rast_physeq_ra = transform_sample_counts(rast_physeq, function(x){x / sum(x)}) 
rast_physeq_ra

otu_table(rast_physeq_ra)

plot_bar(rast_physeq, fill="gene", x="ID1", title = "") +
theme_glab() + labs(x="") +
theme(plot.title = element_text(size=8),legend.position = "right",axis.text.x = element_text(angle = 90, 
vjust = 0.6, hjust=1.2),axis.text=element_text(size=12),
legend.key.size = unit(6, 'mm'),legend.text = element_text(size=12))
#ggsave("plot/methane_cicle_absolute_ab.svg", width=20,height=4)

plot_bar(rast_physeq_ra, fill="gene", x="ID1", title = "") +
theme_glab() + labs(x="") +
theme(plot.title = element_text(size=8),legend.position = "right",axis.text.x = element_text(angle = 90, 
vjust = 0.6, hjust=1.2),axis.text=element_text(size=12),
legend.key.size = unit(6, 'mm'),legend.text = element_text(size=12))
#ggsave("plot/methane_cicle_relative_ab.svg", width=20,height=10)

save.image()
