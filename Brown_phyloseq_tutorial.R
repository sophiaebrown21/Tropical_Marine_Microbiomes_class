library(tidyverse) ; packageVersion("tidyverse") # 2.0.0
BiocManager::install("phyloseq")
library(devtools)
library(phyloseq) ; packageVersion("phyloseq") # 1.42.0
library(vegan) ; packageVersion("vegan") # 2.6.4
install.packages("dendextend")
library(dendextend) ; packageVersion("dendextend") # 1.17.1
library(viridis) ; packageVersion("viridis") # 0.6.2
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library("pairwiseAdonis"); packageVersion("pairwiseAdonis") # 0.4.1
library(ggpubr)
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC", force = TRUE)
library(ANCOMBC)
BiocManager::install("microbiome")
library(microbiome)
install.packages("eulerr")
library(eulerr)
library(RColorBrewer)

###libs from LUTZ code
library(ggplot2)
library(readr)
library(data.table) ; packageVersion("data.table") #1.14.8
library(dplyr) ; packageVersion("dplyr") #1.1.1
install.packages("picante")
library(picante) ; packageVersion("picante") #1.8.2

setwd("/Users/sophie/Dropbox/My Mac (Sophia’s MacBook Air)/Desktop/PRJEB32322")

#uploading your taxonomy and count data
asv_tab = read.delim("ASVs_counts.tsv", sep = "\t")
asv_tax = read.delim("ASVs_taxonomy.txt", sep="\t")

#need to construct around your data using your metadata table from Run Selector
sample_info_tab = read.csv("cuttlefish_mircobiome_mbl_COPY_1.csv") #had to clean metadata so the row.names could be read

#need to make columns the rownames all of them
#need to figure out what the row names are for all your data, likely if you followed the tutorial the tax and tab are "X")
#for your sample data file the first row needs to be the SRR that is listed in your ASV files
asv_tab = asv_tab %>% tibble::column_to_rownames("X")
asv_tax = asv_tax %>% tibble::column_to_rownames("X")
sample_df = sample_info_tab %>% tibble::column_to_rownames("Run")

#transforming your taxonomy data into a matrix for entry into phyloseq
asv_mat = as.matrix(asv_tab)
tax_mat = as.matrix(asv_tax)

#transform into phyloseq objects
asv = otu_table(asv_mat, taxa_are_rows=TRUE)
tax=tax_table(tax_mat)
samples=sample_data(sample_df)
#create subsets to work on later; taking out fecal AND just looking at fecal, esophagus and gill
samples_wo_fecal <- subset_samples(samples, Sample_location!="Fecal")

#make a phyloseq object
phy <-phyloseq(asv, tax, samples)

#subsets for later
phy_wofecal <-phyloseq(asv, tax, samples_wo_fecal)
fec = subset_samples(phy, Sample_location=="Fecal")
eso = subset_samples(phy, Sample_location=="Esophagus")
gil = subset_samples(phy, Sample_location=="Gill")

#remove chloroplast and mitochondria from your phyloseq object
phy <- phy %>% subset_taxa( Family!= "Mitochondria" & Class!="Chloroplast" )
phy_wofecal <- phy_wofecal %>% subset_taxa( Family!= "Mitochondria" & Class!="Chloroplast" )
fec<- fec %>% subset_taxa( Family!= "Mitochondria" & Class!="Chloroplast" )
eso <- eso %>% subset_taxa( Family!= "Mitochondria" & Class!="Chloroplast" )
gil <- gil %>% subset_taxa( Family!= "Mitochondria" & Class!="Chloroplast" )

#number of taxa
ntaxa(phy) #16595 removing Mito and Chlor #9428

#number of samples
nsamples(phy) #455
nsamples(phy_wofecal) #172
nsamples(fec) #283
nsamples(eso) #22
nsamples(gil) #22

#number of variables
sample_variables(phy)

###[1] "cf"                "empo_3"            "env_material"      "Host_sample_ID"    "host_subject_id"  
###[6] "Organism"          "Sample_plate"      "Sample_location"   "Sample_well"       "host_body_habitat"
###[11] "host_body_product" "host_body_site"    "Tank"              "Treatment_period"  "Treatment_type"   

##setting up code for rel abundance
# Create a df with all classes and flag the top 20
# (annotations her are from github) 
#DECIDE how many taxa you want to see
topN = 10 # all others will be grouped in "Others"
topTaxa <- phy %>%  # CHOOSE YOUR DATASET !
  psmelt %>% group_by(Family) %>% # melt the table and group by tax level
  summarise(Abundance=sum(Abundance)) %>% # find most overall abundant taxa
  arrange(desc(Abundance)) %>% # order them 
  mutate(aggTaxo=as.factor(case_when( # aggTaxo will become the plot legend
    row_number()<=topN~Family, row_number()>topN~'Others'))) %>% 
  dplyr::select(-Abundance) %>% #flush abundance value
  head(n=topN+1) # +1 to include the Others section! 

#Family                   #aggTaxo
#1Vibrionaceae            Vibrionaceae
#2Methylophagaceae        Methylophagaceae
#3Colwelliaceae           Colwelliaceae
#4Flavobacteriaceae       Flavobacteriaceae
#5Rhodobacteraceae.       Rhodobacteraceae
#6Pseudoalteromonadaceae  Pseudoalteromonadaceae
#7Rubritaleaceae          Rubritaleaceae
#8Enterobacteriaceae      Enterobacteriaceae
#9Crocinitomicaceae       Crocinitomicaceae
#10Pirellulaceae          Pirellulaceae
#11Saccharospirillaceae   Others

### CREATE an object that will be fed into ggplot

plotData_cuttle <- phy %>% # CHOOSE YOUR DATASET !
  transform_sample_counts(function(x) x/sum(x)) %>% psmelt %>% # "melt" the phyloseq object into a table
  inner_join(.,topTaxa, by="Family", copy=TRUE) %>% #use topTaxa to join aggTaxo variable
  aggregate(Abundance~Sample_location+aggTaxo, #aggregate variables of interest
            data=., FUN=sum) #sum "Other" taxa

mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(topN+1) 

##rel abun plot
ggplot(plotData_cuttle, aes(x=Sample_location, y=Abundance, fill=aggTaxo)) +
  geom_bar(stat="identity",position="stack") + # change to "stack" to plot total counts
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~Sample_location, scales = "free") + # RHS of the tilde is columns; add scales="free" to hide samples without counts
  scale_fill_manual(values = mycolors,breaks=c(pull(topTaxa, Family),"Others")) +
  theme(panel.background = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 7),
        legend.title = aggTaxo="Family", element_text(size = 12))

##try with subsets
#fecal subset
topN = 10 
topTaxa_fec <- fec %>%  # CHOOSE YOUR DATASET !
  psmelt %>% group_by(Family) %>% 
  summarise(Abundance=sum(Abundance)) %>% # find most overall abundant taxa
  arrange(desc(Abundance)) %>% # order them 
  mutate(aggTaxo_fec=as.factor(case_when( # aggTaxo will become the plot legend
    row_number()<=topN~Family, row_number()>topN~'Others'))) %>% 
  dplyr::select(-Abundance) %>% 
  head(n=topN+1) 

plotData_fec <- fec %>% # CHOOSE YOUR DATASET !
  transform_sample_counts(function(x) x/sum(x)) %>% psmelt %>% # "melt" the phyloseq object into a table
  inner_join(.,topTaxa_fec, by="Family", copy=TRUE) %>% #use topTaxa to join aggTaxo variable
  aggregate(Abundance~Treatment_type+aggTaxo_fec, #aggregate variables of interest
            data=., FUN=sum) #sum "Other" taxa

ggplot(plotData_fec, aes(x=Treatment_type, y=Abundance, fill=aggTaxo_fec)) +
  geom_bar(stat="identity",position="stack") + # change to "stack" to plot total counts
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~Treatment_type, scales = "free") + # RHS of the tilde is columns; add scales="free" to hide samples without counts
  scale_fill_manual(values = mycolors,breaks=c(pull(topTaxa_fec, Family),"Others")) +
  theme(panel.background = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 12))

topN = 10 
topTaxa_eso <- eso %>%  # CHOOSE YOUR DATASET !
  psmelt %>% group_by(Family) %>% 
  summarise(Abundance=sum(Abundance)) %>% # find most overall abundant taxa
  arrange(desc(Abundance)) %>% # order them 
  mutate(aggTaxo_eso=as.factor(case_when( # aggTaxo will become the plot legend
    row_number()<=topN~Family, row_number()>topN~'Others'))) %>% 
  dplyr::select(-Abundance) %>%
  head(n=topN+1)  

plotData_eso <- eso %>% # CHOOSE YOUR DATASET !
  transform_sample_counts(function(x) x/sum(x)) %>% psmelt %>% # "melt" the phyloseq object into a table
  inner_join(.,topTaxa_eso, by="Family", copy=TRUE) %>% #use topTaxa to join aggTaxo variable
  aggregate(Abundance~Treatment_type+aggTaxo_eso, #aggregate variables of interest
            data=., FUN=sum) #sum "Other" taxa

ggplot(plotData_eso, aes(x=Treatment_type, y=Abundance, fill=aggTaxo_eso)) +
  geom_bar(stat="identity",position="stack") + # change to "stack" to plot total counts
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~Treatment_type, scales = "free") + # RHS of the tilde is columns; add scales="free" to hide samples without counts
  scale_fill_manual(values = mycolors,breaks=c(pull(topTaxa_eso, Family),"Others")) +
  theme(panel.background = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 12))

topN = 10 
topTaxa_gil <- gil %>%  # CHOOSE YOUR DATASET !
  psmelt %>% group_by(Family) %>% 
  summarise(Abundance=sum(Abundance)) %>% # find most overall abundant taxa
  arrange(desc(Abundance)) %>% # order them 
  mutate(aggTaxo_gil=as.factor(case_when( # aggTaxo will become the plot legend
    row_number()<=topN~Family, row_number()>topN~'Others'))) %>% 
  dplyr::select(-Abundance) %>%
  head(n=topN+1) 

plotData_gil <- gil %>% # CHOOSE YOUR DATASET !
  transform_sample_counts(function(x) x/sum(x)) %>% psmelt %>% # "melt" the phyloseq object into a table
  inner_join(.,topTaxa_gil, by="Family", copy=TRUE) %>% #use topTaxa to join aggTaxo variable
  aggregate(Abundance~Treatment_type+aggTaxo_gil, #aggregate variables of interest
            data=., FUN=sum) #sum "Other" taxa

ggplot(plotData_gil, aes(x=Treatment_type, y=Abundance, fill=aggTaxo_gil)) +
  geom_bar(stat="identity",position="stack") + # change to "stack" to plot total counts
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~Treatment_type, scales = "free") + # RHS of the tilde is columns; add scales="free" to hide samples without counts
  scale_fill_manual(values = mycolors,breaks=c(pull(topTaxa_gil, Family),"Others")) +
  theme(panel.background = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 12))


#TMM CODE AGAIN
#make a stacked barplot of the data at the phylum level
phy_phylum = tax_glom(phy, taxrank = "Phylum",NArm=FALSE)
rel_abund_phylum = phyloseq::transform_sample_counts(phy_phylum,
                                                     function(x){x / sum(x)})

phyloseq::plot_bar(rel_abund_phylum, fill = "Phylum")+
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Sample_location, scales = "free") +
  theme(panel.background = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 12))


#find the core members of the microbiome
#convert to relative abundance
pseq_rel = microbiome::transform(phy, "compositional")

#make a variable for all of the conditions you want to compare
stuff = unique(as.character(meta(phy)$Treatment_period))
#[1] "No. of core taxa in Treatment2 : 153"
#[1] "No. of core taxa in Treatment : 146"
#[1] "No. of core taxa in Pre-treatment : 143"
#[1] "No. of core taxa in Post-treatment : 149"

#make a for loop to go through each bit of "stuff" one by one and combine identified core taxa into a list
list_core <- c() # an empty object to store information
for (n in stuff){ # for each variable n in stuff
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq_rel, stuff== n) # Choose sample from Treatment by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.10) #prevelence really matters--how much of the sample can it make upt 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  print(list_core)
}

#variable for Treatment v control 
stuff_2 = unique(as.character(meta(phy)$Treatment_type))
#[1] "No. of core taxa in Treatment : 136"
#[1] "No. of core taxa in Control : 128"

#loop for treatment v control
list_core2 <- c() # an empty object to store information
for (n in stuff_2){ # for each variable n in stuff
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq_rel, stuff_2== n) # Choose sample from Treatment by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.10) #prevelence really matters--how much of the sample can it make upt 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core2[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

plot(venn(list_core2))

#can also make a list of your factors and the color to give to each
#list of all your factors with the color to asign to each (replace nonCRC, CRC, H and add whatever other factors needed)
mycols = c(nonCRC="#d6e2e9", CRC="#cbf3f0", H="#fcf5c7") 
plot(venn(list_core2),
     fills = mycols)

##The network diagram is a representation of the relationship between microbiome samples in an experiment
##Creates a network display of the “connectedness” of samples according to some user-provided ecological similarity.
ig = make_network(phy, "samples", max.dist=0.1, dist.fun = "euclidean")
plot_network(ig, phy, color="Sample_location", shape="Treatment_type", line_weight=0.3, label=NULL)

#make a pcoa to evaluate each factor
phy_clr = microbiome::transform(phy, 'clr')
phy_fec = microbiome::transform(fec, 'clr')
phy_eso = microbiome::transform(eso, 'clr')
phy_gil = microbiome::transform(gil, 'clr')
phy_ord = ordinate(phy_clr, "RDA", "euclidean")
plot_ordination(phy,phy_ord, type = "samples", color="Sample_location", shape="Treatment_type")

#Test whether the sample locations differ significantly from each other using the permutational ANOVA 
dist.uf <- phyloseq::distance(phy_clr, method = "euclidean")
pairwise.adonis(t(otu_table(phy_clr)), sample_data(phy_clr)$Sample_location, sim.method = "euclidean",
                p.adjust.m = "bonferroni")

#                           pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
#1                Fecal vs Cecum  1 20155.537 15.181623 0.04756421   0.001      0.036   .
#2      Fecal vs Digestive_gland  1 10424.905  7.761953 0.02598039   0.001      0.036   .
#3            Fecal vs Esophagus  1 15139.002 11.940913 0.03791477   0.001      0.036   .
#4                 Fecal vs Gill  1 17617.972 13.697173 0.04325006   0.001      0.036   .
#5            Fecal vs Intestine  1 16571.661 12.396291 0.03905619   0.001      0.036   .
#6               Fecal vs Kidney  1 19175.062 14.550477 0.04567715   0.001      0.036   .
#7                 Fecal vs Skin  1 22361.009 16.780934 0.05215018   0.001      0.036   .
#8              Fecal vs Stomach  1 17477.362 13.081728 0.04112694   0.001      0.036   .
#9      Cecum vs Digestive_gland  1  1441.816  1.334882 0.04128302   0.001      0.036   .
#10           Cecum vs Esophagus  1  1992.586  3.197167 0.06920698   0.001      0.036   .
#11                Cecum vs Gill  1  2596.081  3.447422 0.07422204   0.001      0.036   .
#12           Cecum vs Intestine  1  1501.500  1.341180 0.02894143   0.001      0.036   .
#13              Cecum vs Kidney  1  1053.130  1.070941 0.02376123   0.132      1.000    
#14                Cecum vs Skin  1  1241.916  1.138959 0.02468540   0.019      0.684    
#15             Cecum vs Stomach  1  1269.117  1.139202 0.02469054   0.020      0.720    
#16 Digestive_gland vs Esophagus  1  1915.504  4.094135 0.12008327   0.001      0.036   .
#17      Digestive_gland vs Gill  1  2007.893  3.070486 0.09284672   0.001      0.036   .
#18 Digestive_gland vs Intestine  1  1922.312  1.635315 0.04861900   0.001      0.036   .
#19    Digestive_gland vs Kidney  1  1292.846  1.313820 0.04065816   0.003      0.108    
#20      Digestive_gland vs Skin  1  1623.538  1.431034 0.04280557   0.001      0.036   .
#21   Digestive_gland vs Stomach  1  1562.510  1.338028 0.04013518   0.008      0.288    
#22            Esophagus vs Gill  1  1813.149  5.887861 0.12295101   0.001      0.036   .
#23       Esophagus vs Intestine  1  1927.640  2.742042 0.05866329   0.001      0.036   .
#24          Esophagus vs Kidney  1  1358.014  2.451405 0.05393464   0.001      0.036   .
#25            Esophagus vs Skin  1  2660.469  3.952023 0.08241620   0.001      0.036   .
#26         Esophagus vs Stomach  1  1498.157  2.148283 0.04655174   0.001      0.036   .
#27            Gill vs Intestine  1  2873.004  3.462039 0.07294332   0.001      0.036   .
#28               Gill vs Kidney  1  1953.670  2.857126 0.06230495   0.001      0.036   .
#29                 Gill vs Skin  1  2759.090  3.448620 0.07268114   0.001      0.036   .
#30              Gill vs Stomach  1  2404.192  2.916864 0.06217091   0.001      0.036   .
#31          Intestine vs Kidney  1  1606.970  1.525576 0.03279004   0.001      0.036   .
#32            Intestine vs Skin  1  1888.226  1.632663 0.03427612   0.001      0.036   .
#33         Intestine vs Stomach  1  1227.233  1.040326 0.02211561   0.227      1.000    
#34               Kidney vs Skin  1  1340.389  1.308702 0.02826040   0.001      0.036   .
#35            Kidney vs Stomach  1  1181.669  1.127700 0.02444734   0.026      0.936    
#36              Skin vs Stomach  1  1614.234  1.402272 0.02958239   0.001      0.036   .

#test between treatment periods 
pairwise.adonis(t(otu_table(phy_clr)), sample_data(phy_clr)$Treatment_period, sim.method = "euclidean",
                p.adjust.m = "bonferroni")
                           
                          # pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
#1         Treatment2 vs Treatment  1  7947.733  6.139961 0.03212285   0.001  0.006   *
#2     Treatment2 vs Pre-treatment  1  7910.244  5.399287 0.02820955   0.001  0.006   *
#3    Treatment2 vs Post-treatment  1 54851.276 47.094059 0.15236158   0.001  0.006   *
#4      Treatment vs Pre-treatment  1  5906.494  5.074303 0.02614619   0.001  0.006   *
#5     Treatment vs Post-treatment  1 53794.396 56.423467 0.17554246   0.001  0.006   *
#6 Pre-treatment vs Post-treatment  1 46785.817 43.562350 0.14072238   0.001  0.006   *

#test between treatment groups 
pairwise.adonis(t(otu_table(phy_clr)), sample_data(phy_clr)$Treatment_type, sim.method = "euclidean",
                p.adjust.m = "bonferroni")


#                 pairs Df SumsOfSqs  F.Model          R2 p.value p.adjusted sig
#1 Treatment vs Control  1  5033.592 3.677276 0.008052241   0.001      0.001  **


#pairwise for fecal, esophagus and gill samples 
pairwise.adonis(t(otu_table(phy_fec)), sample_data(phy_fec)$Treatment_type, sim.method = "euclidean",
                p.adjust.m = "bonferroni")
#                 pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
#1 Treatment vs Control  1  6441.355 4.839187 0.01692975   0.001      0.001  **

pairwise.adonis(t(otu_table(phy_eso)), sample_data(phy_eso)$Treatment_type, sim.method = "euclidean",
                p.adjust.m = "bonferroni")
#                 pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
#1 Control vs Treatment  1   182.166 1.120701 0.05306175   0.304      0.304    


pairwise.adonis(t(otu_table(phy_gil)), sample_data(phy_gil)$Treatment_type, sim.method = "euclidean",
                p.adjust.m = "bonferroni")

#                 pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
#1 Control vs Treatment  1  351.9603 0.8960015 0.04287909   0.831      0.831   
