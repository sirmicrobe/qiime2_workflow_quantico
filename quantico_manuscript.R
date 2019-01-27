# R analysis following QIIME2 processing for Quantico Manuscript

library("tidyr")
library("vegan");packageVersion("vegan") #2.5.1
library("phyloseq"); packageVersion("phyloseq") #1.22.3
library("ggplot2"); packageVersion("ggplot2") #2.2.1
library(dplyr) #data frame manipulating
library("FSA")

theme_set(theme_bw())

writeLines(capture.output(sessionInfo()), "RsessionInfo.txt")

#import data from QIIME
otu <- otu_table(features,taxa_are_rows=T)
tax <- tax_table(as.matrix(tax_filtered))
stax <- tax_table(as.matrix(silva_filter))
mapfile = read.csv("~/quantico_map_allR.csv",row.names=1)
mymap <- sample_data(mapfile)
root_tree <- read_tree("root_tree.nwk")
unroot_tree <- read_tree("unroot_tree.nwk")

#check that OTU names are consistent
taxa_names(tax) == sort(taxa_names(otu))
taxa_names(stax) == sort(taxa_names(otu))
#make sure files have same sample names
sort(sample_names(otu)) == sort(sample_names(mymap))


myDatas <- phyloseq(otu,stax,mymap,root_tree)
myDatas_orig <- myData

#use only samples that have >3K reads
pcbs <- prune_samples(sample_sums(myDatas)>=3000, myDatas)
pcbs <- prune_taxa(taxa_sums(pcbs) > 0, pcbs)
sort(sample_sums(pcbs))

#Rarefy
set.seed(04242018)
pcbRs = rarefy_even_depth(pcbs, sample.size=5750) # 9140 taxa and 118 samples
sum(taxa_sums(pcbRs) == 0) #0
pcbRs_rab <- transform_sample_counts(pcbRs, function(x) x / sum(x) )
pcbRs_pdd <- merge_samples(pcbRs, "plot_day_depth")
pcbRs_pdd_rab <- transform_sample_counts(pcbRs_pdd, function(x) x / sum(x) )
sample_data(pcbRs_pdd_rab)$samples <- rownames(sample_data(pcbRs_pdd_rab))

#data frame of samples
d_pcb <- as(sample_data(pcb), "data.frame")
head(d_pcb)
d_pcbR <- as(sample_data(pcbR), "data.frame")



########### Alpha Diversity ###########

#reorder levels to plot clearer
sample_data(pcb)$plot_day_depth <- factor(sample_data(pcb)$plot_day_depth, levels= c("1_0_top", "1_0_bottom","1_140_top","1_140_bottom","1_409_top","1_409_bottom","2_0_top", "2_0_bottom","2_140_top","2_140_bottom","2_409_top","2_409_bottom","3_0_top", "3_0_bottom","3_140_top","3_140_bottom","3_409_top","3_409_bottom","4_0_top", "4_0_bottom","4_140_top","4_140_bottom","4_409_top","4_409_bottom","5_140_top","5_409_top"))
#plot shannon diversity
p_rich_shan_collapse <- plot_richness(subset_samples(pcbs,plot_num %in% c(1,2,3,4)), x = "plot_day_depth", color = "Depth", measures="Shannon") + geom_boxplot()
p_rich_shan_collapse +theme(text = element_text(size=20)) +facet_grid(~plot_num,scale="free") + ylab("Shannon Diversity Index")+ xlab("Sample") 

###quantico stats
d_q_rich <- estimate_richness(pcbs, measures = "Shannon")
d_pcb #samples metadata data frame
d_q_rich <- merge(d_pcb,d_q_rich, by="row.names")

#all by all
kruskal.test(d_q_rich$Shannon~as.factor(plot_day_depth),data=d_q_rich) #Kruskal-Wallis chi-squared = 70.688, df = 25, p-value = 3.039e-06
dunn_all <- dunnTest(Shannon ~ as.factor(plot_day_depth),
         data=d_q_rich,
         method="bh") #multiple comparison test 
print(dunn_all,dunn.test.results=TRUE)
write.csv(dunn_all$res,file="kw_dunn_plot_day_depth_Shannon.csv")

#by site
kruskal.test(d_q_rich$Shannon~as.factor(plot_num),data=d_q_rich) #Kruskal-Wallis chi-squared = 20.207, df = 4, p-value = 0.0004545
dunnTest(Shannon ~ as.factor(plot_num),
         data=d_q_rich,
         method="bh") #multiple comparison test 
# 1 - 4 -2.9724007 2.954808e-03 0.0098493587
# 2 - 4 -3.9610677 7.461537e-05 0.0007461537
# 3 - 4 -3.2674245 1.085308e-03 0.0054265404


#subset within site
d_q_rich_p1 <- subset(d_q_rich, plot_num == "1")
nrow(d_q_rich_p1)
kruskal.test(d_q_rich_p1$Shannon~as.factor(plot_day_depth),data=d_q_rich_p1) #Kruskal-Wallis chi-squared = 10.702, df = 5, p-value = 0.05763
dunnTest(Shannon ~ as.factor(plot_day_depth),
         data=d_q_rich_p1,
         method="bh") #multiple comparison test
#2
d_q_rich_p2 <- subset(d_q_rich, plot_num == "2")
nrow(d_q_rich_p2)
kruskal.test(d_q_rich_p2$Shannon~as.factor(plot_day_depth),data=d_q_rich_p2) #Kruskal-Wallis chi-squared = 19.395, df = 5, p-value = 0.001622
dunnTest(Shannon ~ as.factor(plot_day_depth),
         data=d_q_rich_p2,
         method="bh") #multiple comparison test
# 2_0_top - 2_409_bottom -3.0892112 0.0020068872 0.01505165
# 2_140_bottom - 2_409_bottom -2.8018427 0.0050811640 0.02540582
# 2_140_top - 2_409_bottom -3.3406586 0.0008357992 0.01253699
# 2_0_top - 2_409_top -2.5144742 0.0119210015 0.03576300
# 2_140_top - 2_409_top -2.7659217 0.0056762171 0.02128581
#3
d_q_rich_p3 <- subset(d_q_rich, plot_num == "3")
nrow(d_q_rich_p3)
kruskal.test(d_q_rich_p3$Shannon~as.factor(plot_day_depth),data=d_q_rich_p3) #Kruskal-Wallis chi-squared = 15.978, df = 5, p-value = 0.006907
dunnTest(Shannon ~ as.factor(plot_day_depth),
         data=d_q_rich_p3,
         method="bh") #multiple comparison test
# 3_140_top - 3_409_top -2.5863163 0.009700785 0.04850393
# 3_0_bottom - 3_409_top -2.9814480 0.002868888 0.04303331
#4
d_q_rich_p4 <- subset(d_q_rich, plot_num == "4")
nrow(d_q_rich_p4)
kruskal.test(d_q_rich_p4$Shannon~as.factor(plot_day_depth),data=d_q_rich_p4) #Kruskal-Wallis chi-squared = 4.3445, df = 5, p-value = 0.501
dunnTest(Shannon ~ as.factor(plot_day_depth),
         data=d_q_rich_p4,
         method="bh") #multiple comparison test

#subset by day 0 and top
d_q_rich_d0 <- subset(d_q_rich, Day == "0")
d_q_rich_d0t <- subset(d_q_rich_d0, Depth == "top")
nrow(d_q_rich_d0)
kruskal.test(d_q_rich_d0t$Shannon~as.factor(plot_day_depth),data=d_q_rich_d0t) #Kruskal-Wallis chi-squared = 9.3111, df = 3, p-value = 0.02543
dunnTest(Shannon ~ as.factor(plot_day_depth),
         data=d_q_rich_d0t,
         method="bh") #multiple comparison test
# day 0 bottom
d_q_rich_d0b <- subset(d_q_rich_d0, Depth == "bottom")
kruskal.test(d_q_rich_d0b$Shannon~as.factor(plot_day_depth),data=d_q_rich_d0b) #Kruskal-Wallis chi-squared = 9.3111, df = 3, p-value = 0.02543
dunnTest(Shannon ~ as.factor(plot_day_depth),
         data=d_q_rich_d0b,
         method="bh") #multiple comparison test
# day 409 bottom
d_q_rich_d409 <- subset(d_q_rich, Day == "409")
d_q_rich_d409t <- subset(d_q_rich_d409, Depth == "top")
kruskal.test(d_q_rich_d409t$Shannon~as.factor(plot_day_depth),data=d_q_rich_d409t) 
dunnTest(Shannon ~ as.factor(plot_day_depth),
         data=d_q_rich_d409t,
         method="bh") #multiple comparison test
         

########### NMDS ##############

pcbR_bottom <- subset_samples(pcbRs, Depth == "bottom")
pcbR_bottom = prune_taxa(taxa_sums(pcbR_bottom) > 0, pcbR_bottom) #[ 873 taxa and 134 samples ]
pcbR_bottom <- subset_samples(pcbR_bottom, plot_num != "5")

pcbR_top<- subset_samples(pcbRs, Depth == "top")
pcbR_top = prune_taxa(taxa_sums(pcbR_top) > 0, pcbR_top) #[ 873 taxa and 134 samples ]
pcbR_top <- subset_samples(pcbR_top, plot_num != "5")


pcbR.top.ord.pcoa.wunif = ordinate(pcbR_top, "PCoA", "wunifrac")
plot_ordination(pcbR_top, pcbR.top.ord.pcoa.wunif, color="plot_num", shape="Day")+ scale_color_manual(values=c("blue","darkorange","purple","darkgreen","grey"))+geom_point(size=9) + scale_shape_manual(values=c(20,17,8))

pcbR.bottom.ord.pcoa.wunif = ordinate(pcbR_bottom, "PCoA", "wunifrac")
plot_ordination(pcbR_bottom, pcbR.bottom.ord.pcoa.wunif, color="plot_num", shape="Day")+ scale_color_manual(values=c("blue","darkorange","purple","darkgreen","grey"))+geom_point(size=9) + scale_shape_manual(values=c(20,17,8))

d_top = distance(pcbR_top, "wunifrac")
df_top = as(sample_data(pcbR_top), "data.frame")
adonis(d_top~plot_day, data=df_top)
library("pairwiseAdonis")
pairwise.adonis(d_top,factors=df_top$plot_day,p.adjust.m='bonferroni',perm=999)
#manually reduce dataset for multiple hypothesis testing
df_top_adonis <- read.csv("pairwise_adonis.csv",header=T)
p.adjust(df_top_adonis$p.value,method="BH")
df_top_adonis_day <- read.csv("pairwise_adonis_day.csv",header=T)
p.adjust(df_top_adonis_day$p.value,method="BH")
cbind(df_top_adonis_day, p.adjust(df_top_adonis_day$p.value,method="BH"))

################ BARPLOTS ###################
pcbRs_dehal <- subset_taxa(pcbRs_pdd_rab, Family=="Dehalococcoidaceae")
p_pcbRs_dehal <- plot_bar(pcbRs_dehal, x="samples", y="Abundance", fill="Genus", title="Dehalococcoidaceae-only")
p_pcbRs_dehal + scale_fill_manual(values = c("grey26","red", "chartreuse3","darkorange","cyan2", "darkgreen","blue4", "yellow1","deepskyblue", "mediumorchid3","royalblue","orangered3","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861","brown2","darkcyan","orchid1","orange1","mediumpurple1","skyblue4","salmon1","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4")) + theme(text = element_text(size=14, face="bold")) + coord_flip() + facet_grid(plot_num~.,scales="free",space="free")

pcbRs_burk <- subset_taxa(pcbRs_pdd_rab, Family=="Burkholderiaceae")
p_pcbRs_burk <- plot_bar(pcbRs_burk, x="samples", y="Abundance", fill="Genus", title="Burkholderiaceae-only")
p_pcbRs_burk + scale_fill_manual(values = c("grey26","red", "chartreuse3","darkorange","cyan2", "darkgreen","blue4", "yellow1","deepskyblue", "mediumorchid3","royalblue","orangered3","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861","brown2","darkcyan","orchid1","orange1","mediumpurple1","skyblue4","salmon1","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4")) + theme(text = element_text(size=14, face="bold")) + coord_flip() + facet_grid(plot_num~.,scales="free",space="free") 
pcbRs_burks <- subset_taxa(pcbRs_pdd_rab, Species=="Bacteria;Proteobacteria;Gammaproteobacteria;Betaproteobacteriales;Burkholderiaceae;Burkholderia-Caballeronia-Paraburkholderia")
p_pcbRs_burks <- plot_bar(pcbRs_burks, x="samples", y="Abundance", fill="Species", title="Burkholderia-only")
p_pcbRs_burks + scale_fill_manual(values = c("grey26","red", "chartreuse3","darkorange","cyan2", "darkgreen","blue4", "yellow1","deepskyblue", "mediumorchid3","royalblue","orangered3","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861","brown2","darkcyan","orchid1","orange1","mediumpurple1","skyblue4","salmon1","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4")) + theme(text = element_text(size=14, face="bold")) + coord_flip() + facet_grid(plot_num~.,scales="free",space="free") # + geom_bar(stat="identity", position="stack")

top50otus = names(sort(taxa_sums(pcbRs_pdd_rab), TRUE)[1:50])
taxtab50 = cbind(tax_table(pcbRs_pdd_rab), family50 = NA)
taxtab50[top50otus, "family50"] <- as(tax_table(pcbRs_pdd_rab)[top50otus, "Family"], "character")
tax_table(pcbRs_pdd_rab) <- tax_table(taxtab50)
pcbRs_top50 = prune_taxa(top50otus, pcbRs_pdd_rab)
taxa_sums(pcbRs_top50)
sample_data(pcbRs_top50)$samples <- rownames(sample_data(pcbRs_top50))
title = "Top 50 OTUs-Family"
plot_bar(pcbRs_top50, "samples", fill = "family50", title = title) + coord_flip() + ylab("Percentage of Sequences")  + scale_fill_manual(values = c("grey","darkorange","chartreuse3","darkgreen","red","blue4","yellow1","deepskyblue","#89C5DA","#DA5724","#74D944","#C84248","#673770","cyan2","#D3D93E","royalblue","#38333E","#508578","#D7C1B1","#689030","#AD6F3B","#CD9BCD","#D14285","#6DDE88","#652926","#7FDCC0","#8569D5","#5E738F","#D1A33D","#8A7C64","#599861","brown2","darkcyan","orchid1","orange1","mediumpurple1","skyblue4","steelblue2","yellowgreen","lightslategrey","peachpuff","turquoise4","grey26","#CE50CA","#3F4921","#C0717C","#CBD588","#5F7FC7","#ABC6DF","#272617","#D6D676","#8C378E","#56133D","#C9734D","#6E7424")) + guides(fill=guide_legend(nrow=24,title.position="top",title="Family"))+ theme(text = element_text(size=14, face="bold"))#,legend.position="none")
# top to most abundant OTUs - GENUS
top50otus = names(sort(taxa_sums(pcbRs_pdd_rab), TRUE)[1:50])
taxtab50 = cbind(tax_table(pcbRs_pdd_rab), genus50 = NA)
taxtab50[top50otus, "genus50"] <- as(tax_table(pcbRs_pdd_rab)[top50otus, "Genus"], "character")
tax_table(pcbRs_pdd_rab) <- tax_table(taxtab50)
pcbRs_pdd_rab_top50 = prune_taxa(top50otus, pcbRs_pdd_rab)
taxa_sums(pcbRs_pdd_rab_top50)
title = "Top 50 OTUs-Genus"
sample_data(pcbRs_pdd_rab_top50)$samples <- rownames(sample_data(pcbRs_pdd_rab_top50 ))
plot_bar(pcbRs_pdd_rab_top50, "samples", fill = "genus50", title = title) + coord_flip() + ylab("Percentage of Sequences")  + scale_fill_manual(values = c("grey","darkorange","chartreuse3","darkgreen","red","blue4","yellow1","deepskyblue","violetred","#89C5DA","#DA5724","#74D944","#C84248","#673770","cyan2","#D3D93E","royalblue","#38333E","#508578","#D7C1B1","#689030","#AD6F3B","#CD9BCD","#D14285","#6DDE88","#652926","#7FDCC0","#8569D5","#5E738F","#D1A33D","#8A7C64","#599861","brown2","darkcyan","orchid1","orange1","mediumpurple1","skyblue4","steelblue2","yellowgreen","lightslategrey","peachpuff","turquoise4","grey26","#CE50CA","#3F4921","#C0717C","#CBD588","#5F7FC7","#ABC6DF","#272617","#D6D676","#8C378E","#56133D","#C9734D","#6E7424"))  + guides(fill=guide_legend(nrow=24,title.position="top",title="Genus"))+ theme(text = element_text(size=14, face="bold"),legend.position="none")
# all levels
top50otus = names(sort(taxa_sums(pcbRs_pdd_rab), TRUE)[1:50])
taxtab50 = cbind(tax_table(pcbRs_pdd_rab), taxon50 = NA)
taxtab50[top50otus, "taxon50"] <- as(tax_table(pcbRs_pdd_rab)[top50otus, "Taxon"], "character")
tax_table(pcbRs_pdd_rab) <- tax_table(taxtab50)
pcbRs_pdd_rab_top50 = prune_taxa(top50otus, pcbRs_pdd_rab)
taxa_sums(pcbRs_pdd_rab_top50)
title = "Top 50 OTUs"
sample_data(pcbRs_pdd_rab_top50)$samples <- rownames(sample_data(pcbRs_pdd_rab_top50 ))
plot_bar(pcbRs_pdd_rab_top50, "samples", fill = "taxon50", title = title) + coord_flip() + ylab("Percentage of Sequences")  + scale_fill_manual(values = c("grey","darkorange","chartreuse3","darkgreen","red","blue4","yellow1","deepskyblue","violetred","#89C5DA","#DA5724","#74D944","#C84248","#673770","cyan2","#D3D93E","royalblue","#38333E","#508578","#D7C1B1","#689030","#AD6F3B","#CD9BCD","#D14285","#6DDE88","#652926","#7FDCC0","#8569D5","#5E738F","#D1A33D","#8A7C64","#599861","brown2","darkcyan","orchid1","orange1","mediumpurple1","skyblue4","steelblue2","yellowgreen","lightslategrey","peachpuff","turquoise4","grey26","#CE50CA","#3F4921","#C0717C","#CBD588","#5F7FC7","#ABC6DF","#272617","#D6D676","#8C378E","#56133D","#C9734D","#6E7424"))  + guides(fill=guide_legend(nrow=24,title.position="top",title="Taxonomy"))+ theme(text = element_text(size=14, face="bold"),legend.position="none")
