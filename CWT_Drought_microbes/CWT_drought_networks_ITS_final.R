#####Create environment and load necessary packages#####

library(dplyr)
library(phyloseq)
library(tidyr)
library(igraph)
library(seqtime)
library(SpiecEasi)
library(biomformat)
library(phyloseq)
library(ggplot2)
library(vegan)
library(car)
library(lme4)
library(funrar)
library(reshape)
library(emmeans)
library(cowplot)
library(RANN)
#library(MASS)

setwd('C:/Users/ernie/Desktop/Coweeta_drought')
drought_asvs <- "C:/Users/ernie/Desktop/Coweeta_drought/ITS-table-with-taxonomy.biom"
s1 <- read_biom(drought_asvs)
drought_meta <- "C:/Users/ernie/Desktop/QIIME/16S_metadata_CWT.txt"
drought_meta2 <- import_qiime_sample_data(drought_meta)
drought_otus2 <- import_biom(s1)
drought_otus2 <- merge_phyloseq(drought_otus2, drought_meta2)
drought_otus2 <- subset_samples(drought_otus2, id != "154")
drought_16s.2 <- rarefy_even_depth(drought_otus2, rngseed=TRUE)

#Change Phyla names
drought_16s.3 <- psmelt(drought_16s.2)

drought_16s.3$Class1=drought_16s.3$Rank3

unique(drought_16s.3$Class1)

drought_16s.3[drought_16s.3=="c__Agaricomycetes" ] <- "Agaricomycetes"
drought_16s.3[drought_16s.3=="c__Leotiomycetes"] <-"Leotiomycetes"
drought_16s.3[drought_16s.3=="c__Sordariomycetes"] <-"Sordariomycetes"
drought_16s.3[drought_16s.3=="c__Mortierellomycetes" ] <-"Mortierellomycetes"
drought_16s.3[drought_16s.3=="c__Geminibasidiomycetes" ] <-"Geminibasidiomycetes"
drought_16s.3[drought_16s.3=="c__Tremellomycetes"   ] <-"Tremellomycetes"
drought_16s.3[drought_16s.3=="c__Eurotiomycetes" ] <-"Eurotiomycetes"
drought_16s.3[drought_16s.3=="c__Mucoromycotina_cls_Incertae_sedis"] <-"Mucoromycotina"
drought_16s.3[drought_16s.3=="c__Umbelopsidomycetes"] <-"Umbelopsidomycetes"
drought_16s.3[drought_16s.3=="c__Pezizomycetes"  ] <-"Other"
drought_16s.3[drought_16s.3=="c__Dothideomycetes" ] <-"Other"
drought_16s.3[drought_16s.3=="c__unidentified"  ] <-"Other"
drought_16s.3[drought_16s.3=="c__Xylonomycetes" ] <-"Other"
drought_16s.3[drought_16s.3=="c__Microbotryomycetes"] <-"Other"
drought_16s.3[drought_16s.3=="c__Orbiliomycetes" ] <-"Other"
drought_16s.3[drought_16s.3=="c__Geoglossomycetes" ] <-"Other"
drought_16s.3[drought_16s.3=="c__Saccharomycetes"] <-"Other"
drought_16s.3[drought_16s.3=="c__Laboulbeniomycetes" ] <-"Other"
drought_16s.3[drought_16s.3=="c__Lecanoromycetes" ] <-"Other"
drought_16s.3[drought_16s.3=="c__Zoopagomycetes"  ] <-"Other"
drought_16s.3[drought_16s.3=="c__Calcarisporiellomycetes" ] <-"Other"
drought_16s.3[drought_16s.3=="c__Glomeromycetes" ] <-"Other"
drought_16s.3[drought_16s.3=="c__GS35"] <-"Other"
drought_16s.3[drought_16s.3=="c__Kickxellomycetes"] <-"Other"
drought_16s.3[drought_16s.3=="c__Tritirachiomycetes" ] <-"Other"
drought_16s.3[drought_16s.3=="c__GS19"] <-"Other"
drought_16s.3[drought_16s.3=="c__Basidiobolomycetes" ] <-"Other"
drought_16s.3[drought_16s.3=="c__Endogonomycetes"] <-"Other"
drought_16s.3[drought_16s.3=="c__Mucoromycetes" ] <-"Other"
drought_16s.3[drought_16s.3=="c__Chytridiomycetes"] <-"Other"
drought_16s.3[drought_16s.3=="c__GS18"] <-"Other"
drought_16s.3[drought_16s.3=="c__Taphrinomycetes" ] <-"Other"
drought_16s.3[drought_16s.3=="c__Rhizophydiomycetes"  ] <-"Other"
drought_16s.3[drought_16s.3=="c__Schizosaccharomycetes" ] <-"Other"
drought_16s.3$Class1[is.na(drought_16s.3$Class1)] <- "Other"

unique(drought_16s.3$Class1)

#Assign colors for Phyla
drought_16s.3$ClassColor=drought_16s.3$Class1

library("RColorBrewer")

brewer.pal(9, "Set1")

drought_16s.3$ClassColor[drought_16s.3$ClassColor=="Other"] <-"black"
drought_16s.3$ClassColor[drought_16s.3$ClassColor=="Agaricomycetes"] <-"#E41A1C"
drought_16s.3$ClassColor[drought_16s.3$ClassColor=="Geminibasidiomycetes"] <-"#4DAF4A"
drought_16s.3$ClassColor[drought_16s.3$ClassColor=="Leotiomycetes"] <-"#984EA3"
drought_16s.3$ClassColor[drought_16s.3$ClassColor=="Sordariomycetes" ] <-"#A65628"
drought_16s.3$ClassColor[drought_16s.3$ClassColor=="Tremellomycetes"] <-"#F781BF"
drought_16s.3$ClassColor[drought_16s.3$ClassColor=="Mortierellomycetes"  ] <-"#FF7F00" 
drought_16s.3$ClassColor[drought_16s.3$ClassColor=="Eurotiomycetes" ] <-"#377EB8"
drought_16s.3$ClassColor[drought_16s.3$ClassColor=="Mucoromycotina"] <-"#FFFF33"
drought_16s.3$ClassColor[drought_16s.3$ClassColor=="Umbelopsidomycetes"  ] <-"#999999"
unique(drought_16s.3$ClassColor)


con_ref <- subset_samples(drought_16s.2, Day != "1")
con_ref <- subset_samples(con_ref, Day != "42")
con_ref <- subset_samples(con_ref, LandUse == "Reference")
con_ref <- subset_samples(con_ref, Drought == "Control")


otus5=otu_table(con_ref)
taxa5=tax_table(con_ref)
filterobj5=filterTaxonMatrix(otus5,minocc=10,keepSum = TRUE, return.filtered.indices = TRUE)
otus.f5=filterobj5$mat
taxa.f5=taxa5[setdiff(1:nrow(taxa5),filterobj5$filtered.indices),]
dummyTaxonomy5=c("k__dummy","p__","c__","o__","f__","g__","s__")
taxa.f5=rbind(taxa.f5,dummyTaxonomy5)
rownames(taxa.f5)[nrow(taxa.f5)]="0"
rownames(otus.f5)[nrow(otus.f5)]="0"
updatedotus5=otu_table(otus.f5, taxa_are_rows = TRUE)
updatedtaxa5=tax_table(taxa.f5)
e_con_ref2 =phyloseq(updatedotus5, updatedtaxa5)
a1 <- t(data.frame(otu_table(e_con_ref2)))

Corr <- cor(a1,  method="spearman")

library(corrplot)

res1 <- cor.mtest(a1, conf.level = .99)

res1 <- as.matrix(res1[["p"]])

res1.1 <- as.matrix(p.adjust(res1, method="fdr"))

library(reshape2)
cor1 <- setNames(melt(Corr), c('taxa1', 'taxa2', 'cor'))

p1 <- setNames(melt(res1.1), c('taxa1', 'taxa2', 'p'))

p1.1 <- p1[p1$p < .01,]

cor1.1 <- cor1[rownames(cor1) %in% rownames(p1.1),]

g1<-simplify(graph.edgelist(as.matrix(cor1.1[,c(1,2)]),directed=FALSE))

#Add Phyla information to networks
V(g1)$Class = as.character(drought_16s.3$Class1[match(V(g1)$name, drought_16s.3$OTU)])#Adds taxanomic information to the different verticies
V(g1)$Class[is.na(V(g1)$Phylum)] <- "Other"
V(g1)$classcolor = as.character(drought_16s.3$ClassColor[match(V(g1)$name, drought_16s.3$OTU)])
V(g1)$classcolor[is.na(V(g1)$classcolor)] <- "black"

length(E(g1))
length(V(g1))

jpeg(filename="nfung1.jpeg",res=4000,units = "in", height=5, width=5) 


par(mar = rep(0, 4), xaxs='i', yaxs='i')
plot.igraph(g1, vertex.label=NA, vertex.size=5, edge.color=ifelse(cor1.1$cor > 0, "blue","red"),
            edge.width=c(2), vertex.frame.color=NA,vertex.color=V(g1)$classcolor)


dev.off()


drought_ref <- subset_samples(drought_16s.2, Day != "1")
drought_ref <- subset_samples(drought_ref, Day != "42")
drought_ref <- subset_samples(drought_ref, LandUse == "Reference")
drought_ref <- subset_samples(drought_ref, Drought == "Drought")

otus6=otu_table(drought_ref)
taxa6=tax_table(drought_ref)

filterobj6=filterTaxonMatrix(otus6,minocc=10,keepSum = TRUE, return.filtered.indices = TRUE)
otus.f6=filterobj6$mat
taxa.f6=taxa6[setdiff(1:nrow(taxa6),filterobj6$filtered.indices),]
dummyTaxonomy6=c("k__dummy","p__","c__","o__","f__","g__","s__")
taxa.f6=rbind(taxa.f6,dummyTaxonomy6)
rownames(taxa.f6)[nrow(taxa.f6)]="0"
rownames(otus.f6)[nrow(otus.f6)]="0"

updatedotus6=otu_table(otus.f6, taxa_are_rows = TRUE)
updatedtaxa6=tax_table(taxa.f6)
e_drought_ref2 =phyloseq(updatedotus6, updatedtaxa6)
a2 <- t(data.frame(otu_table(e_drought_ref2)))

Corr2 <- cor(a2,  method="spearman")

library(corrplot)

res2 <- cor.mtest(a2, conf.level = .99)

res2 <- as.matrix(res2[["p"]])

res2.1 <- as.matrix(p.adjust(res2, method="fdr"))

library(reshape2)
cor2 <- setNames(melt(Corr2), c('taxa1', 'taxa2', 'cor'))

p2 <- setNames(melt(res2.1), c('taxa1', 'taxa2', 'p'))

p2.1 <- p2[p2$p < .01,]

cor2.1 <- cor2[rownames(cor2) %in% rownames(p2.1),]

g2<-simplify(graph.edgelist(as.matrix(cor2.1[,c(1,2)]),directed=FALSE))

#Add Phyla information to networks
V(g2)$Class= as.character(drought_16s.3$Class1[match(V(g2)$name, drought_16s.3$OTU)])#Adds taxanomic information to the different verticies
V(g2)$Class[is.na(V(g2)$Class)] <- "Other"
V(g2)$classcolor = as.character(drought_16s.3$ClassColor[match(V(g2)$name, drought_16s.3$OTU)])
V(g2)$classcolor[is.na(V(g2)$classcolor)] <- "black"

length(E(g2))
length(V(g2))


jpeg(filename="nfung2.jpeg", bg="transparent", res=3000, units = "in", height=5, width=5) 

par(mar = rep(0, 4), xaxs='i', yaxs='i')
plot.igraph(g2, vertex.label=NA, vertex.size=5,edge.color=ifelse(cor2.1$cor > 0, "blue","red"),
            edge.width=c(2), vertex.frame.color=NA,vertex.color=V(g2)$classcolor)

dev.off()


con_dis <- subset_samples(drought_16s.2, Day != "1")
con_dis <- subset_samples(con_dis, Day != "42")
con_dis <- subset_samples(con_dis, LandUse == "Disturbed")
con_dis <- subset_samples(con_dis, Drought == "Control")

otus7=otu_table(con_dis)
taxa7=tax_table(con_dis)

filterobj7=filterTaxonMatrix(otus7,minocc=10,keepSum = TRUE, return.filtered.indices = TRUE)
otus.f7=filterobj7$mat
taxa.f7=taxa7[setdiff(1:nrow(taxa7),filterobj7$filtered.indices),]
dummyTaxonomy7=c("k__dummy","p__","c__","o__","f__","g__","s__")
taxa.f7=rbind(taxa.f7,dummyTaxonomy7)
rownames(taxa.f7)[nrow(taxa.f7)]="0"
rownames(otus.f7)[nrow(otus.f7)]="0"

updatedotus7=otu_table(otus.f7, taxa_are_rows = TRUE)
updatedtaxa7=tax_table(taxa.f7)
e_con_dis2 =phyloseq(updatedotus7, updatedtaxa7)

a3 <- t(data.frame(otu_table(e_con_dis2)))

Corr3 <- cor(a3,  method="spearman")

library(corrplot)

res3 <- cor.mtest(a3, conf.level = .99)

res3 <- as.matrix(res3[["p"]])

res3.1 <- as.matrix(p.adjust(res3, method="fdr"))

library(reshape2)
cor3 <- setNames(melt(Corr3), c('taxa1', 'taxa2', 'cor'))

p3 <- setNames(melt(res3.1), c('taxa1', 'taxa2', 'p'))

p3.1 <- p3[p3$p < .01,]

cor3.1 <- cor3[rownames(cor3) %in% rownames(p3.1),]

g3<-simplify(graph.edgelist(as.matrix(cor3.1[,c(1,2)]),directed=FALSE))

#Add Phyla information to networks
V(g3)$Class = as.character(drought_16s.3$Class1[match(V(g3)$name, drought_16s.3$OTU)])#Adds taxanomic information to the different verticies
V(g3)$Class[is.na(V(g3)$Class)] <- "Other"
V(g3)$classcolor = as.character(drought_16s.3$ClassColor[match(V(g3)$name, drought_16s.3$OTU)])
V(g3)$classcolor[is.na(V(g3)$classcolor)] <- "black"

length(E(g3))
length(V(g3))


jpeg(filename="nfung3.jpeg", bg="transparent", res=3000, units = "in", height=5, width=5) 

par(mar = rep(0, 4), xaxs='i', yaxs='i')
plot.igraph(g3, vertex.label=NA,vertex.size=5,edge.color=ifelse(cor3.1$cor > 0, "blue","red"),
            edge.width=c(2),  vertex.frame.color=NA,vertex.color=V(g3)$classcolor)

dev.off()


drought_dis <- subset_samples(drought_16s.2, Day != "1")
drought_dis <- subset_samples(drought_dis, Day != "42")
drought_dis <- subset_samples(drought_dis, LandUse == "Disturbed")
drought_dis <- subset_samples(drought_dis, Drought == "Drought")

otus8=otu_table(drought_dis)
taxa8=tax_table(drought_dis)

filterobj8=filterTaxonMatrix(otus8,minocc=10,keepSum = TRUE, return.filtered.indices = TRUE)
otus.f8=filterobj8$mat
taxa.f8=taxa8[setdiff(1:nrow(taxa8),filterobj8$filtered.indices),]
dummyTaxonomy8=c("k__dummy","p__","c__","o__","f__","g__","s__")
taxa.f8=rbind(taxa.f8,dummyTaxonomy8)
rownames(taxa.f8)[nrow(taxa.f8)]="0"
rownames(otus.f8)[nrow(otus.f8)]="0"

updatedotus8=otu_table(otus.f8, taxa_are_rows = TRUE)
updatedtaxa8=tax_table(taxa.f8)
e_drought_dis2 =phyloseq(updatedotus8, updatedtaxa8)

a4 <- t(data.frame(otu_table(e_drought_dis2)))

Corr4 <- cor(a4,  method="spearman")

library(corrplot)

res4 <- cor.mtest(a4, conf.level = .99)

res4 <- as.matrix(res4[["p"]])

res4.1 <- as.matrix(p.adjust(res4, method="fdr"))

library(reshape2)
cor4 <- setNames(melt(Corr4), c('taxa1', 'taxa2', 'cor'))

p4 <- setNames(melt(res4.1), c('taxa1', 'taxa2', 'p'))

p4.1 <- p4[p4$p < .01,]

cor4.1 <- cor4[rownames(cor4) %in% rownames(p4.1),]

g4<-simplify(graph.edgelist(as.matrix(cor4.1[,c(1,2)]),directed=FALSE))

#Add Phyla information to networks
V(g4)$Class = as.character(drought_16s.3$Class1[match(V(g4)$name, drought_16s.3$OTU)])#Adds taxanomic information to the different verticies
V(g4)$Class[is.na(V(g4)$Class)] <- "Other"
V(g4)$classcolor = as.character(drought_16s.3$ClassColor[match(V(g4)$name, drought_16s.3$OTU)])
V(g4)$classcolor[is.na(V(g4)$classcolor)] <- "black"

length(E(g4))
length(V(g4))


jpeg(filename="nfung4.jpeg", bg="transparent", res=3000, units = "in", height=5, width=5) 

par(mar = rep(0, 4), xaxs='i', yaxs='i')
plot.igraph(g4, vertex.label=NA, vertex.size=5, edge.color=ifelse(cor4.1$cor > 0, "blue","red"),
            edge.width=c(2), vertex.frame.color=NA,vertex.color=V(g4)$classcolor)

dev.off()

jpeg(filename="fungleg.jpeg",res=3000,units = "in", height=5, width=3.25) 
par(mar = rep(0, 4), xaxs='i', yaxs='i')
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend('center',legend=unique(V(g4)$Class),cex=1.5,pt.cex=2.5,pch=21,pt.bg=unique(V(g4)$classcolor))
dev.off()


spiec.graph5 <- g1
spiec.graph6 <- g2
spiec.graph7 <- g3
spiec.graph8 <- g4

#Normalized Degree

degree5=data.frame(degree(spiec.graph5,normalized=TRUE))
colnames(degree5)
names(degree5)[names(degree5) == "degree.spiec.graph5..normalized...TRUE."] <- "clustering"
degree5 <- data.frame(degree = degree5$clustering, Treatment = rep("Reference-Control",length(degree5$clustering)),Drought = rep("Control",length(degree5$clustering)), LandUse=rep("Reference",length(degree5$clustering)),time=rep("Early",length(degree5$clustering)))

degree6=data.frame(degree(spiec.graph6,normalized=TRUE))
colnames(degree6)
names(degree6)[names(degree6) == "degree.spiec.graph6..normalized...TRUE."] <- "clustering"
degree6 <- data.frame(degree = degree6$clustering,Treatment = rep("Reference-Drought",length(degree6$clustering)), Drought = rep("Drought",length(degree6$clustering)), LandUse=rep("Reference",length(degree6$clustering)),time=rep("Early",length(degree6$clustering)))

degree7=data.frame(degree(spiec.graph7,normalized=TRUE))
colnames(degree7)
names(degree7)[names(degree7) == "degree.spiec.graph7..normalized...TRUE."] <- "clustering"
degree7 <- data.frame(degree = degree7$clustering,Treatment = rep("Disturbed-Control",length(degree7$clustering)), Drought = rep("Control",length(degree7$clustering)), LandUse=rep("Disturbed",length(degree7$clustering)),time=rep("Early",length(degree7$clustering)))

degree8=data.frame(degree(spiec.graph8,normalized=TRUE))
colnames(degree8)
names(degree8)[names(degree8) == "degree.spiec.graph8..normalized...TRUE."] <- "clustering"
degree8 <- data.frame(degree = degree8$clustering,Treatment = rep("Disturbed-Drought",length(degree8$clustering)),  Drought = rep("Drought",length(degree8$clustering)), LandUse=rep("Disturbed",length(degree8$clustering)),time=rep("Early",length(degree8$clustering)))


degree <- rbind(degree5,degree6,degree7,degree8)

degree$LandUse <- as.factor(degree$LandUse)

degree$Drought <- as.factor(degree$Drought)


hist(degree$degree)
glm2 <- glm(degree+.001~LandUse*Drought, data=degree,family=Gamma(link=log))
Anova(glm2)

deg <- aggregate(degree~Treatment+time, data=degree, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
deg

deg <- do.call(data.frame, deg)
deg

deg$se <- deg$degree.sd / sqrt(deg$degree.n)
head(deg)

deg$Treatment <- factor(deg$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
deg$Treatment


jpeg(filename="fung_deg.jpeg", bg="transparent", res=600, units = "in", height=4.5, width=6) 


library(ggplot2)
bar_plot_degree<- ggplot(deg, aes(x=time, y=degree.mean, fill=Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=degree.mean-se, ymax=degree.mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate('text', .7, .07, label="Drought: ~italic(P) < 0.001", size=4, parse=TRUE)+
  annotate('text', .7, .065, label="Disturbance: italic(P) < 0.001", size = 4, parse=TRUE)+
  annotate('text', .7
           , .06, label="Disturbance %*% Drought: italic(P) == 0.512", size=4, parse=TRUE)+
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(legend.position="none")+
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  ylab('Normalized Degree') 
#coord_cartesian(ylim=c(.0,.0065)) 
plot(bar_plot_degree)

dev.off()

#Betweenness

between5=data.frame(betweenness(spiec.graph5,normalized=TRUE))
colnames(between5)
names(between5)[names(between5) == "betweenness.spiec.graph5..normalized...TRUE."] <- "clustering"
between5 <- data.frame(between = between5$clustering, Treatment = rep("Reference-Control",length(between5$clustering)),Drought = rep("Control",length(between5$clustering)), LandUse=rep("Reference",length(between5$clustering)),time=rep("Early",length(between5$clustering)))


between6=data.frame(betweenness(spiec.graph6,normalized=TRUE))
colnames(between6)
names(between6)[names(between6) == "betweenness.spiec.graph6..normalized...TRUE."] <- "clustering"
between6 <- data.frame(between = between6$clustering,Treatment = rep("Reference-Drought",length(between6$clustering)), Drought = rep("Drought",length(between6$clustering)), LandUse=rep("Reference",length(between6$clustering)),time=rep("Early",length(between6$clustering)))

between7=data.frame(betweenness(spiec.graph7,normalized=TRUE))
colnames(between7)
names(between7)[names(between7) == "betweenness.spiec.graph7..normalized...TRUE."] <- "clustering"
between7 <- data.frame(between = between7$clustering,Treatment = rep("Disturbed-Control",length(between7$clustering)), Drought = rep("Control",length(between7$clustering)), LandUse=rep("Disturbed",length(between7$clustering)),time=rep("Early",length(between7$clustering)))

between8=data.frame(betweenness(spiec.graph8,normalized=TRUE))
colnames(between8)
names(between8)[names(between8) == "betweenness.spiec.graph8..normalized...TRUE."] <- "clustering"
between8 <- data.frame(between = between8$clustering,Treatment = rep("Disturbed-Drought",length(between8$clustering)),  Drought = rep("Drought",length(between8$clustering)), LandUse=rep("Disturbed",length(between8$clustering)),time=rep("Early",length(between8$clustering)))

between <- rbind(between5,between6,between7,between8)

between$LandUse <- as.factor(between$LandUse)

between$Drought <- as.factor(between$Drought)

hist(between$between)
glm2 <- glm(between+.001~LandUse*Drought, data=between,family=Gamma(link=log))
Anova(glm2)


bet <- aggregate(between~Treatment+time, data=between, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
bet

bet <- do.call(data.frame, bet)
bet

bet$se <- bet$between.sd / sqrt(bet$between.n)
head(bet)

bet$Treatment <- factor(bet$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
bet$Treatment



jpeg(filename="fung_bet.jpeg", bg="transparent", res=600, units = "in", height=4.5, width=6) 


library(ggplot2)
bar_plot_between<- ggplot(bet, aes(x=time, y=between.mean, fill=Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=between.mean-se, ymax=between.mean+se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate('text', 1.1, .0165, label="Drought: ~italic(P) < 0.001", size=4, parse=TRUE)+
  annotate('text', 1.1, .0155, label="Disturbance: italic(P) == 0.147", size = 4, parse=TRUE)+
  annotate('text', 1.1, .0145, label="Disturbance %*% Drought: italic(P) == 0.488", size=4, parse=TRUE)+
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position = c(.6,.68)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  ylab('Normalized Betweenness') 
plot(bar_plot_between)

dev.off()

#z test for negative edges
n_e_cr=sum(cor1.1$cor < 0) #total number of negative edges for the disturbed
t_e_cr=nrow(cor1.1)#total edges for the disturbed
n_e_dr=sum(cor2.1$cor < 0) #total number of positive edges for the disturbed
t_e_dr=nrow(cor2.1)#total edges for the reference

prop.test(c(n_e_cr,n_e_dr),c(t_e_cr,t_e_dr))

n_e_cd=sum(cor3.1$cor < 0) #total number of negative edges for the disturbed
t_e_cd=nrow(cor3.1)#total edges for the disturbed
n_e_dd=sum(cor4.1$cor < 0) #total number of positive edges for the disturbed
t_e_dd=nrow(cor4.1)#total edges for the reference

prop.test(c(n_e_cd,n_e_dd),c(t_e_cd,t_e_dd))


