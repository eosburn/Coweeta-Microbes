library(biomformat)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)
library(car)
library(lme4)
library(reshape)
library(emmeans)
library(cowplot)

soil <- read.csv("C:/Users/ernie/OneDrive/Desktop/Rhodo Sequences/RhodoEnzymesMaster.csv")

m1 <- lm(ph~ Treatment*Date, data=soil)
shapiro.test(resid(m1))
Anova(m1)
m1.1 <- aggregate(ph ~ Treatment, data=soil, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m1.1

m1.2 <- do.call(data.frame, m1.1)

m1.2$se <- m1.2$ph.sd / sqrt(m1.2$ph.n)
m1.2

m2 <- lm(NO3~ Treatment*Date, data=soil)
shapiro.test(resid(m2))
m2 <- glm(NO3~ Treatment*Date, data=soil,family=Gamma(link=log))
Anova(m2)
m2.1 <- aggregate(NO3 ~ Treatment, data=soil, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m2.1

m2.2 <- do.call(data.frame, m2.1)

m2.2$se <- m2.2$NO3.sd / sqrt(m2.2$NO3.n)
m2.2

m3 <- lm(NH4~ Treatment*Date, data=soil)
shapiro.test(resid(m3))
Anova(m3)
m3.1 <- aggregate(NH4 ~ Treatment, data=soil, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m3.1

m3.2 <- do.call(data.frame, m3.1)

m3.2$se <- m3.2$NH4.sd / sqrt(m3.2$NH4.n)
m3.2

m4 <- lm(DOC~ Treatment*Date, data=soil)
shapiro.test(resid(m4))
Anova(m4)
m4.1 <- aggregate(DOC ~ Treatment, data=soil, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m4.1

m4.2 <- do.call(data.frame, m4.1)

m4.2$se <- m4.2$DOC.sd / sqrt(m4.2$DOC.n)
m4.2

m5 <- lm(MBC~ Treatment*Date, data=soil)
shapiro.test(resid(m5))
m5 <- glm(MBC~ Treatment*Date, data=soil,family=Gamma(link=log))
Anova(m5)
m5.1 <- aggregate(MBC ~ Treatment, data=soil, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m5.1

m5.2 <- do.call(data.frame, m5.1)

m5.2$se <- m5.2$MBC.sd / sqrt(m5.2$MBC.n)
m5.2

m6 <- lm(TDN~ Treatment*Date, data=soil)
shapiro.test(resid(m6))
m6 <- glm(TDN~ Treatment*Date, data=soil,family=Gamma(link=log))
Anova(m6)
m6.1 <- aggregate(TDN ~ Treatment, data=soil, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m6.1

m6.2 <- do.call(data.frame, m6.1)

m6.2$se <- m6.2$TDN.sd / sqrt(m6.2$TDN.n)
m6.2


Bac_otus <- "C:/Users/ernie/OneDrive/Desktop/Rhodo sequences/rhodo-16S-table-with-taxonomy.biom"

Fung_otus <- "C:/Users/ernie/OneDrive/Desktop/Rhodo sequences/rhodo-ITS-table-with-taxonomy.biom"

x1 <- read_biom(Bac_otus)

x2 <- read_biom(Fung_otus)

OTU1 <- import_biom(x1)

OTU2 <- import_biom(x2)

meta <- "C:/Users/ernie/OneDrive/Desktop/Rhodo sequences/16S_metadata_rhodo.txt"

meta2 <- import_qiime_sample_data(meta)

OTU1 <- merge_phyloseq(OTU1, meta2)

OTU2 <- merge_phyloseq(OTU2, meta2)

bac <- rarefy_even_depth(OTU1, rngseed=TRUE)

bac1 <- subset_samples(bac, Year == "2014")

bac1 <- subset_samples(bac1, Month == "July")

bac2 <- subset_samples(bac, Year == "2017")

ITS <- rarefy_even_depth(OTU2, rngseed=TRUE)

ITS1 <- subset_samples(ITS, Year == "2014")

ITS1 <- subset_samples(ITS1, Month == "July")

ITS2 <- subset_samples(ITS, Year == "2017")

a1.1 <- t(data.frame(otu_table(bac1)))

a1.2 <- t(data.frame(otu_table(bac2)))

a2.1 <- t(data.frame(otu_table(ITS1)))

a2.2 <- t(data.frame(otu_table(ITS2)))

b1.1 <- vegdist(a1.1, method = "bray")

b1.2 <- vegdist(a1.2, method = "bray")

b2.1 <- vegdist(a2.1, method = "bray")

b2.2 <- vegdist(a2.2, method = "bray")

c1.1 <- data.frame(sample_data(bac1))

c1.2 <- data.frame(sample_data(bac2))

c2.1 <- data.frame(sample_data(ITS1))

c2.2 <- data.frame(sample_data(ITS2))

#Bacterial Beta Diversity

#Before Treatment

set.seed(101)
adonis2(b1.1~Treatment, data=c1.1, perm=999)

set.seed(101)
d1.1 <- metaMDS(b1.1, k=2, trymax=1000)

data.scores1 <- as.data.frame(d1.1$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores1$Treatment <- c1.1$Treatment 
head(data.scores1) 

data.scores1$Treatment <- factor(data.scores1$Treatment, levels = c("REF", "CR", "FF", "CFFR"))
data.scores1$Treatment
head(data.scores1)

ord1 <- ggplot(data.scores1, aes(x=MDS1, y=MDS2)) + 
  geom_point(data=data.scores1,aes(x=MDS1,y=MDS2,shape=data.scores1$Treatment,color=data.scores1$Treatment),size=4) + # add the point markers
  stat_ellipse(data = data.scores1, level = 0.95, geom = "polygon", alpha = 0.25, aes(fill = Treatment),type = "norm") +
  scale_color_manual(name = "Treatment",labels=c("REF","CR","FF","CFFR"),values=c('#D55E00','#56B4E9', '#F0E442', 'darkgreen')) +
  scale_shape(name = "Treatment",labels=c("REF","CR","FF","CFFR"))+
  scale_fill_manual(values=c('#D55E00','#56B4E9', '#F0E442', 'darkgreen'),guide="none") +
  annotate("text", x=-.6,y=.75,label="Treatment: ~italic(P) == 0.29", parse=TRUE, size=4, fontface="bold")+
  labs(color='Treatment') +
  theme(legend.text=element_text(size=12)) +
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  theme_classic() +
  theme(axis.text=element_text(size=12)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=12)) +
  ggtitle(expression("16S ASVs")) 
ord1


#After treatment

set.seed(101)
adonis2(b1.2~Treatment*Month, data=c1.2, perm=999)

library(RVAideMemoire)

set.seed(101)
pairwise.perm.manova(b1.2,c1.2$Treatment,nperm=999,p.method = "fdr")

set.seed(101)
d1.2 <- metaMDS(b1.2, k=2, trymax=1000)

soil$DIN <- soil$NO3+soil$NH4

set.seed(101)
vec1<-envfit(d1.2, soil[,c(33,16,18,34)], perm=1000, na.rm=TRUE)
vec1

soil.scrs1 <- as.data.frame(scores(vec1, display = "vectors"))
soil.scrs1 <- cbind(soil.scrs1, soil = rownames(soil.scrs1))
soil.scrs1

soil.scrs1[soil.scrs1=="ph"] <- "pH**"
soil.scrs1[soil.scrs1=="DIN"] <- "DIN???"

data.scores1.1 <- as.data.frame(d1.2$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores1.1$Treatment <- c1.2$Treatment 
data.scores1.1 

data.scores1.1$Treatment <- factor(data.scores1.1$Treatment, levels = c("REF", "CR", "FF", "CFFR"))
data.scores1.1$Treatment
head(data.scores1.1)

ord1.1 <- ggplot(data.scores1.1, aes(x=MDS1, y=MDS2)) + 
  geom_point(data=data.scores1.1,aes(x=MDS1,y=MDS2,shape=data.scores1.1$Treatment,fill=data.scores1.1$Treatment),size=4) + # add the point markers
  #stat_ellipse(data = data.scores1, level = 0.95, geom = "polygon", alpha = 0.25, aes(fill = Treatment),type = "norm") +
  scale_fill_manual(name = "Treatment",labels=c("REF","CR","FF","CFFR"),values=c('#D55E00','#56B4E9', '#F0E442', 'darkgreen')) +
  scale_shape_manual(name = "Treatment",labels=c("REF","CR","FF","CFFR"),values=c(21,22,23,24))+
  #scale_fill_manual(values=c('#D55E00','#56B4E9', '#F0E442', 'darkgreen'),guide="none") +
  annotate("text", x=-.2,y=.28,label="Treatment: ~italic(P) == 0.066", parse=TRUE, size=5, fontface="bold")+
  annotate("text", x=-.2,y=.25,label="Stress = 0.12", size=4.5)+
  labs(color='Treatment') +
  theme(legend.text=element_text(size=13)) +
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  theme_classic() +
  theme(axis.text=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=13)) +
  ggtitle(expression("16S ASVs")) +
  geom_segment(data=soil.scrs1,aes(x=0,xend=.5*NMDS1,y=0,yend=.5*NMDS2), arrow = arrow(length = unit(0.25, "cm")),colour="black",size=.5,inherit.aes=FALSE) +
  geom_text(data = soil.scrs1, aes(x = .6*NMDS1, y = .59*NMDS2, label = soil),size = 4)
ord1.1


#Fungal Beta Diversity

#Before treatment

set.seed(101)
adonis2(b2.1~Treatment, data=c2.1, perm=999)

set.seed(101)
d2.1 <- metaMDS(b2.1, k=2, trymax=1000)


data.scores2 <- as.data.frame(d2.1$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores2$Treatment <- c2.1$Treatment 
head(data.scores2) 

data.scores2$Treatment <- factor(data.scores2$Treatment, levels = c("REF", "CR", "FF", "CFFR"))
data.scores2$Treatment
head(data.scores2)

ord2 <- ggplot(data.scores2, aes(x=MDS1, y=MDS2)) + 
  geom_point(data=data.scores2,aes(x=MDS1,y=MDS2,shape=data.scores2$Treatment,color=data.scores2$Treatment),size=4) + # add the point markers
  stat_ellipse(data = data.scores2, level = 0.95, geom = "polygon", alpha = 0.25, aes(fill = Treatment),type = "norm") +
  scale_color_manual(name = "Treatment",labels=c("REF","CR","FF","CFFR"),values=c('#D55E00','#56B4E9', '#F0E442', 'darkgreen')) +
  scale_shape(name = "Treatment",labels=c("REF","CR","FF","CFFR"))+
  scale_fill_manual(values=c('#D55E00','#56B4E9', '#F0E442', 'darkgreen'),guide="none") +
  annotate("text", x=-1,y=1.5,label="Treatment: ~italic(P) == 0.62", parse=TRUE, size=4, fontface="bold")+
  labs(color='Treatment') +
  theme(legend.text=element_text(size=12)) +
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  theme_classic() +
  theme(axis.text=element_text(size=12)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=12)) +
  ggtitle(expression("ITS ASVs")) 
ord2

#After treatment

#All

set.seed(101)
adonis2(b2.2~Treatment*Month, data=c2.2, perm=999)

set.seed(101)
d2.2 <- metaMDS(b2.2, k=2, trymax=10000)

set.seed(101)
vec2<-envfit(d2.2, soil[,c(33,16,18,34)], perm=1000, na.rm=TRUE)
vec2

soil.scrs2 <- as.data.frame(scores(vec2, display = "vectors"))
soil.scrs2 <- cbind(soil.scrs2, soil = rownames(soil.scrs2))
soil.scrs2

soil.scrs2[soil.scrs2=="ph"] <- "pH"
soil.scrs2[soil.scrs2=="DIN"] <- "DIN*"

data.scores2.1 <- as.data.frame(d2.2$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores2.1$Treatment <- c2.2$Treatment 
data.scores2.1 

data.scores2.1$Treatment <- factor(data.scores2.1$Treatment, levels = c("REF", "CR", "FF", "CFFR"))
data.scores2.1$Treatment
head(data.scores2.1)

ord2.1 <- ggplot(data.scores2.1, aes(x=MDS1, y=MDS2)) + 
  geom_point(data=data.scores2.1,aes(x=MDS1,y=MDS2,shape=data.scores2.1$Treatment,fill=data.scores2.1$Treatment),size=4) + # add the point markers
  #stat_ellipse(data = data.scores1, level = 0.95, geom = "polygon", alpha = 0.25, aes(fill = Treatment),type = "norm") +
  scale_fill_manual(name = "Treatment",labels=c("REF","CR","FF","CFFR"),values=c('#D55E00','#56B4E9', '#F0E442', 'darkgreen')) +
  scale_shape_manual(name = "Treatment",labels=c("REF","CR","FF","CFFR"),values=c(21,22,23,24))+
  #scale_fill_manual(values=c('#D55E00','#56B4E9', '#F0E442', 'darkgreen'),guide="none") +
  annotate("text", x=-.1,y=.6,label="Treatment: ~italic(P) == 0.064", parse=TRUE, size=5, fontface="bold")+
  annotate("text", x=-.1,y=.53,label="Stress = 0.24", size=4.5)+
  labs(color='Treatment') +
  theme(legend.text=element_text(size=13)) +
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  theme_classic() +
  theme(axis.text=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=13)) +
  ggtitle(expression("ITS ASVs")) +
  geom_segment(data=soil.scrs2,aes(x=0,xend=1.3*NMDS1,y=0,yend=1.3*NMDS2), arrow = arrow(length = unit(0.25, "cm")),colour="black",size=.5,inherit.aes=FALSE) +
  geom_text(data = soil.scrs2, aes(x = 1.45*NMDS1, y = 1.45*NMDS2, label = soil),size = 4)
ord2.1


#Ordination two panel figure

setwd("C:/Users/ernie/OneDrive/Desktop/Rhodo Sequences")

jpeg(filename="fig1.jpeg", bg="transparent", res=500, units = "in", height=4, width=10.5) 

f1 <- plot_grid(ord1.1, ord2.1, ncol = 2, labels=c('A', 'B'),align="hv", label_size=20)

f1

dev.off()


#Bacterial Alpha Diversity


bac_alpha2 <- estimate_richness(bac2, split = TRUE)

bac_alpha2 <- cbind(bac_alpha2, c1.2)

glm2 = lm(Observed~ Treatment*Month,data=bac_alpha2)
shapiro.test(resid(glm2))
glm2 = glm(Observed~ Treatment*Month,data=bac_alpha2,family=Gamma(link=log))
Anova(glm2)
rich<- aggregate(Observed ~ Treatment, data=bac_alpha2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

rich<- do.call(data.frame, rich)
rich

rich$se <- rich$Observed.sd / sqrt(rich$Observed.n)
rich


glm2 = lm(Simpson~ Treatment*Month,data=bac_alpha2)
shapiro.test(resid(glm2))
Anova(glm2)
simp<- aggregate(Simpson ~ Treatment, data=bac_alpha2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

simp<- do.call(data.frame, simp)
simp

simp$se <- simp$Simpson.sd / sqrt(simp$Simpson.n)
simp


glm2 = lm(Shannon~ Treatment*Month,data=bac_alpha2)
shapiro.test(resid(glm2))
Anova(glm2)


bac_alpha2$Treatment <- factor(bac_alpha2$Treatment, levels = c("REF", "CR", "FF", "CFFR"))
bac_alpha2$Treatment

bac_box <-  ggplot(bac_alpha2, aes(x=bac_alpha2$Treatment, y=bac_alpha2$Shannon, fill=bac_alpha2$Treatment)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.shape=NA) +
  annotate('text', 2, 6.7, label="Treatment: italic(P) == 0.59", size = 5, parse=TRUE)+
  theme_classic() +
  geom_point(position=position_jitter(.25),size=2) +
  stat_summary(fun.y=mean,geom="point",shape=23,size=3,color="black",fill="white")+
  scale_fill_manual(values=c('#D55E00','#56B4E9', '#F0E442', 'darkgreen')) +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=18)) +
  theme(legend.text = element_text(color = "black", size = 18)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Bacteria") +
  scale_y_continuous(expression("Shannon Diversity (H')")) +
  theme(axis.title.x=element_blank())
bac_box

#Fungal Alpha Diversity

fung_alpha2 <- estimate_richness(ITS2, split = TRUE)

fung_alpha2 <- cbind(fung_alpha2, c2.2)

glm2 = lm(Observed~ Treatment*Month,data=fung_alpha2)
shapiro.test(resid(glm2))
Anova(glm2)
rich<- aggregate(Observed ~ Treatment, data=fung_alpha2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

rich<- do.call(data.frame, rich)
rich

rich$se <- rich$Observed.sd / sqrt(rich$Observed.n)
rich


glm2 = lm(Simpson~ Treatment*Month,data=fung_alpha2)
shapiro.test(resid(glm2))
glm2 = glm(Simpson~ Treatment*Month,data=fung_alpha2,family=Gamma(link=log))
Anova(glm2)
simp<- aggregate(Simpson ~ Treatment, data=fung_alpha2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

simp<- do.call(data.frame, simp)
simp

simp$se <- simp$Simpson.sd / sqrt(simp$Simpson.n)
simp

glm2 = lm(Shannon~ Treatment*Month,data=fung_alpha2)
shapiro.test(resid(glm2))
Anova(glm2)

fung_alpha2$Treatment <- factor(fung_alpha2$Treatment, levels = c("REF", "CR", "FF", "CFFR"))
fung_alpha2$Treatment

fung_box <-  ggplot(fung_alpha2, aes(x=fung_alpha2$Treatment, y=fung_alpha2$Shannon, fill=fung_alpha2$Treatment)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.shape=NA) +
  annotate('text', 2, 4.25, label="Treatment: italic(P) == 0.65", size = 5, parse=TRUE)+
  theme_classic() +
  geom_point(position=position_jitter(.25),size=2) +
  stat_summary(fun.y=mean,geom="point",shape=23,size=3,color="black",fill="white")+
  scale_fill_manual(values=c('#D55E00','#56B4E9', '#F0E442', 'darkgreen')) +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=18)) +
  theme(legend.text = element_text(color = "black", size = 18)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Fungi") +
  scale_y_continuous(expression("Shannon Diversity (H')")) +
  theme(axis.title.x=element_blank())
fung_box

#Multipanel shannon diversity figure

jpeg(filename="fig2.jpeg", bg="transparent", res=500, units = "in", height=4, width=7) 

f2 <- plot_grid(bac_box, fung_box, ncol = 2, labels=c('A', 'B'),align="hv", label_size=20)

f2

dev.off()

#bacterial phyla stuff

colnames(tax_table(bac2)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

rk <- transform_sample_counts(bac2, function(x) x/sum(x))

# agglomerate taxa
glom_rk <- tax_glom(rk, taxrank = 'Phylum')

# create dataframe from phyloseq object
dat_rk <- psmelt(glom_rk)
head(dat_rk)

# convert Phylum to a character vector from a factor because R
dat_rk$Phylum <- as.character(dat_rk$Phylum)

library(plyr)
# group dataframe by Phylum, calculate median rel. abundance
means <- ddply(dat_rk, ~Phylum, function(x) c(mean=mean(x$Abundance)))
means

# find Phyla whose rel. abund. is less than 1%
remainder <- means[means$mean <= 0.01,]$Phylum

# change their name to "Other"
dat_rk[dat_rk$Phylum %in% remainder,]$Phylum <- 'Other'

Coweeta_phyla <- aggregate(Abundance~Phylum, dat_rk, FUN=mean)

Coweeta_phyla$Abundance <- Coweeta_phyla$Abundance * 100

is.num <- sapply(Coweeta_phyla, is.numeric)
Coweeta_phyla[is.num] <- lapply(Coweeta_phyla[is.num], round, 1)

Coweeta_phyla

rk_phyla <- aggregate(Abundance~Sample+Phylum+Treatment+Month, dat_rk, FUN=sum)

rk_phyla1 <- cast(rk_phyla, Sample+Treatment+Month ~ Phylum, value="Abundance")

p1 <- lm(p__Acidobacteria~ Treatment*Month, data=rk_phyla1)
shapiro.test(resid(p1))
Anova(p1)

p2 <- lm(p__Proteobacteria~Treatment*Month, data=rk_phyla1)
shapiro.test(resid(p2))
Anova(p2)

p3 <- lm(p__Verrucomicrobia~ Treatment*Month, data=rk_phyla1)
shapiro.test(resid(p3))
Anova(p3)
aggregate(p__Verrucomicrobia~Treatment,data=rk_phyla1,FUN=mean)

p4 <- lm(p__Actinobacteria~ Treatment*Month, data=rk_phyla1)
shapiro.test(resid(p4))
Anova(p4)

p5 <- lm(p__Bacteroidetes~ Treatment*Month, data=rk_phyla1)
shapiro.test(resid(p5))
Anova(p5)

p7 <- lm(p__Planctomycetes~ Treatment*Month, data=rk_phyla1)
shapiro.test(resid(p7))
p7 <- glm(p__Planctomycetes~ Treatment*Month, data=rk_phyla1,family=Gamma(link=log))
Anova(p7)
aggregate(p__Planctomycetes~Treatment,data=rk_phyla1,FUN=mean)

p8 <- lm(p__Chloroflexi~ Treatment*Month, data=rk_phyla1)
shapiro.test(resid(p8))
p8 <- glm(p__Chloroflexi~ Treatment*Month, data=rk_phyla1,family=Gamma(link=log))
Anova(p8)


rk_phyla$Treatment <- factor(rk_phyla$Treatment, levels = c("REF", "CR", "FF", "CFFR"))
rk_phyla$Treatment

rk_phyla[rk_phyla=="p__Acidobacteria"] <- "Acidobacteria"
rk_phyla[rk_phyla=="p__Actinobacteria"] <- "Actinobacteria"
rk_phyla[rk_phyla=="p__Bacteroidetes"] <- "Bacteroidetes"
rk_phyla[rk_phyla=="p__Proteobacteria"] <- "Proteobacteria"
rk_phyla[rk_phyla=="p__Verrucomicrobia"] <- "Verrucomicrobia???"
rk_phyla[rk_phyla=="p__Chloroflexi"] <- "Chloroflexi"
rk_phyla[rk_phyla=="p__Planctomycetes"] <- "Planctomycetes???"


bac_bar <-ggplot(rk_phyla, aes(fill=Phylum, y=Abundance, x=Treatment)) + 
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous("Relative Abundance") +
  scale_fill_brewer(palette="Set1") +
  theme_classic() +
  theme(legend.box.margin=margin(-10,-10,-10,-10))+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Bacterial Phyla") +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=10)) 
bac_bar

#Bacterial functional groups

dat_rk2 <- psmelt(glom_rk)
head(dat_rk2)

# convert Phylum to a character vector from a factor because R
dat_rk2$Phylum <- as.character(dat_rk2$Phylum)

Coweeta_phyla2 <- aggregate(Abundance~Phylum, dat_rk2, FUN=mean)

Coweeta_phyla2$Abundance <- Coweeta_phyla2$Abundance * 100

is.num <- sapply(Coweeta_phyla2, is.numeric)
Coweeta_phyla2[is.num] <- lapply(Coweeta_phyla2[is.num], round, 1)

Coweeta_phyla2

rk_phyla2 <- aggregate(Abundance~Sample+Phylum+Treatment+Month, dat_rk2, FUN=sum)

rk_phyla2.1 <- cast(rk_phyla2, Sample+Treatment+Month ~ Phylum, value="Abundance")

p9 <- lm(p__Nitrospirae~ Treatment*Month, data=rk_phyla2.1)
shapiro.test(resid(p9))
p9 <- glm(p__Nitrospirae+.001~ Treatment*Month, data=rk_phyla2.1,family=Gamma(link=log))
Anova(p9)
nitro <- aggregate(p__Nitrospirae~Treatment,data=rk_phyla2.1,FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
nitro
nitro <- do.call(data.frame, nitro)
nitro
nitro$se <- nitro$p__Nitrospirae.sd/ sqrt(nitro$p__Nitrospirae.n)
nitro

rk_phyla2.1$rk <- ((rk_phyla2.1$p__Proteobacteria + rk_phyla2.1$p__Bacteroidetes)/(rk_phyla2.1$p__Acidobacteria+rk_phyla2.1$p__Verrucomicrobia))

p10 <- lm(rk~ Treatment*Month, data=rk_phyla2.1)
shapiro.test(resid(p10))
p10 <- glm(rk~ Treatment*Month, data=rk_phyla2.1,family=Gamma(link=log))
Anova(p10)
rk <- aggregate(rk~Treatment,data=rk_phyla2.1,FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
rk
rk <- do.call(data.frame, rk)
rk
rk$se <- rk$rk.sd/ sqrt(rk$rk.n)
rk

#fungal classes

colnames(tax_table(ITS2)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

fung_c <- transform_sample_counts(ITS2, function(x) x/sum(x))

# agglomerate taxa
glom_fung <- tax_glom(fung_c, taxrank = 'Class')

# create dataframe from phyloseq object
dat_fung <- psmelt(glom_fung)
head(dat_fung)

# convert Phylum to a character vector from a factor because R
dat_fung$Class <- as.character(dat_fung$Class)

library(plyr)
# group dataframe by Phylum, calculate median rel. abundance
means_f <- ddply(dat_fung, ~Class, function(x) c(mean=mean(x$Abundance)))
means_f

# find Phyla whose rel. abund. is less than 1%
remainder2 <- means_f[means_f$mean <= 0.01,]$Class

# change their name to "Other"
dat_fung[dat_fung$Class %in% remainder2,]$Class <- 'Other'

Coweeta_fung <- aggregate(Abundance~Class, dat_fung, FUN=mean)

Coweeta_fung$Abundance <- Coweeta_fung$Abundance * 100

is.num <- sapply(Coweeta_fung, is.numeric)
Coweeta_fung[is.num] <- lapply(Coweeta_fung[is.num], round, 1)

Coweeta_fung

fung_classes <- aggregate(Abundance~Sample+Class+Treatment+Month, dat_fung, FUN=sum)

fung_classes1 <- cast(fung_classes, Sample+Treatment+Month ~ Class, value="Abundance")


p1 <- lm(c__Agaricomycetes~Treatment*Month, data=fung_classes1)
shapiro.test(resid(p1))
Anova(p1)

p2 <- lm(c__Eurotiomycetes~ Treatment*Month, data=fung_classes1)
shapiro.test(resid(p2))
p2 <- glm(c__Eurotiomycetes+.001~ Treatment*Month, data=fung_classes1,family=Gamma(link=log))
Anova(p2)

p3 <- lm(c__Geminibasidiomycetes ~Treatment*Month, data=fung_classes1)
shapiro.test(resid(p3))
p2 <- glm(c__Geminibasidiomycetes+.001~ Treatment*Month, data=fung_classes1,family=Gamma(link=log))
Anova(p3)

p4 <- lm(c__Leotiomycetes~ Treatment*Month, data=fung_classes1)
shapiro.test(resid(p4))
p4 <- glm(c__Leotiomycetes~ Treatment*Month, data=fung_classes1,family=Gamma(link=log))
Anova(p4)

p5 <- lm(c__Mortierellomycetes~ Treatment*Month, data=fung_classes1)
shapiro.test(resid(p5))
p5 <- glm(c__Mortierellomycetes~ Treatment*Month, data=fung_classes1, family=Gamma(link=log))
Anova(p5)

p6 <- lm(c__Archaeorhizomycetes~ Treatment*Month, data=fung_classes1)
shapiro.test(resid(p6))
p6 <- glm(c__Archaeorhizomycetes+.001~ Treatment*Month, data=fung_classes1,family=Gamma(link=log))
Anova(p6)

p7 <- lm(c__Sordariomycetes~ Treatment*Month, data=fung_classes1)
shapiro.test(resid(p7))
Anova(p7)

p8 <- lm(c__Tremellomycetes~ Treatment*Month, data=fung_classes1)
shapiro.test(resid(p8))
p8 <- glm(c__Tremellomycetes~ Treatment*Month, data=fung_classes1,family=Gamma(link=log))
Anova(p8)

p9 <- lm(c__Saccharomycetes~ Treatment*Month, data=fung_classes1)
shapiro.test(resid(p9))
p9 <- glm(c__Saccharomycetes+.001~ Treatment*Month, data=fung_classes1,family=Gamma(link=log))
Anova(p9)

fung_classes$Treatment <- factor(fung_classes$Treatment, levels = c("REF", "CR", "FF", "CFFR"))
fung_classes$Treatment

fung_classes[fung_classes=="c__Agaricomycetes"] <- "Agaricomycetes"
fung_classes[fung_classes=="c__Eurotiomycetes"] <- "Eurotiomycetes"
fung_classes[fung_classes=="c__Geminibasidiomycetes"] <- "Geminibasidiomycetes"
fung_classes[fung_classes=="c__Leotiomycetes"] <- "Leotiomycetes"
fung_classes[fung_classes=="c__Mortierellomycetes"] <- "Mortierellomycetes"
fung_classes[fung_classes=="c__Archaeorhizomycetes"] <- "Archaeorhizomycetes"
fung_classes[fung_classes=="c__Sordariomycetes"] <- "Sordariomycetes"
fung_classes[fung_classes=="c__Tremellomycetes"] <- "Tremellomycetes"
fung_classes[fung_classes=="c__Saccharomycetes"] <- "Saccharomycetes"

brewer.pal(9,"Set1")

fung_bar <-ggplot(fung_classes, aes(fill=Class, y=Abundance, x=Treatment)) + 
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous("Relative Abundance") +
  scale_fill_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
                      "#999999","black")) +
  theme_classic() +
  theme(legend.box.margin=margin(-10,-10,-10,-10))+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Fungal Classes") +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=10)) 
fung_bar

jpeg(filename="fig3.jpeg", bg="transparent", res=500, units = "in", height=4, width=8) 

f3 <- plot_grid(bac_bar, fung_bar, ncol = 2, labels=c('A', 'B'),align="hv", label_size=20)

f3

dev.off()

#Import ITS taxonomy

tax <- read.csv("C:/Users/ernie/Desktop/Rhodo sequences/its-taxonomy.csv")

its_asvs <- data.frame(otu_table(ITS2))

tax2 <- tax[tax$X.OTUID %in% rownames(its_asvs),]

tax2 <- tax2[order(tax2$X.OTUID),]

its_asvs <- cbind(rownames(its_asvs), data.frame(its_asvs, row.names=NULL))

colnames(its_asvs)[1] <- "OTUID"

its_asvs2 <- its_asvs[order(its_asvs$OTUID),]

its_asvs2$taxonomy <- tax2$taxonomy

write.csv(its_asvs2,"C:/Users/ernie/Desktop/Rhodo sequences/its-asvs-with-taxonomy.csv", row.names = FALSE)

#Funguild Analysis

guilds <- read.csv("C:/Users/ernie/OneDrive/Desktop/Rhodo sequences/rhodo-funguild.csv")

library(plyr)
guilds <- ddply(guilds,"Guild",numcolwise(sum))

library(funrar)

guilds.matrix2 <- t(guilds[,c(2:33)])

guilds.rel <- data.frame(make_relative(guilds.matrix2))

colnames(guilds.rel) <- guilds[,1]

guilds <- guilds.rel


guilds3 <- cbind(guilds, c1.2)


library(tidyr)
guilds4 <- gather(guilds3, Guild, Abundance, `Animal Pathogen`:`Wood Saprotroph`, factor_key=TRUE)

guilds4 <- guilds4[,-c(2:4)]

guilds5 <- aggregate(Abundance~id+Guild+Treatment+Month, guilds4, FUN=sum)

aggregate(Abundance*100~Guild, guilds4, FUN=mean)



guilds6 <- cast(guilds5, id+Treatment+Month ~ Guild, value="Abundance")

guilds6$Treatment <- as.factor(guilds6$Treatment)

guilds6$Month <- as.factor(guilds6$Month)

glm40 <- lm(`Animal Pathogen`~ Treatment*Month,data=guilds6)
shapiro.test(resid(glm40))
glm40 <- glm(`Animal Pathogen`+.001~ Treatment*Month,data=guilds6,family=Gamma(link=log))
Anova(glm40)
aggregate(`Animal Pathogen`~Treatment,data=guilds6, FUN=mean)

glm41 <- lm(Ectomycorrhizal~ Treatment*Month,data=guilds6)
shapiro.test(resid(glm41))
Anova(glm41)

glm42 <- lm(Endophyte~ Treatment*Month,data=guilds6)
shapiro.test(resid(glm42))
glm42 <- glm(Endophyte+.001~ Treatment*Month,data=guilds6,family=Gamma(link=log))
Anova(glm42)

glm43 <- lm(`Ericoid Mycorrhizal`~ Treatment*Month,data=guilds6)
shapiro.test(resid(glm43))
glm43 <- glm(`Ericoid Mycorrhizal`+.001~ Treatment*Month,data=guilds6,family=Gamma(link=log))
Anova(glm43)
aggregate(`Ericoid Mycorrhizal`*100~Treatment,data=guilds6, FUN=mean)

glm44 <- lm(`Plant Pathogen`~ Treatment*Month,data=guilds6)
shapiro.test(resid(glm44))
glm44 <- glm(`Plant Pathogen`+.001~ Treatment*Month,data=guilds6,family=Gamma(link=log))
Anova(glm44)
aggregate(`Plant Pathogen`~Treatment,data=guilds6, FUN=mean)

glm45 <- lm(`Soil Saprotroph`~ Treatment*Month,data=guilds6)
shapiro.test(resid(glm45))
Anova(glm45)

glm46 <- lm(`Undefined Saprotroph`~ Treatment*Month,data=guilds6)
shapiro.test(resid(glm46))
glm46 <- glm(`Undefined Saprotroph`+.001~ Treatment*Month,data=guilds6,family=Gamma(link=log))
Anova(glm46)

glm47 <- lm(`Fungal Parasite`~ Treatment*Month,data=guilds6)
shapiro.test(resid(glm47))
glm47 <- glm(`Fungal Parasite`+.01~ Treatment*Month,data=guilds6,family=Gamma(link=log))
Anova(glm47)
aggregate(`Fungal Parasite`~Treatment,data=guilds6, FUN=mean)

glm48 <- lm(`Wood Saprotroph`~ Treatment*Month,data=guilds6)
shapiro.test(resid(glm48))
glm48 <- glm(`Wood Saprotroph`+.001~ Treatment*Month,data=guilds6,family=Gamma(link=log))
Anova(glm48)
aggregate(`Wood Saprotroph`*100~Treatment,data=guilds6, FUN=mean)

glm49 <- lm(`Orchid Mycorrhizal`~ Treatment*Month,data=guilds6)
shapiro.test(resid(glm49))
glm49 <- glm(`Orchid Mycorrhizal`+.001~ Treatment*Month,data=guilds6,family=Gamma(link=log))
Anova(glm49)
aggregate(`Orchid Mycorrhizal`*100~Treatment,data=guilds6, FUN=mean)

glm50 <- lm(`Arbuscular Mycorrhizal`~ Treatment*Month,data=guilds6)
shapiro.test(resid(glm50))
glm50 <- glm(`Arbuscular Mycorrhizal`+.001~ Treatment*Month,data=guilds6,family=Gamma(link=log))
Anova(glm50)
aggregate(`Arbuscular Mycorrhizal`*100~Treatment,data=guilds6, FUN=mean)

FD_funguild3 <- aggregate(Abundance~Treatment+Guild, guilds5, FUN=mean)

head(FD_funguild3)

FD_funguild3$Guild <- as.character(FD_funguild3$Guild)

FD_funguild3[FD_funguild3=="Fungal Parasite"] <- "Fungal Parasite"
FD_funguild3[FD_funguild3=="Arbuscular Mycorrhizal"] <- "Arbuscular Mycorrhizal***"
FD_funguild3[FD_funguild3=="Endophyte"] <- "Endophyte"
FD_funguild3[FD_funguild3=="Animal Pathogen"] <- "Animal Pathogen"
FD_funguild3[FD_funguild3=="Ericoid Mycorrhizal"] <- "Ericoid Mycorrhizal*"
FD_funguild3[FD_funguild3=="Plant Pathogen"] <- "Plant Pathogen"
FD_funguild3[FD_funguild3=="Wood Saprotroph"] <- "Wood Saprotroph*"
FD_funguild3[FD_funguild3=="Orchid Mycorrhizal"] <- "Orchid Mycorrhizal***"

head(FD_funguild3)

library(tidyr)

FD_funguild3 <- spread(FD_funguild3, Guild, Abundance)
FD_funguild3

FD_funguild4 <- as.matrix(FD_funguild3[,c(2:13),])
FD_funguild4 <- FD_funguild4[c(4,2,3,1),]

rownames(FD_funguild4) <- c("REF", "CR","FF", "CFFR")
FD_funguild4 <- t(scale(FD_funguild4)) 

#Make funguild heatmap

my_palette <- colorRampPalette(c("white", "#009E73"))(n = 299)

library(grid)

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("white", "#377EB8"))(n = 299)

library(pheatmap)

jpeg(filename="fig4.jpeg", bg="transparent", res=500, units = "in", height=4, width=5) 

heat1 <- pheatmap(FD_funguild4,color=my_palette,cluster_cols = F, cellwidth=25,cluster_rows = T, main = "Fungal Guilds", cex=1.1)
heat1

dev.off()

#Differentially Abundant ASVs

library(edgeR)

counts1 <- data.frame(otu_table(bac2))

e1 <- DGEList(counts=counts1, group = c1.2$Treatment)

e1 <- calcNormFactors(e1)

design1 <- model.matrix(~0+c1.2$Treatment)
e1 <- estimateDisp(e1,design1)

et1 <- exactTest(e1, pair=c("FF","REF"))
et2 <- exactTest(e1, pair=c("CR","REF"))
et3 <- exactTest(e1, pair=c("CFFR","REF"))

t1.1 <- as.data.frame(et1$table)
pvalue.1 <- subset(t1.1, PValue<=0.05)

t1.2 <- as.data.frame(et2$table)
pvalue.2 <- subset(t1.2, PValue<=0.05)

t1.3 <- as.data.frame(et3$table)
pvalue.3 <- subset(t1.3, PValue<=0.05)

pvalue <- rbind(pvalue.1,pvalue.2,pvalue.3)

pvalue2 <- data.frame(unique(row.names(pvalue)))

unique(row.names(pvalue))

#Stats on phyla of differentially abundant 16s OTUs

bac_diffabun <- psmelt(bac2)

bac_diffabun2 <- bac_diffabun[bac_diffabun$OTU %in% rownames(pvalue),]

bac_diffabun2.11 <- aggregate(Abundance ~ OTU + Phylum, data=bac_diffabun2, FUN=sum)

bac_diffabun2.1 <- aggregate(Abundance ~ Treatment + Month + Sample + Phylum, data=bac_diffabun2, FUN=sum)

bac_diffabun2.1$Abundance <- bac_diffabun2.1$Abundance / sum(bac_diffabun2.1$Abundance)

bac_diffabun2.1$Phylum <- as.character(bac_diffabun2.1$Phylum)

# group dataframe by Phylum, calculate mean rel. abundance
means3 <-aggregate(Abundance~ Phylum, data=bac_diffabun2.1, FUN=sum)

means3

means3$Abundance <- means3$Abundance * 100

is.num <- sapply(means3, is.numeric)
means3[is.num] <- lapply(means3[is.num], round, 1)

means3


# find Phyla whose rel. abund. is less than 1%
#remainder3 <- means3[means3$Abundance <= 0.01,]$Phylum

bac_diffabun2$Phylum <- as.character(bac_diffabun2$Phylum)

bac_diffabun2 <- aggregate(Abundance ~ Treatment + Month + Sample + Phylum, data=bac_diffabun2, FUN=sum)

#bac_diffabun2[bac_diffabun2$Phylum %in% remainder3,]$Phylum <- 'Other'

bac_diffabun2 <-aggregate(Abundance~ Phylum + Treatment + Sample + Month, data=bac_diffabun2, FUN=sum)

library(reshape)

bac_diffabun3 <- cast(bac_diffabun2, Sample+Treatment+Month ~ Phylum, value="Abundance")

mean(rowSums(bac_diffabun3))

relab1 <- as.matrix(bac_diffabun3[,c(4:13)])

relab2 <- apply(relab1, 1, function(i) i/sum(i))

relab3 <- as.data.frame(t(relab2))

bac_diffabun4 <- cbind(bac_diffabun3[,c(1:3)],relab3)

colnames(bac_diffabun4) <- colnames(bac_diffabun3)

p1 <- lm(p__Acidobacteria~ Treatment*Month, data=bac_diffabun4)
shapiro.test(resid(p1))
Anova(p1)

p2 <- lm(p__Proteobacteria~Treatment*Month, data=bac_diffabun4)
shapiro.test(resid(p2))
Anova(p2)
aggregate(p__Proteobacteria~Treatment,data=bac_diffabun4,FUN=mean)

p3 <- lm(p__Verrucomicrobia~ Treatment*Month, data=bac_diffabun4)
shapiro.test(resid(p3))
p3 <- glm(p__Verrucomicrobia+.001~ Treatment*Month, data=bac_diffabun4,family=Gamma(link=log))
Anova(p3)

p4 <- lm(p__Actinobacteria~ Treatment*Month, data=bac_diffabun4)
shapiro.test(resid(p4))
p4 <- glm(p__Actinobacteria+.001~ Treatment*Month, data=bac_diffabun4,family=Gamma(link=log))
Anova(p4)
aggregate(p__Actinobacteria*100~Treatment,data=bac_diffabun4,FUN=mean)

p5 <- lm(p__Bacteroidetes~ Treatment*Month, data=bac_diffabun4)
shapiro.test(resid(p5))
Anova(p5)
aggregate(p__Bacteroidetes~Treatment,data=bac_diffabun4,FUN=mean)

p7 <- lm(p__Planctomycetes~ Treatment*Month, data=bac_diffabun4)
shapiro.test(resid(p7))
p7 <- glm(p__Planctomycetes+.001~ Treatment*Month, data=bac_diffabun4,family=Gamma(link=log))
Anova(p7)
aggregate(p__Planctomycetes*100~Treatment,data=bac_diffabun4,FUN=mean)

p8 <- lm(p__FCPU426~ Treatment*Month, data=bac_diffabun4)
shapiro.test(resid(p8))
p8 <- glm(p__FCPU426+.001~ Treatment*Month, data=bac_diffabun4,family=Gamma(link=log))
Anova(p8)
aggregate(p__FCPU426~Treatment,data=bac_diffabun4,FUN=mean)

p9 <- lm(p__Elusimicrobia~ Treatment*Month, data=bac_diffabun4)
shapiro.test(resid(p9))
p9 <- glm(p__Elusimicrobia+.001~ Treatment*Month, data=bac_diffabun4,family=Gamma(link=log))
Anova(p9)
aggregate(p__Elusimicrobia*100~Treatment,data=bac_diffabun4,FUN=mean)

p10 <- lm(p__Cyanobacteria~ Treatment*Month, data=bac_diffabun4)
shapiro.test(resid(p10))
p10 <- glm(p__Cyanobacteria+.001~ Treatment*Month, data=bac_diffabun4,family=Gamma(link=log))
Anova(p10)
aggregate(p__Cyanobacteria*100~Treatment,data=bac_diffabun4,FUN=mean)

p11 <- lm(p__Chlamydiae~ Treatment*Month, data=bac_diffabun4)
shapiro.test(resid(p11))
p11 <- glm(p__Chlamydiae+.001~ Treatment*Month, data=bac_diffabun4,family=Gamma(link=log))
Anova(p11)
aggregate(p__Chlamydiae*100~Treatment,data=bac_diffabun4,FUN=mean)


bac_diffabun5 <- aggregate(Abundance~Treatment+Phylum, bac_diffabun2, FUN=mean)

bac_diffabun5[bac_diffabun5=="p__Acidobacteria"] <- "Acidobacteria"
bac_diffabun5[bac_diffabun5=="p__Actinobacteria"] <- "Actinobacteria***"
bac_diffabun5[bac_diffabun5=="p__FCPU426"] <- "FCPU426**"
bac_diffabun5[bac_diffabun5=="p__Bacteroidetes"] <- "Bacteroidetes???"
bac_diffabun5[bac_diffabun5=="p__Proteobacteria"] <- "Proteobacteria???"
bac_diffabun5[bac_diffabun5=="p__Verrucomicrobia"] <- "Verrucomicrobia"
bac_diffabun5[bac_diffabun5=="p__Planctomycetes"] <- "Planctomycetes???"
bac_diffabun5[bac_diffabun5=="p__Elusimicrobia"] <- "Elusimocrobia**"
bac_diffabun5[bac_diffabun5=="p__Cyanobacteria"] <- "Cyanobacteria***"
bac_diffabun5[bac_diffabun5=="p__Chlamydiae"] <- "Chlamydiae***"

head(bac_diffabun5)

library(tidyr)

bac_diffabun5a <- spread(bac_diffabun5, Phylum, Abundance)
bac_diffabun5a

bac_diffabun6 <- as.matrix(bac_diffabun5a[,c(2:11)])
bac_diffabun6

bac_diffabun7 <- bac_diffabun6[c(4,2,3,1),]
bac_diffabun7

rownames(bac_diffabun7) <- c("REF","CR","FF","CFFR")

#Differentially Abundant ASVs - fungi

library(edgeR)

counts1 <- data.frame(otu_table(ITS2))

e1 <- DGEList(counts=counts1, group = c2.2$Treatment)

e1 <- calcNormFactors(e1)

design1 <- model.matrix(~0+c2.2$Treatment)
e1 <- estimateDisp(e1,design1)

et1 <- exactTest(e1, pair=c("FF","REF"))
et2 <- exactTest(e1, pair=c("CR","REF"))
et3 <- exactTest(e1, pair=c("CFFR","REF"))

t1.1 <- as.data.frame(et1$table)
pvalue.1 <- subset(t1.1, PValue<=0.05)

t1.2 <- as.data.frame(et2$table)
pvalue.2 <- subset(t1.2, PValue<=0.05)

t1.3 <- as.data.frame(et3$table)
pvalue.3 <- subset(t1.3, PValue<=0.05)

pvalue <- rbind(pvalue.1,pvalue.2,pvalue.3)

pvalue2 <- data.frame(unique(row.names(pvalue)))

unique(row.names(pvalue))

#Stats on phyla of differentially abundant  ITS

fung_diffabun <- psmelt(ITS2)

fung_diffabun2 <- fung_diffabun[fung_diffabun$OTU %in% rownames(pvalue),]

fung_diffabun2.1 <- aggregate(Abundance ~ Treatment + Month + Sample + Class, data=fung_diffabun2, FUN=sum)

fung_diffabun2.1$Abundance <- fung_diffabun2.1$Abundance / sum(fung_diffabun2.1$Abundance)

fung_diffabun2.1$Class <- as.character(fung_diffabun2.1$Class)

# group dataframe by Phylum, calculate mean rel. abundance
means3 <-aggregate(Abundance~ Class, data=fung_diffabun2.1, FUN=sum)

means3

# find Phyla whose rel. abund. is less than 1%
#remainder3 <- means3[means3$Abundance <= 0.01,]$Phylum

fung_diffabun2$Class <- as.character(fung_diffabun2$Class)

fung_diffabun2 <- aggregate(Abundance ~ Treatment + Month + Sample + Class, data=fung_diffabun2, FUN=sum)

#fung_diffabun2[fung_diffabun2$Phylum %in% remainder3,]$Phylum <- 'Other'

fung_diffabun2 <-aggregate(Abundance~ Class + Treatment + Sample + Month, data=fung_diffabun2, FUN=sum)

library(reshape)

fung_diffabun3 <- cast(fung_diffabun2, Sample+Treatment+Month ~ Class, value="Abundance")

mean(rowMeans(fung_diffabun3))

relab1 <- as.matrix(fung_diffabun3[,c(4:5)])

relab2 <- apply(relab1, 1, function(i) i/sum(i))

relab3 <- as.data.frame(t(relab2))

fung_diffabun4 <- cbind(fung_diffabun3[,c(1:3)],relab3)

colnames(fung_diffabun4) <- colnames(fung_diffabun3)

p7 <- lm(c__Archaeorhizomycetes~ Treatment*Month, data=fung_diffabun4)
shapiro.test(resid(p7))
p7 <- glm(c__Archaeorhizomycetes+.01~ Treatment*Month, data=fung_diffabun4,family=Gamma(link=log))
Anova(p7)
aggregate(c__Archaeorhizomycetes~Treatment,data=fung_diffabun4,FUN=mean)


p8 <- lm(c__Agaricomycetes~ Treatment*Month, data=fung_diffabun4)
shapiro.test(resid(p8))
p8 <- glm(c__Agaricomycetes+.001~ Treatment*Month, data=fung_diffabun4,family=Gamma(link=log))
Anova(p8)
aggregate(c__Agaricomycetes~Treatment,data=fung_diffabun4,FUN=mean)

its_diffabun5 <- aggregate(Abundance~Treatment+Class, fung_diffabun2, FUN=mean)

its_diffabun5[its_diffabun5 =="c__Agaricomycetes"] <- "Agaricomycetes"
its_diffabun5[its_diffabun5 =="c__Archaeorhizomycetes"] <- "Archaeorhizomycetes*"

head(its_diffabun5)

library(tidyr)

its_diffabun5a <- spread(its_diffabun5, Class, Abundance)
its_diffabun5a

its_diffabun6 <- as.matrix(its_diffabun5a[,c(2:3)])
its_diffabun6

its_diffabun7 <- its_diffabun6[c(4,2,3,1),]
its_diffabun7

rownames(its_diffabun7) <- c("REF","CR","FF","CFFR")

diffabun <- cbind(bac_diffabun7, its_diffabun7)

diffabun2 <- t(scale(diffabun)) 

#Make funguild heatmap

my_palette <- colorRampPalette(c("white", "#009E73"))(n = 299)

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

library(pheatmap)

its_diffabun8 <- data.frame(t(its_diffabun7))

its_diffabun8$Taxon <- "Fungi"

bac_diffabun8 <- data.frame(t(bac_diffabun7))

bac_diffabun8$Taxon <- "Bacteria"

diffabun3 <- rbind(its_diffabun8, bac_diffabun8)

aka2 = data.frame(Taxon = factor(c("Fungi","Fungi","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria")))
rownames(aka2)<-rownames(diffabun3)

aka3 = list(Taxon = c(Bacteria = "#4DAF4A" , Fungi= "#E41A1C"))

my_palette <- colorRampPalette(c("white", "#377EB8"))(n = 299)


jpeg(filename="fig4.jpeg", bg="transparent", res=500, units = "in", height=4, width=7.5) 

heat2 <- pheatmap(diffabun2,color=my_palette,cluster_cols = F,annotation_names_row=F,annotation_colors = aka3[1], cellwidth=25,annotation_row = aka2,cluster_rows = T, main = "Differentially abundant ASVs", cex=1.1)
heat2

dev.off()
