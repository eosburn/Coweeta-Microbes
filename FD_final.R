library(biomformat)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)
library(car)
library(lme4)
library(funrar)
library(reshape)
library(emmeans)
library(cowplot)
library(caret)
library(RANN)
library(MASS)

#Import and prepare data 

setwd("C:/Users/ernie/Desktop/")

soil <- read.csv("C:/Users/ernie/Desktop/Data/Chapter 2 Data/FD_soil.csv")

names(soil)[names(soil) == "C.N"] <- "C:N"

soil2 <- soil[-37,]

veg <- read.csv("C:/Users/ernie/Desktop/Data/Chapter 2 Data/FD_veg.csv")

Bac_otus <- "C:/Users/ernie/Desktop/Data/Chapter 2 Data/FD-16s-table-with-taxonomy-final.biom"

Fung_otus <- "C:/Users/ernie/Desktop/Data/Chapter 2 Data/FD-ITS-table-with-taxonomy-final.biom"

funguild <- read.csv("C:/Users/ernie/Desktop/Data/Chapter 2 Data/funguild_final_relabun.csv")

funguild = setNames(data.frame(t(funguild[,-1])), funguild[,1])

x1 <- read_biom(Bac_otus)

x2 <- read_biom(Fung_otus)

meta <- "C:/Users/ernie/Desktop/Data/Chapter 2 Data/FD_16s_metadata.txt"

meta2 <- import_qiime_sample_data(meta)

OTU1 <- import_biom(x1)

OTU2 <- import_biom(x2)

OTU1 <- merge_phyloseq(OTU1, meta2)

OTU2 <- merge_phyloseq(OTU2, meta2)

OTU2 <- subset_samples(OTU2, id != "FD37")

FD_ITS <- rarefy_even_depth(OTU2, rngseed=TRUE)

Pine_ITS <- subset_samples(FD_ITS, Watershed.Pair == "Pine Conv.")

Clear_ITS <- subset_samples(FD_ITS, Watershed.Pair == "Clear Cut")

Pasture_ITS <- subset_samples(FD_ITS, Watershed.Pair == "Pasture Conv.")

Cable_ITS <- subset_samples(FD_ITS, Watershed.Pair == "Cable Logged")

FD_16s <- rarefy_even_depth(OTU1, rngseed=TRUE)

Pine_16s <- subset_samples(FD_16s, Watershed.Pair == "Pine Conv.")

Clear_16s <- subset_samples(FD_16s, Watershed.Pair == "Clear Cut")

Pasture_16s <- subset_samples(FD_16s, Watershed.Pair == "Pasture Conv.")

Cable_16s <- subset_samples(FD_16s, Watershed.Pair == "Cable Logged")

a1 <- t(data.frame(otu_table(FD_16s)))

a1Pine <- t(data.frame(otu_table(Pine_16s)))

a1Clear<- t(data.frame(otu_table(Clear_16s)))

a1Pasture<- t(data.frame(otu_table(Pasture_16s)))

a1Cable <- t(data.frame(otu_table(Cable_16s)))

a2 <- t(data.frame(otu_table(FD_ITS)))

a2Pine <- t(data.frame(otu_table(Pine_ITS)))

a2Clear<- t(data.frame(otu_table(Clear_ITS)))

a2Pasture<- t(data.frame(otu_table(Pasture_ITS)))

a2Cable <- t(data.frame(otu_table(Cable_ITS)))

b1 <- vegdist(a1, method = "bray")

b2 <- vegdist(a2, method = "bray")

c1 <- data.frame(sample_data(FD_16s))

c1Pine <- data.frame(sample_data(Pine_16s))

c1Clear <- data.frame(sample_data(Clear_16s))

c1Pasture <- data.frame(sample_data(Pasture_16s))

c1Cable <- data.frame(sample_data(Cable_16s))

c2 <- data.frame(sample_data(FD_ITS))

c2Pine <- data.frame(sample_data(Pine_ITS))

c2Clear <- data.frame(sample_data(Clear_ITS))

c2Pasture <- data.frame(sample_data(Pasture_ITS))

c2Cable <- data.frame(sample_data(Cable_ITS))

d1 <- cbind(a1, c1)

d1Pine <- cbind(a1Pine, c1Pine)

d1Clear <- cbind(a1Clear, c1Clear)

d1Pasture <- cbind(a1Pasture, c1Pasture)

d1Cable <- cbind(a1Cable, c1Cable)

d2 <- cbind(a2, c2)

d2Pine <- cbind(a2Pine, c2Pine)

d2Clear <- cbind(a2Clear, c2Clear)

d2Pasture <- cbind(a2Pasture, c2Pasture)

d2Cable <- cbind(a2Cable, c2Cable)

#Soil Chemistry 

sc1 <- aov(pH~Use*Pair, data=soil)
shapiro.test(resid(sc1))
summary(sc1)
sc1.2 <- aggregate(pH~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc1.2

sc2 <- aov(log10(Moisture)~Use*Pair, data=soil)
shapiro.test(resid(sc2))
summary(sc2)
sc2.2 <- aggregate(Moisture~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc2.2

sc3 <- aov(log10(NO3)~Use*Pair, data=soil)
shapiro.test(resid(sc3))
summary(sc3)
hist(soil$NO3)
glm1 = glm(NO3 ~ Use*Pair,data=soil, family=Gamma(link=log))
Anova(glm1)
sc3.2 <- aggregate(NO3~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc3.2

sc4 <- aov(NH4~Use*Pair, data=soil)
shapiro.test(resid(sc4))
summary(sc4)
sc4.2 <- aggregate(NH4~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc4.2

sc5 <- aov(log(TDN)~Use*Pair, data=soil)
shapiro.test(resid(sc5))
summary(sc5)
hist(soil$TDN)
glm2 = glm(TDN ~ Use*Pair,data=soil, family=Gamma(link=log));
Anova(glm2)
sc5.2 <- aggregate(TDN~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc5.2

sc6 <- aov(log10(DON)~Use*Pair, data=soil)
shapiro.test(resid(sc6))
summary(sc6)
hist(soil$DON)
glm6 = glm(DON ~ Use*Pair,data=soil, family=Gamma(link=log));
Anova(glm6)
sc6.2 <- aggregate(DON~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc6.2

sc7 <- aov(log10(TN) ~ Use*Pair, data=soil)
shapiro.test(resid(sc7))
summary(sc7)
sc9.2 <- aggregate(TN~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc9.2

sc8 <- aov(log10(TC) ~ Use*Pair, data=soil)
shapiro.test(resid(sc8))
summary(sc8)
sc8.2 <- aggregate(TC~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc8.2

sc9 <- aov(`C:N` ~ Use*Pair, data=soil)
shapiro.test(resid(sc9))
summary(sc9)
sc9.2 <- aggregate(`C:N`~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc9.2

sc10 <- aov(DOC ~ Use*Pair, data=soil)
shapiro.test(resid(sc10))
summary(sc10)
sc10.2 <- aggregate(DOC~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc10.2

sc11 <- aov(DOC.TDN ~ Use*Pair, data=soil)
shapiro.test(resid(sc11))
summary(sc11)
sc11.2 <- aggregate(DOC.TDN~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc11.2

sc12 <- aov(SIR ~ Use*Pair, data=soil)
shapiro.test(resid(sc12))
summary(sc12)
sc12.2 <- aggregate(SIR~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc12.2

sc14 <- aov(MBC ~ Use*Pair, data=soil)
shapiro.test(resid(sc14))
summary(sc14)
sc14.2 <- aggregate(MBC~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc14.2

sc15 <- aov(MBN ~ Use*Pair, data=soil)
shapiro.test(resid(sc15))
summary(sc15)
sc15.2 <- aggregate(MBN~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc15.2

sc16 <- aov(log10(MBC.N) ~ Use*Pair, data=soil)
shapiro.test(resid(sc16))
summary(sc16)
sc16.2 <- aggregate(MBC.N~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc16.2

#Bacterial and Fungal abundance analysis

bac_abun <- aov(bac ~ Use * Pair, data=soil)
shapiro.test(resid(bac_abun))
summary(bac_abun)
ag1 <- aggregate(bac~Use, data=soil, FUN=mean)
ag1
#Disturbed is 17% higher
em1 <- emmeans(bac_abun, ~ Use | Pair)
contrast(em1, method = "pairwise")

fung_abun <- aov(ITS ~ Use * Pair, data = soil)
shapiro.test(resid(fung_abun))
Anova(fung_abun)

FB_ratio <- aov(log10(soil$FB) ~ Use*Pair, data=soil)
shapiro.test(resid(FB_ratio))
Anova(FB_ratio)
ag2 <- aggregate(FB~Use, data=soil, FUN=mean)
ag2
#Disturbed is 23% higher
ag2.1 <- aggregate(FB~Use*Pair, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
ag2.1
em2 <- emmeans(FB_ratio, ~ Use | Pair)
contrast(em2, method = "pairwise")

#Bacterial abundance (for figure S1)

bacbox <- ggplot(soil, aes(x=Pair, y=bac, fill=Use)) + 
  stat_boxplot(geom = "errorbar", size=.25,width = 0.25, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.25, outlier.size=.25) +
  scale_fill_manual(values=c("#D55E00", "#0072B2")) +
  labs(list(x ="")) +
  theme_classic() +
  theme(axis.title=element_text(size=8)) +
  theme(legend.text=element_text(size=6)) +
  theme(text=element_text(size=8)) +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  theme(axis.title.x=element_blank()) +
  annotate("text", 1.5, 10, label="Disturbance: ~italic(P) == 0.07", parse=TRUE, size=1.75, fontface="bold")+  
  annotate("text", 1.5, 9.95, label="`Watershed Pair`: ~italic(P) == 0.72", parse=TRUE, size=1.75, fontface="bold")+           
  annotate("text", 1.5, 9.9, label = "Disturbance %*% Pair: ~italic(P) == 0.56", parse=TRUE, size=1.75, fontface="bold")+   
  annotate("text", 3, 9.95, label="*", size=3, fontface="bold")+
  theme(axis.line.x=element_line(colour="black", size=.25)) + 
  theme(axis.line.y=element_line(colour="black", size=.25)) + 
  theme(legend.title=element_text(colour="white")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  ggtitle("Bacterial Abundance") +
  theme(legend.position = "none") +
  coord_fixed(ratio=4.5) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=7)) +
  scale_y_continuous(bquote(~'Log'[10]~'(16s gene copies'~~gdw^-1*')'))
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

bacbox <- bacbox + 
  scale_x_discrete(breaks=unique(soil$Pair), 
                   labels=addline_format(c("Cable Logged", 
                                           "Pasture Conv.", "Pine Conv.",  "Clear Cut")))
plot(bacbox)

#Fungal abundance (For figure S1)

fungbox <- ggplot(soil, aes(x=Pair, y=ITS, fill=Use)) + 
  stat_boxplot(geom = "errorbar", width = 0.25, size=.25, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE, size=.25, outlier.size=.25) +
  scale_fill_manual(values=c("#D55E00", "#0072B2")) +
  labs(list(x ="")) +
  theme_classic() +
  annotate("text", 3.4, 8.6, label="Disturbance: ~italic(P) == 0.54", parse=TRUE, size=1.75, fontface="bold")+  
  annotate("text", 3.4, 8.55, label="`Watershed Pair`: ~italic(P) == 0.38", parse=TRUE, size=1.75, fontface="bold")+           
  annotate("text", 3.4, 8.5, label = "Disturbance %*% Pair: ~italic(P) == 0.31", parse=TRUE, size=1.75, fontface="bold")+   
  theme(axis.title=element_text(size=8)) +
  theme(text=element_text(size=8)) +
  theme(legend.position=c(.8, 1)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x=element_line(colour="black", size=.25)) + 
  theme(axis.line.y=element_line(colour="black", size=.25)) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.title=element_text(colour="white")) +
  theme(legend.position = "none") +
  ggtitle("Fungal Abundance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=7)) +
  coord_fixed(ratio=4.5) +
  scale_y_continuous(bquote(~'Log'[10]~'(ITS gene copies'~~gdw^-1*')'))
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

fungbox <- fungbox + 
  scale_x_discrete(breaks=unique(soil$Pair), 
                   labels=addline_format(c("Cable Logged", 
                                           "Pasture Conv.", "Pine Conv.",  "Clear Cut")))
plot(fungbox)

#F:B, watersheds separated (for FigureS1)

fbbox <- ggplot(soil, aes(x=Pair, y=FB, fill=Use)) + 
  stat_boxplot(geom = "errorbar", width = 0.25,size=.25, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.25, outlier.size=.25) +
  scale_fill_manual(values=c("#D55E00", "#0072B2")) +
  labs(list(x ="")) +
  theme_classic() +
  theme(axis.title=element_text(size=8)) +
  theme(text=element_text(size=8)) +
  theme(legend.position=c(.8, 1)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x=element_line(colour="black", size=.25)) + 
  theme(axis.line.y=element_line(colour="black", size=.25)) + 
  theme(legend.title=element_text(colour="white")) +
  annotate("text", 1.6, .095, label="Disturbance: ~italic(P) == 0.001", parse=TRUE, size=1.75, fontface="bold")+  
  annotate("text", 1.6, .09, label="`Watershed Pair`: ~italic(P) == 0.017", parse=TRUE, size=1.75, fontface="bold")+
  annotate("text", 1.6, .085, label="Disturbance %*% Pair: ~italic(P) == 0.006", parse=TRUE, size=1.75, fontface="bold")+ 
  ggtitle("Fungal:Bacterial Ratios") +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(axis.line.x=element_line(colour="black", size=.25)) + 
  theme(axis.line.y=element_line(colour="black", size=.25)) + 
  annotate("text", 3, .085, label="***", size=3, fontface="bold")+
  annotate("text", 4, .085, label="**", size=3, fontface="bold")+
  coord_fixed(ratio=42) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=7)) +
  scale_y_continuous(bquote(~'ITS:16S Gene copies')) 
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

fbbox <- fbbox + 
  scale_x_discrete(breaks=unique(soil$Pair), 
                   labels=addline_format(c("Cable Logged", 
                                           "Pasture Conv.", "Pine Conv.",  "Clear Cut")))
plot(fbbox)

#Aggregated F:B ratios (for Figure 2)

fbbox2 <- ggplot(soil, aes(x=Use, y=FB, fill=Use)) + 
  stat_boxplot(geom = "errorbar", width = 0.25,size=.25, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.25,outlier.size=-1) +
  geom_point(position=position_jitter(.25),aes(shape=Pair),size=0.9) +
  scale_fill_manual(values=c("#D55E00", "#0072B2")) +
  labs(list(x ="")) +
  theme_classic() +
  theme(axis.title=element_text(size=8)) +
  theme(text=element_text(size=9)) +
  theme(legend.position=c(.8, 1)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x=element_line(colour="black", size=.25)) + 
  theme(axis.line.y=element_line(colour="black", size=.25)) + 
  theme(legend.title=element_text(colour="white")) +
  annotate("text", 1.6, .085, label="Disturbance: ~italic(P) < 0.001", parse=TRUE, size=2.5, fontface="bold")+  
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(axis.line.x=element_line(colour="black", size=.25)) + 
  theme(axis.line.y=element_line(colour="black", size=.25)) + 
  coord_fixed(ratio=42) +
  scale_y_continuous(bquote(~'ITS:16S Gene copies')) 
plot(fbbox2)

#Make Panel Figure S1

jpeg(filename="FigureS1.jpeg", bg="transparent", res=300, units = "in", height=2, width=7.5) 

S1 <- plot_grid(bacbox, fungbox, fbbox, ncol = 3, labels=c('A', 'B','C'), align="hv", label_size=8)

S1

dev.off()

#16s Alpha diversity analysis

bac_alpha <- estimate_richness(FD_16s, split = TRUE, measures="shannon")

bac_alpha <- cbind(bac_alpha, c1)

shannon <- aov(Shannon ~ Treatment*Watershed.Pair, data = bac_alpha)
shapiro.test(resid(shannon))
Anova(shannon)
hist(bac_alpha$Shannon)
shannon1 <- glm(Shannon ~ Treatment*Watershed.Pair, data=bac_alpha, family=Gamma(link=log))
Anova(shannon1)
ag3 <- aggregate(Shannon~Treatment, data=bac_alpha, FUN=mean)
ag3
ag3.1 <- aggregate(Shannon~Treatment*Watershed.Pair, data=bac_alpha, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
ag3.1
em_shan <- emmeans(shannon1, ~ Treatment | Watershed.Pair)
contrast(em_shan, method = "pairwise")

#16s shannon watersheds separated (for figure S2)

BacShanBox <- ggplot(bac_alpha, aes(x=Watershed.Pair, y=Shannon, fill=Treatment)) + 
  stat_boxplot(geom = "errorbar", width = 0.25, size=.25,position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE, size=.25,outlier.size=.25) +
  scale_fill_manual(values=c("#D55E00", "#0072B2")) +
  labs(list(x ="")) +
  theme_classic() +
  theme(axis.title=element_text(size=8)) +
  annotate("text", 1.5, 6, label="Disturbance: ~italic(P) < 0.001", parse=TRUE, size=1.5, fontface="bold")+  
  annotate("text", 1.5, 5.9, label="`Watershed Pair`: ~italic(P) == 0.006", parse=TRUE, size=1.5, fontface="bold")+
  annotate("text", 1.5, 5.8, label="Disturbance %*% Pair: ~italic(P) == 0.25", parse=TRUE, size=1.5, fontface="bold")+ 
  theme(text=element_text(size=8)) +
  theme(legend.position=c(.2,1)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x=element_line(colour="black", size=.25)) + 
  theme(axis.line.y=element_line(colour="black", size=.25)) + 
  theme(legend.title=element_text(colour="white")) +
  ggtitle("Bacterial Shannon Diversity") +
  theme(plot.title = element_text(size=7)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  coord_fixed(ratio=2.5) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", 3, 5.9, label="**", size=3, fontface="bold")+
  theme(legend.position = "none") +
  scale_y_continuous(bquote(~'H')) 
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

BacShanBox <- BacShanBox + 
  scale_x_discrete(breaks=unique(bac_alpha$Watershed.Pair), 
                   labels=addline_format(c("Cable Logged", 
                                           "Pasture Conv.", "Pine Conv.",  "Clear Cut")))
plot(BacShanBox)

#16s shannon aggregated (For figure 2)

BacShanBox2 <- ggplot(bac_alpha, aes(x=Treatment, y=Shannon, fill=Treatment)) + 
  stat_boxplot(geom = "errorbar", width = 0.25, size=.25,position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE, size=.25, outlier.size=-1) +
  geom_point(position=position_jitter(.25), aes(shape=Watershed.Pair),size=.9) +
  scale_fill_manual(values=c("#D55E00", "#0072B2")) +
  labs(list(x ="")) +
  theme_classic() +
  theme(axis.title=element_text(size=8)) +
  annotate("text", 1.5, 5.8, label="Disturbance: ~italic(P) < 0.001", parse=TRUE, size=2.5, fontface="bold")+  
  theme(text=element_text(size=9)) +
  theme(legend.position=c(.2,1)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x=element_line(colour="black", size=.25)) + 
  theme(axis.line.y=element_line(colour="black", size=.25)) + 
  theme(legend.title=element_text(colour="white")) +
  theme(plot.title = element_text(size=6)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  coord_fixed(ratio=2.5) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  scale_y_continuous(expression("16S Shannon Diversity (H')"))
BacShanBox2

fung_alpha <- estimate_richness(FD_ITS, split = TRUE, measures=c("Shannon"))

fung_alpha <- cbind(fung_alpha, c2)

shannon.2 <- aov(Shannon ~ Treatment*Watershed.Pair, data = fung_alpha)
shapiro.test(resid(shannon.2))
Anova(shannon.2)

#ITS shannon watersheds separated (for figure S2)

FungShanBox <- ggplot(fung_alpha, aes(x=Watershed.Pair, y=Shannon, fill=Treatment)) + 
  stat_boxplot(geom = "errorbar", size=.25, width = 0.25, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.25, outlier.size=.25) +
  scale_fill_manual(values=c("#D55E00", "#0072B2")) +
  labs(list(x ="")) +
  theme_classic() +
  theme(axis.title=element_text(size=8)) +
  theme(text=element_text(size=8)) +
  theme(legend.position=c(.2,1)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x=element_line(colour="black", size=.25)) + 
  theme(axis.line.y=element_line(colour="black", size=.25)) + 
  theme(legend.title=element_text(colour="white")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  ggtitle("Fungal Shannon Diversity") +
  annotate("text", 3.25, 4.85, label="Disturbance: ~italic(P) == 0.24", parse=TRUE, size=1.5, fontface="bold")+  
  annotate("text", 3.25, 4.65, label="`Watershed Pair`: ~italic(P) == 0.43", parse=TRUE, size=1.5, fontface="bold")+
  annotate("text", 3.25, 4.45, label="Disturbance %*% Pair: ~italic(P) == 0.67", parse=TRUE, size=1.5, fontface="bold")+ 
  theme(plot.title = element_text(size=7)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.line.x=element_line(colour="black", size=.25)) + 
  theme(axis.line.y=element_line(colour="black", size=.25)) +
  coord_fixed(ratio=1.1) +
  theme(legend.position = "none") +
  scale_y_continuous(bquote(~'H')) 
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

FungShanBox <- FungShanBox + 
  scale_x_discrete(breaks=unique(fung_alpha$Watershed.Pair), 
                   labels=addline_format(c("Cable Logged", 
                                           "Pasture Conv.", "Pine Conv.",  "Clear Cut")))
plot(FungShanBox)

#ITS shannon aggregated (for figure 2)

FungShanBox2 <- ggplot(fung_alpha, aes(x=Treatment, y=Shannon, fill=Treatment)) + 
  stat_boxplot(geom = "errorbar", size=.25, width = 0.25, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.25, outlier.size=-1) +
  geom_point(position=position_jitter(.25),aes(shape=Watershed.Pair),size=.9) +
  scale_fill_manual(values=c("#D55E00", "#0072B2")) +
  labs(list(x ="")) +
  theme_classic() +
  theme(axis.title=element_text(size=8)) +
  theme(text=element_text(size=9)) +
  theme(legend.position=c(.2,1)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x=element_line(colour="black", size=.25)) + 
  theme(axis.line.y=element_line(colour="black", size=.25)) + 
  theme(legend.title=element_text(colour="white")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  annotate("text", 1.5, 4.6, label="Disturbance: ~italic(P) == 0.24", parse=TRUE, size=2.5, fontface="bold")+  
  theme(axis.line.x=element_line(colour="black", size=.25)) + 
  theme(axis.line.y=element_line(colour="black", size=.25)) +
  coord_fixed(ratio=1.1) +
  theme(legend.position = "none") +
  scale_y_continuous(expression("ITS Shannon Diversity (H')"))
FungShanBox2

#Make Panel Figure S2

jpeg(filename="FigureS2.jpeg", bg="transparent", res=300, units = "in", height=2, width=4) 

S2 <- plot_grid(BacShanBox, FungShanBox, ncol = 2, labels=c('A', 'B'), align="hv", label_size=6)
S2

dev.off()

#16s Beta diversity analysis

adonis2(a1 ~Treatment*Watershed.Pair, data = d1, permutations = 999, method = "bray")

adonis2(a1Pine ~Treatment, data = d1Pine, permutations = 9999, method = "bray")

adonis2(a1Clear ~Treatment, data = d1Clear, permutations = 9999, method = "bray")

adonis2(a1Pasture~Treatment, data = d1Pasture, permutations = 9999, method = "bray")

adonis2(a1Cable ~Treatment, data = d1Cable, permutations = 9999, method = "bray")

p1 <- c(.002, .002, .004, .025)

p.adjust(p1, method="hochberg")

e1 <- metaMDS(b1, k=2, trymax=1000)

soil <- soil[order(soil$id),]
head(soil)

vec<-envfit(e1, soil[,c(6,8,9,10,19)], perm=1000, na.rm=TRUE)
vec

soil.scrs <- as.data.frame(scores(vec, display = "vectors"))
soil.scrs <- cbind(soil.scrs, soil = rownames(soil.scrs))
soil.scrs

data.scores1 <- as.data.frame(e1$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores1$site <- rownames(data.scores1)  # create a column of site names, from the rownames of data.scores
data.scores1$treat <- c1$Treatment 
data.scores1$pair <- c1$Watershed.Pair
head(data.scores1) 

#16s ordination for figure 3

ord1 <- ggplot() + 
  geom_point(data=data.scores1,aes(x=MDS1,y=MDS2,shape=data.scores1$pair,color=data.scores1$treat),size=3) + # add the point markers
  scale_color_manual(values=c("#D55E00", "#0072B2")) +
  annotate("text", -.1, -.2, label="Disturbance: ~italic(P) == 0.001", parse=TRUE, size=4.5, fontface="bold")+
  #annotate("text", 0.027, -.022, label="`Watershed Pair`: ~italic(P) == 0.001", parse=TRUE, size=3, fontface="bold")+
  #annotate("text", 0.027, -.026, label="`Treatment x Pair`: ~italic(P) == 0.004", parse=TRUE, size=3, fontface="bold")+
  labs(color='Treatment') +
  labs(shape='Watershed Pair') +
  theme(legend.text=element_text(size=14)) +
  theme_classic() +
  theme(axis.text=element_text(size=14)) +
  ggtitle("16S OTUs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=14)) +
  coord_fixed() +
  geom_segment(data=soil.scrs,aes(x=0,xend=.21*NMDS1,y=0,yend=.21*NMDS2), arrow = arrow(length = unit(0.25, "cm")),colour="black",size=.25,inherit.aes=FALSE) +
  geom_text(data = soil.scrs, aes(x = .25*NMDS1, y = .25*NMDS2, label = soil),size = 5)
ord1

#legend for figure 2

legend2 <- ggplot() + 
  geom_point(data=data.scores1,aes(x=MDS1,y=MDS2,shape=data.scores1$pair),size=3) + # add the point markers
  scale_color_manual(values=c("#D55E00", "#0072B2")) +
  annotate("text", -.1, -.2, label="Disturbance: ~italic(P) == 0.001", parse=TRUE, size=4, fontface="bold")+
  #annotate("text", 0.027, -.022, label="`Watershed Pair`: ~italic(P) == 0.001", parse=TRUE, size=3, fontface="bold")+
  #annotate("text", 0.027, -.026, label="`Treatment x Pair`: ~italic(P) == 0.004", parse=TRUE, size=3, fontface="bold")+
  labs(color='Treatment') +
  labs(shape='Watershed Pair') +
  theme(legend.text=element_text(size=12)) +
  theme_classic() +
  theme(axis.text=element_text(size=12)) +
  ggtitle("16S OTUs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=12)) +
  coord_fixed() +
  geom_segment(data=soil.scrs,aes(x=0,xend=.3*NMDS1,y=0,yend=.3*NMDS2), arrow = arrow(length = unit(0.25, "cm")),colour="black",size=.25,inherit.aes=FALSE) 
legend2

#ITS Beta diversity analysis

adonis2(a2 ~Treatment*Watershed.Pair, data = d2, permutations = 999, method = "bray")

adonis2(a2Pine ~Treatment, data = d2Pine, permutations = 9999, method = "bray")

adonis2(a2Clear ~Treatment, data = d2Clear, permutations = 9999, method = "bray")

adonis2(a2Pasture~Treatment, data = d2Pasture, permutations = 9999, method = "bray")

adonis2(a2Cable ~Treatment, data = d2Cable, permutations = 9999, method = "bray")

p2 <- c(.0025, .002, .009, .065)

p.adjust(p2, method="hochberg")

e2 <- metaMDS(b2, k=2, trymax=1000)
e2

soil2 <- soil2[order(soil2$id),]
head(soil2)

vec2<-envfit(e2, soil2[,c(6,8,9,10,19)], perm=1000, na.rm=TRUE)
vec2

soil.scrs2 <- as.data.frame(scores(vec2, display = "vectors"))
soil.scrs2 <- cbind(soil.scrs2, soil2 = rownames(soil.scrs2))
soil.scrs2

data.scores2 <- as.data.frame(e2$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores2$site <- rownames(data.scores2)  # create a column of site names, from the rownames of data.scores
data.scores2$treat <- c2$Treatment 
data.scores2$pair <- c2$Watershed.Pair
head(data.scores2) 

#ITS ordination for figure 3

ord2 <- ggplot() + 
  geom_point(data=data.scores2,aes(x=MDS1,y=MDS2,shape=data.scores2$pair,color=data.scores2$treat),size=3) + # add the point markers
  scale_color_manual(values=c("#D55E00", "#0072B2")) +
  annotate("text", -.25, -.4, label="Disturbance: ~italic(P) == 0.001", parse=TRUE, size=4.5, fontface="bold")+
  #annotate("text", -.4, .52, label="`Watershed Pair`: ~italic(P) == 0.001", parse=TRUE, size=3, fontface="bold")+
  #annotate("text", -.4, .42, label="`Treatment x Pair`: ~italic(P) == 0.001", parse=TRUE, size=3, fontface="bold")+
  labs(color='Treatment') +
  labs(shape='Watershed Pair') +
  theme(legend.text=element_text(size=14)) +
  theme_classic() +
  theme(axis.text=element_text(size=14)) +
  ggtitle("ITS OTUs") +
  theme(text=element_text(size=14)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed() +
  geom_segment(data=soil.scrs2,aes(x=0,xend=.6*NMDS1,y=0,yend=.6*NMDS2), arrow = arrow(length = unit(0.25, "cm")),colour="black",size=.25,inherit.aes=FALSE) +
  geom_text(data = soil.scrs2, aes(x = .67*NMDS1, y = .67*NMDS2, label = soil2),size = 5)
ord2

#Make panel figure 3

jpeg(filename="Figure3.jpeg", bg="transparent", res=300, units = "in", height=9, width=12) 

f2 <- plot_grid(ord1, ord2,NULL,NULL, ncol = 2, labels=c('A', 'B', 'C', 'D'), align="hv", label_size=15)
f2

dev.off()

#Variance partitioning

soil <- soil[order(soil$id),]

soil_vp <- soil[,c(6:19)]

veg <- veg[order(veg$id),]

veg_vp <- veg[,c(7:24)]

soil.impute <- preProcess(soil_vp, method = c("knnImpute"))

soil_vp2 <- predict(soil.impute, soil_vp)

soil_ITS1 <- soil_vp2[-37,]

veg_ITS1 <- veg_vp[-37,]

varpart_16s <- varpart(b1, as.matrix(soil_vp2), as.matrix(veg_vp))

varpart_16s

plot(varpart_16s)

anova(dbrda(b1~as.matrix(soil_vp2)))

anova(dbrda(b1~as.matrix(veg_vp)))

varpart_ITS <- varpart(b2, as.matrix(soil_ITS1), as.matrix(veg_ITS1))

varpart_ITS

plot(varpart_ITS)

anova(dbrda(b2~as.matrix(soil_ITS1)))

anova(dbrda(b2~as.matrix(veg_ITS1)))

#16s Phylum analysis

colnames(tax_table(FD_16s)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

phy <- transform_sample_counts(FD_16s, function(x) x/sum(x))

# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Phylum')

# create dataframe from phyloseq object
dat <- psmelt(glom)

# convert Phylum to a character vector from a factor because R
dat$Phylum <- as.character(dat$Phylum)

library(plyr)
# group dataframe by Phylum, calculate median rel. abundance
means <- ddply(dat, ~Phylum, function(x) c(mean=mean(x$Abundance)))

# find Phyla whose rel. abund. is less than 1%
remainder <- means[means$mean <= 0.01,]$Phylum

# change their name to "Other"
dat[dat$Phylum %in% remainder,]$Phylum <- 'Other'

FD_16s_phyla <- aggregate(Abundance~Sample+Phylum+Treatment+Watershed.Pair+Description, dat, FUN=sum)

Coweeta_phyla <- aggregate(Abundance~Phylum, FD_16s_phyla, FUN=mean)

Coweeta_phyla[Coweeta_phyla=="p__Acidobacteria"] <- "Acidobacteria"
Coweeta_phyla[Coweeta_phyla=="p__Actinobacteria"] <- "Actinobacteria"
Coweeta_phyla[Coweeta_phyla=="p__Nitrospirae"] <- "Nitrospirae"
Coweeta_phyla[Coweeta_phyla=="p__Bacteroidetes"] <- "Bacteroidetes"
Coweeta_phyla[Coweeta_phyla=="p__Proteobacteria"] <- "Proteobacteria"
Coweeta_phyla[Coweeta_phyla=="p__Verrucomicrobia"] <- "Verrucomicrobia"
Coweeta_phyla[Coweeta_phyla=="p__Chloroflexi"] <- "Chloroflexi"
Coweeta_phyla[Coweeta_phyla=="p__Planctomycetes"] <- "Planctomycetes"

Coweeta_phyla

Coweeta_phyla$Abundance <- Coweeta_phyla$Abundance * 100

is.num <- sapply(Coweeta_phyla, is.numeric)
Coweeta_phyla[is.num] <- lapply(Coweeta_phyla[is.num], round, 1)

Coweeta_phyla

cbPalette1 <- c("#CC79A7","#E69F00",  "#009E73", "#F0E442","black", "#0072B2", "#D55E00", "#56B4E9", "#999999")

Coweeta_phyla$Phylum <- factor(Coweeta_phyla$Phylum, levels = c("Proteobacteria", "Acidobacteria","Nitrospirae", "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Planctomycetes", "Other", "Verrucomicrobia"))
Coweeta_phyla

#16s Phyla Pie Chart

jpeg(filename="16spie.jpeg", bg="transparent", res=500, units = "in", height=3, width=4) 

with(Coweeta_phyla,pie(Abundance, labels=paste0(as.character(Phylum), " ", Abundance, "%"), cex = .4, col=rev(cbPalette1), radius=1))

dev.off()

library(reshape)
library(car)

FD_16s_phyla1 <- cast(FD_16s_phyla, Sample+Treatment+Watershed.Pair ~ Phylum, value="Abundance")
FD_16s_phyla1

p1 <- aov(p__Acidobacteria ~ Treatment * Watershed.Pair, data = FD_16s_phyla1)
shapiro.test(resid(p1))
p1.1 <- aggregate(p__Acidobacteria ~ Treatment, data=FD_16s_phyla1, FUN=mean)
p1.1
p1.2 <- aggregate(p__Acidobacteria~Treatment*Watershed.Pair, data=FD_16s_phyla1, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
p1.2
#Reference is ~20% higher
Anova(p1)
em_p1 <- emmeans(p1, ~ Treatment | Watershed.Pair)
contrast(em_p1, method = "pairwise")

p2 <- aov(p__Proteobacteria ~ Treatment * Watershed.Pair, data = FD_16s_phyla1)
shapiro.test(resid(p2))
Anova(p2)
p2.1 <- aggregate(p__Proteobacteria ~ Treatment, data=FD_16s_phyla1, FUN=mean)
p2.1
p2.2 <- aggregate(p__Proteobacteria~Treatment*Watershed.Pair, data=FD_16s_phyla1, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
p2.2
#Disturbed is ~11% higher
em_p2 <- emmeans(p2, ~ Treatment | Watershed.Pair)
contrast(em_p2, method = "pairwise")

p3 <- aov(p__Verrucomicrobia ~ Treatment * Watershed.Pair, data = FD_16s_phyla1)
shapiro.test(resid(p3))
Anova(p3)

p4 <- aov(p__Chloroflexi ~ Treatment * Watershed.Pair, data = FD_16s_phyla1)
shapiro.test(resid(p4))
Anova(p4)
p4.1 <- aggregate(p__Chloroflexi ~ Treatment, data=FD_16s_phyla1, FUN=mean)
p4.1
p4.2 <- aggregate(p__Chloroflexi~Treatment*Watershed.Pair, data=FD_16s_phyla1, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
p4.2
#Disturbed is ~23% higher
em_p4 <- emmeans(p4, ~ Treatment | Watershed.Pair)
contrast(em_p4, method = "pairwise")

p5 <- aov(log10(p__Planctomycetes) ~ Treatment * Watershed.Pair, data = FD_16s_phyla1)
shapiro.test(resid(p5))
Anova(p5)
p5.1 <- aggregate(p__Planctomycetes ~ Treatment, data=FD_16s_phyla1, FUN=mean)
p5.1
p5.2 <- aggregate(p__Planctomycetes~Treatment*Watershed.Pair, data=FD_16s_phyla1, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
p5.2
#Reference is ~12% higher
em_p5 <- emmeans(p5, ~ Treatment | Watershed.Pair)
contrast(em_p5, method = "pairwise")

p6 <- aov(log10(p__Actinobacteria) ~ Treatment * Watershed.Pair, data = FD_16s_phyla1)
shapiro.test(resid(p6))
hist(FD_16s_phyla1$p__Actinobacteria)
p6.2 <- aggregate(p__Actinobacteria ~ Treatment, data=FD_16s_phyla1, FUN=mean)
p6.2
p6.3 <- aggregate(p__Actinobacteria~Treatment*Watershed.Pair, data=FD_16s_phyla1, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
p6.3
#Disturbed is ~39% higher
p6.1 <- glm(p__Actinobacteria ~ Treatment*Watershed.Pair,data=FD_16s_phyla1, family=Gamma(link=log))
Anova(p6.1)
em_p6.1 <- emmeans(p6.1, ~ Treatment | Watershed.Pair)
contrast(em_p6.1, method = "pairwise")

p7 <- aov(log10(p__Nitrospirae+1) ~ Treatment * Watershed.Pair, data = FD_16s_phyla1)
shapiro.test(resid(p7))
p7.2 <- aggregate(p__Nitrospirae ~ Treatment, data=FD_16s_phyla1, FUN=mean)
p7.2
p7.3 <- aggregate(p__Nitrospirae~Treatment*Watershed.Pair, data=FD_16s_phyla1, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
p7.3
#disturbed is 280% higher
p7.1 <- glm(p__Nitrospirae+1 ~ Treatment*Watershed.Pair,data=FD_16s_phyla1, family=Gamma(link=log))
Anova(p7)
em_p7 <- emmeans(p7, ~ Treatment | Watershed.Pair)
contrast(em_p7, method = "pairwise")

p8 <- aov(p__Bacteroidetes ~ Treatment * Watershed.Pair, data = FD_16s_phyla1)
shapiro.test(resid(p8))
Anova(p8)

#Copiotroph:Oligotroph analysis

FD_16s_phyla1$OligoCopio <- (FD_16s_phyla1$p__Bacteroidetes+FD_16s_phyla1$p__Proteobacteria)/(FD_16s_phyla1$p__Actinobacteria+FD_16s_phyla1$p__Acidobacteria)

p9 <- aov(OligoCopio ~ Treatment * Watershed.Pair, data = FD_16s_phyla1)
shapiro.test(resid(p9))
Anova(p9)
p9.1 <- aggregate(OligoCopio~Treatment, data=FD_16s_phyla1, FUN=mean)
p9.1
ag9.2 <- aggregate(OligoCopio~Treatment*Watershed.Pair, data=FD_16s_phyla1, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
ag9.2
#Disturbed is 29% higher
emoc <- emmeans(p9, ~ Treatment | Watershed.Pair)
contrast(emoc, method = "pairwise")

#Copiotrophs to Oligotrophs watersheds separated (for figure s3)

jpeg(filename="FigureS3.jpeg", bg="transparent", res=500, units = "in", height=3.75, width=5.5)

box_plotOC <- ggplot(FD_16s_phyla1, aes(x=Watershed.Pair, y=OligoCopio, fill=Treatment)) + 
  geom_boxplot() +
  stat_boxplot(geom = "errorbar", width = 0.25, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE) +
  annotate('text', 2.5, .78, label="Disturbance: italic(P) < 0.001", size=4.5, fontface="bold", parse=TRUE)+
  annotate('text', 2.5, .68, label="`Watershed Pair`: italic(P) < 0.001", size=4.5, fontface="bold", parse=TRUE)+
  annotate('text', 2.5, .58, label="`Disturbance x Pair`: italic(P) == 0.440", size=4.5, fontface="bold", parse=TRUE)+
  annotate('text', 2, 1.8, label="***", size=4.5, fontface="bold")+
  annotate('text', 3, 1.65, label="**", size=4.5, fontface="bold")+
  annotate('text', 4, 1.4, label="*", size=4.5, fontface="bold")+
  annotate('text', 1, 1.3, label="*", size=4.5, fontface="bold")+
  scale_fill_manual(values=c("#D55E00", "#0072B2")) +
  labs(list(x ="")) +
  theme_classic() +
  theme(axis.title=element_text(size=19)) +
  theme(text=element_text(size=19)) +
  theme(legend.position=c(.8, 1)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  theme(legend.title=element_text(colour="white")) +
  scale_y_continuous(bquote(~'Copiotrophs:Oligotrophs')) +
  coord_cartesian(ylim=c(.5,1.8)) 
plot(box_plotOC)

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

box_plotOC + 
  scale_x_discrete(breaks=unique(FD_16s_phyla1$Watershed.Pair), 
                   labels=addline_format(c("Cable Logged", 
                                           "Pasture Conv.", "Pine Conv.",  "Clear Cut")))

dev.off()

#Copiotrophs to oligotrophs aggregated (for figure 2)

box_plotOC2 <- ggplot(FD_16s_phyla1, aes(x=Treatment, y=OligoCopio, fill=Treatment)) + 
  stat_boxplot(geom = "errorbar", size=.25, width = 0.25, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.25, outlier.size=-1) +
  geom_point(position=position_jitter(.25),aes(shape=Watershed.Pair),size=0.9) +
  scale_fill_manual(values=c("#D55E00", "#0072B2")) +
  labs(list(x ="")) +
  theme_classic() +
  theme(axis.title=element_text(size=8)) +
  theme(text=element_text(size=9)) +
  theme(legend.position=c(.2,1)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x=element_line(colour="black", size=.25)) + 
  theme(axis.line.y=element_line(colour="black", size=.25)) + 
  theme(legend.title=element_text(colour="white")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  annotate("text", 1.5, 2, label="Disturbance: ~italic(P) < 0.001", parse=TRUE, size=2.5, fontface="bold")+  
  theme(plot.title = element_text(size=6)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.line.x=element_line(colour="black", size=.25)) + 
  theme(axis.line.y=element_line(colour="black", size=.25)) +
  coord_fixed(ratio=2) +
  theme(legend.position = "none") +
  scale_y_continuous(bquote(~'Copiotrophs:Oligotrophs'))
box_plotOC2

#Correlogram, plus panel figure 2 stuff

bac_alpha <- bac_alpha[order(bac_alpha$Number),]

bac_alpha2 <- bac_alpha[-37,]

bac_alpha2 <- bac_alpha2[order(bac_alpha2$id),]

fung_alpha <- fung_alpha[order(fung_alpha$id),]

FD_16s_phyla1.2 <- FD_16s_phyla1[-31,]

soil2 <- soil2[order(soil2$id),]

soil2$`ITS (H')` <- fung_alpha$Shannon

soil2$`16s (H')` <- bac_alpha2$Shannon

soil2$`Copio:Oligo` <- FD_16s_phyla1.2$OligoCopio

soil2$`ITS:16s` <- soil2$FB

library(corrplot)

Corr <- cor(soil2[,c(6,8,9,10,19,26,27,23,24,25)], use="complete.obs", method="spearman")

res1 <- cor.mtest(soil2[,c(6,8,9,10,19,26,27,23,24,25)], conf.level = .95)

Corr

jpeg(filename="correlogram.jpeg", bg="transparent", res=500, units = "in", height=4, width=4) 

corrplot(Corr, method="circle", type="upper", tl.col="black", cl.ratio=0.2, tl.cex=1,p.mat = res1$p, insig = "blank")

dev.off()

#New Panel Figure 2

legend3 <- get_legend(legend2)

jpeg(filename="Figure2_new.jpeg", bg="transparent", res=500, units = "in", height=4, width=6) 

fig1 <- plot_grid(fbbox2, box_plotOC2, NULL,BacShanBox2, FungShanBox2,NULL , ncol = 3, labels=c('A', 'B', '','C', 'D', 'E'), align="hv", label_size=10)
fig1 + draw_grob(legend3, x = .3, y = .3)

dev.off()

#16s phyla Heatmap stuff

FD_16s_phyla2 <- aggregate(Abundance~Description+Treatment+Phylum, FD_16s_phyla, FUN=mean)

FD_16s_phyla2[FD_16s_phyla2=="p__Acidobacteria"] <- "Acidobacteria***"
FD_16s_phyla2[FD_16s_phyla2=="p__Actinobacteria"] <- "Actinobacteria**"
FD_16s_phyla2[FD_16s_phyla2=="p__Nitrospirae"] <- "Nitrospirae***"
FD_16s_phyla2[FD_16s_phyla2=="p__Bacteroidetes"] <- "Bacteroidetes"
FD_16s_phyla2[FD_16s_phyla2=="p__Proteobacteria"] <- "Proteobacteria***"
FD_16s_phyla2[FD_16s_phyla2=="p__Verrucomicrobia"] <- "Verrucomicrobia"
FD_16s_phyla2[FD_16s_phyla2=="p__Chloroflexi"] <- "Chloroflexi*"
FD_16s_phyla2[FD_16s_phyla2=="p__Planctomycetes"] <- "Planctomycetes**"

head(FD_16s_phyla2)

library(tidyr)

FD_16s_phyla2a <- spread(FD_16s_phyla2, Phylum, Abundance)
FD_16s_phyla2a

FD_16s_phyla3 <- as.matrix(FD_16s_phyla2a[,c(3:7,9:11)])
FD_16s_phyla3

FD_16s_phyla4 <- FD_16s_phyla3[c(4,1,3,5,8,7,2,6),]
FD_16s_phyla4

rownames(FD_16s_phyla4) <- c("CL (Ref)", "Pas (Ref)","Pine (Ref)", "CC (Ref)","CL", "Pas","Pine", "CC")

aka2 = data.frame(ID = factor(rep(c("Reference","Disturbed"), each=4)))
rownames(aka2)<-rownames(FD_16s_phyla4)
aka3 = list(ID = c(Reference = "#0072B2", Disturbed="#D55E00"))

FD_16s_phyla4 <- t(scale(FD_16s_phyla4)) 

my_palette <- colorRampPalette(c("white", "#800080"))(n = 299)

library(pheatmap)

#16s phyla heatmap

heat1 <- pheatmap(FD_16s_phyla4,color=my_palette,cluster_cols = F, cluster_rows = T, annotation_col = aka2, 
                  annotation_colors = aka3[1],annotation_names_col =F,legend=T,annotation_legend = FALSE, border_color = NA, main = "Bacterial Phyla", gaps_col = c(4,8), cex=1)

#ITS Class analysis

colnames(tax_table(FD_ITS)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

FD_ITS2 <- transform_sample_counts(FD_ITS, function(x) x/sum(x))

# agglomerate taxa
glom2 <- tax_glom(FD_ITS2, taxrank = 'Class')

# create dataframe from phyloseq object
dat2 <- psmelt(glom2)

head(dat2)

# convert Class to a character vector from a factor because R
dat2$Class <- as.character(dat2$Class)

library(plyr)

# group dataframe by Class, calculate mean rel. abundance
means2 <- ddply(dat2, ~Class, function(x) c(mean=mean(x$Abundance)))

# find Phyla whose rel. abund. is less than 1%
remainder2 <- means2[means2$mean <= 0.01,]$Class

# change their name to "Other"
dat2[dat2$Class %in% remainder2,]$Class <- 'Other'

FD_ITS_class <- aggregate(Abundance~Sample+Class+Treatment+Watershed.Pair+Description, dat2, FUN=sum)

Coweeta_fungal_classes <- aggregate(Abundance~Class, FD_ITS_class, FUN=mean)

Coweeta_fungal_classes

Coweeta_fungal_classes[Coweeta_fungal_classes=="c__Agaricomycetes"] <- "Agaricomycetes"
Coweeta_fungal_classes[Coweeta_fungal_classes=="c__Archaeorhizomycetes"] <- "Archaeorhizomycetes"
Coweeta_fungal_classes[Coweeta_fungal_classes=="c__Leotiomycetes"] <- "Leotiomycetes"
Coweeta_fungal_classes[Coweeta_fungal_classes=="c__Mortierellomycetes"] <- "Mortierellomycetes"
Coweeta_fungal_classes[Coweeta_fungal_classes=="c__Sordariomycetes"] <- "Sordariomycetes"
Coweeta_fungal_classes[Coweeta_fungal_classes=="c__Tremellomycetes"] <- "Tremellomycetes"
Coweeta_fungal_classes[Coweeta_fungal_classes=="c__Geminibasidiomycetes"] <- "Geminibasidiomycetes"
Coweeta_fungal_classes[Coweeta_fungal_classes=="c__Mucoromycotina_cls_Incertae_sedis"] <- "Mucoromycotina incertae sedis"
Coweeta_fungal_classes[Coweeta_fungal_classes=="c__Eurotiomycetes"] <- "Eurotiomycetes"

Coweeta_fungal_classes$Abundance <- Coweeta_fungal_classes$Abundance * 100

is.num <- sapply(Coweeta_fungal_classes, is.numeric)
Coweeta_fungal_classes[is.num] <- lapply(Coweeta_fungal_classes[is.num], round, 1)

Coweeta_fungal_classes <- Coweeta_fungal_classes[order(Coweeta_fungal_classes$Abundance),]

cbPalette2 <- c("#D55E00","black","#0072B2","#999999","#CC79A7","#E69F00","#F0E442",  "#009E73",    "purple", "#56B4E9")

#Fungal Classes pie chart (for figure 3)

jpeg(filename="ITSpie.jpeg", bg="transparent", res=500, units = "in", height=3, width=4) 

with(Coweeta_fungal_classes,pie(Abundance, labels=paste0(as.character(Class), " ", Abundance, "%"), cex = .4, col=cbPalette2, radius=1))

dev.off()

library(reshape)

FD_ITS_class2 <- cast(FD_ITS_class, Sample+Treatment+Watershed.Pair ~ Class, value="Abundance")
FD_ITS_class2

p11 <- aov(c__Agaricomycetes ~ Treatment * Watershed.Pair, data = FD_ITS_class2)
shapiro.test(resid(p11))
Anova(p11)
p11.1 <- aggregate(c__Agaricomycetes ~ Treatment, data = FD_ITS_class2, FUN=mean)
p11.1
#Reference is 17% higher
em_p11 <- emmeans(p11, ~ Treatment | Watershed.Pair)
contrast(em_p11, method = "pairwise")

p12 <- aov(log10(c__Archaeorhizomycetes) ~ Treatment * Watershed.Pair, data = FD_ITS_class2)
shapiro.test(resid(p12))
Anova(p12)

p13 <- aov(c__Leotiomycetes ~ Treatment * Watershed.Pair, data = FD_ITS_class2)
shapiro.test(resid(p13))
Anova(p13)

p14 <- aov(log10(c__Sordariomycetes) ~ Treatment * Watershed.Pair, data = FD_ITS_class2)
shapiro.test(resid(p14))
Anova(p14)
p14.2 <- aggregate(c__Sordariomycetes ~ Treatment, data=FD_ITS_class2, FUN=mean)
p14.2
#Disturbed is 67% higher
em_p14 <- emmeans(p14, ~ Treatment | Watershed.Pair)
contrast(em_p14, method = "pairwise")

p15 <- aov(log10(c__Tremellomycetes) ~ Treatment * Watershed.Pair, data = FD_ITS_class2)
shapiro.test(resid(p15))
hist(FD_ITS_class2$c__Tremellomycetes)
p15.1 <- glm(c__Tremellomycetes ~ Treatment * Watershed.Pair, data = FD_ITS_class2, family=Gamma(link=log))
Anova(p15.1)

p16 <- aov(c__Geminibasidiomycetes ~ Treatment * Watershed.Pair, data = FD_ITS_class2)
shapiro.test(resid(p16))
p16.2 <- aggregate(c__Geminibasidiomycetes ~ Treatment, data=FD_ITS_class2, FUN=mean)
p16.2
#Reference is 73% higher
p16.1 <- glm((c__Geminibasidiomycetes+1) ~ Treatment * Watershed.Pair, data = FD_ITS_class2, family=Gamma(link=log))
Anova(p16.1)
em_p16.1 <- emmeans(p16, ~ Treatment | Watershed.Pair)
contrast(em_p16.1, method = "pairwise")

p17 <- aov(c__Mucoromycotina_cls_Incertae_sedis ~ Treatment * Watershed.Pair, data = FD_ITS_class2)
shapiro.test(resid(p17))
hist(FD_ITS_class2$c__Mucoromycotina_cls_Incertae_sedis)
p17.2 <- aggregate(c__Mucoromycotina_cls_Incertae_sedis ~ Treatment, data=FD_ITS_class2, FUN=mean)
p17.2
#Reference is 85% higher
p17.1 <- glm((c__Mucoromycotina_cls_Incertae_sedis+1) ~ Treatment * Watershed.Pair, data = FD_ITS_class2, family=Gamma(link=log))
Anova(p17.1)
em_p17.1 <- emmeans(p17, ~ Treatment | Watershed.Pair)
contrast(em_p17.1, method = "pairwise")

p18 <- aov(log10(c__Mortierellomycetes) ~ Treatment * Watershed.Pair, data = FD_ITS_class2)
shapiro.test(resid(p18))
Anova(p18)

p19 <- aov(log10(c__Eurotiomycetes+1) ~ Treatment * Watershed.Pair, data = FD_ITS_class2)
shapiro.test(resid(p19))
hist(FD_ITS_class2$c__Eurotiomycetes)
p19.2 <- aggregate(c__Eurotiomycetes ~ Treatment, data=FD_ITS_class2, FUN=mean)
p19.2
#Disturbed is 150% higher
p19.1 <- glm((c__Eurotiomycetes+1) ~ Treatment * Watershed.Pair, data = FD_ITS_class2, family=Gamma(link=log))
Anova(p19.1)
em_p19.1 <- emmeans(p19, ~ Treatment | Watershed.Pair)
contrast(em_p19.1, method = "pairwise")

#ITS classes heatmap

FD_ITS_class3 <- aggregate(Abundance~Description+Treatment+Class, FD_ITS_class, FUN=mean)

head(FD_ITS_class3)

FD_ITS_class3[FD_ITS_class3=="c__Agaricomycetes"] <- "Agaricomycetes"
FD_ITS_class3[FD_ITS_class3=="c__Archaeorhizomycetes"] <- "Archaeorhizomycetes"
FD_ITS_class3[FD_ITS_class3=="c__Leotiomycetes"] <- "Leotiomycetes"
FD_ITS_class3[FD_ITS_class3=="c__Mortierellomycetes"] <- "Mortierellomycetes"
FD_ITS_class3[FD_ITS_class3=="c__Sordariomycetes"] <- "Sordariomycetes*"
FD_ITS_class3[FD_ITS_class3=="c__Tremellomycetes"] <- "Tremellomycetes"
FD_ITS_class3[FD_ITS_class3=="c__Geminibasidiomycetes"] <- "Geminibasidiomycetes*"
FD_ITS_class3[FD_ITS_class3=="c__Mucoromycotina_cls_Incertae_sedis"] <- "Mucoromycotina incertae sedis*"
FD_ITS_class3[FD_ITS_class3=="c__Eurotiomycetes"] <- "Eurotiomycetes*"

head(FD_ITS_class3)

library(tidyr)

FD_ITS_class3 <- spread(FD_ITS_class3, Class, Abundance)
FD_ITS_class3

FD_ITS_class4 <- as.matrix(FD_ITS_class3[,c(3:9,11:12)])
FD_ITS_class4 <- FD_ITS_class4[c(4,1,3,5,8,7,2,6),]

rownames(FD_ITS_class4) <- c("CL (Ref)", "Pas (Ref)","Pine (Ref)", "CC (Ref)","CL", "Pas","Pine", "CC")

aka2 = data.frame(ID = factor(rep(c("Reference","Disturbed"), each=4)))
rownames(aka2)<-rownames(FD_ITS_class4)
aka3 = list(ID = c(Reference = "#0072B2", Disturbed="#D55E00"))

FD_ITS_class4 <- t(scale(FD_ITS_class4)) 

my_palette <- colorRampPalette(c("white", "#800080"))(n = 299)

library(pheatmap)

#ITS classes heatmap

heat2 <- pheatmap(FD_ITS_class4,color=my_palette,cluster_cols = F, cluster_rows = T, annotation_col = aka2, 
                  annotation_colors = aka3[1],annotation_names_col =F,legend=T,annotation_legend = FALSE, border_color = NA, main = "Fungal Classes", gaps_col = c(4,8), cex=1)

#Funguild analysis

funguild3 <- cbind(funguild, c2)

head(funguild3)

library(tidyr)

funguild4 <- gather(funguild3, Guild, Abundance, `Animal Pathogen`:`Wood Saprotroph`, factor_key=TRUE)
head(funguild4)

FD_funguild <- aggregate(Abundance~id+Guild+Treatment+Watershed.Pair+Description, funguild4, FUN=sum)

Coweeta_funguild <- aggregate(Abundance~Guild, FD_funguild, FUN=mean)

Coweeta_funguild$Abundance <- Coweeta_funguild$Abundance * 100

is.num <- sapply(Coweeta_funguild, is.numeric)
Coweeta_funguild[is.num] <- lapply(Coweeta_funguild[is.num], round, 1)

Coweeta_funguild <- Coweeta_funguild[order(Coweeta_funguild$Abundance),]

Coweeta_funguild

Coweeta_funguild2 <- Coweeta_funguild[-c(1:9),]

other <- data.frame("Other", 0.9)

names(other)<-c("Guild","Abundance")

Coweeta_funguild3 <- rbind(Coweeta_funguild2, other)

Coweeta_funguild3

cbPalette3 <- c("#F0E442","#800080", "#D55E00", "#0072B2", "#999999","#CC79A7","#E69F00",  "#009E73","#56B4E9")

#Funguild pie chart

jpeg(filename="Guildpie.jpeg", bg="transparent", res=500, units = "in", height=3, width=4) 

with(Coweeta_funguild3,pie(Abundance, labels=paste0(as.character(Guild), " ", Abundance, "%"), cex = .4, col=cbPalette3, radius=1))

dev.off()

FD_funguild2 <- cast(FD_funguild, id+Treatment+Watershed.Pair ~ Guild, value="Abundance")
FD_funguild2

p21 <- aov(Ectomycorrhizal ~ Treatment * Watershed.Pair, data = FD_funguild2)
shapiro.test(resid(p21))
Anova(p21)
p21.1 <- aggregate(Ectomycorrhizal ~ Treatment, data=FD_funguild2, FUN=mean)
p21.1
#Reference is 52% higher
em_p21 <- emmeans(p21, ~ Treatment | Watershed.Pair)
contrast(em_p21, method = "pairwise")

p22 <- aov(log10(`Soil Saprotroph`) ~ Treatment * Watershed.Pair, data = FD_funguild2)
shapiro.test(resid(p22))
Anova(p22)

p23 <- aov(log10(`Undefined Saprotroph`) ~ Treatment * Watershed.Pair, data = FD_funguild2)
shapiro.test(resid(p23))
Anova(p23)

p24 <- aov(log10(`Arbuscular Mycorrhizal`+1) ~ Treatment * Watershed.Pair, data = FD_funguild2)
shapiro.test(resid(p24))
p24.1 <- glm((`Arbuscular Mycorrhizal`+1) ~ Treatment * Watershed.Pair,data=FD_funguild2, family=Gamma(link=log))
Anova(p24)
p24.1 <- aggregate(`Arbuscular Mycorrhizal` ~ Treatment, data=FD_funguild2, FUN=mean)
p24.1
#Disturbed is 83% higher
em_p24.1 <- emmeans(p24, ~ Treatment | Watershed.Pair)
contrast(em_p24.1, method = "pairwise")

p25 <- aov(log10(`Wood Saprotroph`+1) ~ Treatment * Watershed.Pair, data = FD_funguild2)
shapiro.test(resid(p25))
p25.1 <- glm((`Wood Saprotroph`+1) ~ Treatment * Watershed.Pair,data=FD_funguild2, family=Gamma(link=log))
Anova(p25.1)

p26 <- aov(log10(`Ericoid Mycorrhizal`+1) ~ Treatment * Watershed.Pair, data = FD_funguild2)
shapiro.test(resid(p26))
p26.1 <- glm((`Ericoid Mycorrhizal`+1) ~ Treatment * Watershed.Pair,data=FD_funguild2, family=Gamma(link=log))
Anova(p26.1)

p27 <- aov(log10(`Orchid Mycorrhizal`+1) ~ Treatment * Watershed.Pair, data = FD_funguild2)
shapiro.test(resid(p27))
p27.1 <- glm((`Orchid Mycorrhizal`+1) ~ Treatment * Watershed.Pair,data=FD_funguild2, family=Gamma(link=log))
Anova(p27.1)

p28 <- aov(log10(`Dung Saprotroph`+1) ~ Treatment * Watershed.Pair, data = FD_funguild2)
shapiro.test(resid(p28))
p28.1 <- glm((`Dung Saprotroph`+1) ~ Treatment * Watershed.Pair,data=FD_funguild2, family=Gamma(link=log))
Anova(p28.1)

p29 <- aov(`Animal Pathogen` ~ Treatment * Watershed.Pair, data = FD_funguild2)
shapiro.test(resid(p29))
Anova(p29)
p29.1 <- glm((`Animal Pathogen`+1) ~ Treatment * Watershed.Pair,data=FD_funguild2, family=Gamma(link=log))
Anova(p29.1)
p29.1 <- aggregate(`Animal Pathogen` ~ Treatment, data=FD_funguild2, FUN=mean)
p29.1
#Disturbed is 3,620% higher . . . 
em_p29.1 <- emmeans(p29, ~ Treatment | Watershed.Pair)
contrast(em_p29.1, method = "pairwise")

p30 <- aov(`Fungal Parasite` ~ Treatment * Watershed.Pair, data = FD_funguild2)
shapiro.test(resid(p30))
Anova(p30)
p30.1 <- glm((`Fungal Parasite`+1) ~ Treatment * Watershed.Pair,data=FD_funguild2, family=Gamma(link=log))
Anova(p30.1)

p31 <- aov(Lichenized ~ Treatment * Watershed.Pair, data = FD_funguild2)
shapiro.test(resid(p31))
Anova(p31)
p31.1 <- glm((Lichenized+1) ~ Treatment * Watershed.Pair,data=FD_funguild2, family=Gamma(link=log))
Anova(p31.1)

p32 <- aov(`Plant Pathogen` ~ Treatment * Watershed.Pair, data = FD_funguild2)
shapiro.test(resid(p32))
Anova(p32)
p32.1 <- glm((`Plant Pathogen`+1) ~ Treatment * Watershed.Pair,data=FD_funguild2, family=Gamma(link=log))
Anova(p32.1)

p33 <- aov(Endophyte ~ Treatment * Watershed.Pair, data = FD_funguild2)
shapiro.test(resid(p33))
Anova(p33)
p33.1 <- glm((Endophyte+1) ~ Treatment * Watershed.Pair,data=FD_funguild2, family=Gamma(link=log))
Anova(p33.1)
p33.2 <- aggregate(Endophyte ~ Treatment, data=FD_funguild2, FUN=mean)
p33.2
#Disturbed is 150% higher
em_p33.1 <- emmeans(p33, ~ Treatment | Watershed.Pair)
contrast(em_p33.1, method = "pairwise")

#Funguild heatmap

FD_funguild3 <- aggregate(Abundance~Description+Treatment+Guild, FD_funguild, FUN=mean)

head(FD_funguild3)

FD_funguild3$Guild <- as.character(FD_funguild3$Guild)

FD_funguild3[FD_funguild3=="Ectomycorrhizal"] <- "Ectomycorrhizal*"
FD_funguild3[FD_funguild3=="Arbuscular Mycorrhizal"] <- "Arbuscular Mycorrhizal*"
FD_funguild3[FD_funguild3=="Endophyte"] <- "Endophyte*"
FD_funguild3[FD_funguild3=="Animal Pathogen"] <- "Animal Pathogen*"

head(FD_funguild3)

library(tidyr)

FD_funguild3 <- spread(FD_funguild3, Guild, Abundance)
FD_funguild3

FD_funguild4 <- as.matrix(FD_funguild3[,c(3:10,13:15,17:18),])
FD_funguild4 <- FD_funguild4[c(4,1,3,5,8,7,2,6),]

rownames(FD_funguild4) <- c("CL (Ref)", "Pas (Ref)","Pine (Ref)", "CC (Ref)","CL", "Pas","Pine", "CC")

aka2 = data.frame(ID = factor(rep(c("Reference","Disturbed"), each=4)))
rownames(aka2)<-rownames(FD_funguild4)
aka3 = list(ID = c(Reference = "#0072B2", Disturbed="#D55E00"))

FD_funguild4 <- t(scale(FD_funguild4)) 

my_palette <- colorRampPalette(c("white", "#800080"))(n = 299)

library(pheatmap)

#Funguild heatmap

heat3 <- pheatmap(FD_funguild4,color=my_palette,cluster_cols = F, cluster_rows = T, annotation_col = aka2, 
                  annotation_colors = aka3[1],annotation_names_col =F,legend=T,annotation_legend = FALSE, border_color = NA, main = "Fungal Guilds", gaps_col = c(4,8), cex=1)

#Make panel figure 4

jpeg(filename="Figure4.jpeg", bg="transparent", res=500, units = "in", height=10, width=8.5) 

fig4 <- plot_grid(NULL, heat1[[4]],NULL, heat2[[4]], NULL, heat3[[4]], ncol = 2, axis="tblr",labels=c('A', 'B', 'C', 'D', 'E', 'F'), align="hv", label_size=18)
fig4

dev.off()

#16s Differentially Abundant OTUs

library(edgeR)

counts1 <- data.frame(otu_table(FD_16s))

e1 <- DGEList(counts=counts1, group = d1$Treatment)

e1 <- calcNormFactors(e1)
design1 <- model.matrix(~d1$Treatment)
e1 <- estimateDisp(e1,design1)

et <- exactTest(e1, pair=c("Reference","Disturbed"))

t1.1 <- as.data.frame(et$table)

pvalue.1 <- subset(t1.1, PValue<=0.05)

sigtab_bac1 <- psmelt(FD_16s)

sigtab_bac <- sigtab_bac1[sigtab_bac1$OTU %in% rownames(pvalue.1),]

#16s edgeR results

sigtab_bac2 <- aggregate(Abundance~Phylum, sigtab_bac, FUN=sum)

sigtab_bac2$Abundance <- sigtab_bac2$Abundance / sum(sigtab_bac2$Abundance)

sigtab_bac2$Phylum <- as.character(sigtab_bac2$Phylum)

remainder3 <- sigtab_bac2[c(sigtab_bac2$Abundance <= 0.01),]$Phylum

# change their name to "Other"
sigtab_bac2[sigtab_bac2$Phylum %in% remainder3,]$Phylum <- 'Other'

sigtab_bac3 <- aggregate(Abundance~Phylum, sigtab_bac2, FUN=sum)

sigtab_bac3[sigtab_bac3=="p__Acidobacteria"] <- "Acidobacteria"
sigtab_bac3[sigtab_bac3=="p__Actinobacteria"] <- "Actinobacteria"
sigtab_bac3[sigtab_bac3=="p__Nitrospirae"] <- "Nitrospirae"
sigtab_bac3[sigtab_bac3=="p__Bacteroidetes"] <- "Bacteroidetes"
sigtab_bac3[sigtab_bac3=="p__Proteobacteria"] <- "Proteobacteria"
sigtab_bac3[sigtab_bac3=="p__Verrucomicrobia"] <- "Verrucomicrobia"
sigtab_bac3[sigtab_bac3=="p__Chloroflexi"] <- "Chloroflexi"
sigtab_bac3[sigtab_bac3=="p__Planctomycetes"] <- "Planctomycetes"

sigtab_bac3$Abundance <- sigtab_bac3$Abundance * 100

is.num <- sapply(sigtab_bac3, is.numeric)
sigtab_bac3[is.num] <- lapply(sigtab_bac3[is.num], round, 1)

cbPalette3 <- c("#CC79A7","#E69F00",  "#009E73", "#F0E442","black", "#0072B2", "#D55E00", "#56B4E9","#999999")

sigtab_bac3$Phylum <- factor(sigtab_bac3$Phylum, levels = c("Proteobacteria", "Acidobacteria","Nitrospirae", "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Planctomycetes", "Other", "Verrucomicrobia"))
sigtab_bac3

#16s Phyla Pie Chart (differentially abundant otus)

jpeg(filename="16spie2.jpeg", bg="transparent", res=500, units = "in", height=4, width=6.5) 

with(sigtab_bac3,pie(Abundance, labels=paste0(as.character(Phylum), " ", Abundance, "%"), cex = .6, col=rev(cbPalette3), radius=1))

dev.off()

#Stats on phyla of differentially abundant 16s OTUs

bac_diffabun <- psmelt(FD_16s)

bac_diffabun2 <- bac_diffabun[bac_diffabun$OTU %in% rownames(pvalue.1),]

bac_diffabun2.1 <- aggregate(Abundance ~ Treatment + Watershed.Pair + Sample + Phylum, data=bac_diffabun2, FUN=sum)

bac_diffabun2.1$Abundance <- bac_diffabun2.1$Abundance / sum(bac_diffabun2.1$Abundance)

bac_diffabun2.1$Phylum <- as.character(bac_diffabun2.1$Phylum)

# group dataframe by Phylum, calculate mean rel. abundance
means3 <-aggregate(Abundance~ Phylum, data=bac_diffabun2.1, FUN=sum)

# find Phyla whose rel. abund. is less than 1%
remainder3 <- means3[means3$Abundance <= 0.01,]$Phylum

bac_diffabun2$Phylum <- as.character(bac_diffabun2$Phylum)

bac_diffabun2 <- aggregate(Abundance ~ Treatment + Watershed.Pair + Sample + Phylum, data=bac_diffabun2, FUN=sum)

bac_diffabun2[bac_diffabun2$Phylum %in% remainder3,]$Phylum <- 'Other'

bac_diffabun2 <-aggregate(Abundance~ Phylum + Treatment + Sample + Watershed.Pair, data=bac_diffabun2, FUN=sum)

library(reshape)

bac_diffabun3 <- cast(bac_diffabun2, Sample+Treatment+Watershed.Pair ~ Phylum, value="Abundance")

relab1 <- as.matrix(bac_diffabun3[,c(4:12)])

relab2 <- apply(relab1, 1, function(i) i/sum(i))

relab3 <- as.data.frame(t(relab2))

bac_diffabun4 <- cbind(bac_diffabun3[,c(1:3)],relab3)

colnames(bac_diffabun4) <- colnames(bac_diffabun3)

da1 <- aov(p__Acidobacteria ~ Treatment * Watershed.Pair, data=bac_diffabun4)
shapiro.test(resid(da1))
Anova(da1)
da1.1 <- aggregate(p__Acidobacteria~Treatment*Watershed.Pair, data=bac_diffabun4, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
da1.1
em_da1 <- emmeans(da1, ~ Treatment | Watershed.Pair)
contrast(em_da1, method = "pairwise")

da2 <- aov(p__Proteobacteria ~ Treatment * Watershed.Pair, data=bac_diffabun4)
shapiro.test(resid(da2))
Anova(da2)
da2.1 <- aggregate(p__Proteobacteria~Treatment*Watershed.Pair, data=bac_diffabun4, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
da2.1
da2.2 <- glm(p__Proteobacteria~ Treatment*Watershed.Pair, data=bac_diffabun3, family=Gamma(link=log))
Anova(da2.2)
em_da2 <- emmeans(da2, ~ Treatment | Watershed.Pair)
contrast(em_da2, method = "pairwise")

da3 <- aov(log10(p__Actinobacteria+1) ~ Treatment * Watershed.Pair, data=bac_diffabun4)
shapiro.test(resid(da3))
da3.1 <- glm(p__Actinobacteria+1 ~ Treatment * Watershed.Pair, data=bac_diffabun4, family=Gamma(link=log))
Anova(da3)
da3.2 <- aggregate(p__Actinobacteria~Treatment*Watershed.Pair, data=bac_diffabun3, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
da3.2
em_da3 <- emmeans(da3, ~ Treatment | Watershed.Pair)
contrast(em_da3, method = "pairwise")

da4 <- aov(log10(p__Nitrospirae+1) ~ Treatment * Watershed.Pair, data=bac_diffabun4)
shapiro.test(resid(da4))
da4.1 <- glm(p__Nitrospirae+1 ~ Treatment * Watershed.Pair, data=bac_diffabun4, family=Gamma(link=log))
Anova(da4)
da4.2 <- aggregate(p__Nitrospirae~Treatment*Watershed.Pair, data=bac_diffabun4, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
da4.2
em_da4 <- emmeans(da4, ~ Treatment | Watershed.Pair)
contrast(em_da4, method = "pairwise")

da5 <- aov(p__Planctomycetes ~ Treatment * Watershed.Pair, data=bac_diffabun4)
shapiro.test(resid(da5))
Anova(da5)
da5.1 <- aggregate(p__Planctomycetes~Treatment*Watershed.Pair, data=bac_diffabun4, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
da5.1
da5.2 <- glm(p__Planctomycetes ~ Treatment * Watershed.Pair, data=bac_diffabun4, family=Gamma(link=log))
Anova(da5.2)
em_da5 <- emmeans(da5, ~ Treatment | Watershed.Pair)
contrast(em_da5, method = "pairwise")

da6 <- aov(p__Chloroflexi ~ Treatment * Watershed.Pair, data=bac_diffabun4)
shapiro.test(resid(da6))
Anova(da6)
da6.1 <- aggregate(p__Chloroflexi~Treatment*Watershed.Pair, data=bac_diffabun4, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
da6.1
em_da6 <- emmeans(da6, ~ Treatment | Watershed.Pair)
contrast(em_da6, method = "pairwise")

da7 <- aov(log10(p__Bacteroidetes+1) ~ Treatment * Watershed.Pair, data=bac_diffabun4)
shapiro.test(resid(da7))
da7.1 <- glm(p__Bacteroidetes+1 ~ Treatment * Watershed.Pair, data=bac_diffabun4, family=Gamma(link=log))
Anova(da7)
da7.1 <- aggregate(p__Bacteroidetes~Treatment*Watershed.Pair, data=bac_diffabun4, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
da7.1
em_da7 <- emmeans(da7, ~ Treatment | Watershed.Pair)
contrast(em_da7, method = "pairwise")

da8 <- aov(log10(p__Verrucomicrobia) ~ Treatment * Watershed.Pair, data=bac_diffabun4)
shapiro.test(resid(da8))
Anova(da8)
da8.1 <- aggregate(p__Verrucomicrobia~Treatment*Watershed.Pair, data=bac_diffabun4, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
da8.1
em_da8 <- emmeans(da8, ~ Treatment | Watershed.Pair)
contrast(em_da8, method = "pairwise")

#Bac diff abun heatmap

bac_diffabun5 <- aggregate(Abundance~Watershed.Pair+Treatment+Phylum, bac_diffabun2, FUN=mean)

bac_diffabun5[bac_diffabun5=="p__Acidobacteria"] <- "Acidobacteria***"
bac_diffabun5[bac_diffabun5=="p__Actinobacteria"] <- "Actinobacteria***"
bac_diffabun5[bac_diffabun5=="p__Nitrospirae"] <- "Nitrospirae***"
bac_diffabun5[bac_diffabun5=="p__Bacteroidetes"] <- "Bacteroidetes***"
bac_diffabun5[bac_diffabun5=="p__Proteobacteria"] <- "Proteobacteria***"
bac_diffabun5[bac_diffabun5=="p__Verrucomicrobia"] <- "Verrucomicrobia***"
bac_diffabun5[bac_diffabun5=="p__Chloroflexi"] <- "Chloroflexi***"
bac_diffabun5[bac_diffabun5=="p__Planctomycetes"] <- "Planctomycetes**"

head(bac_diffabun5)

library(tidyr)

bac_diffabun5a <- spread(bac_diffabun5, Phylum, Abundance)
bac_diffabun5a

bac_diffabun6 <- as.matrix(bac_diffabun5a[,c(3:7,9:11)])
bac_diffabun6

bac_diffabun7 <- bac_diffabun6[c(2,6,8,4,1,5,7,3),]
bac_diffabun7

rownames(bac_diffabun7) <- c("CL (Ref)", "Pas (Ref)","Pine (Ref)", "CC (Ref)","CL", "Pas","Pine", "CC")

aka2 = data.frame(ID = factor(rep(c("Reference","Disturbed"), each=4)))
rownames(aka2)<-rownames(bac_diffabun7)
aka3 = list(ID = c(Reference = "#0072B2", Disturbed="#D55E00"))

bac_diffabun7 <- t(scale(bac_diffabun7)) 

my_palette <- colorRampPalette(c("white", "#800080"))(n = 299)

library(pheatmap)

#16s phyla of differentially abundant OTUs heatmap

heat3.1 <- pheatmap(bac_diffabun7,color=my_palette,cluster_cols = F, cluster_rows = T, annotation_col = aka2, 
                    annotation_colors = aka3[1],annotation_names_col =F,legend=F,annotation_legend = FALSE, border_color = NA, main = "", gaps_col = c(4,8), cex=1)

#ITS Differentially Abundant OTUs

counts2 <- data.frame(otu_table(FD_ITS))

e2 <- DGEList(counts=counts2, group = d2$Treatment)

e2 <- calcNormFactors(e2)
design2 <- model.matrix(~d2$Treatment)
e2 <- estimateDisp(e2,design2)

et2 <- exactTest(e2, pair=c("Reference","Disturbed"))

t2.1 <- as.data.frame(et2$table)

pvalue2.1 <- subset(t2.1, PValue<=0.05)

FD_ITS.1 <- transform_sample_counts(FD_ITS, function(x) x/sum(x))

colnames(tax_table(FD_ITS.1)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

sigtab_its1 <- psmelt(FD_ITS.1)

sigtab_its <- sigtab_its1[sigtab_its1$OTU %in% rownames(pvalue2.1),]

#ITS edgeR results

sigtab_its2 <- aggregate(Abundance~Class, sigtab_its, FUN=sum)

sigtab_its2$Abundance <- sigtab_its2$Abundance / sum(sigtab_its2$Abundance)

sigtab_its2$Class <- as.character(sigtab_its2$Class)

remainder4 <- sigtab_its2[c(sigtab_its2$Abundance <= 0.01),]$Class

# change their name to "Other"
sigtab_its2[sigtab_its2$Class %in% remainder4,]$Class <- 'Other'

sigtab_its3 <- aggregate(Abundance~Class, sigtab_its2, FUN=sum)

sigtab_its3[sigtab_its3=="c__Agaricomycetes"] <- "Agaricomycetes"
sigtab_its3[sigtab_its3=="c__Archaeorhizomycetes"] <- "Archaeorhizomycetes"
sigtab_its3[sigtab_its3=="c__Eurotiomycetes"] <- "Eurotiomycetes"
sigtab_its3[sigtab_its3=="c__Sordariomycetes"] <- "Sordariomycetes"

sigtab_its3$Abundance <- sigtab_its3$Abundance * 100

is.num <- sapply(sigtab_its3, is.numeric)
sigtab_its3[is.num] <- lapply(sigtab_its3[is.num], round, 1)

cbPalette4 <- c("#CC79A7","#E69F00",  "#009E73","#999999","#D55E00",  "#F0E442", "#56B4E9")

#ITS Phyla Pie Chart (differentially abundant otus)

jpeg(filename="ITSpie2.jpeg", bg="transparent", res=500, units = "in", height=4, width=6.5) 

with(sigtab_its3,pie(Abundance, labels=paste0(as.character(Class), " ", Abundance, "%"), cex = .6, col=rev(cbPalette4), radius=1))

dev.off()

#Stats on differentially abundant ITS OTUs

its_diffabun <- psmelt(FD_ITS)

its_diffabun2 <- its_diffabun[its_diffabun$OTU %in% rownames(pvalue2.1),]

its_diffabun2.1 <- aggregate(Abundance ~ Treatment + Watershed.Pair + Sample + Class, data=its_diffabun2, FUN=sum)

its_diffabun2.1$Abundance <- its_diffabun2.1$Abundance / sum(its_diffabun2.1$Abundance)

its_diffabun2.1$Class <- as.character(its_diffabun2.1$Class)

remainder5 <- its_diffabun2.1[its_diffabun2.1$Class == c("c__Geoglossomycetes","c__Microbotryomycetes","c__Orbiliomycetes","c__Pezizomycetes","c__Umbelopsidomycetes","c__unidentified"),]$Class

# change their name to "Other"

its_diffabun2$Class <- as.character(its_diffabun2$Class)

its_diffabun2[its_diffabun2$Class %in% remainder5,]$Class <- 'Other'

its_diffabun2 <-aggregate(Abundance~ Class + Treatment + Sample + Watershed.Pair, data=its_diffabun2, FUN=sum)

library(reshape)

its_diffabun2.1 <- cast(its_diffabun2, Sample+Treatment+Watershed.Pair ~ Class, value="Abundance")

relab4 <- as.matrix(its_diffabun2.1[,c(4:10)])

relab5 <- apply(relab4, 1, function(i) i/sum(i))

relab6 <- as.data.frame(t(relab5))

its_diffabun3 <- cbind(its_diffabun2.1[,c(1:3)],relab6)

colnames(its_diffabun3) <- colnames(its_diffabun2.1)

da9 <- aov(c__Agaricomycetes ~ Treatment * Watershed.Pair, data=its_diffabun3)
shapiro.test(resid(da9))
Anova(da9)
da9.1 <- glm(c__Agaricomycetes ~ Treatment * Watershed.Pair, data=its_diffabun3, family=Gamma(link=log))
Anova(da9.1)
da9.2 <- aggregate(c__Agaricomycetes~Treatment*Watershed.Pair, data=its_diffabun3, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
da9.2
em_da9 <- emmeans(da9, ~ Treatment | Watershed.Pair)
contrast(em_da9, method = "pairwise")

da10 <- aov(log10(c__Archaeorhizomycetes+1) ~ Treatment * Watershed.Pair, data=its_diffabun3)
shapiro.test(resid(da10))
da10.1 <- glm(c__Archaeorhizomycetes+1 ~ Treatment * Watershed.Pair, data= its_diffabun3, family=Gamma(link=log))
Anova(da10)

da11 <- aov(log10(c__Sordariomycetes+1) ~ Treatment * Watershed.Pair, data=its_diffabun3)
shapiro.test(resid(da11))
da11.1 <- glm(c__Sordariomycetes+1 ~ Treatment * Watershed.Pair, data= its_diffabun3, family=Gamma(link=log))
Anova(da11)
da11.2 <- aggregate(c__Sordariomycetes~Treatment*Watershed.Pair, data=its_diffabun3, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
da11.2
em_da11 <- emmeans(da11, ~ Treatment | Watershed.Pair)
contrast(em_da11, method = "pairwise")

da12 <- aov(log10(c__Eurotiomycetes+1) ~ Treatment * Watershed.Pair, data=its_diffabun3)
shapiro.test(resid(da12))
da12.1 <- glm(c__Eurotiomycetes+1 ~ Treatment * Watershed.Pair, data= its_diffabun3, family=Gamma(link=log))
Anova(da12.1)
da12.2 <- aggregate(c__Eurotiomycetes~Treatment*Watershed.Pair, data=its_diffabun3, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
da12.2
em_da12 <- emmeans(da12, ~ Treatment | Watershed.Pair)
contrast(em_da12, method = "pairwise")

da13 <- aov(log10(c__Mortierellomycetes+1) ~ Treatment * Watershed.Pair, data=its_diffabun3)
shapiro.test(resid(da13))
da13.1 <- glm(c__Mortierellomycetes+1 ~ Treatment * Watershed.Pair, data= its_diffabun3, family=Gamma(link=log))
Anova(da13)
da13.2 <- aggregate(c__Mortierellomycetes~Treatment*Watershed.Pair, data=its_diffabun3, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
da13.2
em_da13 <- emmeans(da13, ~ Treatment | Watershed.Pair)
contrast(em_da13, method = "pairwise")

da14 <- aov(log10(c__Leotiomycetes+1) ~ Treatment * Watershed.Pair, data=its_diffabun3)
shapiro.test(resid(da14))
da14.1 <- glm(c__Leotiomycetes+1 ~ Treatment * Watershed.Pair, data= its_diffabun3, family=Gamma(link=log))
Anova(da14.1)

#Heat map for phyla of differentially abundant ITS OTUs

its_diffabun5 <- aggregate(Abundance~Watershed.Pair+Treatment+Class, its_diffabun2, FUN=mean)

its_diffabun5[its_diffabun5 =="c__Agaricomycetes"] <- "Agaricomycetes*"
its_diffabun5[its_diffabun5 =="c__Archaeorhizomycetes"] <- "Archaeorhizomycetes"
its_diffabun5[its_diffabun5 =="c__Eurotiomycetes"] <- "Eurotiomycetes***"
its_diffabun5[its_diffabun5 =="c__Sordariomycetes"] <- "Sordariomycetes**"
its_diffabun5[its_diffabun5 =="c__Leotiomycetes"] <- "Leotiomycetes"
its_diffabun5[its_diffabun5 =="c__Mortierellomycetes"] <- "Mortierellomycetes**"

head(bac_diffabun5)

library(tidyr)

its_diffabun5a <- spread(its_diffabun5, Class, Abundance)
its_diffabun5a

its_diffabun6 <- as.matrix(its_diffabun5a[,c(3:7,9)])
its_diffabun6

its_diffabun7 <- its_diffabun6[c(2,6,8,4,1,5,7,3),]
its_diffabun7

rownames(its_diffabun7) <- c("CL (Ref)", "Pas (Ref)","Pine (Ref)", "CC (Ref)","CL", "Pas","Pine", "CC")

aka2 = data.frame(ID = factor(rep(c("Reference","Disturbed"), each=4)))
rownames(aka2)<-rownames(its_diffabun7)
aka3 = list(ID = c(Reference = "#0072B2", Disturbed="#D55E00"))

its_diffabun7 <- t(scale(its_diffabun7)) 

my_palette <- colorRampPalette(c("white", "#800080"))(n = 299)

library(pheatmap)

#ITS phyla heatmap

heat4 <- pheatmap(its_diffabun7,color=my_palette,cluster_cols = F, cluster_rows = T, annotation_col = aka2, 
                  annotation_colors = aka3[1],annotation_names_col =F,legend=F,annotation_legend = FALSE, border_color = NA, main = "", gaps_col = c(4,8), cex=1)

#Figure 5 panel code

jpeg(filename="Figure5.jpeg", bg="transparent", res=500, units = "in", height=9, width=10.5) 

q5 <- plot_grid( NULL, NULL,heat3.1[[4]], heat4[[4]], ncol = 2,labels=c('A', 'B','C', 'D'), align="hv", label_size=18)
q5

dev.off()

#Vegetation analysis

library(vegan)

adonis2(veg[,c(7:24)]~ Use*Pair, data=veg,permutations = 999, method = "bray")

vegbray <- vegdist(veg[,c(7:24)], method="bray")

vegbray <- as.matrix(vegbray)



vegMDS <- metaMDS(veg[,c(7:24)], k=2, trymax=1000,distance="bray")

vegMDS <- metaMDS(veg[,c(7:24)], distance="bray", k=2, trymax=1000, autotransform=FALSE)

veg.scrs <- scores(vegMDS, display = "species") 
veg.scrs2 <- veg.scrs[c(1,4,7,9,13),]
veg.scrs2 <- as.data.frame(veg.scrs2)
veg.scrs2 <- cbind(veg.scrs2, vegspecies = rownames(veg.scrs2))
veg.scrs2

veg.scores <- as.data.frame(vegMDS$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
veg.scores$treat <- veg$Use 
veg.scores$pair <- veg$Pair

library(ggplot2)

vegord <- ggplot() + 
  geom_point(data=veg.scores,aes(x=MDS1,y=MDS2,shape=veg.scores$pair,color=veg.scores$treat),size=2) + # add the point markers
  scale_color_manual(values=c("#D55E00", "#0072B2")) +
  annotate("text", -.45, 1.4, label="Disturbance: ~italic(P) == 0.001", parse=TRUE, size=3, fontface="bold")+
  annotate("text", -.45, 1.2, label="`Watershed Pair`: ~italic(P) == 0.001", parse=TRUE, size=3, fontface="bold")+
  annotate("text", -.45, 1, label="`Treatment x Pair`: ~italic(P) == 0.001", parse=TRUE, size=3, fontface="bold")+
  labs(color='Treatment') +
  labs(shape='Watershed Pair') +
  theme(legend.text=element_text(size=12)) +
  theme_classic() +
  theme(axis.text=element_text(size=12)) +
  ggtitle("Woody vegetation communities") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=12)) +
  geom_segment(data=veg.scrs2,aes(x=0,xend=1*NMDS1,y=0,yend=1*NMDS2), arrow = arrow(length = unit(0.25, "cm")),colour="black",size=.25,inherit.aes=FALSE) +
  geom_text(data = veg.scrs2, aes(x = 1.25*NMDS1, y = 1.25*NMDS2, label = vegspecies),size = 3)

jpeg(filename="FigureS4.jpeg", bg="transparent", res=500, units = "in", height=4.5, width=6) 

vegord

dev.off()
