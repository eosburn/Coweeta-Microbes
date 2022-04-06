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

Bac_otus <- "C:/Users/ernie/OneDrive/Desktop/Completed Projects/Coweeta_drought/16S-table-with-taxonomy.biom"

Fung_otus <- "C:/Users/ernie/OneDrive/Desktop/Completed Projects/Coweeta_drought/ITS-table-with-taxonomy.biom"

x1 <- read_biom(Bac_otus)

x2 <- read_biom(Fung_otus)

OTU1 <- import_biom(x1)

OTU2 <- import_biom(x2)

meta <- "C:/Users/ernie/OneDrive/Desktop/Completed Projects/Coweeta_drought/16S_metadata_CWT.txt"

meta2 <- import_qiime_sample_data(meta)

OTU1 <- merge_phyloseq(OTU1, meta2)

OTU2 <- merge_phyloseq(OTU2, meta2)

meta3 <- data.frame(sample_data(OTU1))

OTU1 <- subset_samples(OTU1, id != "94")

bac <- rarefy_even_depth(OTU1, rngseed=TRUE)

OTU2 <- subset_samples(OTU2, id != "154")

ITS <- rarefy_even_depth(OTU2, rngseed=TRUE)

a1 <- t(data.frame(otu_table(bac)))

a2 <- t(data.frame(otu_table(ITS)))

b1 <- vegdist(a1, method = "bray")

b2 <- vegdist(a2, method = "bray")

c1 <- data.frame(sample_data(bac))

c1$Day <- as.factor(c1$Day)

c2 <- data.frame(sample_data(ITS))

c2$Day <- as.factor(c2$Day)

d1 <- cbind(a1, c1)

d2 <- cbind(a2, c2)

#Beta diversity stuff 16S

adonis2(a1 ~Drought*LandUse*Day, data = d1, permutations = 999, method = "bray")

bac.1 <- subset_samples(bac, Drought == "Drought")

a1.1 <- t(data.frame(otu_table(bac.1)))

b1.1 <- vegdist(a1.1, method = "bray")

c1.1 <- data.frame(sample_data(bac.1))

LandUse <- c1.1$LandUse

disp <- betadisper(b1.1,LandUse,type="centroid")
anova(disp)

e1 <- metaMDS(b1, k=2, trymax=1000)

data.scores1 <- as.data.frame(e1$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores1$site <- rownames(data.scores1)  # create a column of site names, from the rownames of data.scores
data.scores1$Treatment <- c1$Treatment 
data.scores1$Day <- c1$Day
head(data.scores1) 

data.scores1$Treatment <- factor(data.scores1$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
data.scores1$Treatment
head(data.scores1)

multi_ord <- aggregate(MDS1 ~ Treatment+Day, data.scores1,FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x)))
multi_ord

multi_ord2 <- aggregate(MDS2 ~ Treatment+Day, data.scores1,FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x)))
multi_ord2

multi_ord <- do.call(data.frame, multi_ord)
multi_ord

multi_ord2 <- do.call(data.frame, multi_ord2)
multi_ord2

multi_ord3 <- cbind(multi_ord, multi_ord2[,c(3:5)])
multi_ord3

multi_ord3$MDS1se <- multi_ord3$MDS1.sd / sqrt(multi_ord3$MDS1.n)
multi_ord3

multi_ord3$MDS2se <- multi_ord3$MDS2.sd / sqrt(multi_ord3$MDS2.n)
multi_ord3

ord_1 <- ggplot(multi_ord3, aes(x=MDS1.mean, y=MDS2.mean,color=Treatment,shape=Day)) + 
  geom_errorbar(aes(ymin=multi_ord3$MDS2.mean-multi_ord3$MDS2se, ymax=multi_ord3$MDS2.mean+multi_ord3$MDS2se), width=0.01, size=.75, position=position_dodge(0.9)) +
  geom_errorbarh(aes(xmin=multi_ord3$MDS1.mean-multi_ord3$MDS1se, xmax=multi_ord3$MDS1.mean+multi_ord3$MDS1se), height=0.01, size=.75) +
  geom_point(data=multi_ord3, aes(x=MDS1.mean,y=MDS2.mean,color=Treatment, shape=Day),size=6) + # add the point markers
  annotate("text", -.09, .13, label="Land~Use: ~italic(P) == 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.12, -.2, label="'2D Stress' == 0.12", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.09, .15, label="Drought: ~italic(P) == 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.09, .11, label="`Land Use x Drought`: ~italic(P) == 0.08", parse=TRUE, size=4.5, fontface="bold")+
  scale_color_brewer(palette="Set2") +
  scale_shape_manual(values = c(15, 16, 17,18,7,4)) +
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  theme_classic() +
  theme(axis.text=element_text(size=15)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color='Treatment') +
  
  labs(shape='Day') +
  theme(text=element_text(size=15)) +
  theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"))+
  theme(legend.margin=margin(c(-1,-1,-1,-1))) +
  ggtitle("16S Sequences") +
  theme(legend.text=element_text(size=14)) 
ord_1


#PLFA Bacterial composition

plfa_multi1 <- read.csv("C:/Users/ernie/OneDrive/Desktop/Current Projects/PLFA_drought/PLFA1.csv")
plfa_multi2 <- read.csv("C:/Users/ernie/OneDrive/Desktop/Current Projects/PLFA_drought/PLFA2.csv")

#Merging the two datasets 

plfa_multi <- merge(plfa_multi1, plfa_multi2,all.x=TRUE, all.y=TRUE)

plfa_multi[is.na(plfa_multi)] <- 0

plfa_multi <- plfa_multi[order(plfa_multi$ID_Only),]

plfa_multi <- plfa_multi[-c(1,14,27,40,53,66,79,92,105,118,131,144),]

meta3 <- meta3[order(meta3$id),]

plfa_multi <- cbind(plfa_multi, meta3)

plfa_matrix <- plfa_multi[,c(7:110)]

plfa_matrix2 <- decostand(plfa_matrix, "total")

rowSums(plfa_matrix2)

plfa_multi2 <- decostand(plfa_multi[,c(7:110)], "total")

a <- vegdist(plfa_multi2, method="bray")

adonis2(a~LandUse * Drought, data=plfa_multi, permutations=999)

e2 <- metaMDS(a, k=2, trymax=1000)

data.scores2 <- as.data.frame(e2$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores2$site <- rownames(data.scores2)  # create a column of site names, from the rownames of data.scores
data.scores2$Treatment <- plfa_multi$Treatment 
data.scores2$Day <- plfa_multi$Day
head(data.scores2) 

data.scores2$Treatment <- factor(data.scores2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
data.scores2$Treatment
head(data.scores2)

multi_ord4 <- aggregate(MDS1 ~ Treatment+Day, data.scores2,FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x)))
multi_ord4

multi_ord5 <- aggregate(MDS2 ~ Treatment+Day, data.scores2,FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x)))
multi_ord5

multi_ord4 <- do.call(data.frame, multi_ord4)
multi_ord4

multi_ord5 <- do.call(data.frame, multi_ord5)
multi_ord5

multi_ord6 <- cbind(multi_ord4, multi_ord5[,c(3:5)])
multi_ord6

multi_ord6$MDS1se <- multi_ord6$MDS1.sd / sqrt(multi_ord6$MDS1.n)
multi_ord6

multi_ord6$MDS2se <- multi_ord6$MDS2.sd / sqrt(multi_ord6$MDS2.n)
multi_ord6

multi_ord6$Day <- as.factor(multi_ord6$Day)

ord_2 <- ggplot(multi_ord6, aes(x=MDS1.mean, y=MDS2.mean,color=Treatment,shape=Day)) + 
  geom_errorbar(aes(ymin=multi_ord6$MDS2.mean-multi_ord6$MDS2se, ymax=multi_ord6$MDS2.mean+multi_ord6$MDS2se), width=0.01, size=.75, position=position_dodge(0.9)) +
  geom_errorbarh(aes(xmin=multi_ord6$MDS1.mean-multi_ord6$MDS1se, xmax=multi_ord6$MDS1.mean+multi_ord6$MDS1se), height=0.01, size=.75) +
  geom_point(data=multi_ord6, aes(x=MDS1.mean,y=MDS2.mean,color=Treatment, shape=Day),size=6) + # add the point markers
  annotate("text", -.045, -.03, label="Land~Use: ~italic(P) == 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.045, -.06, label="'2D Stress' == 0.13", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.045, -.04, label="Drought: ~italic(P) == 0.005", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.045, -.05, label="`Land Use x Drought`: ~italic(P) == 0.540", parse=TRUE, size=4.5, fontface="bold")+
  scale_color_brewer(palette="Set2") +
  scale_shape_manual(values = c(15, 16, 17,18,7,4)) +
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  theme_classic() +
  theme(axis.text=element_text(size=15)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color='Treatment') +
  labs(shape='Day') +
  theme(legend.margin=margin(c(-1,-1,-1,-1))) +
  theme(text=element_text(size=15)) +
  theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"))+
  ggtitle("Bacterial PLFAs") +
  theme(legend.text=element_text(size=15)) 
ord_2

setwd("C:/Users/ernie/OneDrive/Desktop/")

jpeg(filename="beta.jpeg", bg="transparent", res=500, units = "in", height=5.3, width=15) 

f4 <- plot_grid(ord_1, ord_2, ncol = 2, labels=c('A', 'B'),align="hv", label_size=25)

f4

dev.off()

#Correlations between the two ordinations


plfa_multi3 <- plfa_multi2[-58,]

a1 <- vegdist(plfa_multi3, method="bray")

mantel(b1, a1)

e3 <- metaMDS(a1,k=2, trymax=1000)

protest(e3, e1)

summary(procrustes(e1,e3))

p <- procrustes(e1,e3)

plot(p)



#Fungal:Bacterial stuff

din <- read.csv("C:/Users/ernie/OneDrive/Desktop/Completed Projects/Coweeta_drought/drought_chem.csv")

din$Day <- as.factor(din$Day)

plfa1 <- read.csv("C:/Users/ernie/OneDrive/Desktop/Current Projects/PLFA_drought/PLFA_groups.csv")

plfa1$Day <- as.factor(plfa1$Day)

plfa1$Total_Bac <- plfa1$Gram_neg + plfa1$Gram_pos

plfa1$FB <- plfa1$Total_Fungi/plfa1$Total_Bac

plfa1$ITS.16S <- din$ITS.16S

plfa1$ITS.16S2 <- din$ITS.16S2

glm1 <- lmer(FB ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm1))

glm1 <- glmer(FB~Drought*LandUse*Day + (1|Plot), data=plfa1 , family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Anova(glm1)

aggregate(FB~Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))


FB <- aggregate(FB~Treatment+LandUse+Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(FB)

FB2 <- do.call(data.frame, FB)
FB2

FB2$se <- FB2$FB.sd / sqrt(FB2$FB.n)
head(FB2)

FB2$Treatment <- factor(FB2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
FB2$Treatment
head(FB2)

line_plot_FB <- ggplot(FB2, aes(x=FB2$Day, y=FB2$FB.mean, group=FB2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=FB2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=FB2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=FB2$FB.mean-FB2$se, ymax=FB2$FB.mean+FB2$se), width=.25, size=.75) +
  annotate('text', 4, .07, label="Drought: ~italic(P) == 0.007", size=5, parse=TRUE)+
  annotate('text', 4, .065, label="Land~Use: italic(P) == 0.265", size = 5, parse=TRUE)+
  annotate('text', 4, .06, label="Land~Use %*% Drought: italic(P) == 0.97", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set2") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.6,.5)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("Fungal PLFA:Bacterial PLFA")) 
plot(line_plot_FB)

glm2 <- lmer(ITS.16S ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm2))

glm2 <- glmer(ITS.16S ~ Drought*LandUse*Day + (1|Plot), data=plfa1,family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Anova(glm2)

aggregate(ITS.16S~LandUse, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

aggregate(ITS.16S~Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

ITS.16S <- aggregate(ITS.16S~Treatment+LandUse+Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(ITS.16S)

ITS.16S2 <- do.call(data.frame, ITS.16S)
ITS.16S2

ITS.16S2$se <- ITS.16S2$ITS.16S.sd / sqrt(ITS.16S2$ITS.16S.n)
head(ITS.16S2)

ITS.16S2$Treatment <- factor(ITS.16S2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
ITS.16S2$Treatment
head(ITS.16S2)

line_plot_ITS.16S <- ggplot(ITS.16S2, aes(x=ITS.16S2$Day, y=ITS.16S2$ITS.16S.mean, group=ITS.16S2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=ITS.16S2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=ITS.16S2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=ITS.16S2$ITS.16S.mean-ITS.16S2$se, ymax=ITS.16S2$ITS.16S.mean+ITS.16S2$se), width=.25, size=.75) +
  annotate('text', 3, .17, label="Drought: ~italic(P) == 0.067", size=5, parse=TRUE)+
  annotate('text', 3, .16, label="Land~Use: italic(P) == .021", size = 5, parse=TRUE)+
  annotate('text', 3, .15, label="Land~Use %*% Drought: italic(P) == 0.703", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set2") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position="none") +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("ITS:16S Gene Copies (qPCR)")) 
plot(line_plot_ITS.16S)


cor.test(plfa1$ITS.16S, plfa1$FB, method="spearman")


setwd("C:/Users/ernie/OneDrive/Desktop/")

jpeg(filename="FB.jpeg", bg="transparent", res=500, units = "in", height=5.3, width=15) 

f4 <- plot_grid(line_plot_ITS.16S,line_plot_FB,  ncol = 2, labels=c('A', 'B'),align="hv", label_size=25)

f4

dev.off()

glm2 <- lmer(ITS.16S2 ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm2))
Anova(glm2)

glm2 <- glmer(ITS.16S2 ~ Drought*LandUse*Day + (1|Plot), data=plfa1,family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Anova(glm2)

aggregate(ITS.16S2~LandUse, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

aggregate(ITS.16S2~Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

ITS.16S2 <- aggregate(ITS.16S2~Treatment+LandUse+Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(ITS.16S2)

ITS.16S21 <- do.call(data.frame, ITS.16S2)
ITS.16S21

ITS.16S21$se <- ITS.16S21$ITS.16S2.sd / sqrt(ITS.16S21$ITS.16S2.n)
head(ITS.16S21)

ITS.16S21$Treatment <- factor(ITS.16S21$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
ITS.16S21$Treatment
head(ITS.16S21)

line_plot_ITS.16S2 <- ggplot(ITS.16S21, aes(x=ITS.16S21$Day, y=ITS.16S21$ITS.16S2.mean, group=ITS.16S21$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=ITS.16S21$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=ITS.16S21$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=ITS.16S21$ITS.16S2.mean-ITS.16S21$se, ymax=ITS.16S21$ITS.16S2.mean+ITS.16S21$se), width=.25, size=.75) +
  annotate('text', 3.5, .25, label="Drought: ~italic(P) == 0.56", size=5, parse=TRUE)+
  annotate('text', 3.5, .23, label="Land~Use: italic(P) == .012", size = 5, parse=TRUE)+
  annotate('text', 3.5, .21, label="Land~Use %*% Drought: italic(P) == 0.574", size=5, parse=TRUE)+
  labs(list(x =""),title="ITS:16S gene copies (corrected for operon number)") +
  theme_classic() +
  scale_color_brewer(palette="Set2") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("ITS:16S Gene Copies (qPCR)")) 
plot(line_plot_ITS.16S2)


cor.test(plfa1$ITS.16S, plfa1$FB, method="spearman")

#Bacterial Abundance

plfa1$bac <- din$X16S

glm5 <- lmer(bac ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm5))
Anova(glm5)

bacagg <- aggregate(bac~Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

bacagg <- do.call(data.frame,bacagg)

bacagg$bac.mean2 <- 10^(bacagg$bac.mean)

aggregate(bac.mean2~Drought+Day, data=bacagg, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

bac <- aggregate(bac~Treatment+LandUse+Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(bac)

bac2 <- do.call(data.frame, bac)
bac2

bac2$se <- bac2$bac.sd / sqrt(bac2$bac.n)
head(bac2)

bac2$Treatment <- factor(bac2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
bac2$Treatment
head(bac2)

line_plot_bac <- ggplot(bac2, aes(x=bac2$Day, y=bac2$bac.mean, group=bac2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=bac2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=bac2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=bac2$bac.mean-bac2$se, ymax=bac2$bac.mean+bac2$se), width=.25, size=.75) +
  annotate('text', 5, 9.05, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 5, 9, label="Land~Use: italic(P) == 0.833", size = 5, parse=TRUE)+
  annotate('text', 5, 8.95, label="Land~Use %*% Drought: italic(P) == 0.837", size=5, parse=TRUE)+
  labs(list(x =""),title="16S gene copies") +
  theme_classic() +
  scale_color_brewer(palette="Set2") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.2,.2)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(bquote(~'Log'[10]~'('~gene~copies~~gdw^-1*')'))
plot(line_plot_bac)

plfa1$bac2 <- din$X16S2

glm5 <- lmer(bac2 ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm5))
Anova(glm5)

bacagg <- aggregate(bac2~Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

bacagg <- do.call(data.frame,bacagg)

bacagg$bac2.mean2 <- 10^(bacagg$bac2.mean)

aggregate(bac2.mean2~Drought+Day, data=bacagg, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

bac <- aggregate(bac2~Treatment+LandUse+Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(bac)

bac2 <- do.call(data.frame, bac)
bac2

bac2$se <- bac2$bac2.sd / sqrt(bac2$bac2.n)
head(bac2)

bac2$Treatment <- factor(bac2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
bac2$Treatment
head(bac2)

line_plot_bac <- ggplot(bac2, aes(x=bac2$Day, y=bac2$bac2.mean, group=bac2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=bac2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=bac2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=bac2$bac2.mean-bac2$se, ymax=bac2$bac2.mean+bac2$se), width=.25, size=.75) +
  annotate('text', 5, 8.85, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 5, 8.8, label="Land~Use: italic(P) == 0.921", size = 5, parse=TRUE)+
  annotate('text', 5, 8.75, label="Land~Use %*% Drought: italic(P) == 0.869", size=5, parse=TRUE)+
  labs(list(x =""),title="16S gene copies (corrected for operon number)") +
  theme_classic() +
  scale_color_brewer(palette="Set2") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.2,.2)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(bquote(~'Log'[10]~'('~gene~copies~~gdw^-1*')'))
plot(line_plot_bac)

glm5 <- lmer(Total_Bac ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm5))
Anova(glm5)

bac <- aggregate(Total_Bac~Treatment+LandUse+Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(bac)

bac2 <- do.call(data.frame, bac)
bac2

bac2$se <- bac2$Total_Bac.sd / sqrt(bac2$Total_Bac.n)
head(bac2)

bac2$Treatment <- factor(bac2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
bac2$Treatment
head(bac2)

line_plot_bac <- ggplot(bac2, aes(x=bac2$Day, y=bac2$Total_Bac.mean, group=bac2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=bac2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=bac2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=bac2$Total_Bac.mean-bac2$se, ymax=bac2$Total_Bac.mean+bac2$se), width=.25, size=.75) +
  annotate('text', 5, 100000, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 5, 95000, label="Land~Use: italic(P) < 0.001", size = 5, parse=TRUE)+
  annotate('text', 5, 90000, label="Land~Use %*% Drought: italic(P) == 0.532", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set2") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.4,.95)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(bquote(~'Bacterial PLFA'~'('~nmol~~gdw^-1*')'))
plot(line_plot_bac)

cor.test(plfa1$Total_Bac, plfa1$bac, method="spearman")

#Total Fungi

plfa1$ITS <- din$ITS

plfa2 <- plfa1

plfa2 <- plfa2[-30,]


glm6 <- lmer(Total_Fungi ~ Drought*LandUse*Day + (1|Plot), data=plfa2)
qqPlot(resid(glm6))
Anova(glm6)

aggregate(Total_Fungi~Drought+Day, data=plfa2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

fung <- aggregate(Total_Fungi~Treatment+LandUse+Drought+Day, data=plfa2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(fung)

fung2 <- do.call(data.frame, fung)
fung2

fung2$se <- fung2$Total_Fungi.sd / sqrt(fung2$Total_Fungi.n)
head(fung2)

fung2$Treatment <- factor(fung2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
fung2$Treatment
head(fung2)

line_plot_fung <- ggplot(fung2, aes(x=fung2$Day, y=fung2$Total_Fungi.mean, group=fung2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=fung2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=fung2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=fung2$Total_Fungi.mean-fung2$se, ymax=fung2$Total_Fungi.mean+fung2$se), width=.25, size=.75) +
  annotate('text', 5.5, 5000, label="Drought: ~italic(P) == 0.740", size=5, parse=TRUE)+
  annotate('text', 5.5, 4800, label="Land~Use: italic(P) < 0.001", size = 5, parse=TRUE)+
  annotate('text', 5.5, 4600, label="Land~Use %*% Drought: italic(P) == 0.807", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set2") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.5,.65)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(bquote(~'Fungal PLFA'~'('~nmol~~gdw^-1*')')) 
plot(line_plot_fung)

glm7 <- lmer(ITS ~ Drought*LandUse*Day + (1|Plot), data=din)
qqPlot(resid(glm7))
Anova(glm7)

ITS <- aggregate(ITS~Treatment+LandUse+Drought+Day, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(ITS)

ITS2 <- do.call(data.frame, ITS)
ITS2

ITS2$se <- ITS2$ITS.sd / sqrt(ITS2$ITS.n)
head(ITS2)

ITS2$Treatment <- factor(ITS2$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
ITS2$Treatment
head(ITS2)

line_plot_ITS <- ggplot(ITS2, aes(x=ITS2$Day, y=ITS2$ITS.mean, group=ITS2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=ITS2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=ITS2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=ITS2$ITS.mean-ITS2$se, ymax=ITS2$ITS.mean+ITS2$se), width=.25, size=.75) +
  annotate('text', 5, 7.8, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 5, 7.7, label="Land~Use: italic(P) < 0.001", size = 5, parse=TRUE)+
  annotate('text', 5, 7.6, label="Land~Use %*% Drought: italic(P) == 0.567", size=5, parse=TRUE)+
  labs(list(x =""),title="ITS gene copies") +
  theme_classic() +
  scale_color_brewer(palette="Set2") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.2,.15)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(bquote(~'Log'[10]~'('~gene~copies~~gdw^-1*')'))
plot(line_plot_ITS)

#Bacterial phyla stuff

colnames(tax_table(bac)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

rk <- transform_sample_counts(bac, function(x) x/sum(x))

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

rk_phyla <- aggregate(Abundance~Sample+Phylum+Treatment+Day+LandUse+Plot+Drought, dat_rk, FUN=sum)

rk_phyla1 <- cast(rk_phyla, Sample+Treatment+Day+LandUse+Drought+Plot ~ Phylum, value="Abundance")

rk_phyla1$rk <- (rk_phyla1$p__Bacteroidetes+rk_phyla1$p__Proteobacteria)/(rk_phyla1$p__Verrucomicrobia+rk_phyla1$p__Acidobacteria)

rk_phyla1$Day <- as.factor(rk_phyla1$Day)

cor.test(plfa2$ITS, plfa2$Total_Fungi, method="spearman")

rk_phyla1$Sample <- as.numeric(rk_phyla1$Sample)

rk_phyla1 <- rk_phyla1[order(rk_phyla1$Sample),]

plfa2 <- plfa2[order(plfa2$ID_Only),]

cor.test(rk_phyla1$rk, plfa2$FB, method="spearman")

#16S copy number stuff

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

Bac_otus2 <- "C:/Users/ernie/OneDrive/Desktop/Current Projects/PLFA_Drought/16S-table.biom"

x1 <- read_biom(Bac_otus2)

x2 <- as.matrix(biom_data(x1))

head(x2)

colSums(x2)

copies <- read.csv("C:/Users/ernie/OneDrive/Desktop/Current Projects/PLFA_Drought/16S-copynumbers.csv")

x3 <- x2[rownames(x2) %in% copies$sequence,]

colSums(x3)

min(colSums(x3))

x4 <- t(x3)

x5 <- rrarefy(x4, 14265)

rowSums(x5)

library(funrar)

x6 <- make_relative(x5)

rowSums(x6)

copies <- copies[order(copies$sequence),]

x7 <- t(x6)

x8 <- cbind(rownames(x7), data.frame(x7, row.names=NULL,check.names=FALSE))

colnames(x8)[1] <- c("sequence")

x8 <- x8[order(x8$sequence),]

x9 = sapply(x8[,c(2:144)], '*', copies$X16S_rRNA_Count) 

colSums(x9)

table <- data.frame(t(x9))

table$copynumber <- rowSums(table)

table2 <- cbind(rownames(table), data.frame(table, row.names=NULL,check.names=FALSE))

colnames(table2)[1] <- c("id")

table2$id <- as.numeric(table2$id)

table2 <- table2[order(table2$id),]

meta <- read.csv("C:/Users/ernie/OneDrive/Desktop/Completed Projects/Coweeta_drought/16S_metadata_CWT2.csv")

meta <- meta[order(meta$id),]

table2 <- cbind(table2, meta)

table3 <- table2[,c(1,11509:11521)]

write.csv(table3,"C:/Users/ernie/OneDrive/Desktop/Current Projects/PLFA_drought/copies.csv", row.names = FALSE)

aggregate(copynumber~Treatment+Day,data=table2, FUN=mean)

table2$Day <- as.factor(table2$Day)

a <- lmer(copynumber ~ Drought * LandUse*Day + (1|Plot), data=table2)
qqPlot(resid(a))
a2 <- glmer(copynumber ~ Drought*LandUse*Day + (1|Plot), data=table2, family=Gamma(link=log))
AIC(a,a2)
Anova(a2)

aggregate(copynumber~Drought+Day, data=table2, FUN=mean)

aggregate(copynumber~Drought+Day+LandUse, data=table2, FUN=mean)

copy <- aggregate(copynumber~Treatment+LandUse+Drought+Day, data=table2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
copy

copy2 <- do.call(data.frame, copy)
copy2

copy2$se <- copy2$copynumber.sd / sqrt(copy2$copynumber.n)
head(copy2)

copy2$Treatment <- factor(copy2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
copy2$Treatment
head(copy2)

setwd("C:/Users/ernie/OneDrive/Desktop/")

jpeg(filename="fig3.jpeg", bg="transparent", res=1200, units = "in", height=4.5, width=6.6)

line_plot_copynumber <- ggplot(copy2, aes(x=copy2$Day, y=copy2$copynumber.mean, group=copy2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=copy2$Treatment), size=2) +
  geom_line(position=position_dodge(width=0.25),aes(color=copy2$Treatment),size=2) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=copy2$copynumber.mean-copy2$se, ymax=copy2$copynumber.mean+copy2$se), width=.5, size=.75) +
  annotate('text', 1.8, 1.85, label="Drought: ~italic(P) < 0.001", size=4, parse=TRUE)+
  annotate('text',1.8, 1.8, label="Land~Use: italic(P) == 0.434", size = 4, parse=TRUE)+
  annotate('text', 1.8, 1.75, label="Land~Use %*% Drought: italic(P) == 0.068", size=4, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set2") +
  theme(axis.title=element_text(size=17)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.85,.85)) +
  theme(legend.text = element_text(color = "black", size = 12)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("16S operon number per organism")) 
plot(line_plot_copynumber)

dev.off()

cor.test(table2$copynumber, plfa2$ITS.16S)

table2.45 <- table2[table2$Day!="1",]

table2.45 <- table2.45[table2.45$Day!="42",]

plfa2.45 <- plfa2[plfa2$Day!="1",]

plfa2.45 <- plfa2.45[plfa2.45$Day!="42",]

plot(x=plfa2.45$ITS.16S,y=table2.45$copynumber)

cor.test(table2.45$copynumber, plfa2.45$bac,method="spearman")

cor.test(table2.45$copynumber, plfa2.45$ITS.16S,method="spearman")

plot(x=rk_phyla1$rk, y=table2$copynumber)

cor.test(table2$copynumber, rk_phyla1$rk)

