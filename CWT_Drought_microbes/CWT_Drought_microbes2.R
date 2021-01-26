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

Bac_otus <- "C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/16S-table-with-taxonomy.biom"

Fung_otus <- "C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/ITS-table-with-taxonomy.biom"

x1 <- read_biom(Bac_otus)

x2 <- read_biom(Fung_otus)

OTU1 <- import_biom(x1)

OTU2 <- import_biom(x2)

meta <- "C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/16S_metadata_CWT.txt"

meta2 <- import_qiime_sample_data(meta)

OTU1 <- merge_phyloseq(OTU1, meta2)

OTU2 <- merge_phyloseq(OTU2, meta2)

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

#Alpha diversity stuff

bac_alpha <- estimate_richness(bac, split = TRUE)

bac_alpha <- cbind(bac_alpha, c1)

bac_alpha <- bac_alpha[order(bac_alpha$id),]

glm1 = glmer(Shannon~ Drought*LandUse*Day+(1|Plot), family=Gamma(link=log), data=bac_alpha,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Anova(glm1)
em_glm1 <- emmeans(glm1, ~ LandUse * Drought | Day)
contrast(em_glm1, method = "pairwise")

aggregate(Shannon~LandUse,data=bac_alpha, FUN=mean)
aggregate(Shannon~Drought+Day,data=bac_alpha, FUN=mean)


fung_alpha <- estimate_richness(ITS, split = TRUE)

fung_alpha <- cbind(fung_alpha, c2)

fung_alpha$Day <- as.factor(fung_alpha$Day)

aggregate(Shannon~Drought+Day,data=fung_alpha, FUN=mean)

glm2 = glmer(Shannon~ Drought*LandUse*Day+(1|Plot), family=Gamma(link=log), data=fung_alpha,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Anova(glm2)
em_glm2 <- emmeans(glm2, ~ LandUse * Drought | Day)
contrast(em_glm2, method = "pairwise")

bac_shan <- aggregate(Shannon~Treatment+LandUse+Drought+Day, data=bac_alpha, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(bac_shan)

bac_shan.2 <- aggregate(Observed~Treatment+Day, data=bac_alpha, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
bac_shan.2

bac_shan.3 <- aggregate(Simpson~Treatment+Day, data=bac_alpha, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
bac_shan.3

bac_shan2 <- do.call(data.frame, bac_shan)
bac_shan2

bac_shan2$se <- bac_shan2$Shannon.sd / sqrt(bac_shan2$Shannon.n)
head(bac_shan2)

bac_shan2$Treatment <- factor(bac_shan2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
bac_shan2$Treatment
head(bac_shan2)

library(grid)
text_1 <- textGrob("Day 1:\nPre-drought", gp=gpar(fontsize=15))
text_42 <- textGrob("Day 42:\nPost-drought", gp=gpar(fontsize=15))
text_43 <- textGrob("Day 43: 1 Day\nafter re-wet", gp=gpar(fontsize=15))
text_45 <- textGrob("Day 45: 3 Days\nafter re-wet", gp=gpar(fontsize=15))
text_56 <- textGrob("Day 56: 14 Days\nafter re-wet", gp=gpar(fontsize=15))
text_84 <- textGrob("Day 84: 42 Days\nafter re-wet", gp=gpar(fontsize=15))

line_plot_bac_alpha <- ggplot(bac_shan2, aes(x=bac_shan2$Day, y=bac_shan2$Shannon.mean, group=bac_shan2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=bac_shan2$Treatment), size=2) +
  geom_line(position=position_dodge(width=0.25),aes(color=bac_shan2$Treatment),size=1) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=bac_shan2$Shannon.mean-bac_shan2$se, ymax=bac_shan2$Shannon.mean+bac_shan2$se), width=.5, size=.75) +
  annotate('text', 2.15, 5.4, label="Drought: ~italic(P) < 0.001", size=4, parse=TRUE)+
  annotate('text', 2.15, 5.33, label="Disturbance: italic(P) < 0.001", size = 4, parse=TRUE)+
  annotate('text', 2.15, 5.26, label="Disturbance %*% Drought: italic(P) == 0.664", size=4, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=20)) +
  theme(legend.position="none") +
  theme(legend.text = element_text(color = "black", size = 10)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("16S Shannon Diversity (H')")) 
#theme(axis.title.x=element_blank(),
#axis.text.x=element_blank()) +
#theme(plot.margin = unit(c(1,1,4,1), "lines")) +
#annotation_custom(text_1,xmin=1,xmax=1,ymin=5.5,ymax=4.6) + 
#annotation_custom(text_42,xmin=2,xmax=2,ymin=5.5,ymax=4.6) +
#annotation_custom(text_43,xmin=3,xmax=3,ymin=5.5,ymax=4.6) + 
#annotation_custom(text_45,xmin=4,xmax=4,ymin=5.5,ymax=4.6) +
#annotation_custom(text_56,xmin=5,xmax=5,ymin=5.5,ymax=4.6) + 
#annotation_custom(text_84,xmin=6,xmax=6,ymin=5.5,ymax=4.6) +
#coord_cartesian(clip="off") 
plot(line_plot_bac_alpha)

fung_shan <- aggregate(Shannon~Treatment+LandUse+Drought+Day, data=fung_alpha, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(fung_shan)

fung_shan.2 <- aggregate(Observed~Treatment+Day, data=fung_alpha, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
fung_shan.2

fung_shan.3 <- aggregate(Simpson~Treatment+Day, data=fung_alpha, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
fung_shan.3

fung_shan2 <- do.call(data.frame, fung_shan)
fung_shan2

fung_shan2$se <- fung_shan2$Shannon.sd / sqrt(fung_shan2$Shannon.n)
head(fung_shan2)

fung_shan2$Treatment <- factor(fung_shan2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
fung_shan2$Treatment
head(fung_shan2)

line_plot_fung_alpha <- ggplot(fung_shan2, aes(x=fung_shan2$Day, y=fung_shan2$Shannon.mean, group=fung_shan2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=fung_shan2$Treatment), size=2) +
  geom_line(position=position_dodge(width=0.25),aes(color=fung_shan2$Treatment),size=1) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=fung_shan2$Shannon.mean-fung_shan2$se, ymax=fung_shan2$Shannon.mean+fung_shan2$se), width=.5, size=.75) +
  annotate('text', 4.8, 2.0, label="Drought: ~italic(P) == 0.024", size=4, parse=TRUE)+
  annotate('text', 4.8, 1.88, label="Disturbance: italic(P) == 0.67", size = 4, parse=TRUE)+
  annotate('text', 4.8, 1.76, label="Disturbance %*% Drought: italic(P) == 0.684", size=4, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.25,.88)) +
  theme(legend.text = element_text(color = "black", size = 10)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("ITS Shannon Diversity (H')")) 
#theme(axis.title.x=element_blank(),
#axis.text.x=element_blank()) +
#theme(plot.margin = unit(c(1,1,4,1), "lines")) +
#annotation_custom(text_1,xmin=1,xmax=1,ymin=5.5,ymax=-2.25) + 
#annotation_custom(text_42,xmin=2,xmax=2,ymin=5.5,ymax=-2.25) +
#annotation_custom(text_43,xmin=3,xmax=3,ymin=5.5,ymax=-2.25) + 
#annotation_custom(text_45,xmin=4,xmax=4,ymin=5.5,ymax=-2.25) +
#annotation_custom(text_56,xmin=5,xmax=5,ymin=5.5,ymax=-2.25) + 
#annotation_custom(text_84,xmin=6,xmax=6,ymin=5.5,ymax=-2.25) +
#coord_cartesian(clip="off") 
plot(line_plot_fung_alpha)
addSmallLegend <- function(myPlot, pointSize = 1.5, textSize = 12, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_blank(), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

line_plot_fung_alpha <- addSmallLegend(line_plot_fung_alpha)

line_plot_fung_alpha

#Resistance of Alpha Diversity

bac_alpha2 <- bac_alpha[bac_alpha$Day == '45',] 

bac_alpha3 <- bac_alpha2[bac_alpha2$Drought == 'Control',] 

bac_alpha4 <- bac_alpha2[bac_alpha2$Drought == 'Drought',] 

bac_alpha3 <- bac_alpha3[order(bac_alpha3$Plot),]

bac_alpha4 <- bac_alpha4[order(bac_alpha4$Plot),]

colnames(bac_alpha4)[6] <- "Shannon2"

bac_alpha5 <- cbind(bac_alpha3, bac_alpha4)

bac_alpha5$Taxon <- "Bacteria"

#Calculate Resistance

bac_alpha5$resist <- (1 - (2*abs(bac_alpha5$Shannon-bac_alpha5$Shannon2))/(bac_alpha5$Shannon+abs(bac_alpha5$Shannon-bac_alpha5$Shannon2)))

r1 <- lm(resist~LandUse,data=bac_alpha5)
anova(r1)


fung_alpha2 <- fung_alpha[fung_alpha$Day == '45',] 

fung_alpha3 <- fung_alpha2[fung_alpha2$Drought == 'Control',] 

fung_alpha4 <- fung_alpha2[fung_alpha2$Drought == 'Drought',] 

fung_alpha3 <- fung_alpha3[order(fung_alpha3$Plot),]

fung_alpha4 <- fung_alpha4[order(fung_alpha4$Plot),]

colnames(fung_alpha4)[6] <- "Shannon2"

fung_alpha5 <- cbind(fung_alpha3, fung_alpha4)

fung_alpha5$Taxon <- "Fungi"

#Calculate Resistance

fung_alpha5$resist <- (1 - (2*abs(fung_alpha5$Shannon-fung_alpha5$Shannon2))/(fung_alpha5$Shannon+abs(fung_alpha5$Shannon-fung_alpha5$Shannon2)))

r1 <- lm(resist~LandUse,data=fung_alpha5)
anova(r1)

resist1 <- rbind(bac_alpha5, fung_alpha5)

resist1 <- resist1[,c(45,53,54)]

r1 <- aov(resist ~ Taxon, data=resist1)
summary(r1)
hist(resist1$resist)
shapiro.test(resid(r1))

kruskal.test(resist ~ Taxon, data=resist1)

aggregate(resist~Taxon, data=resist1, FUN=mean)

resist1[resist1=="Reference"] <- "Reference-Drought"
resist1[resist1=="Disturbed"] <- "Disturbed-Drought"

resist1$LandUse <- factor(resist1$LandUse, levels = c("Reference-Drought", "Disturbed-Drought"))

box_plotResist_alpha<- ggplot(resist1, aes(x=resist1$Taxon, y=resist1$resist, fill=resist1$LandUse)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.size=.25) +
  annotate('text', 1, 0, label="Taxon: italic(P) == 0.017", size = 5, parse=TRUE)+
  theme_classic() +
  scale_fill_manual(values=c("#377EB8", "#984EA3")) +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  scale_y_continuous(expression("Resistance")) +
  theme(legend.position=c(.32,.2)) +
  theme(axis.title.x=element_blank())
plot(box_plotResist_alpha)

#Resilience of Alpha Diversity

bac_alpha6 <- bac_alpha[bac_alpha$Day == '84',] 

bac_alpha7 <- bac_alpha6[bac_alpha6$Drought == 'Control',] 

bac_alpha8 <- bac_alpha6[bac_alpha6$Drought == 'Drought',] 

bac_alpha7 <- bac_alpha7[order(bac_alpha7$Plot),]

bac_alpha8 <- bac_alpha8[order(bac_alpha8$Plot),]

colnames(bac_alpha7)[6] <- "Shannon3"

colnames(bac_alpha8)[6] <- "Shannon4"

bac_alpha5$Shannon3 <- bac_alpha7$Shannon3 

bac_alpha5$Shannon4 <- bac_alpha8$Shannon4

#Calculate Resilience

bac_alpha5$resil <- ((2*abs(bac_alpha5$Shannon-bac_alpha5$Shannon2)/(abs(bac_alpha5$Shannon-bac_alpha5$Shannon2)+abs(bac_alpha5$Shannon3-bac_alpha5$Shannon4)))-1)

r2 <- lm(resil ~ LandUse, data=bac_alpha5)
shapiro.test(resid(r2))
anova(r2)

bac_r <- bac_alpha5[,c(19,6,32,55,56,57,53)]

fung_alpha6 <- fung_alpha[fung_alpha$Day == '84',] 

fung_alpha7 <- fung_alpha6[fung_alpha6$Drought == 'Control',] 

fung_alpha8 <- fung_alpha6[fung_alpha6$Drought == 'Drought',] 

fung_alpha7 <- fung_alpha7[order(fung_alpha7$Plot),]

fung_alpha8 <- fung_alpha8[order(fung_alpha8$Plot),]

colnames(fung_alpha7)[6] <- "Shannon3"

colnames(fung_alpha8)[6] <- "Shannon4"

fung_alpha5$Shannon3 <- fung_alpha7$Shannon3 

fung_alpha5$Shannon4 <- fung_alpha8$Shannon4

#Calculate Resilience

fung_alpha5$resil <- ((2*abs(fung_alpha5$Shannon-fung_alpha5$Shannon2)/(abs(fung_alpha5$Shannon-fung_alpha5$Shannon2)+abs(fung_alpha5$Shannon3-fung_alpha5$Shannon4)))-1)

r2 <- lm(resil~LandUse,data=fung_alpha5)
shapiro.test(resid(r2))
anova(r2)


fung_r <- fung_alpha5[,c(19,6,32,55,56,57,53)]

resil1 <- rbind(bac_r, fung_r)

r2 <- lm(resil ~ Taxon, data=resil1)
anova(r2)
hist(resil1$resil)
shapiro.test(resid(r2))

aggregate(resil~Taxon, data=resil1, FUN=mean)


resil1[resil1=="Reference"] <- "Reference-Drought"
resil1[resil1=="Disturbed"] <- "Disturbed-Drought"

resil1$LandUse <- factor(resil1$LandUse, levels = c("Reference-Drought", "Disturbed-Drought"))


box_plotResil_alpha<- ggplot(resil1, aes(x=resist1$Taxon, y=resil1$resil, fill=resil1$LandUse)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.size=1) +
  annotate('text', 1.1, -.65, label="Taxon: italic(P) == 0.014", size = 5, parse=TRUE)+
  theme_classic() +
  scale_fill_manual(values=c("#377EB8", "#984EA3")) +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  scale_y_continuous(expression("Resiliance")) +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank())
plot(box_plotResil_alpha)

#Alpha Diversity multipanel figure

jpeg(filename="alpha.jpeg", bg="transparent", res=600, units = "in", height=8, width=10.25) 

f5 <- plot_grid(line_plot_bac_alpha, line_plot_fung_alpha, box_plotResist_alpha,box_plotResil_alpha, ncol = 2, labels=c('A', 'B','C','D'),align="hv", label_size=20)

f5

dev.off()

#Beta diversity stuff 16S

adonis2(a1 ~Drought*LandUse, data = d1, permutations = 999, method = "bray")

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
  geom_point(data=multi_ord3, aes(x=MDS1.mean,y=MDS2.mean,color=Treatment, shape=Day),size=5) + # add the point markers
  annotate("text", -.09, .15, label="Disturbance: ~italic(P) == 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.12, -.2, label="'2D Stress' == 0.12", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.09, .13, label="Drought: ~italic(P) == 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.09, .11, label="`Disturbance x Drought`: ~italic(P) == 0.08", parse=TRUE, size=4.5, fontface="bold")+
  scale_color_brewer(palette="Set1") +
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
  ggtitle("16S ASVs") +
  theme(legend.text=element_text(size=14)) 
ord_1

#Bray-Curtis distance from t1

bac_dist <- read.csv("C:/Users/ernie/Desktop/Coweeta_drought/16S-distances.csv")

bac_dist$Day <- as.factor(bac_dist$Day)

glm3 <- lmer(Distance ~ Drought*LandUse*Day + (1|Plot), data=bac_dist)
qqPlot(resid(glm3))
Anova(glm3)
em_glm3 <- emmeans(glm3, ~ LandUse * Drought | Day)
contrast(em_glm3, method = "pairwise")

bac_beta <- aggregate(Distance~Treatment+LandUse+Drought+Day, data=bac_dist, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(bac_beta)

bac_beta2 <- do.call(data.frame, bac_beta)
bac_beta2

bac_beta2$se <- bac_beta2$Distance.sd / sqrt(bac_beta2$Distance.n)
head(bac_beta2)

bac_beta2$Treatment <- factor(bac_beta2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
bac_beta2$Treatment
head(bac_beta2)

line_plot_bac_beta <- ggplot(bac_beta2, aes(x=bac_beta2$Day, y=bac_beta2$Distance.mean, group=bac_beta2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=bac_beta2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=bac_beta2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=bac_beta2$Distance.mean-bac_beta2$se, ymax=bac_beta2$Distance.mean+bac_beta2$se), width=.25, size=.75) +
  annotate('text', 1.5, .68, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 1.5, .66, label="Disturbance: italic(P) < 0.001", size = 5, parse=TRUE)+
  annotate('text', 1.5, .64, label="Disturbance %*% Drought: italic(P) == 0.722", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.68,.11)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("16S Bray-Curtis Distances")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  annotation_custom(text_42,xmin=1,xmax=1,ymin=.5,ymax=.15) +
  annotation_custom(text_43,xmin=2,xmax=2,ymin=.5,ymax=.15) + 
  annotation_custom(text_45,xmin=3,xmax=3,ymin=.5,ymax=.15) +
  annotation_custom(text_56,xmin=4,xmax=4,ymin=.5,ymax=.15) + 
  annotation_custom(text_84,xmin=5,xmax=5,ymin=.5,ymax=.15) +
  coord_cartesian(clip="off") 
plot(line_plot_bac_beta)

#Bray-curtis distance from control

bac_pairs <- read.csv("C:/Users/ernie/Desktop/Coweeta_drought/16S-pairs.csv")

bac_pairs$Day <- as.factor(bac_pairs$Day)

glm3.1 <- lmer(Distance ~ LandUse*Day + (1|Plot), data=bac_pairs)
qqPlot(resid(glm3.1))
Anova(glm3.1)
em_glm3.1 <- emmeans(glm3.1, ~ LandUse | Day)
contrast(em_glm3.1, method = "pairwise")

aggregate(Distance ~ LandUse, data=bac_pairs, FUN=mean)

aggregate(Distance ~ LandUse+Day, data=bac_pairs, FUN=mean)


bac_pairs <- aggregate(Distance~LandUse+Day, data=bac_pairs, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(bac_pairs)

bac_pairs2 <- do.call(data.frame, bac_pairs)
bac_pairs2

bac_pairs2$se <- bac_pairs2$Distance.sd / sqrt(bac_pairs2$Distance.n)
head(bac_pairs2)

bac_pairs2[bac_pairs2=="Reference"] <- "Reference-Drought"
bac_pairs2[bac_pairs2=="Disturbed"] <- "Disturbed-Drought"

bac_pairs2$LandUse <- as.factor(bac_pairs2$LandUse)

bac_pairs2$LandUse <- factor(bac_pairs2$LandUse, levels = c("Reference-Drought", "Disturbed-Drought"))
bac_pairs2$LandUse
head(bac_pairs2)

line_plot_bac_pairs <- ggplot(bac_pairs2, aes(x=bac_pairs2$Day, y=bac_pairs2$Distance.mean, group=bac_pairs2$LandUse)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=bac_pairs2$LandUse), size=3) +
  geom_line(position=position_dodge(width=0.25),aes(color=bac_pairs2$LandUse),size=2) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=bac_pairs2$Distance.mean-bac_pairs2$se, ymax=bac_pairs2$Distance.mean+bac_pairs2$se), width=.25, size=.75) +
  annotate('text', 2.2, .65, label="Disturbance: italic(P) == 0.004", size = 6, parse=TRUE)+
  annotate('text', 2.2, .62, label="Day: italic(P) < 0.001", size = 6, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_manual(values=c("#377EB8", "#984EA3")) +
  theme(axis.title=element_text(size=17)) +
  theme(text=element_text(size=24)) +
  theme(legend.position=c(.65,.15)) +
  theme(legend.text = element_text(color = "black", size = 18)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("16S Bray-Curtis Distance from Control")) 
#theme(axis.title.x=element_blank(),
#axis.text.x=element_blank()) +
#theme(plot.margin = unit(c(1,1,4,1), "lines")) +
#annotation_custom(text_1,xmin=1,xmax=1,ymin=.57,ymax=.15) +
#annotation_custom(text_42,xmin=2,xmax=2,ymin=.57,ymax=.15) +
#annotation_custom(text_43,xmin=3,xmax=3,ymin=.57,ymax=.15) + 
#annotation_custom(text_45,xmin=4,xmax=4,ymin=.57,ymax=.15) +
#annotation_custom(text_56,xmin=5,xmax=5,ymin=.57,ymax=.15) + 
#annotation_custom(text_84,xmin=6,xmax=6,ymin=.57,ymax=.15) +
#coord_cartesian(clip="off") 
plot(line_plot_bac_pairs)

#Beta diversity ITS

adonis2(a2 ~Drought*LandUse, data = d2, permutations = 999, method = "bray")

e2 <- metaMDS(b2, k=2, trymax=1000)

data.scores2 <- as.data.frame(e2$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores2$site <- rownames(data.scores2)  # create a column of site names, from the rownames of data.scores
data.scores2$Treatment <- c2$Treatment 
data.scores2$Day <- c2$Day
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

ord_2 <- ggplot(multi_ord6, aes(x=MDS1.mean, y=MDS2.mean,color=Treatment,shape=Day)) + 
  geom_errorbar(aes(ymin=multi_ord6$MDS2.mean-multi_ord6$MDS2se, ymax=multi_ord6$MDS2.mean+multi_ord6$MDS2se), width=0.01, size=.75, position=position_dodge(0.9)) +
  geom_errorbarh(aes(xmin=multi_ord6$MDS1.mean-multi_ord6$MDS1se, xmax=multi_ord6$MDS1.mean+multi_ord6$MDS1se), height=0.01, size=.75) +
  geom_point(data=multi_ord6, aes(x=MDS1.mean,y=MDS2.mean,color=Treatment, shape=Day),size=5) + # add the point markers
  annotate("text", .1, -.12, label="Disturbance: ~italic(P) == 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.15, .11, label="'2D Stress' == 0.20", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", .1, -.14, label="Drought: ~italic(P) == 0.09", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", .1, -.16, label="`Disturbance x Drought`: ~italic(P) == 0.851", parse=TRUE, size=4.5, fontface="bold")+
  scale_color_brewer(palette="Set1") +
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
  ggtitle("ITS ASVs") +
  theme(legend.text=element_text(size=15)) 
ord_2

#Bray-Curtis distance from t1

fung_dist <- read.csv("C:/Users/ernie/Desktop/Coweeta_drought/ITS-distances.csv")

fung_dist$Day <- as.factor(fung_dist$Day)

glm4 <- lmer(Distance ~ Drought*LandUse*Day + (1|Plot), data=fung_dist)
qqPlot(resid(glm4))
Anova(glm4)
em_glm4 <- emmeans(glm4, ~ LandUse * Drought | Day)
contrast(em_glm4, method = "pairwise")

fung_beta <- aggregate(Distance~Treatment+LandUse+Drought+Day, data=fung_dist, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(fung_beta)

fung_beta2 <- do.call(data.frame, fung_beta)
fung_beta2

fung_beta2$se <- fung_beta2$Distance.sd / sqrt(fung_beta2$Distance.n)
head(fung_beta2)

fung_beta2$Treatment <- factor(fung_beta2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
fung_beta2$Treatment
head(fung_beta2)

line_plot_fung_beta <- ggplot(fung_beta2, aes(x=fung_beta2$Day, y=fung_beta2$Distance.mean, group=fung_beta2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=fung_beta2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=fung_beta2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=fung_beta2$Distance.mean-fung_beta2$se, ymax=fung_beta2$Distance.mean+fung_beta2$se), width=.25, size=.75) +
  annotate('text', 2, .89, label="Drought: ~italic(P) == 0.045", size=5, parse=TRUE)+
  annotate('text', 2, .87, label="Disturbance: italic(P) == 0.284", size = 5, parse=TRUE)+
  annotate('text', 2, .85, label="Disturbance %*% Drought: italic(P) == 0.227", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.68,.1)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("ITS Bray-Curtis Distances")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  annotation_custom(text_42,xmin=1,xmax=1,ymin=.9,ymax=.15) +
  annotation_custom(text_43,xmin=2,xmax=2,ymin=.9,ymax=.15) + 
  annotation_custom(text_45,xmin=3,xmax=3,ymin=.9,ymax=.15) +
  annotation_custom(text_56,xmin=4,xmax=4,ymin=.9,ymax=.15) + 
  annotation_custom(text_84,xmin=5,xmax=5,ymin=.9,ymax=.15) +
  coord_cartesian(clip="off") 
plot(line_plot_fung_beta)

#Bray-Curtis distance from control

fung_pairs <- read.csv("C:/Users/ernie/Desktop/Coweeta_drought/ITS-pairs.csv")

fung_pairs$Day <- as.factor(fung_pairs$Day)

glm4.1 <- lmer(Distance ~ LandUse*Day + (1|Plot), data=fung_pairs)
qqPlot(resid(glm4.1))
Anova(glm4.1)
em_glm4.1 <- emmeans(glm4.1, ~ LandUse | Day)
contrast(em_glm4.1, method = "pairwise")

fung_pairs <- aggregate(Distance~LandUse+Day, data=fung_pairs, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(fung_pairs)

fung_pairs2 <- do.call(data.frame, fung_pairs)
fung_pairs2

fung_pairs2$se <- fung_pairs2$Distance.sd / sqrt(fung_pairs2$Distance.n)
head(fung_pairs2)

fung_pairs2$LandUse <- factor(fung_pairs2$LandUse, levels = c("Reference", "Disturbed"))
fung_pairs2$LandUse
head(fung_pairs2)

line_plot_fung_pairs <- ggplot(fung_pairs2, aes(x=fung_pairs2$Day, y=fung_pairs2$Distance.mean, group=fung_pairs2$LandUse)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=fung_pairs2$LandUse), size=3) +
  geom_line(position=position_dodge(width=0.25),aes(color=fung_pairs2$LandUse),size=2) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=fung_pairs2$Distance.mean-fung_pairs2$se, ymax=fung_pairs2$Distance.mean+fung_pairs2$se), width=.25, size=.75) +
  annotate('text', 4.5, .54, label="Disturbance: italic(P) == 0.598", size = 6, parse=TRUE)+
  annotate('text', 4.5, .5, label="Day: italic(P) == 0.14", size = 6, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_manual(values=c("#377EB8", "#984EA3")) +
  theme(axis.title=element_text(size=17)) +
  theme(text=element_text(size=24)) +
  theme(legend.position="none") +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("ITS Bray-Curtis Distance from Control")) 
#theme(axis.title.x=element_blank(),
#axis.text.x=element_blank()) +
#theme(plot.margin = unit(c(1,1,4,1), "lines")) +
#annotation_custom(text_1,xmin=1,xmax=1,ymin=.72,ymax=.15) +
#annotation_custom(text_42,xmin=2,xmax=2,ymin=.72,ymax=.15) +
#annotation_custom(text_43,xmin=3,xmax=3,ymin=.72,ymax=.15) + 
#annotation_custom(text_45,xmin=4,xmax=4,ymin=.72,ymax=.15) +
#annotation_custom(text_56,xmin=5,xmax=5,ymin=.72,ymax=.15) + 
#annotation_custom(text_84,xmin=6,xmax=6,ymin=.72,ymax=.15) +
#coord_cartesian(clip="off") 
plot(line_plot_fung_pairs)

#Resistance of Beta Diversity

bac_beta2 <- bac_dist[bac_dist$Day == '44',] 

###Day 44 is a mistake, it's the same as day 45 everywhere else

bac93 <- subset_samples(bac, id == "93")

bac100 <- subset_samples(bac, id == "100")

bac93.1 <- t(data.frame(otu_table(bac93)))

bac100.1 <- t(data.frame(otu_table(bac100)))

bac100.2 <- rbind(bac93.1, bac100.1)

vegdist(bac100.2, method="bray")

bac_beta3 <- bac_beta2[bac_beta2$Drought == 'Control',] 

bac_beta3 <- bac_beta3[,c(7,8,11,12,14)]

bac_beta4 <- bac_beta2[bac_beta2$Drought == 'Drought',] 

bac_beta4 <- bac_beta4[,c(7,8,11,12,14)]

bac100.4 <- c("Disturbed", "W6P2", "Drought", "44","0.578")

bac_beta4 <- rbind(bac_beta4, bac100.4)

bac_beta3 <- bac_beta3[order(bac_beta3$Plot),]

bac_beta4 <- bac_beta4[order(bac_beta4$Plot),]

colnames(bac_beta4)[5] <- "Distance2"

bac_beta5 <- cbind(bac_beta3, bac_beta4)

bac_beta5 <- bac_beta5[,c(1,5,10)]

bac_beta5$Taxon <- "Bacteria"

bac_beta5$Distance2 <- as.numeric(bac_beta5$Distance2)

#Calculate Resistance

bac_beta5$resist <- (1 - (2*abs(bac_beta5$Distance-bac_beta5$Distance2))/(bac_beta5$Distance+abs(bac_beta5$Distance-bac_beta5$Distance2)))

r1 <- lm(resist~LandUse, data=bac_beta5)
hist(bac_beta5$resist)
shapiro.test(resid(r1))
anova(r1)

fung_beta2 <- fung_dist[fung_dist$Day == '44',] 

###Day 44 is a mistake, it's the same as day 45 everywhere else

fung_beta3 <- fung_beta2[fung_beta2$Drought == 'Control',] 

fung_beta4 <- fung_beta2[fung_beta2$Drought == 'Drought',] 

fung_beta3 <- fung_beta3[order(fung_beta3$Plot),]

fung_beta4 <- fung_beta4[order(fung_beta4$Plot),]

colnames(fung_beta4)[14] <- "Distance2"

fung_beta5 <- cbind(fung_beta3, fung_beta4)

fung_beta5 <- fung_beta5[,c(7,14,28)]

fung_beta5$Taxon <- "Fungi"

#Calculate Resistance

fung_beta5$resist <- (1 - (2*abs(fung_beta5$Distance-fung_beta5$Distance2))/(fung_beta5$Distance+abs(fung_beta5$Distance-fung_beta5$Distance2)))

r2 <- lm(resist~LandUse, data=fung_beta5)
hist(fung_beta5$resist)
shapiro.test(resid(r2))
anova(r2)

resist2 <- rbind(bac_beta5, fung_beta5)

aggregate(resist~Taxon,data=resist2,FUN=mean)

r3 <- lm(resist~Taxon, data=resist2)
hist(resist2$resist)
shapiro.test(resid(r3))
anova(r3)

resist2$LandUse <- factor(resist2$LandUse, levels = c("Reference", "Disturbed"))

box_plotResist_beta<- ggplot(resist2, aes(x=resist2$Taxon, y=resist2$resist, fill=resist2$LandUse)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.size=.25) +
  annotate('text', 1.5, .95, label="Taxon: italic(P) == 0.027", size = 6, parse=TRUE)+
  theme_classic() +
  scale_fill_manual(values=c("#377EB8", "#984EA3")) +
  theme(axis.title=element_text(size=25)) +
  theme(text=element_text(size=25)) +
  theme(legend.text = element_text(color = "black", size = 18)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  theme(legend.position="none") +
  ylab('Resistance') +
  theme(axis.title.x=element_blank())
plot(box_plotResist_beta)

#Resilience of Beta Diversity

bac_beta6 <- bac_dist[bac_dist$Day == '84',] 

bac_beta7 <- bac_beta6[bac_beta6$Drought == 'Control',] 

bac_beta8 <- bac_beta6[bac_beta6$Drought == 'Drought',] 

bac104 <- subset_samples(bac, id == "104")

bac104.1 <- t(data.frame(otu_table(bac104)))

bac104.2 <- rbind(bac93.1, bac104.1)

vegdist(bac104.2, method="bray")

bac_beta8 <- bac_beta8[,c(7,8,11,12,14)]

bac100.4 <- c("Disturbed", "W6P2", "Drought", "84","0.6245")

bac_beta8 <- rbind(bac_beta8, bac100.4)

bac_beta7 <- bac_beta7[order(bac_beta7$Plot),]

bac_beta8 <- bac_beta8[order(bac_beta8$Plot),]

bac_beta7 <- bac_beta7[,c(7,8,11,12,14)]

colnames(bac_beta7)[5] <- "Distance3"

colnames(bac_beta8)[5] <- "Distance4"

bac_beta5$Distance3 <- bac_beta7$Distance3

bac_beta5$Distance4<- bac_beta8$Distance4

bac_beta5$Distance4 <- as.numeric(bac_beta5$Distance4)

#Calculate Resilience

bac_beta5$resil <- ((2*abs(bac_beta5$Distance-bac_beta5$Distance2)/(abs(bac_beta5$Distance-bac_beta5$Distance2)+abs(bac_beta5$Distance3-bac_beta5$Distance4)))-1)

r4 <- aov(resil ~ LandUse, data=bac_beta5)
shapiro.test(resid(r4))
summary(r4)
anova(r4)

aggregate(resil~LandUse,data=bac_beta5,FUN=mean)

fung_beta6 <- fung_dist[fung_dist$Day == '84',] 

fung_beta7 <- fung_beta6[fung_beta6$Drought == 'Control',] 

fung_beta8 <- fung_beta6[fung_beta6$Drought == 'Drought',] 

fung_beta7 <- fung_beta7[order(fung_beta7$Plot),]

fung_beta8 <- fung_beta8[order(fung_beta8$Plot),]

colnames(fung_beta7)[14] <- "Distance3"

colnames(fung_beta8)[14] <- "Distance4"

fung_beta5$Distance3 <- fung_beta7$Distance3

fung_beta5$Distance4<- fung_beta8$Distance4

fung_beta5 <- fung_beta5[-10,]

#Calculate Resilience

fung_beta5$resil <- ((2*abs(fung_beta5$Distance-fung_beta5$Distance2)/(abs(fung_beta5$Distance-fung_beta5$Distance2)+abs(fung_beta5$Distance3-fung_beta5$Distance4)))-1)

r4 <- aov(resil ~ LandUse, data=fung_beta5)
anova(r4)

resil2 <- rbind(bac_beta5, fung_beta5)

r4 <- aov(resil ~ Taxon, data=resil2)
hist(resil2$resil)
shapiro.test(resid(r4))
anova(r4)

resil2[resil2=="Reference"] <- "Reference-Drought"
resil2[resil2=="Disturbed"] <- "Disturbed-Drought"

resil2$LandUse <- factor(resil2$LandUse, levels = c("Reference-Drought", "Disturbed-Drought"))

box_plotResil_beta<- ggplot(resil2, aes(x=resil2$Taxon, y=resil2$resil, fill=resil2$LandUse)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.size=.25) +
  annotate('text', 1,1, label="Taxon: italic(P) == 0.834", size = 6, parse=TRUE)+
  annotate('text', 1, .7, label="???", size = 8)+
  theme_classic() +
  scale_fill_manual(values=c("#377EB8", "#984EA3")) +
  theme(axis.title=element_text(size=25)) +
  theme(text=element_text(size=25)) +
  theme(legend.position=c(.25,.2)) +
  theme(legend.text = element_text(color = "black", size = 18)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  ylab('Resilience') +
  theme(axis.title.x=element_blank())
plot(box_plotResil_beta)

#Beta Diversity multi-panel figure

jpeg(filename="beta.jpeg", bg="transparent", res=500, units = "in", height=16, width=14) 

f4 <- plot_grid(ord_1, ord_2,line_plot_bac_pairs,line_plot_fung_pairs, box_plotResist_beta, box_plotResil_beta, ncol = 2, labels=c('A', 'B','C','D','E','F'),align="hv", label_size=25)

f4

dev.off()

#bacterial phyla stuff

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

p1 <- lmer(p__Acidobacteria~ LandUse*Drought*Day+(1|Plot), data=rk_phyla1)
qqPlot(resid(p1))
Anova(p1)

aggregate(p__Acidobacteria ~ Drought+Day, data=rk_phyla1,FUN=mean)

aggregate(p__Acidobacteria ~ Drought+Day+LandUse, data=rk_phyla1,FUN=mean)

p2 <- lmer(p__Proteobacteria~ LandUse*Drought*Day+(1|Plot), data=rk_phyla1)
qqPlot(resid(p2))
Anova(p2)

aggregate(p__Proteobacteria ~ Drought+Day, data=rk_phyla1,FUN=mean)


p3 <- lmer(p__Verrucomicrobia~ LandUse*Drought*Day+(1|Plot), data=rk_phyla1)
qqPlot(resid(p3))
Anova(p3)

p4 <- lmer(p__Actinobacteria~ LandUse*Drought*Day+(1|Plot), data=rk_phyla1)
qqPlot(resid(p4))
Anova(p4)

aggregate(p__Actinobacteria~ Drought+Day, data=rk_phyla1,FUN=mean)

p5 <- glmer(p__Bacteroidetes~ LandUse*Drought*Day+(1|Plot), data=rk_phyla1,family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(p5))
Anova(p5)

p6 <- glmer(p__Firmicutes+.0001~ LandUse*Drought*Day+(1|Plot), data=rk_phyla1,family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(p6))
Anova(p6)

aggregate(p__Firmicutes~ Drought+Day, data=rk_phyla1,FUN=mean)

p7 <- lmer(p__Planctomycetes~ LandUse*Drought*Day+(1|Plot), data=rk_phyla1)
qqPlot(resid(p7))
Anova(p7)

p8 <- lmer(p__Chloroflexi~ LandUse*Drought*Day+(1|Plot), data=rk_phyla1)
qqPlot(resid(p8))
Anova(p8)

p9 <- glmer(p__Nitrospirae~ LandUse*Drought*Day+(1|Plot), data=rk_phyla1,family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(p9))
Anova(p9)

aggregate(p__Nitrospirae~ Drought+Day, data=rk_phyla1,FUN=mean)

glm0 = glmer(rk~ LandUse*Drought*Day+(1|Plot), data=rk_phyla1,family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(glm0))
Anova(glm0)
em_glm0 <- emmeans(glm0, ~ LandUse * Drought | Day)
contrast(em_glm0, method = "pairwise")

#Bacterial phyla line plots

bac_ag <- aggregate(Abundance~ Treatment + Drought + Day + LandUse + Phylum, data=rk_phyla, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

bac_ag <- bac_ag[bac_ag$Phylum != 'Other',]

bac_ag <- bac_ag[bac_ag$Drought == 'Drought',]    

bac_ag <- bac_ag[bac_ag$LandUse == 'Reference',]   

bac_ag <- do.call(data.frame, bac_ag)
bac_ag

bac_ag$Phylum <- as.character(bac_ag$Phylum)

bac_ag[bac_ag=="p__Acidobacteria"] <- "Acidobacteria*???"
bac_ag[bac_ag=="p__Actinobacteria"] <- "Actinobacteria*???"
bac_ag[bac_ag=="p__Nitrospirae"] <- "Nitrospirae*???"
bac_ag[bac_ag=="p__Bacteroidetes"] <- "Bacteroidetes"
bac_ag[bac_ag=="p__Proteobacteria"] <- "Proteobacteria*"
bac_ag[bac_ag=="p__Verrucomicrobia"] <- "Verrucomicrobia*"
bac_ag[bac_ag=="p__Chloroflexi"] <- "Chloroflexi"
bac_ag[bac_ag=="p__Planctomycetes"] <- "Planctomycetes*"
bac_ag[bac_ag=="p__Firmicutes"] <- "Firmicutes*"


bac_ag$Phylum <- as.factor(bac_ag$Phylum)

bac_ag$Day <- as.factor(bac_ag$Day)

bac_ag1 <- bac_ag

legend_title1 <- "Bacterial Phyla"

bac1 <- ggplot(bac_ag1, aes(x=bac_ag1$Day, y=bac_ag1$Abundance.mean, group=bac_ag1$Phylum)) + 
  geom_point(position=position_dodge(width=0.1),aes(color=bac_ag1$Phylum), size=2) +
  geom_line(position=position_dodge(width=0.1),aes(color=bac_ag1$Phylum),size=1) +
  #geom_errorbar(position=position_dodge(width=0.9),aes(ymin=resp2.2$Cmin.mean-resp2.2$se, ymax=resp2.2$Cmin.mean+resp2.2$se), width=3, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(legend_title1,palette="Set1",direction=-1) +
  labs(color='Bacterial Phyla') +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 12)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  #theme(
        #legend.margin=margin(c(1,5,5,5))) +
  ylab("Relative Abundance") +
  scale_y_continuous(limits=c(0, .47)) 
#theme(axis.title.x=element_blank(),
#axis.text.x=element_blank()) +
#theme(plot.margin = unit(c(1,1,4,1), "lines")) +
#annotation_custom(text_1,xmin=1,xmax=1,ymin=.1,ymax=-.2) + 
#annotation_custom(text_42,xmin=2,xmax=2,ymin=.1,ymax=-.2) +
# annotation_custom(text_43,xmin=3,xmax=3,ymin=.1,ymax=-.2) + 
#annotation_custom(text_45,xmin=4,xmax=4,ymin=.1,ymax=-.2) +
#annotation_custom(text_56,xmin=5,xmax=5,ymin=.1,ymax=-.2) + 
#annotation_custom(text_84,xmin=6,xmax=6,ymin=.1,ymax=-.2) +
#coord_cartesian(clip="off") 
bac1

bac_ag <- aggregate(Abundance~ Treatment + Drought + Day + LandUse + Phylum, data=rk_phyla, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

bac_ag <- bac_ag[bac_ag$Phylum != 'Other',]

bac_ag <- bac_ag[bac_ag$Drought == 'Drought',]    

bac_ag <- bac_ag[bac_ag$LandUse == 'Disturbed',]   

bac_ag <- do.call(data.frame, bac_ag)
bac_ag

bac_ag$Phylum <- as.character(bac_ag$Phylum)

bac_ag[bac_ag=="p__Acidobacteria"] <- "Acidobacteria*???"
bac_ag[bac_ag=="p__Actinobacteria"] <- "Actinobacteria*???"
bac_ag[bac_ag=="p__Nitrospirae"] <- "Nitrospirae*???"
bac_ag[bac_ag=="p__Bacteroidetes"] <- "Bacteroidetes"
bac_ag[bac_ag=="p__Proteobacteria"] <- "Proteobacteria*"
bac_ag[bac_ag=="p__Verrucomicrobia"] <- "Verrucomicrobia*"
bac_ag[bac_ag=="p__Chloroflexi"] <- "Chloroflexi"
bac_ag[bac_ag=="p__Planctomycetes"] <- "Planctomycetes*"
bac_ag[bac_ag=="p__Firmicutes"] <- "Firmicutes*"


bac_ag$Phylum <- as.factor(bac_ag$Phylum)

bac_ag$Day <- as.factor(bac_ag$Day)

bac_ag2 <- bac_ag

bac2 <- ggplot(bac_ag2, aes(x=bac_ag2$Day, y=bac_ag2$Abundance.mean, group=bac_ag2$Phylum)) + 
  geom_point(position=position_dodge(width=0.1),aes(color=bac_ag2$Phylum), size=2) +
  geom_line(position=position_dodge(width=0.1),aes(color=bac_ag2$Phylum),size=1) +
  #geom_errorbar(position=position_dodge(width=0.9),aes(ymin=resp2.2$Cmin.mean-resp2.2$se, ymax=resp2.2$Cmin.mean+resp2.2$se), width=3, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1",direction=-1) +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 12)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  ylab("Relative Abundance") +
  scale_y_continuous(limits=c(0, .47)) 
#theme(axis.title.x=element_blank(),
#axis.text.x=element_blank()) +
#theme(plot.margin = unit(c(1,1,4,1), "lines")) +
#annotation_custom(text_1,xmin=1,xmax=1,ymin=.1,ymax=-.2) + 
#annotation_custom(text_42,xmin=2,xmax=2,ymin=.1,ymax=-.2) +
# annotation_custom(text_43,xmin=3,xmax=3,ymin=.1,ymax=-.2) + 
#annotation_custom(text_45,xmin=4,xmax=4,ymin=.1,ymax=-.2) +
#annotation_custom(text_56,xmin=5,xmax=5,ymin=.1,ymax=-.2) + 
#annotation_custom(text_84,xmin=6,xmax=6,ymin=.1,ymax=-.2) +
#coord_cartesian(clip="off") 
bac2


bac_ag <- aggregate(Abundance~ Treatment + Drought + Day + LandUse + Phylum, data=rk_phyla, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

bac_ag <- bac_ag[bac_ag$Phylum != 'Other',]

bac_ag <- bac_ag[bac_ag$Drought == 'Control',]    

bac_ag <- bac_ag[bac_ag$LandUse == 'Reference',]   

bac_ag <- do.call(data.frame, bac_ag)
bac_ag

bac_ag$Phylum <- as.character(bac_ag$Phylum)

bac_ag[bac_ag=="p__Acidobacteria"] <- "Acidobacteria"
bac_ag[bac_ag=="p__Actinobacteria"] <- "Actinobacteria"
bac_ag[bac_ag=="p__Nitrospirae"] <- "Nitrospirae"
bac_ag[bac_ag=="p__Bacteroidetes"] <- "Bacteroidetes"
bac_ag[bac_ag=="p__Proteobacteria"] <- "Proteobacteria"
bac_ag[bac_ag=="p__Verrucomicrobia"] <- "Verrucomicrobia"
bac_ag[bac_ag=="p__Chloroflexi"] <- "Chloroflexi"
bac_ag[bac_ag=="p__Planctomycetes"] <- "Planctomycetes"
bac_ag[bac_ag=="p__Firmicutes"] <- "Firmicutes"


bac_ag$Phylum <- as.factor(bac_ag$Phylum)

bac_ag$Day <- as.factor(bac_ag$Day)

bac_ag3 <- bac_ag

bac3<- ggplot(bac_ag, aes(x=bac_ag3$Day, y=bac_ag3$Abundance.mean, group=bac_ag3$Phylum)) + 
  geom_point(position=position_dodge(width=0.1),aes(color=bac_ag3$Phylum), size=2) +
  geom_line(position=position_dodge(width=0.1),aes(color=bac_ag3$Phylum),size=1) +
  #geom_errorbar(position=position_dodge(width=0.9),aes(ymin=resp2.2$Cmin.mean-resp2.2$se, ymax=resp2.2$Cmin.mean+resp2.2$se), width=3, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(legend_title1, palette="Set1", direction=-1) +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 12)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(
        legend.margin=margin(c(1,5,5,5))) +
  ylab("Relative Abundance") +
  scale_y_continuous(limits=c(0, .47)) 
#theme(axis.title.x=element_blank(),
#axis.text.x=element_blank()) +
#theme(plot.margin = unit(c(1,1,4,1), "lines")) +
#annotation_custom(text_1,xmin=1,xmax=1,ymin=.1,ymax=-.2) + 
#annotation_custom(text_42,xmin=2,xmax=2,ymin=.1,ymax=-.2) +
# annotation_custom(text_43,xmin=3,xmax=3,ymin=.1,ymax=-.2) + 
#annotation_custom(text_45,xmin=4,xmax=4,ymin=.1,ymax=-.2) +
#annotation_custom(text_56,xmin=5,xmax=5,ymin=.1,ymax=-.2) + 
#annotation_custom(text_84,xmin=6,xmax=6,ymin=.1,ymax=-.2) +
#coord_cartesian(clip="off") 
bac3


bac_ag <- aggregate(Abundance~ Treatment + Drought + Day + LandUse + Phylum, data=rk_phyla, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

bac_ag <- bac_ag[bac_ag$Phylum != 'Other',]

bac_ag <- bac_ag[bac_ag$Drought == 'Control',]    

bac_ag <- bac_ag[bac_ag$LandUse == 'Disturbed',]   

bac_ag <- do.call(data.frame, bac_ag)
bac_ag

bac_ag$Phylum <- as.character(bac_ag$Phylum)

bac_ag[bac_ag=="p__Acidobacteria"] <- "Acidobacteria"
bac_ag[bac_ag=="p__Actinobacteria"] <- "Actinobacteria"
bac_ag[bac_ag=="p__Nitrospirae"] <- "Nitrospirae"
bac_ag[bac_ag=="p__Bacteroidetes"] <- "Bacteroidetes"
bac_ag[bac_ag=="p__Proteobacteria"] <- "Proteobacteria"
bac_ag[bac_ag=="p__Verrucomicrobia"] <- "Verrucomicrobia"
bac_ag[bac_ag=="p__Chloroflexi"] <- "Chloroflexi"
bac_ag[bac_ag=="p__Planctomycetes"] <- "Planctomycetes"
bac_ag[bac_ag=="p__Firmicutes"] <- "Firmicutes"


bac_ag$Phylum <- as.factor(bac_ag$Phylum)

bac_ag$Day <- as.factor(bac_ag$Day)

bac_ag4 <- bac_ag

bac4 <- ggplot(bac_ag4, aes(x=bac_ag4$Day, y=bac_ag4$Abundance.mean, group=bac_ag4$Phylum)) + 
  geom_point(position=position_dodge(width=0.1),aes(color=bac_ag4$Phylum), size=2) +
  geom_line(position=position_dodge(width=0.1),aes(color=bac_ag4$Phylum),size=1) +
  #geom_errorbar(position=position_dodge(width=0.9),aes(ymin=resp2.2$Cmin.mean-resp2.2$se, ymax=resp2.2$Cmin.mean+resp2.2$se), width=3, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1",direction=-1) +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 12)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  ylab("Relative Abundance") +
  scale_y_continuous(limits=c(0, .47)) 
#theme(axis.title.x=element_blank(),
#axis.text.x=element_blank()) +
#theme(plot.margin = unit(c(1,1,4,1), "lines")) +
#annotation_custom(text_1,xmin=1,xmax=1,ymin=.1,ymax=-.2) + 
#annotation_custom(text_42,xmin=2,xmax=2,ymin=.1,ymax=-.2) +
# annotation_custom(text_43,xmin=3,xmax=3,ymin=.1,ymax=-.2) + 
#annotation_custom(text_45,xmin=4,xmax=4,ymin=.1,ymax=-.2) +
#annotation_custom(text_56,xmin=5,xmax=5,ymin=.1,ymax=-.2) + 
#annotation_custom(text_84,xmin=6,xmax=6,ymin=.1,ymax=-.2) +
#coord_cartesian(clip="off") 
bac4

aggregate(rk~Drought+Day,data=rk_phyla1, FUN=mean)
aggregate(rk~Drought+Day+LandUse,data=rk_phyla1, FUN=mean)

rk1 <- aggregate(rk~Treatment+LandUse+Drought+Day, data=rk_phyla1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(rk1)

rk1.1 <- do.call(data.frame, rk1)
rk1.1

rk1.1$se <- rk1.1$rk.sd / sqrt(rk1.1$rk.n)
head(rk1.1)

rk1.1$Treatment <- factor(rk1.1$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
rk1.1$Treatment
head(rk1.1)

library(grid)
text_1 <- textGrob("Day 1:\nPre-drought", gp=gpar(fontsize=15))
text_42 <- textGrob("Day 42:\nPost-drought", gp=gpar(fontsize=15))
text_43 <- textGrob("Day 43: 1 Day\nafter re-wet", gp=gpar(fontsize=15))
text_45 <- textGrob("Day 45: 3 Days\nafter re-wet", gp=gpar(fontsize=15))
text_56 <- textGrob("Day 56: 14 Days\nafter re-wet", gp=gpar(fontsize=15))
text_84 <- textGrob("Day 84: 42 Days\nafter re-wet", gp=gpar(fontsize=15))

jpeg(filename="drought_rk.jpeg", bg="transparent", res=600, units = "in", height=4.5, width=5.5) 


line_plot_rk <- ggplot(rk1.1, aes(x=rk1.1$Day, y=rk1.1$rk.mean, group=rk1.1$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=rk1.1$Treatment), size=3) +
  geom_line(position=position_dodge(width=0.25),aes(color=rk1.1$Treatment),size=2) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=rk1.1$rk.mean-rk1.1$se, ymax=rk1.1$rk.mean+rk1.1$se), width=.5, size=.75) +
  annotate('text', 2, 2.8, label="Drought: ~italic(P) < 0.001", size=4, parse=TRUE)+
  annotate('text', 2, 2.6, label="Disturbance: italic(P) == 0.014", size = 4, parse=TRUE)+
  annotate('text', 2, 2.4, label="Disturbance %*% Drought: italic(P) == 0.278", size=4, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.85,.85)) +
  theme(legend.text = element_text(color = "black", size = 10)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("16S  r:K")) +
  #theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank()) +
  #theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  #annotation_custom(text_1,xmin=1,xmax=1,ymin=-3.85,ymax=4.6) + 
  #annotation_custom(text_42,xmin=2,xmax=2,ymin=-3.85,ymax=4.6) +
  #annotation_custom(text_43,xmin=3,xmax=3,ymin=-3.85,ymax=4.6) + 
  #annotation_custom(text_45,xmin=4,xmax=4,ymin=-3.85,ymax=4.6) +
  #annotation_custom(text_56,xmin=5,xmax=5,ymin=-3.85,ymax=4.6) + 
  #annotation_custom(text_84,xmin=6,xmax=6,ymin=-3.85,ymax=4.6) +
  coord_cartesian(clip="off") 
plot(line_plot_rk)

dev.off()

#fungal classes

colnames(tax_table(ITS)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

fung_c <- transform_sample_counts(ITS, function(x) x/sum(x))

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

fung_classes <- aggregate(Abundance~Sample+Class+Treatment+Day+LandUse+Plot+Drought, dat_fung, FUN=sum)

fung_classes1 <- cast(fung_classes, Sample+Treatment+Day+LandUse+Drought+Plot ~ Class, value="Abundance")

fung_classes1$Day <- as.factor(fung_classes1$Day)

p1 <- lmer(c__Agaricomycetes~ Drought*LandUse*Day+(1|Plot), data=fung_classes1)
qqPlot(resid(p1))
Anova(p1)

aggregate(c__Agaricomycetes~Drought+Day, data=fung_classes1, FUN=mean)

p2 <- lmer(c__Eurotiomycetes~ Drought*LandUse*Day+(1|Plot), data=fung_classes1)
qqPlot(resid(p2))
Anova(p2)

p3 <- glmer(c__Geminibasidiomycetes ~Drought*LandUse*Day+(1|Plot), fung_classes1,family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(p3))
Anova(p3)

p4 <- glmer(c__Leotiomycetes~ Drought*LandUse*Day+(1|Plot), data=fung_classes1,family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(p4))
Anova(p4)

p5 <- glmer(c__Mortierellomycetes~ Drought*LandUse*Day+(1|Plot), data=fung_classes1,family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(p5))
Anova(p5)

p6 <- glmer(c__Mucoromycotina_cls_Incertae_sedis+.0001~ Drought*LandUse*Day+(1|Plot), data=fung_classes1,family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(p6))
Anova(p6)

p7 <- glmer(c__Sordariomycetes~ Drought*LandUse*Day+(1|Plot), data=fung_classes1,family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(p7))
Anova(p7)

p8 <- glmer(c__Tremellomycetes+.0001~ Drought*LandUse*Day+(1|Plot), data=fung_classes1,family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(p8))
Anova(p8)

p9 <- glmer(c__Umbelopsidomycetes+.00001~ Drought*LandUse*Day+(1|Plot), data=fung_classes1,family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(p9))
Anova(p9)

#fungal classes line plots

fung_ag <- aggregate(Abundance~ Treatment + Drought + Day + LandUse + Class, data=fung_classes, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

fung_ag <- fung_ag[fung_ag$Class != 'Other',]

fung_ag <- fung_ag[fung_ag$Drought == 'Drought',]    

fung_ag <- fung_ag[fung_ag$LandUse == 'Reference',]   

fung_ag <- do.call(data.frame, fung_ag)
fung_ag

fung_ag$Class <- as.character(fung_ag$Class)

fung_ag[fung_ag=="c__Agaricomycetes"] <- "Agaricomycetes*???"
fung_ag[fung_ag=="c__Eurotiomycetes"] <- "Eurotiomycetes*"
fung_ag[fung_ag=="c__Geminibasidiomycetes"] <- "Geminibasidiomycetes*"
fung_ag[fung_ag=="c__Leotiomycetes"] <- "Leotiomycetes*"
fung_ag[fung_ag=="c__Mortierellomycetes"] <- "Mortierellomycetes*"
fung_ag[fung_ag=="c__Mucoromycotina_cls_Incertae_sedis"] <- "Mucoromycotina*"
fung_ag[fung_ag=="c__Sordariomycetes"] <- "Sordariomycetes*???"
fung_ag[fung_ag=="c__Tremellomycetes"] <- "Tremellomycetes*"
fung_ag[fung_ag=="c__Umbelopsidomycetes"] <- "Umbelopsidomycetes*???"

fung_ag$Class <- as.factor(fung_ag$Class)

fung_ag$Day <- as.factor(fung_ag$Day)

fung_ag1 <- fung_ag

legend_fung <- 'Fungal Classes'

fung1 <- ggplot(fung_ag1, aes(x=fung_ag1$Day, y=fung_ag1$Abundance.mean, group=fung_ag1$Class)) + 
  geom_point(position=position_dodge(width=0.1),aes(color=fung_ag1$Class), size=2) +
  geom_line(position=position_dodge(width=0.1),aes(color=fung_ag1$Class),size=1) +
  #geom_errorbar(position=position_dodge(width=0.9),aes(ymin=resp2.2$Cmin.mean-resp2.2$se, ymax=resp2.2$Cmin.mean+resp2.2$se), width=3, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(legend_fung, palette="Set1") +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 12)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  ylab("Relative Abundance") +
scale_y_continuous(limits=c(0, .76)) 
#theme(axis.title.x=element_blank(),
#axis.text.x=element_blank()) +
#theme(plot.margin = unit(c(1,1,4,1), "lines")) +
#annotation_custom(text_1,xmin=1,xmax=1,ymin=.1,ymax=-.2) + 
#annotation_custom(text_42,xmin=2,xmax=2,ymin=.1,ymax=-.2) +
# annotation_custom(text_43,xmin=3,xmax=3,ymin=.1,ymax=-.2) + 
#annotation_custom(text_45,xmin=4,xmax=4,ymin=.1,ymax=-.2) +
#annotation_custom(text_56,xmin=5,xmax=5,ymin=.1,ymax=-.2) + 
#annotation_custom(text_84,xmin=6,xmax=6,ymin=.1,ymax=-.2) +
#coord_cartesian(clip="off") 
fung1

fung_ag <- aggregate(Abundance~ Treatment + Drought + Day + LandUse + Class, data=fung_classes, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

fung_ag <- fung_ag[fung_ag$Class != 'Other',]

fung_ag <- fung_ag[fung_ag$Drought == 'Drought',]    

fung_ag <- fung_ag[fung_ag$LandUse == 'Disturbed',]   

fung_ag <- do.call(data.frame, fung_ag)
fung_ag

fung_ag$Class <- as.character(fung_ag$Class)

fung_ag[fung_ag=="c__Agaricomycetes"] <- "Agaricomycetes*???"
fung_ag[fung_ag=="c__Eurotiomycetes"] <- "Eurotiomycetes*"
fung_ag[fung_ag=="c__Geminibasidiomycetes"] <- "Geminibasidiomycetes*"
fung_ag[fung_ag=="c__Leotiomycetes"] <- "Leotiomycetes*"
fung_ag[fung_ag=="c__Mortierellomycetes"] <- "Mortierellomycetes*"
fung_ag[fung_ag=="c__Mucoromycotina_cls_Incertae_sedis"] <- "Mucoromycotina incertae sedis*"
fung_ag[fung_ag=="c__Sordariomycetes"] <- "Sordariomycetes*???"
fung_ag[fung_ag=="c__Tremellomycetes"] <- "Tremellomycetes*"
fung_ag[fung_ag=="c__Umbelopsidomycetes"] <- "Umbelopsidomycetes*???"

fung_ag$Class <- as.factor(fung_ag$Class)

fung_ag$Day <- as.factor(fung_ag$Day)

fung_ag2 <- fung_ag

fung2<- ggplot(fung_ag2, aes(x=fung_ag2$Day, y=fung_ag2$Abundance.mean, group=fung_ag2$Class)) + 
  geom_point(position=position_dodge(width=0.1),aes(color=fung_ag2$Class), size=2) +
  geom_line(position=position_dodge(width=0.1),aes(color=fung_ag2$Class),size=1) +
  #geom_errorbar(position=position_dodge(width=0.9),aes(ymin=resp2.2$Cmin.mean-resp2.2$se, ymax=resp2.2$Cmin.mean+resp2.2$se), width=3, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 10)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  #theme(legend.margin=margin(c(1,5,5,5))) 
  ylab("Relative Abundance") +
  scale_y_continuous(limits=c(0, .76)) 
#theme(axis.title.x=element_blank(),
#axis.text.x=element_blank()) +
#theme(plot.margin = unit(c(1,1,4,1), "lines")) +
#annotation_custom(text_1,xmin=1,xmax=1,ymin=.1,ymax=-.2) + 
#annotation_custom(text_42,xmin=2,xmax=2,ymin=.1,ymax=-.2) +
# annotation_custom(text_43,xmin=3,xmax=3,ymin=.1,ymax=-.2) + 
#annotation_custom(text_45,xmin=4,xmax=4,ymin=.1,ymax=-.2) +
#annotation_custom(text_56,xmin=5,xmax=5,ymin=.1,ymax=-.2) + 
#annotation_custom(text_84,xmin=6,xmax=6,ymin=.1,ymax=-.2) +
#coord_cartesian(clip="off") 
fung2

fung_ag <- aggregate(Abundance~ Treatment + Drought + Day + LandUse + Class, data=fung_classes, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

fung_ag <- fung_ag[fung_ag$Class != 'Other',]

fung_ag <- fung_ag[fung_ag$Drought == 'Control',]    

fung_ag <- fung_ag[fung_ag$LandUse == 'Reference',]   

fung_ag <- do.call(data.frame, fung_ag)
fung_ag

fung_ag$Class <- as.character(fung_ag$Class)

fung_ag[fung_ag=="c__Agaricomycetes"] <- "Agaricomycetes"
fung_ag[fung_ag=="c__Eurotiomycetes"] <- "Eurotiomycetes"
fung_ag[fung_ag=="c__Geminibasidiomycetes"] <- "Geminibasidiomycetes"
fung_ag[fung_ag=="c__Leotiomycetes"] <- "Leotiomycetes"
fung_ag[fung_ag=="c__Mortierellomycetes"] <- "Mortierellomycetes"
fung_ag[fung_ag=="c__Mucoromycotina_cls_Incertae_sedis"] <- "Mucoromycotina"
fung_ag[fung_ag=="c__Sordariomycetes"] <- "Sordariomycetes"
fung_ag[fung_ag=="c__Tremellomycetes"] <- "Tremellomycetes"
fung_ag[fung_ag=="c__Umbelopsidomycetes"] <- "Umbelopsidomycetes"

fung_ag$Class <- as.factor(fung_ag$Class)

fung_ag$Day <- as.factor(fung_ag$Day)

fung_ag3 <- fung_ag

fung3 <- ggplot(fung_ag3, aes(x=fung_ag$Day, y=fung_ag3$Abundance.mean, group=fung_ag3$Class)) + 
  geom_point(position=position_dodge(width=0.1),aes(color=fung_ag3$Class), size=2) +
  geom_line(position=position_dodge(width=0.1),aes(color=fung_ag3$Class),size=1) +
  #geom_errorbar(position=position_dodge(width=0.9),aes(ymin=resp2.2$Cmin.mean-resp2.2$se, ymax=resp2.2$Cmin.mean+resp2.2$se), width=3, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(legend_fung, palette="Set1") +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 12)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  #theme(legend.margin=margin(c(1,5,5,5))) 
  ylab("Relative Abundance") 
  #scale_y_continuous(limits=c(0, .8)) 
#theme(axis.title.x=element_blank(),
#axis.text.x=element_blank()) +
#theme(plot.margin = unit(c(1,1,4,1), "lines")) +
#annotation_custom(text_1,xmin=1,xmax=1,ymin=.1,ymax=-.2) + 
#annotation_custom(text_42,xmin=2,xmax=2,ymin=.1,ymax=-.2) +
# annotation_custom(text_43,xmin=3,xmax=3,ymin=.1,ymax=-.2) + 
#annotation_custom(text_45,xmin=4,xmax=4,ymin=.1,ymax=-.2) +
#annotation_custom(text_56,xmin=5,xmax=5,ymin=.1,ymax=-.2) + 
#annotation_custom(text_84,xmin=6,xmax=6,ymin=.1,ymax=-.2) +
#coord_cartesian(clip="off") 
fung3

fung_ag <- aggregate(Abundance~ Treatment + Drought + Day + LandUse + Class, data=fung_classes, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

fung_ag <- fung_ag[fung_ag$Class != 'Other',]

fung_ag <- fung_ag[fung_ag$Drought == 'Control',]    

fung_ag <- fung_ag[fung_ag$LandUse == 'Disturbed',]   

fung_ag <- do.call(data.frame, fung_ag)
fung_ag

fung_ag$Class <- as.character(fung_ag$Class)

fung_ag[fung_ag=="c__Agaricomycetes"] <- "Agaricomycetes"
fung_ag[fung_ag=="c__Eurotiomycetes"] <- "Eurotiomycetes"
fung_ag[fung_ag=="c__Geminibasidiomycetes"] <- "Geminibasidiomycetes"
fung_ag[fung_ag=="c__Leotiomycetes"] <- "Leotiomycetes"
fung_ag[fung_ag=="c__Mortierellomycetes"] <- "Mortierellomycetes"
fung_ag[fung_ag=="c__Mucoromycotina_cls_Incertae_sedis"] <- "Mucoromycotina"
fung_ag[fung_ag=="c__Sordariomycetes"] <- "Sordariomycetes"
fung_ag[fung_ag=="c__Tremellomycetes"] <- "Tremellomycetes"
fung_ag[fung_ag=="c__Umbelopsidomycetes"] <- "Umbelopsidomycetes"

fung_ag$Class <- as.factor(fung_ag$Class)

fung_ag$Day <- as.factor(fung_ag$Day)

fung_ag4 <- fung_ag

fung4 <- ggplot(fung_ag4, aes(x=fung_ag$Day, y=fung_ag4$Abundance.mean, group=fung_ag4$Class)) + 
  geom_point(position=position_dodge(width=0.1),aes(color=fung_ag4$Class), size=2) +
  geom_line(position=position_dodge(width=0.1),aes(color=fung_ag4$Class),size=1) +
  #geom_errorbar(position=position_dodge(width=0.9),aes(ymin=resp2.2$Cmin.mean-resp2.2$se, ymax=resp2.2$Cmin.mean+resp2.2$se), width=3, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 12)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  #theme(legend.margin=margin(c(1,5,5,5))) 
  ylab("Relative Abundance") 
#scale_y_continuous(limits=c(0, .8)) 
#theme(axis.title.x=element_blank(),
#axis.text.x=element_blank()) +
#theme(plot.margin = unit(c(1,1,4,1), "lines")) +
#annotation_custom(text_1,xmin=1,xmax=1,ymin=.1,ymax=-.2) + 
#annotation_custom(text_42,xmin=2,xmax=2,ymin=.1,ymax=-.2) +
# annotation_custom(text_43,xmin=3,xmax=3,ymin=.1,ymax=-.2) + 
#annotation_custom(text_45,xmin=4,xmax=4,ymin=.1,ymax=-.2) +
#annotation_custom(text_56,xmin=5,xmax=5,ymin=.1,ymax=-.2) + 
#annotation_custom(text_84,xmin=6,xmax=6,ymin=.1,ymax=-.2) +
#coord_cartesian(clip="off") 
fung4


#Phylum, class, guilds multipanel figures


jpeg(filename="line_plot1.jpeg", bg="transparent", res=600, units = "in", height=10, width=7.5)


prow <- plot_grid(
  bac1 + theme(legend.position="none"),
  bac2 + theme(legend.position="none"),
  fung1 + theme(legend.position="none"),
  fung2+ theme(legend.position="none"),
  ncol = 2,
  align = 'hv',
  labels = c("A             Reference - Drought", "B             Disturbed - Drought", "C             Reference - Drought", "D             Disturbed - Drought"),
  label_size=16,
  hjust=0,
  label_y=1.025, scale = 0.9
)
prow

dev.off()

jpeg(filename="line_plot2.jpeg", bg="transparent", res=600, units = "in", height=10, width=7.5)


prow <- plot_grid(
  bac3 + theme(legend.position="none"),
  bac4 + theme(legend.position="none"),
  fung3 + theme(legend.position="none"),
  fung4+ theme(legend.position="none"),
  ncol = 2,
  align = 'hv',
  labels = c("A             Reference - Control", "B             Disturbed - Control", "C             Reference - Control", "D             Disturbed - Control"),
  label_size=16,
  hjust=0,
  label_y=1.025, scale = 0.9
)
prow

dev.off()

legend1 <- get_legend(bac1 + theme(legend.box.margin = margin(0, 0, 0, 2)))
legend2 <- get_legend(fung1 + theme(legend.box.margin = margin(0, 0, 0, 2)))


jpeg(filename="legends1.jpeg", bg="transparent", res=600, units = "in", height=10, width=4)


leg <- plot_grid(legend1, legend2,
              
                    ncol = 1)
leg

dev.off()

legend4 <- get_legend(bac3 + theme(legend.box.margin = margin(0, 0, 0, 2)))
legend5 <- get_legend(fung3 + theme(legend.box.margin = margin(0, 0, 0, 2)))

jpeg(filename="legends2.jpeg", bg="transparent", res=600, units = "in", height=10, width=4)


leg <- plot_grid(legend4, legend5,
                 
                 ncol = 1)
leg

dev.off()

