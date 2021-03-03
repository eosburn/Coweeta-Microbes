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

crp <- read.csv("C:/Users/ernie/OneDrive/Desktop/FD_CRP/CRP.csv")

Bac_otus <- "C:/Users/ernie/OneDrive/Desktop/Data/Chapter 2 Data/FD-16s-table-with-taxonomy-final.biom"

Fung_otus <- "C:/Users/ernie/OneDrive/Desktop/Data/Chapter 2 Data/FD-ITS-table-with-taxonomy-final.biom"

x1 <- read_biom(Bac_otus)

x2 <- read_biom(Fung_otus)

meta <- "C:/Users/ernie/OneDrive/Desktop/Data/Chapter 2 Data/FD_16s_metadata.txt"

meta2 <- import_qiime_sample_data(meta)

OTU1 <- import_biom(x1)

OTU2 <- import_biom(x2)

OTU1 <- merge_phyloseq(OTU1, meta2)

OTU2 <- merge_phyloseq(OTU2, meta2)

FD_ITS <- rarefy_even_depth(OTU2, rngseed=TRUE)

FD_16s <- rarefy_even_depth(OTU1, rngseed=TRUE)

#FD_16s <- transform_sample_counts(FD_16s, function(x) x/sum(x))

a1 <- t(data.frame(otu_table(FD_16s)))

a2 <- t(data.frame(otu_table(FD_ITS)))

a1 <- a1[row.names(a1) != "FD39",]

a1 <- a1[row.names(a1) != "FD47",]

a2 <- a2[row.names(a2) != "FD39",]

a2 <- a2[row.names(a2) != "FD47",]

b1 <- vegdist(a1, method = "bray")

b2 <- vegdist(a2, method = "bray")

library(vegan)

crp2 <- na.omit(crp)

crp2 <- crp2[order(crp2$id),]

soil_dist <- vegdist(crp2[,c(22:32)], method="euclidean")

crp_dist <- vegdist(crp2[,c(37:43)], method="euclidean")

mantel(soil_dist,crp_dist)

mantel(crp_dist, b1)

mantel(crp_dist, b2)

crp_dist2 <- vegdist(crp2[,c(3:9)], method="euclidean")

mantel(soil_dist,crp_dist2)

mantel(crp_dist2, b1)

mantel(crp_dist2, b2)

crp_dist3 <- vegdist(crp2[,c(11:17)], method="euclidean")

mantel(soil_dist,crp_dist3)

mantel(crp_dist3, b1)

mantel(crp_dist3, b2)

#Multivariate analysis of carbon profiles per mg microbial biomass

names(crp2)[names(crp2) == "DOC.TDN"] <- "DOC:TDN"

crp2.1 <- data.frame(scale(crp2[,c(37:43)]))

crp2.1$Use <- crp2$Use

adonis2(crp2.1[,c(1:7)] ~ Use, data=crp2.1, permutations = 999, method = "euclidean",na.rm=TRUE)

a <- vegdist(crp2.1[,c(1:7)], method="euclidean")

b <- cmdscale(a, k=2, eig=TRUE)

b$points

vec<-envfit(b, crp2[,c(22,23,24,26,27,28)], perm=1000, na.rm=TRUE)
vec

soil.scrs <- as.data.frame(scores(vec, display = "vectors"))
soil.scrs <- cbind(soil.scrs, soil = rownames(soil.scrs))
soil.scrs

data.scores1 <- as.data.frame(b$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores1$site <- rownames(data.scores1)  # create a column of site names, from the rownames of data.scores
data.scores1$Treatment <- crp2$Use 
data.scores1$pair <- crp2$Pair
head(data.scores1) 

var1 <- round(b$eig[1]/sum(b$eig)*100,1)
var2 <- round(b$eig[2]/sum(b$eig)*100,1)

library(ggplot2)
library(tidyverse)

ord1 <- ggplot(data.scores1, aes(x=V1, y=V2)) + 
  geom_point(data=data.scores1,aes(x=V1,y=V2,shape=data.scores1$pair,color=data.scores1$Treatment),size=3) + # add the point markers
  stat_ellipse(data = data.scores1, level = 0.95, geom = "polygon", alpha = 0.1, aes(fill = Treatment),type = "norm") +
  scale_color_manual(values=c("#104e8b","darkgoldenrod1")) +
  scale_fill_manual(values=c("#104e8b","darkgoldenrod1")) +
  annotate("text", x=-3,y=-3,label="Disturbance: ~italic(P) == 0.001", parse=TRUE, size=4, fontface="bold")+
  labs(color='Treatment') +
  labs(shape='Watershed Pair') +
  theme(legend.text=element_text(size=12)) +
  xlab("PCA 1: 60.9%") +
  ylab("PCA 2: 19.7%") +
  theme_classic() +
  theme(axis.text=element_text(size=12)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=12)) +
  ggtitle(expression("CRPs per mg microbial biomass C")) +
  geom_segment(data=soil.scrs,aes(x=0,xend=4*Dim1,y=0,yend=4*Dim2), arrow = arrow(length = unit(0.25, "cm")),colour="black",size=.25,inherit.aes=FALSE) +
  geom_text(data = soil.scrs, aes(x = 4.5*Dim1, y = 4.5*Dim2, label = soil),size = 3.5)
ord1

#Multivariate analysis of carbon profiles per g dry soil

crp2.2 <- data.frame(scale(crp2[,c(3:9)]))

crp2.2$Use <- crp2$Use

adonis2(crp2.2[,c(1:7)] ~ Use, data=crp2.2, permutations = 999, method = "euclidean",na.rm=TRUE)

a2 <- vegdist(crp2.2[,c(1:7)], method="euclidean")

b2 <- cmdscale(a2, k=2, eig=TRUE)

Use <- crp2.2$Use

disp <- betadisper(a2,Use,type="centroid")
anova(disp)

b2$points

vec2<-envfit(b2, crp2[,c(22,23,24,26,27,28)], perm=1000, na.rm=TRUE)
vec2

soil.scrs2 <- as.data.frame(scores(vec2, display = "vectors"))
soil.scrs2 <- cbind(soil.scrs2, soil = rownames(soil.scrs2))
soil.scrs2

data.scores2 <- as.data.frame(b2$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores2$site <- rownames(data.scores2)  # create a column of site names, from the rownames of data.scores
data.scores2$Treatment <- crp2$Use 
data.scores2$pair <- crp2$Pair
head(data.scores2) 

var1 <- round(b2$eig[1]/sum(b2$eig)*100,1)
var2 <- round(b2$eig[2]/sum(b2$eig)*100,1)

library(ggplot2)

ord2 <- ggplot(data.scores2, aes(x=V1, y=V2)) + 
  geom_point(data=data.scores2,aes(x=V1,y=V2,shape=data.scores2$pair,color=data.scores2$Treatment),size=3) + # add the point markers
  stat_ellipse(data = data.scores2, level = 0.95, geom = "polygon", alpha = 0.1, aes(fill = Treatment),type = "norm") +
  scale_color_manual(values=c("#104e8b","darkgoldenrod1")) +
  scale_fill_manual(values=c("#104e8b","darkgoldenrod1")) +
  annotate("text", -4, -3.5, label="Disturbance: ~italic(P) == 0.269", parse=TRUE, size=4, fontface="bold")+
  labs(color='Treatment') +
  labs(shape='Watershed Pair') +
  theme(legend.text=element_text(size=12)) +
  xlab("PCA 1: 62.3%") +
  ylab("PCA 2: 17.3%") +
  theme_classic() +
  theme(axis.text=element_text(size=12)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(expression("CRPs per g dry soil")) +
  theme(text=element_text(size=12)) +
  geom_segment(data=soil.scrs2,aes(x=0,xend=5*Dim1,y=0,yend=5*Dim2), arrow = arrow(length = unit(0.25, "cm")),colour="black",size=.25,inherit.aes=FALSE) +
  geom_text(data = soil.scrs2, aes(x = 5.5*Dim1, y = 5.5*Dim2, label = soil),size = 3.5)
ord2

#Multivariate analysis of carbon profiles as percentages

crp2.3 <- data.frame(scale(crp2[,c(11:17)]))

crp2.3$Use <- crp2$Use

adonis2(crp2.2[,c(1:7)] ~ Use, data=crp2.2, permutations = 999, method = "euclidean",na.rm=TRUE)

a3 <- vegdist(crp2.3[,c(1:7)], method="euclidean")

b3 <- cmdscale(a3, k=2, eig=TRUE)

b3$points

vec3<-envfit(b3, crp2[,c(23,24,26,27,28)], perm=1000, na.rm=TRUE)
vec3

soil.scrs3 <- as.data.frame(scores(vec3, display = "vectors"))
soil.scrs3 <- cbind(soil.scrs3, soil = rownames(soil.scrs3))
soil.scrs3

data.scores3 <- as.data.frame(b3$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores3$site <- rownames(data.scores3)  # create a column of site names, from the rownames of data.scores
data.scores3$Treatment <- crp2$Use 
data.scores3$pair <- crp2$Pair
head(data.scores3) 

var1 <- round(b3$eig[1]/sum(b3$eig)*100,1)
var2 <- round(b3$eig[2]/sum(b3$eig)*100,1)

library(ggplot2)

ord3 <- ggplot(data.scores3, aes(x=V1, y=V2)) + 
  geom_point(data=data.scores3,aes(x=V1,y=V2,shape=data.scores3$pair,color=data.scores3$Treatment),size=3) + # add the point markers
  stat_ellipse(data = data.scores3, level = 0.95, geom = "polygon", alpha = 0.1, aes(fill = Treatment),type = "norm") +
  scale_color_manual(values=c("#104e8b","darkgoldenrod1")) +
  scale_fill_manual(values=c("#104e8b","darkgoldenrod1")) +
  annotate("text", -2.5,-3.5, label="Disturbance: ~italic(P) == 0.283", parse=TRUE, size=4, fontface="bold")+
  labs(color='Treatment') +
  labs(shape='Watershed Pair') +
  theme(legend.text=element_text(size=12)) +
  xlab("PCA 1: 37.6%") +
  ylab("PCA 2: 21.8%") +
  theme_classic() +
  theme(axis.text=element_text(size=12)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=12)) +
  ggtitle(expression("CRPs as percentages")) +
  geom_segment(data=soil.scrs3,aes(x=0,xend=5*Dim1,y=0,yend=5*Dim2), arrow = arrow(length = unit(0.25, "cm")),colour="black",size=.25,inherit.aes=FALSE) +
  geom_text(data = soil.scrs3, aes(x = 5.7*Dim1, y = 5.7*Dim2, label = soil),size = 3.5)
ord3

#Multipanel ordination figure

setwd("C:/Users/ernie/OneDrive/Desktop/FD_CRP")

jpeg(filename="fig1.jpeg", bg="transparent", res=500, units = "in", height=10, width=6) 

f1 <- plot_grid(ord2, ord3, ord1, ncol = 1, labels=c('A', 'B', 'C'),align="hv", label_size=30)

f1

dev.off()


#Univariate analysis of CRPs per gram dry soil

library(lme4)

library(car)

m1 <- lmer(Glucose ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m1))
Anova(m1)
m1.1 <- aggregate(Glucose ~ Use, data=crp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m1.1

m1.2 <- do.call(data.frame, m1.1)

m1.2$se <- m1.2$Glucose.sd / sqrt(m1.2$Glucose.n)
m1.2

m2 <- lmer(Cellulose ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m2))
m2 <- glmer(Cellulose ~ Use + (1|Pair), data=crp,family=Gamma(link=log))
Anova(m2)
m2.1 <- aggregate(Cellulose ~ Use, data=crp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m2.1

m2.2 <- do.call(data.frame, m2.1)

m2.2$se <- m2.2$Cellulose.sd / sqrt(m2.2$Cellulose.n)
m2.2

m3 <- lmer(Lignin ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m3))
Anova(m3)
m3.1 <- aggregate(Lignin ~ Use, data=crp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m3.1

m3.2 <- do.call(data.frame, m3.1)

m3.2$se <- m3.2$Lignin.sd / sqrt(m3.2$Lignin.n)
m3.2

m4 <- lmer(Chitin ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m4))
Anova(m4)
m4.1 <- aggregate(Chitin ~ Use, data=crp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m4.1

m4.2 <- do.call(data.frame, m4.1)

m4.2$se <- m4.2$Chitin.sd / sqrt(m4.2$Chitin.n)
m4.2

m5 <- lmer(Glycine ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m5))
Anova(m5)
m5.1 <- aggregate(Glycine ~ Use, data=crp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m5.1

m5.2 <- do.call(data.frame, m5.1)

m5.2$se <- m5.2$Glycine.sd / sqrt(m5.2$Glycine.n)
m5.2

m6 <- lmer(Oxalic ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m6))
m6 <- glmer(Oxalic ~ Use + (1|Pair), data=crp,family=Gamma(link=log))
Anova(m6)
m6.1 <- aggregate(Oxalic ~ Use, data=crp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m6.1

m6.2 <- do.call(data.frame, m6.1)

m6.2$se <- m6.2$Oxalic.sd / sqrt(m6.2$Oxalic.n)
m6.2

m7 <- lmer(Yeast ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m7))
Anova(m7)
m7.1 <- aggregate(Yeast ~ Use, data=crp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m7.1

m7.2 <- do.call(data.frame, m7.1)

m7.2$se <- m7.2$Yeast.sd / sqrt(m7.2$Yeast.n)
m7.2

m8 <- lmer(Total ~ Use + (1|Pair), data=crp2)
shapiro.test(resid(m8))
Anova(m8)
m8.1 <- aggregate(Total ~ Use, data=crp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m8.1

m8.2 <- do.call(data.frame, m8.1)

m8.2$se <- m8.2$Total.sd / sqrt(m8.2$Total.n)
m8.2



crp.1 <- crp[,c(1:9,21)]

library(reshape2)
crp.2 <- melt(crp.1, id.vars = c("Soil","id", "Use"))

carpag <- aggregate(value ~ variable +Use, data=crp.2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
carpag

carpag <- do.call(data.frame, carpag)

carpag$se <- carpag$value.sd / sqrt(carpag$value.n)
carpag

carpag$variable <- as.character(carpag$variable)

carpag[carpag=="Oxalic"] <- "Oxalic Acid"

carpag[carpag=="Yeast"] <- "Autolyzed Yeast"

carpag$variable <- as.factor((carpag$variable))

levels(carpag$variable) <- gsub(" ", "\n", levels(carpag$variable))

crp_plot1 <- ggplot(carpag, aes(x=carpag$variable, y=carpag$value.mean, fill=carpag$Use)) + 
  geom_errorbar(aes(ymin=carpag$value.mean-carpag$se, ymax=carpag$value.mean+carpag$se), width=0.2, size=1, position=position_dodge(0.5)) +
  stat_summary(fun.y=mean,geom="point",shape=21,size=5,color="black", position=position_dodge(0.5))+
  annotate('text',5,.9 , label="**", size=8)+
  annotate('text', 1,1.8 , label="*", size=8)+
  theme_classic() +
  scale_fill_manual(values=c("#104e8b","darkgoldenrod1")) +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.title=element_blank()) +
  theme(legend.position=c(.85,.9)) +
  xlab("") +
  ggtitle("Catabolic Responses per gram dry soil") +
  scale_y_continuous(bquote(~mu*'g'~'CO'[2]~'-C'~~'gdw'^-1~~h^-1))
plot(crp_plot1)

#Univariate analysis of CRPs as a percentage

m1 <- lmer(X.Glucose ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m1))
Anova(m1)
m1.1 <- aggregate(X.Glucose ~ Use, data=crp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m1.1

m1.2 <- do.call(data.frame, m1.1)

m1.2$se <- m1.2$X.Glucose.sd / sqrt(m1.2$X.Glucose.n)
m1.2

m2 <- lmer(X.Cellulose ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m2))
m2 <- glmer(X.Cellulose ~ Use + (1|Pair), data=crp,family=Gamma(link=log))
Anova(m2)
m2.1 <- aggregate(X.Cellulose ~ Use, data=crp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m2.1

m2.2 <- do.call(data.frame, m2.1)

m2.2$se <- m2.2$X.Cellulose.sd / sqrt(m2.2$X.Cellulose.n)
m2.2

m3 <- lmer(X.lignin ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m3))
Anova(m3)
m3.1 <- aggregate(X.lignin ~ Use, data=crp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m3.1

m3.2 <- do.call(data.frame, m3.1)

m3.2$se <- m3.2$X.lignin.sd / sqrt(m3.2$X.lignin.n)
m3.2

m4 <- lmer(X.Chitin ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m4))
Anova(m4)
m4.1 <- aggregate(X.Chitin ~ Use, data=crp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m4.1

m4.2 <- do.call(data.frame, m4.1)

m4.2$se <- m4.2$X.Chitin.sd / sqrt(m4.2$X.Chitin.n)
m4.2

m5 <- lmer(X.Glycine ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m5))
Anova(m5)
m5.1 <- aggregate(X.Glycine ~ Use, data=crp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m5.1

m5.2 <- do.call(data.frame, m5.1)

m5.2$se <- m5.2$X.Glycine.sd / sqrt(m5.2$X.Glycine.n)
m5.2

m6 <- lmer(X.Oxalic ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m6))
Anova(m6)
m6.1 <- aggregate(X.Oxalic ~ Use, data=crp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m6.1

m6.2 <- do.call(data.frame, m6.1)

m6.2$se <- m6.2$X.Oxalic.sd / sqrt(m6.2$X.Oxalic.n)
m6.2

m7 <- lmer(X.Yeast ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m7))
Anova(m7)
m7.1 <- aggregate(X.Yeast ~ Use, data=crp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
m7.1

m7.2 <- do.call(data.frame, m7.1)

m7.2$se <- m7.2$X.Yeast.sd / sqrt(m7.2$X.Yeast.n)
m7.2

crp.11 <- crp[,c(1:2,11:17,21)]

library(reshape2)
crp.12 <- melt(crp.11, id.vars = c("Soil","id", "Use"))

carpag2 <- aggregate(value ~ variable +Use, data=crp.12, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
carpag2

carpag2 <- do.call(data.frame, carpag2)

carpag2$se <- carpag2$value.sd / sqrt(carpag2$value.n)
carpag2

carpag2$variable <- as.character(carpag2$variable)

carpag2[carpag2=="X.Oxalic"] <- "Oxalic Acid"
carpag2[carpag2=="X.Yeast"] <- "Autolyzed Yeast"
carpag2[carpag2=="X.Glucose"] <- "Glucose"
carpag2[carpag2=="X.Glycine"] <- "Glycine"
carpag2[carpag2=="X.Chitin"] <- "Chitin"
carpag2[carpag2=="X.lignin"] <- "Lignin"
carpag2[carpag2=="X.Cellulose"] <- "Cellulose"

carpag2$variable <- as.factor((carpag2$variable))

levels(carpag2$variable) <- gsub(" ", "\n", levels(carpag2$variable))

crp_plot2 <- ggplot(carpag2, aes(x=carpag2$variable, y=carpag2$value.mean, fill=carpag2$Use)) + 
  geom_errorbar(aes(ymin=carpag2$value.mean-carpag2$se, ymax=carpag2$value.mean+carpag2$se), width=0.2, size=1, position=position_dodge(0.5)) +
  stat_summary(fun.y=mean,geom="point",shape=21,size=5,color="black", position=position_dodge(0.5))+
  annotate('text',5,.125 , label="*", size=8)+
  theme_classic() +
  scale_fill_manual(values=c("#104e8b","darkgoldenrod1")) +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.title=element_blank()) +
  theme(legend.position=c(.85,.9)) +
  xlab("") +
  ggtitle("Catabolic Responses as percentages") +
  scale_y_continuous("% of total C respired")
plot(crp_plot2)

#Univariate analysis of CRPs per mg microbial biomass C

m1 <- lmer(Glucose2~ Use + (1|Pair), data=crp)
shapiro.test(resid(m1))
Anova(m1)
m1.1 <- aggregate(Glucose2 ~ Use, data=crp, FUN=mean)
m1.1

glu_box <-  ggplot(crp2, aes(x=crp2$Use, y=crp2$Glucose2, fill=crp2$Use)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.size=.25) +
  annotate('text', 1.5, 16.5, label="Disturbance: italic(P) == 0.003", size = 4, parse=TRUE)+
  theme_classic() +
  geom_point(position=position_jitter(.25),size=2) +
  stat_summary(fun.y=mean,geom="point",shape=23,size=3,color="black",fill="white")+
  scale_fill_manual(values=c("#104e8b","darkgoldenrod1")) +
  theme(axis.title=element_text(size=12)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 18)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  theme(legend.position="none") +
  ggtitle("Glucose") +
  scale_y_continuous(bquote(~mu*'g'~'CO'[2]~'-C'~~'mg'~'Microbial'~'Biomass'~'C'^-1~~h^-1))+
  theme(axis.title.x=element_blank())
glu_box

m2 <- lmer(Cellulose2 ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m2))
Anova(m2)
m2.1 <- aggregate(Cellulose2 ~ Use, data=crp, FUN=mean)
m2.1

cell_box <-  ggplot(crp2, aes(x=crp2$Use, y=crp2$Cellulose2, fill=crp2$Use)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.size=.25) +
  annotate('text', 1.5, 6.25, label="Disturbance: italic(P) == 0.01", size = 4, parse=TRUE)+
  theme_classic() +
  geom_point(position=position_jitter(.25),size=2) +
  stat_summary(fun.y=mean,geom="point",shape=23,size=3,color="black",fill="white")+
  scale_fill_manual(values=c("#104e8b","darkgoldenrod1")) +
  theme(axis.title=element_text(size=12)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 18)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  theme(legend.position="none") +
  ggtitle("Cellulose") +
  scale_y_continuous(bquote(~mu*'g'~'CO'[2]~'-C'~~'mg'~'Microbial'~'Biomass'~'C'^-1~~h^-1))+
  theme(axis.title.x=element_blank())
cell_box


m3 <- lmer(Lignin2 ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m3))
m3 <- glmer(Lignin2 ~ Use + (1|Pair), data=crp,family=Gamma(link=log))
Anova(m3)
m3.1 <- aggregate(Lignin2 ~ Use, data=crp, FUN=mean)
m3.1

lig_box <-  ggplot(crp2, aes(x=crp2$Use, y=crp2$Lignin2, fill=crp2$Use)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.shape=NA) +
  annotate('text', 1.5, 15.5, label="Disturbance: italic(P) == 0.002", size = 4, parse=TRUE)+
  theme_classic() +
  geom_point(position=position_jitter(.25),size=2) +
  stat_summary(fun.y=mean,geom="point",shape=23,size=3,color="black",fill="white")+
  scale_fill_manual(values=c("#104e8b","darkgoldenrod1")) +
  theme(axis.title=element_text(size=12)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 18)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  theme(legend.position="none") +
  ggtitle("Lignin") +
  scale_y_continuous(bquote(~mu*'g'~'CO'[2]~'-C'~~'mg'~'Microbial'~'Biomass'~'C'^-1~~h^-1))+
  theme(axis.title.x=element_blank())
lig_box

m4 <- lmer(Chitin2 ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m4))
Anova(m4)
m4.1 <- aggregate(Chitin2 ~ Use, data=crp, FUN=mean)
m4.1

chi_box <-  ggplot(crp2, aes(x=crp2$Use, y=crp2$Chitin2, fill=crp2$Use)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.shape=NA) +
  annotate('text', 1.5, 5, label="Disturbance: italic(P) == 0.057", size = 4, parse=TRUE)+
  theme_classic() +
  geom_point(position=position_jitter(.25),size=2) +
  stat_summary(fun.y=mean,geom="point",shape=23,size=3,color="black",fill="white")+
  scale_fill_manual(values=c("#104e8b","darkgoldenrod1")) +
  theme(axis.title=element_text(size=12)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 18)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  theme(legend.position="none") +
  ggtitle("Chitin") +
  scale_y_continuous(bquote(~mu*'g'~'CO'[2]~'-C'~~'mg'~'Microbial'~'Biomass'~'C'^-1~~h^-1))+
  theme(axis.title.x=element_blank())
chi_box

m5 <- lmer(Glycine2 ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m5))
Anova(m5)
m5.1 <- aggregate(Glycine2 ~ Use, data=crp, FUN=mean)
m5.1

gly_box <-  ggplot(crp2, aes(x=crp2$Use, y=crp2$Glycine2, fill=crp2$Use)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.size=.25) +
  annotate('text', 1.5, 7, label="Disturbance: italic(P) < 0.001", size = 4, parse=TRUE)+
  theme_classic() +
  geom_point(position=position_jitter(.25),size=2) +
  stat_summary(fun.y=mean,geom="point",shape=23,size=3,color="black",fill="white")+
  scale_fill_manual(values=c("#104e8b","darkgoldenrod1")) +
  theme(axis.title=element_text(size=12)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 18)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  theme(legend.position="none") +
  ggtitle("Glycine") +
  scale_y_continuous(bquote(~mu*'g'~'CO'[2]~'-C'~~'mg'~'Microbial'~'Biomass'~'C'^-1~~h^-1))+
  theme(axis.title.x=element_blank())
gly_box

m6 <- lmer(Oxalic2 ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m6))
Anova(m6)
m6.1 <- aggregate(Oxalic2 ~ Use, data=crp, FUN=mean)
m6.1

ox_box <-  ggplot(crp2, aes(x=crp2$Use, y=crp2$Oxalic2, fill=crp2$Use)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.shape=NA) +
  annotate('text', 1.5, 13.5, label="Disturbance: italic(P) == 0.19", size = 4, parse=TRUE)+
  theme_classic() +
  geom_point(position=position_jitter(.25),size=2) +
  stat_summary(fun.y=mean,geom="point",shape=23,size=3,color="black",fill="white")+
  scale_fill_manual(values=c("#104e8b","darkgoldenrod1")) +
  theme(axis.title=element_text(size=12)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 18)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  theme(legend.position="none") +
  ggtitle("Oxalic Acid") +
  scale_y_continuous(bquote(~mu*'g'~'CO'[2]~'-C'~~'mg'~'Microbial'~'Biomass'~'C'^-1~~h^-1))+
  theme(axis.title.x=element_blank())
ox_box

m7 <- lmer(Yeast2 ~ Use + (1|Pair), data=crp)
shapiro.test(resid(m7))
Anova(m7)
m7.1 <- aggregate(Yeast2 ~ Use, data=crp, FUN=mean)
m7.1

yeast_box <-  ggplot(crp2, aes(x=crp2$Use, y=crp2$Yeast2, fill=crp2$Use)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.size=.25) +
  annotate('text', 1.5, 16, label="Disturbance: italic(P) < 0.001", size = 4, parse=TRUE)+
  theme_classic() +
  geom_point(position=position_jitter(.25),size=2) +
  stat_summary(fun.y=mean,geom="point",shape=23,size=3,color="black",fill="white")+
  scale_fill_manual(values=c("#104e8b","darkgoldenrod1")) +
  theme(axis.title=element_text(size=12)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 18)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  theme(legend.position="none") +
  ggtitle("Autolyzed Yeast") +
  scale_y_continuous(bquote(~mu*'g'~'CO'[2]~'-C'~~'mg'~'Microbial'~'Biomass'~'C'^-1~~h^-1))+
  theme(axis.title.x=element_blank())
yeast_box

crp2$Total2 <- rowSums(crp2[,c(37:43)])

m8 <- lmer(Total2 ~ Use + (1|Pair), data=crp2)
shapiro.test(resid(m8))
m8 <- glmer(Total2~Use + (1|Pair),data=crp2, family=Gamma(link=log))
Anova(m8)

total_box <-  ggplot(crp2, aes(x=crp2$Use, y=crp2$Total2, fill=crp2$Use)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.shape=NA) +
  annotate('text', 1.5, 63, label="Disturbance: italic(P) < 0.001", size = 4, parse=TRUE)+
  theme_classic() +
  geom_point(position=position_jitter(.25),size=2) +
  stat_summary(fun.y=mean,geom="point",shape=23,size=3,color="black",fill="white")+
  scale_fill_manual(values=c("#104e8b","darkgoldenrod1")) +
  theme(axis.title=element_text(size=12)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 18)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  theme(legend.position="none") +
  ggtitle("Total Response") +
  scale_y_continuous(bquote(~mu*'g'~'CO'[2]~'-C'~~'mg'~'Microbial'~'Biomass'~'C'^-1~~h^-1))+
  theme(axis.title.x=element_blank())
total_box

setwd("C:/Users/ernie/Desktop/Manuscripts/FD_CRP")

jpeg(filename="fig2.jpeg", bg="transparent", res=500, units = "in", height=8, width=9.75) 

f2 <- plot_grid(glu_box, cell_box, lig_box,chi_box,gly_box, ox_box,yeast_box,total_box , ncol = 4, labels=c('A', 'B', 'C','D','E','F','G','H'),align="hv", label_size=20)

f2

dev.off()

crp.21 <- crp[,c(1:2,37:43,21)]

library(reshape2)
crp.22 <- melt(crp.21, id.vars = c("Soil","id", "Use"))

carpag3 <- aggregate(value ~ variable +Use, data=crp.22, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
carpag3

carpag3 <- do.call(data.frame, carpag3)

carpag3$se <- carpag3$value.sd / sqrt(carpag3$value.n)
carpag3

carpag3$variable <- as.character(carpag3$variable)

carpag3[carpag3=="Oxalic2"] <- "Oxalic Acid"
carpag3[carpag3=="Yeast2"] <- "Autolyzed Yeast"
carpag3[carpag3=="Glucose2"] <- "Glucose"
carpag3[carpag3=="Glycine2"] <- "Glycine"
carpag3[carpag3=="Chitin2"] <- "Chitin"
carpag3[carpag3=="Lignin2"] <- "Lignin"
carpag3[carpag3=="Cellulose2"] <- "Cellulose"

carpag3$variable <- as.factor((carpag3$variable))

levels(carpag3$variable) <- gsub(" ", "\n", levels(carpag3$variable))

crp_plot3 <- ggplot(carpag3, aes(x=carpag3$variable, y=carpag3$value.mean, fill=carpag3$Use)) + 
  geom_errorbar(aes(ymin=carpag3$value.mean-carpag3$se, ymax=carpag3$value.mean+carpag3$se), width=0.2, size=1, position=position_dodge(0.5)) +
  stat_summary(fun.y=mean,geom="point",shape=21,size=5,color="black", position=position_dodge(0.5))+
  annotate('text',1,10.5 , label="***", size=8)+
  annotate('text',2,4 , label="*", size=8)+
  annotate('text',3,3 , label="???", size=6)+
  annotate('text',4,11 , label="**", size=8)+
  annotate('text',5,5 , label="***", size=8)+
  annotate('text',6,8.3 , label="**", size=8)+
  theme_classic() +
  scale_fill_manual(values=c("#104e8b","darkgoldenrod1")) +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.title=element_blank()) +
  theme(legend.position=c(.85,.9)) +
  xlab("") +
  ggtitle("Catabolic Responses per mg microbial biomass C") +
  scale_y_continuous(bquote(~mu*'g'~'CO'[2]~'-C'~~'mg'~'MBC'^-1~~h^-1))
plot(crp_plot3)

setwd("C:/Users/ernie/OneDrive/Desktop/FD_CRP")

jpeg(filename="fig2.jpeg", bg="transparent", res=500, units = "in", height=15, width=8.5) 

f2 <- plot_grid(crp_plot1, crp_plot2, crp_plot3, ncol = 1, labels=c('A', 'B', 'C'),align="hv", label_size=35)
f2

dev.off()
