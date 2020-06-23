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

setwd("C:/Users/ernie/Desktop/FD_isotopes")

soil <- read.csv("C:/Users/ernie/Desktop/FD_isotopes/FD_microbesandprocesses.csv")

Bac_otus <- "C:/Users/ernie/Desktop/Data/Chapter 2 Data/FD-16s-table-with-taxonomy-final.biom"

Fung_otus <- "C:/Users/ernie/Desktop/Data/Chapter 2 Data/FD-ITS-table-with-taxonomy-final.biom"

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

FD_16s <- rarefy_even_depth(OTU1, rngseed=TRUE)

a1 <- t(data.frame(otu_table(FD_16s)))

a2 <- t(data.frame(otu_table(FD_ITS)))

b1 <- vegdist(a1, method = "bray")

b2 <- vegdist(a2, method = "bray")

sam_data16s <- data.frame(sample_data(FD_16s))

sam_dataITS <- data.frame(sample_data(FD_ITS))

d1 <- cbind(a1, sam_data16s)

d2 <- cbind(a2, sam_dataITS)

#Soil Chemistry 

sc1 <- lmer(pH~ Use+(1|Pair), data=soil)
shapiro.test(resid(sc1))
Anova(sc1)
sc1.2 <- aggregate(pH~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc1.2

sc2 <- lmer(Moisture~Use+(1|Pair), data=soil)
shapiro.test(resid(sc2))
Anova(sc2)
sc2.2 <- aggregate(Moisture~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc2.2

sc3 <- lmer(log10(NO3)~Use + (1|Pair), data=soil)
shapiro.test(resid(sc3))
Anova(sc3)
hist(soil$NO3)
glm1 = glmer(NO3 ~ Use +(1|Pair),data=soil, family=Gamma(link=log), nAGQ=0, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1000000)))
Anova(glm1)
sc3.2 <- aggregate(NO3~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc3.2

sc4 <- lmer(NH4~Use+(1|Pair),data=soil)
shapiro.test(resid(sc4))
Anova(sc4)
sc4.2 <- aggregate(NH4~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc4.2

sc5 <- lmer(log10(TDN)~Use+(1|Pair),data=soil)
shapiro.test(resid(sc5))
Anova(sc5)
hist(soil$TDN)
glm2 = glmer(TDN ~ Use+(1|Pair),data=soil, family=Gamma(link=log))
Anova(glm2)
sc5.2 <- aggregate(TDN~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc5.2

sc6 <- lmer(log10(DON) ~ Use+(1|Pair), data=soil)
shapiro.test(resid(sc6))
Anova(sc6)
hist(soil$DON)
glm6 = glmer(DON ~ Use+(1|Pair),data=soil, family=Gamma(link=log));
Anova(glm6)
sc6.2 <- aggregate(DON~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc6.2

sc7 <- lmer(log10(TN) ~ Use + (1|Pair), data=soil)
shapiro.test(resid(sc7))
Anova(sc7)
sc9.2 <- aggregate(TN~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc9.2

sc8 <- lmer(log10(TC) ~ Use+(1|Pair), data=soil)
shapiro.test(resid(sc8))
Anova(sc8)
sc8.2 <- aggregate(TC~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc8.2

sc9 <- lmer(C.N~ Use+(1|Pair), data=soil)
shapiro.test(resid(sc9))
Anova(sc9)
sc9.2 <- aggregate(C.N~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc9.2

sc10 <- lmer(DOC~ Use + (1|Pair), data=soil)
shapiro.test(resid(sc10))
Anova(sc10)
sc10.2 <- aggregate(DOC~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc10.2

sc11 <- lmer(DOC.TDN ~ Use+(1|Pair), data=soil)
shapiro.test(resid(sc11))
Anova(sc11)
sc11.2 <- aggregate(DOC.TDN~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc11.2

sc12 <- lmer(SIR ~ Use + (1|Pair), data=soil)
shapiro.test(resid(sc12))
Anova(sc12)
sc12.2 <- aggregate(SIR~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc12.2

sc14 <- lmer(MBC ~ Use+(1|Pair), data=soil)
shapiro.test(resid(sc14))
Anova(sc14)
sc14.2 <- aggregate(MBC~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc14.2

sc15 <- lmer(log10(MBN) ~ Use + (1|Pair), data=soil)
shapiro.test(resid(sc15))
Anova(sc15)
sc15.2 <- aggregate(MBN~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc15.2

sc16 <- lmer(log10(MBC.N) ~ Use+(1|Pair), data=soil)
shapiro.test(resid(sc16))
Anova(sc16)
sc16.2 <- aggregate(MBC.N~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc16.2

#Enzymes

e1 <- lmer(BG ~ Use + (1|Pair), data=soil)
shapiro.test(resid(e1))
e1.1<- glmer(BG ~ Use+(1|Pair), data=soil, family=Gamma(link=log))
Anova(e1.1)
e1.2 <- aggregate(BG~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
e1.2

e2 <- lmer(XYL ~ Use + (1|Pair), data=soil)
shapiro.test(resid(e2))
Anova(e2)
e2.2 <- aggregate(XYL~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
e2.2

e3 <- lmer(log10(CHB) ~ Use + (1|Pair), data=soil)
shapiro.test(resid(e3))
hist(soil$CHB)
e3.1<- glmer(CHB ~ Use+(1|Pair), data=soil, family=Gamma(link=log))
Anova(e3.1)
e3.2 <- aggregate(CHB~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
e3.2

e4 <- lmer(log10(NAG) ~ Use + (1|Pair), data=soil)
shapiro.test(resid(e4))
e4.1<- glmer(NAG ~ Use+(1|Pair), data=soil, family=Gamma(link=log))
Anova(e4.1)
e4.2 <- aggregate(NAG~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
e4.2

e5 <- lmer(LAP ~ Use + (1|Pair), data=soil)
shapiro.test(resid(e5))
Anova(e5)
e5.2 <- aggregate(LAP~Use, data=soil, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
e5.2

#Process Rates

#Net nitrification

r1 <- lmer(log10(Nit+1) ~ Use + (1|Pair), data=soil)
shapiro.test(resid(r1))
r1.3 <- glmer(Nit+1 ~ Use+(1|Pair), data=soil, family=Gamma(link=log))
Anova(r1.3)

#Net N mineralization

r2 <- lmer(Nmin ~ Use + (1|Pair), data=soil)
shapiro.test(resid(r2))
r2.1<- glmer(Nmin+2 ~ Use+(1|Pair), data=soil, family=Gamma(link=log))
Anova(r2.1)

#C mineralization

r3 <- lmer(log10(CO2) ~ Use + (1|Pair), data=soil)
shapiro.test(resid(r3))
Anova(r3)

#Metabolic quotient

r4 <- lmer(qCO2 ~ Use+(1|Pair), data=soil)
shapiro.test(resid(r4))
Anova(r4)

#Gross N mineralization

r5 <- lmer(log10(gNmin)~ Use+(1|Pair),data=soil)
shapiro.test(resid(r5))
r5.1 <- glmer(gNmin ~ Use + (1|Pair), family=Gamma(link=log), data=soil)
Anova(r5.1)

#Gross nitrification

r6 <- lmer(gNit~ Use+(1|Pair),data=soil)
shapiro.test(resid(r6))
r6.1 <- glmer(gNit ~ Use + (1|Pair), family=Gamma(link=log), data=soil)
Anova(r6.1)

#Biomass specific N mineralization

r7 <- lmer(qNH4~ Use+(1|Pair),data=soil)
shapiro.test(resid(r7))
r7.1 <- glmer(qNH4 ~ Use + (1|Pair), family=Gamma(link=log), data=soil)
Anova(r7.1)

#NH4 immobilization

r8 <- lmer(NH4imm~ Use+(1|Pair),data=soil)
shapiro.test(resid(r8))
Anova(r8)

#Nitrate immobilization

r9 <- lmer(log10(NO3imm)~ Use+(1|Pair),data=soil)
shapiro.test(resid(r9))
r9.1 <- glmer(NO3imm+1 ~ Use + (1|Pair), family=Gamma(link=log), data=soil)
Anova(r9)
Anova(r9.1)

#Plots for process rates

soil$CO2 <- soil$CO2*24

cmin <- aggregate(CO2~ Use, soil,FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x)))
cmin

cmin <- do.call(data.frame, cmin)
cmin

cmin$se <- cmin$CO2.sd / sqrt(cmin$CO2.n)
cmin

bar_plotcmin <- ggplot(cmin, aes(x=as.factor(cmin$Use), y=cmin$CO2.mean, fill = cmin$Use)) + 
  geom_bar(stat="identity", color="black", size=.5, position="dodge") +
  geom_errorbar(aes(ymin=cmin$CO2.mean-cmin$se, ymax=cmin$CO2.mean+cmin$se), width=0.2, size=.5, position=position_dodge(0.9)) +
  annotate('text', 1.75, 20.5, label="Disturbance: italic(P) == 0.14", size=5, parse=TRUE)+
  scale_fill_manual(values=c("white", "gray50")) +
  labs(list(x ="")) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=22)) +
  theme(legend.position="none") +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(bquote(~'C mineralization'~~ (mu*'g'~'CO'[2]~'-C'~~gdw^-1~~d^-1))) 
plot(bar_plotcmin)

gNmin <- aggregate(gNmin~ Use, soil,FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x)))
gNmin

gNmin <- do.call(data.frame, gNmin)
gNmin

gNmin$se <- gNmin$gNmin.sd / sqrt(gNmin$gNmin.n)
gNmin

bar_plotgNmin <- ggplot(cmin, aes(x=as.factor(gNmin$Use), y=gNmin$gNmin.mean, fill = gNmin$Use)) + 
  geom_bar(stat="identity", color="black", size=.5, position="dodge") +
  geom_errorbar(aes(ymin=gNmin$gNmin.mean-gNmin$se, ymax=gNmin$gNmin.mean+gNmin$se), width=0.2, size=.5, position=position_dodge(0.9)) +
  annotate('text', 1.75, 3, label="Disturbance: italic(P) == 0.02", size=5, parse=TRUE)+
  scale_fill_manual(values=c("white", "gray50")) +
  labs(list(x ="")) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=22)) +
  theme(legend.position="none") +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(bquote(~'N mineralization'~~ (mu*'g'~'NH'[4]~'-N'~~gdw^-1~~d^-1))) 
plot(bar_plotgNmin)

gNit <- aggregate(gNit~ Use, soil,FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x)))
gNit

gNit <- do.call(data.frame, gNit)
gNit

gNit$se <- gNit$gNit.sd / sqrt(gNit$gNit.n)
gNit

bar_plotgNit <- ggplot(gNit, aes(x=as.factor(gNit$Use), y=gNit$gNit.mean, fill = gNit$Use)) + 
  geom_bar(stat="identity", color="black", size=.5, position="dodge") +
  geom_errorbar(aes(ymin=gNit$gNit.mean-gNit$se, ymax=gNit$gNit.mean+gNit$se), width=0.2, size=.5, position=position_dodge(0.9)) +
  annotate('text', 1.77, .4, label="Disturbance: italic(P) < 0.001", size=5, parse=TRUE)+
  scale_fill_manual(values=c("white", "gray50")) +
  labs(list(x ="")) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=22)) +
  theme(legend.position="none") +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(bquote(~'Nitrification'~~ (mu*'g'~'NO'[3]~'-N'~~gdw^-1~~d^-1))) 
plot(bar_plotgNit)

qnh4 <- aggregate(qNH4~ Use, soil,FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x)))
qnh4

qnh4 <- do.call(data.frame, qnh4)
qnh4

qnh4$se <- qnh4$qNH4.sd / sqrt(qnh4$qNH4.n)
qnh4

bar_plotqNH4 <- ggplot(qnh4, aes(x=as.factor(qnh4$Use), y=qnh4$qNH4.mean, fill = qnh4$Use)) + 
  geom_bar(stat="identity", color="black", size=.5, position="dodge") +
  geom_errorbar(aes(ymin=qnh4$qNH4.mean-qnh4$se, ymax=qnh4$qNH4.mean+qnh4$se), width=0.2, size=.5, position=position_dodge(0.9)) +
  annotate('text', 1.75, 17, label="Disturbance: italic(P) < 0.001", size=5, parse=TRUE)+
  scale_fill_manual(values=c("white", "gray50")) +
  labs(list(x ="")) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=22)) +
  theme(legend.position="none") +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(bquote(~'qNH'[4]~~ (mu*'g'~'NH'[4]~'-N'~~'mg Microbial C'^-1~~d^-1))) 
plot(bar_plotqNH4)

soil$qCO2 <- soil$qCO2*24

qco2 <- aggregate(qCO2~ Use, soil,FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x)))
qco2

qco2 <- do.call(data.frame, qco2)
qco2

qco2$se <- qco2$qCO2.sd / sqrt(qco2$qCO2.n)
qco2

bar_plotqCO2 <- ggplot(qco2, aes(x=as.factor(qco2$Use), y=qco2$qCO2.mean, fill = qco2$Use)) + 
  geom_bar(stat="identity", color="black", size=.5, position="dodge") +
  geom_errorbar(aes(ymin=qco2$qCO2.mean-qco2$se, ymax=qco2$qCO2.mean+qco2$se), width=0.2, size=.5, position=position_dodge(0.9)) +
  annotate('text', 1.75, 120, label="Disturbance: italic(P) == 0.003", size=5, parse=TRUE)+
  scale_fill_manual(values=c("white", "gray50")) +
  labs(list(x ="")) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=22)) +
  theme(legend.position="none") +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(bquote(~'qCO'[2]~~ (mu*'g'~'CO'[2]~'-C'~~'mg Microbial C'^-1~~d^-1))) 
plot(bar_plotqCO2)

#Prepare Microbial variables, add to soil dataframe

bac_alpha <- estimate_richness(FD_16s, split = TRUE, measures=c("shannon","Observed","Simpson"))

bac_alpha <- cbind(bac_alpha, sam_data16s)

bac_alpha <- bac_alpha[order(bac_alpha$Number),]
head(bac_alpha)

bac_alpha <- bac_alpha[-39,]

soil$bac_alpha <- bac_alpha$Shannon

soil2 <- soil[-37,]

fung_alpha <- estimate_richness(FD_ITS, split = TRUE, measures=c("shannon","Observed","Simpson"))

fung_alpha <- cbind(fung_alpha, sam_dataITS)

fung_alpha <- fung_alpha[,c(1,2,3,12)]

fung_alpha <- rbind(fung_alpha, "FD37" = c(NA, NA, NA, 37))

fung_alpha <- fung_alpha[order(fung_alpha$Number),]
head(fung_alpha)

fung_alpha <- fung_alpha[-39,]

soil$fung_alpha <- fung_alpha$Shannon

bac_beta.1<- cmdscale(b1, k=2, eig=TRUE)

bac_beta <- as.data.frame(bac_beta.1$points)
bac_beta$number <- sam_data16s$Number
head(bac_beta)

bac_beta <- bac_beta[order(bac_beta$number),]
head(bac_beta)

bac_beta <- bac_beta[-39,]

soil$bac_beta <- bac_beta$V1

fung_beta.1 <- cmdscale(b2, k=2, eig=TRUE)

fung_beta <- as.data.frame(fung_beta.1$points)
fung_beta$number <- sam_dataITS$Number
head(fung_beta)

fung_beta <- rbind(fung_beta, "FD37" = c(NA, NA, 37))

fung_beta <- fung_beta[order(fung_beta$number),]
head(fung_beta)

fung_beta <- fung_beta[-39,]

soil$fung_beta <- fung_beta$V1

var1 <- round(bac_beta.1$eig[1]/sum(bac_beta.1$eig)*100,1)
var2 <- round(bac_beta.1$eig[2]/sum(bac_beta.1$eig)*100,1)

var1.1 <- round(fung_beta.1$eig[1]/sum(fung_beta.1$eig)*100,1)
var2.1 <- round(fung_beta.1$eig[2]/sum(fung_beta.1$eig)*100,1)

#PCoA plots for supplement

Bacord <- plot_ordination(
  physeq = FD_16s,
  ordination = bac_beta.1,
  color = "Treatment",
  shape = "Treatment",
  title = "") + 
  scale_color_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_classic() +
  annotate("text", -.15, .15, label="", parse=TRUE, size=6, fontface="bold")+
  annotate("text", -.23, -.16, label="", size=5)+
  xlab("PCoA 1 (31.2%)") + 
  ylab("PCoA 2 (9.0%)") +
  theme(axis.text=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  geom_point(aes(color = Treatment), size = 3) +
  geom_point(size = 1.5) +
  coord_fixed()
Bacord

ITSord <- plot_ordination(
  physeq = FD_ITS,
  ordination = fung_beta.1,
  color = "Treatment",
  shape = "Treatment",
  title = "") + 
  scale_color_manual(values = c("#FC8D62", "#8DA0CB")) +
  theme_classic() +
  annotate("text", -.15, .15, label="", parse=TRUE, size=6, fontface="bold")+
  annotate("text", -.23, -.16, label="", size=5)+
  xlab("PCoA 1 (9.5%)") + 
  ylab("PCoA 2 (7.0%)") +
  theme(axis.text=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  geom_point(aes(color = Treatment), size = 3) +
  geom_point(size = 1.5) +
  coord_fixed()
ITSord

#r:K stuff

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

library(reshape)

FD_16s_phyla1 <- cast(FD_16s_phyla, Sample+Treatment+Watershed.Pair ~ Phylum, value="Abundance")
FD_16s_phyla1

FD_16s_phyla1$OligoCopio <- (FD_16s_phyla1$p__Bacteroidetes+FD_16s_phyla1$p__Proteobacteria)/(FD_16s_phyla1$p__Actinobacteria+FD_16s_phyla1$p__Acidobacteria)
head(FD_16s_phyla1)

FD_16s_phyla1$Number <- sam_data16s$Number

FD_16s_phyla1 <- FD_16s_phyla1[order(FD_16s_phyla1$Number),]

FD_16s_phyla1 <- FD_16s_phyla1[-39,]

soil$`r:K` <- FD_16s_phyla1$OligoCopio

#Regression stuff

library(MuMIn)
library(metafor)
library(mice)

eval(metafor:::.MuMIn)

soil.0 <- soil

names(soil.0)[names(soil.0) == "r:K"] <- "rK"

soil.reg <- mice(data = soil.0, seed=100)

soil.2 <- complete(soil.reg,1)

#best subsets regression model for C mineralization - soil chemsitry

library(MASS)
library(leaps)

CO2reg1 <- data.frame(soil.2$Moisture, soil.2$pH, soil.2$NH4, soil.2$NO3, soil.2$DOC, soil.2$DOC.TDN, soil.2$C.N,soil.2$CO2) 
library(bestglm)
CO2.bestglm1 <-
  bestglm(Xy = CO2reg1,
          family = gaussian,
          IC = "AIC",                 # Information criteria 
          method = "exhaustive")
summary(CO2.bestglm1$BestModel)

#best subsets regression model for C mineralization -  community composition variables

CO2reg2 <- data.frame(soil.2$MBC,soil.2$CHB, soil.2$MBC.N,soil.2$BG,soil.2$XYL,soil.2$fung_alpha, soil.2$fung_beta, soil.2$bac_alpha,soil.2$bac_beta,soil.2$FB,soil.2$rK, soil.2$CO2) 
CO2.bestglm2 <-
  bestglm(Xy = CO2reg2,
          family = gaussian,
          IC = "AIC",                 # Information criteria 
          method = "exhaustive")
summary(CO2.bestglm2$BestModel)

#Full best subsets regression model for C mineralization 

library(bestglm)

CO2reg3 <- data.frame(soil.2$BG,soil.2$XYL,soil.2$CHB,soil.2$fung_alpha,soil.2$fung_beta, soil.2$bac_alpha,soil.2$bac_beta,soil.2$FB,soil.2$rK,soil.2$Moisture,soil.2$pH,soil.2$NH4,soil.2$NO3,soil.2$DOC,soil.2$DOC.TDN,soil.2$C.N,soil.2$MBC,soil.2$MBC.N, soil.2$CO2) 
CO2.bestglm3 <-
  bestglm(Xy = CO2reg3,
          family = gaussian,
          IC = "AIC",                 # Information criteria 
          method = "exhaustive")
summary(CO2.bestglm3$BestModel)

#Model selection

library(AICcmodavg)

aictab(cand.set=list(CO2.bestglm1$BestModel, CO2.bestglm2$BestModel, CO2.bestglm3$BestModel),modnames=c("soil","microbes","full"))

#Dredge models for C-mineralization

library(arm)

CO2fullreg <- lm(CO2 ~ Moisture+pH+NH4+NO3+DOC+DOC.TDN+C.N+BG+XYL+CHB+MBC+MBC.N+FB+bac_alpha+fung_alpha+bac_beta+fung_beta+rK, data=soil.2,na.action=na.fail)
stdz.model <- standardize(CO2fullreg,standardize.y = TRUE)
res<- dredge(stdz.model, trace=2)
importance(res)
avg <- model.avg(subset(res, delta <= 4))
avg
#Turn off AICcmodavg package for the importance function to work with MuMIn objects
importance(avg)

#best subsets regression model for N mineralization - soil chemsitry

NH4reg1 <- data.frame(soil.2$Moisture, soil.2$pH, soil.2$NH4, soil.2$NO3, soil.2$DOC, soil.2$DOC.TDN, soil.2$C.N,soil.2$gNmin) 
library(bestglm)
NH4.bestglm1 <-
  bestglm(Xy = NH4reg1,
          family = gaussian,
          IC = "AIC",                 # Information criteria 
          method = "exhaustive")
summary(NH4.bestglm1$BestModel)

#best subsets regression model for N mineralization -  community composition variables

NH4reg2 <- data.frame(soil.2$MBC,soil.2$LAP,soil.2$NAG, soil.2$MBC.N,soil.2$fung_alpha, soil.2$fung_beta, soil.2$bac_alpha,soil.2$bac_beta,soil.2$FB,soil.2$rK, soil.2$gNmin) 
NH4.bestglm2 <-
  bestglm(Xy = NH4reg2,
          family = gaussian,
          IC = "AIC",                 # Information criteria 
          method = "exhaustive")
summary(NH4.bestglm2$BestModel)

#Full best subsets regression model for N mineralization 

NH4reg3 <- data.frame(soil.2$fung_alpha,soil.2$LAP,soil.2$NAG, soil.2$fung_beta, soil.2$bac_alpha,soil.2$bac_beta,soil.2$FB,soil.2$rK,soil.2$Moisture,soil.2$pH,soil.2$NH4,soil.2$NO3,soil.2$DOC,soil.2$DOC.TDN,soil.2$C.N,soil.2$MBC,soil.2$MBC.N, soil.2$gNmin) 
NH4.bestglm3 <-
  bestglm(Xy = NH4reg3,
          family = gaussian,
          IC = "AIC",                 # Information criteria 
          method = "exhaustive")
summary(NH4.bestglm3$BestModel)

#Model selection for N mineralization

aictab(cand.set=list(NH4.bestglm1$BestModel, NH4.bestglm2$BestModel),modnames=c("soil","microbes"))

#Full dredge model for N-mineralization

NH4fullreg <- lm(gNmin ~ Moisture+pH+NH4+NO3+DOC+DOC.TDN+C.N+MBC+MBC.N+LAP+NAG+FB+bac_alpha+fung_alpha+bac_beta+fung_beta+rK, data=soil.2,na.action=na.fail)
stdz.model2 <- standardize(NH4fullreg,standardize.y = TRUE)
res2<- dredge(stdz.model2, trace=2)
avg2 <- model.avg(subset(res2, delta <= 4))
importance(res2)
avg2

#best subsets regression model for nitrification - soil chemsitry

NO3reg1 <- data.frame(soil.2$Moisture, soil.2$pH, soil.2$NH4, soil.2$NO3, soil.2$DOC, soil.2$DOC.TDN, soil.2$C.N,soil.2$gNit) 
library(bestglm)
NO3.bestglm1 <-
  bestglm(Xy = NO3reg1,
          family = gaussian,
          IC = "AIC",                 # Information criteria 
          method = "exhaustive")
summary(NO3.bestglm1$BestModel)

#best subsets regression model for Nitrification -  community composition variables

NO3reg2 <- data.frame(soil.2$MBC,soil.2$MBC.N,soil.2$logAOA,soil.2$logAOB,soil.2$logCAOB,soil.2$fung_alpha, soil.2$fung_beta, soil.2$bac_alpha,soil.2$bac_beta,soil.2$FB,soil.2$rK, soil.2$gNit) 
NO3.bestglm2 <-
  bestglm(Xy = NO3reg2,
          family = gaussian,
          IC = "AIC",                 # Information criteria 
          method = "exhaustive")
summary(NO3.bestglm2$BestModel)

#Full best subsets regression model for Nitrification

NO3reg3 <- data.frame(soil.2$fung_alpha,soil.2$fung_beta,soil.2$MBN, soil.2$bac_alpha,soil.2$bac_beta,soil.2$FB,soil.2$rK,soil.2$Moisture,soil.2$pH,soil.2$NH4,soil.2$NO3,soil.2$DOC,soil.2$DOC.TDN,soil.2$C.N,soil.2$logAOA,soil.2$logAOB,soil.2$logCAOB,soil.2$MBC,soil.2$MBC.N, soil.2$gNit) 
NO3.bestglm3 <-
  bestglm(Xy = NO3reg3,
          family = gaussian,
          IC = "AIC",                 # Information criteria 
          method = "exhaustive")
summary(NO3.bestglm3$BestModel)

#Model selection for Nitrification

aictab(cand.set=list(NO3.bestglm1$BestModel, NO3.bestglm2$BestModel,NO3.bestglm3$BestModel),modnames=c("soil","microbes","full"))

#Full dredge model for Nitrification

NO3fullreg <- lm(gNit ~ Moisture+pH+NH4+NO3+DOC+DOC.TDN+C.N+logAOA+logAOB+logCAOB+MBC+MBC.N+FB+bac_alpha+fung_alpha+bac_beta+fung_beta+rK, data=soil.2,na.action=na.fail)
stdz.model3 <- standardize(NO3fullreg,standardize.y = TRUE)
res3<- dredge(stdz.model3, trace=2)
avg3 <- model.avg(subset(res3, delta <= 4))
importance(res3)
avg3

#Multifunctionality approach

soil_multi <- soil.2[,c(1,2,3,4,5,9,11,26,40,39,31,32,33,34,35)]

soil_multi2 <- soil_multi[,c(6:15)]

soil_multi3 <- scale(soil_multi2, center=T, scale=T)

soil_multi4 <- cbind(soil_multi[,c(1:5)],soil_multi3)

set.seed(101)

adonis2(soil_multi4[,c(6:15)] ~ Use, data = soil_multi4, permutations = 9999, method = "euclidean", na.rm=T, strata=Pair)

multi_beta <- vegdist(soil_multi4[6:15], method="euclidean",na.rm=T)

multi <- cmdscale(multi_beta, k=2, eig=TRUE)

var1 <- round(multi$eig[1]/sum(multi$eig)*100,1)
var2 <- round(multi$eig[2]/sum(multi$eig)*100,1)

P <- princomp(soil_multi4[6:15],cor=FALSE, scores=TRUE)

summary(P)

loadings(P)

data.scores1 <- as.data.frame(multi$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores1$Treatment <- soil_multi4$Use

multi_ord <- aggregate(V1 ~ Treatment, data.scores1,FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x)))
multi_ord

multi_ord2 <- aggregate(V2 ~ Treatment, data.scores1,FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x)))
multi_ord2

multi_ord <- do.call(data.frame, multi_ord)
multi_ord

multi_ord2 <- do.call(data.frame, multi_ord2)
multi_ord2

multi_ord3 <- cbind(multi_ord, multi_ord2[,c(2:4)])
multi_ord3

multi_ord3$V1se <- multi_ord3$V1.sd / sqrt(multi_ord3$V1.n)
multi_ord3

multi_ord3$V2se <- multi_ord3$V2.sd / sqrt(multi_ord3$V2.n)
multi_ord3

ord_multi <- ggplot(multi_ord3, aes(x=V1.mean, y=V2.mean, fill=Treatment)) + 
  geom_errorbar(aes(ymin=multi_ord3$V2.mean-multi_ord3$V2se, ymax=multi_ord3$V2.mean+multi_ord3$V2se), width=0.1, size=.5, position=position_dodge(0.9)) +
  geom_errorbarh(aes(xmin=multi_ord3$V1.mean-multi_ord3$V1se, xmax=multi_ord3$V1.mean+multi_ord3$V1se), height=0.1, size=.5) +
  geom_point(data=multi_ord3,aes(x=V1.mean,y=V2.mean, fill=Treatment),pch=22,size=8) + # add the point markers
  xlab("PCA 1 (49%)" ) +
  ylab("PCA 2 (20%)") +
  annotate("text", -.15, -.15, label="Disturbance: ~italic(P) == 0.059", parse=TRUE, size=5, fontface="bold")+
  annotate("text", -.27, .7, label="Multifunctionality", parse=TRUE, size=5, fontface="bold")+
  scale_fill_manual(values=c("white","gray50")) +
  theme_classic() +
  theme(axis.text=element_text(size=14)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=14)) +
  theme(legend.position = c(0.78, .13)) +
  theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"))+
  theme(legend.title = element_blank(),legend.text=element_text(size=16)) 
ord_multi

#Make panel Figure 1

jpeg(filename="Figure1.jpeg", bg="transparent", res=500, units = "in", height=10, width=13) 

fig1 <- plot_grid(bar_plotcmin, bar_plotgNmin, bar_plotgNit,bar_plotqCO2, bar_plotqNH4,ord_multi , ncol = 3, labels=c('A', 'B', 'C','D', 'E', 'F'), align="hv", label_size=20,label_y=1.02)
fig1

dev.off()

library(corrplot)
library(RColorBrewer)
scalebluered <- colorRampPalette(brewer.pal(8, "RdBu"))(8)

#Correlation figure 

soil_multi6 <- soil_multi4

soil_multi6$Multifunctionality <- rowMeans(soil_multi3, na.rm=T)

soil.3 <- soil.2

soil.3$Multifunctionality <- soil_multi6$Multifunctionality

names(soil.3)[names(soil.3) == "CO2"] <- "C min"

names(soil.3)[names(soil.3)== "gNmin"] <- "N min"

names(soil.3)[names(soil.3) == "gNit"] <- "Nitrification"

names(soil.3)[names(soil.3) == "DOC.TDN"] <- "DOC:TDN"

names(soil.3)[names(soil.3) == "MBC.N"] <- "MBC:MBN"

names(soil.3)[names(soil.3) == "logAOB"] <- "AOB"

names(soil.3)[names(soil.3)== "logAOA"] <- "AOA"

names(soil.3)[names(soil.3)== "logCAOB"] <- "CAOB"

names(soil.3)[names(soil.3) == "FB"] <- "ITS:16S"

names(soil.3)[names(soil.3)== "bac_alpha"] <- "16S (H')"

names(soil.3)[names(soil.3) == "fung_alpha"] <- "ITS (H')"

names(soil.3)[names(soil.3) == "bac_beta"] <- "16S (PCoA 1)"

names(soil.3)[names(soil.3)== "fung_beta"] <- "ITS (PCoA 1)"

names(soil.3)[names(soil.3)== "rK"] <- "16S (r:K)"

names(soil.3)[names(soil.3) == "Nmin"] <- "Net N mineralization"

names(soil.3)[names(soil.3)== "Nit"] <- "Net Nitrification"

names(soil.3)[names(soil.3) == "C.N"] <- "C:N"

library(corrplot)
library(RColorBrewer)
scalebluered <- colorRampPalette(brewer.pal(8, "RdBu"))(8)

Corr <- cor(soil.3[,c(6,7,8,10,13,15,19,9,11,26,39,40,49,32,33,31,34,35,36,37,38,21,24,30,44,45,46,47,48)], use="complete.obs", method="pearson")

res1 <- cor.mtest(soil.3[,c(6,7,8,10,13,15,19,9,11,26,39,40,49,32,33,31,34,35,36,37,38,21,24,30,44,45,46,47,48)], conf.level = .95)

jpeg(filename="correlogram.jpeg", bg="transparent", res=500, units = "in", height=7, width=7) 

corrplot(Corr, method="color", col=scalebluered,number.cex=.35, addCoef.col = "black",  type="upper",cl.cex=1, tl.cex=.75, tl.col="black", tl.srt=45,sig.level=0.05,p.mat = res1$p, insig = "blank")

dev.off()

#C mineralization correlations

Corr2 <- cor(soil.3[,c(26,21,24,32,33,35, 30,44, 45,46,47,48)], use="complete.obs", method="pearson")

res2 <- cor.mtest(soil.3[,c(26,21,24,32,33,35, 30,44, 45,46,47,48)], conf.level = .95)

jpeg(filename="CO2correlogram.jpeg", bg="transparent", res=500, units = "in", height=1.5, width=6) 

corrplot(Corr2[1,1:12, drop=FALSE], method="color", col=scalebluered, diag=FALSE, addCoef.col = "black",  type="upper",  tl.col="black",tl.srt=45, sig.level=0.05,p.mat = res2$p, insig = "blank", cl.pos='n')

dev.off()

#N mineralization correlations

Corr3 <- cor(soil.3[,c(40,21,24, 31,34, 30,44, 45,46,47,48)], use="complete.obs", method="pearson")

res3 <- cor.mtest(soil.3[,c(40,21,24, 31,34, 30,44, 45,46,47,48)], conf.level = .95)

jpeg(filename="NH4correlogram.jpeg", bg="transparent", res=500, units = "in", height=1.5, width=5.5) 

corrplot(Corr3[1,1:11, drop=FALSE], method="color", col=scalebluered, diag=FALSE, addCoef.col = "black",  type="upper",  tl.col="black",tl.srt=45, sig.level=0.05,p.mat = res3$p, insig = "blank", cl.pos='n')

dev.off()

#nitrification correlations

Corr4 <- cor(soil.3[,c(39,21,24, 36,37,38, 30, 44,45,46,47,48)],use = "complete.obs", method="pearson")

res4 <- cor.mtest(soil.3[,c(39,21,24, 36,37,38, 30, 44,45,46,47,48)])

jpeg(filename="NO3correlogram.jpeg", bg="transparent", res=500, units = "in", height=1.5, width=6.7)

corrplot(Corr4[1,1:12, drop=FALSE], method="color", col=scalebluered, diag=FALSE, addCoef.col = "black",  type="upper",  tl.col="black", tl.srt=45,sig.level=0.05,p.mat = res4$p, insig = "blank", cl.pos='n')

dev.off()

#Multifunctionality correlations

Corr5 <- cor(soil.3[,c(49,21,24, 44,45,46,47,48)],use = "complete.obs", method="pearson")

res5<- cor.mtest(soil.3[,c(49,21,24, 44,45,46,47,48)])

jpeg(filename="Mcorrelogram.jpeg", bg="transparent", res=500, units = "in", height=1.5, width=6)

corrplot(Corr5[1,1:8, drop=FALSE], method="color", col=scalebluered, diag=FALSE, addCoef.col = "black",  type="upper",  tl.col="black", tl.srt=45,sig.level=0.05,p.mat = res5$p, insig = "blank", cl.pos='n')

dev.off()

#Additional correlations

soil.3$bac_rich <- bac_alpha$Observed

cor.test(soil.3$bac_rich, soil.3$`C min`)
cor.test(soil.3$bac_rich, soil.3$Nitrification)
cor.test(soil.3$bac_rich, soil.3$`N min`)

soil.3$bac_even <- bac_alpha$Simpson

cor.test(soil.3$bac_even, soil.3$`C min`)
cor.test(soil.3$bac_even, soil.3$Nitrification)
cor.test(soil.3$bac_even, soil.3$`N min`)

soil.3$fung_rich <- fung_alpha$Observed

cor.test(soil.3$fung_rich, soil.3$`C min`)
cor.test(soil.3$fung_rich, soil.3$Nitrification)
cor.test(soil.3$fung_rich, soil.3$`N min`)

soil.3$fung_even <- fung_alpha$Simpson

cor.test(soil.3$fung_even, soil.3$`C min`)
cor.test(soil.3$fung_even, soil.3$Nitrification)
cor.test(soil.3$fung_even, soil.3$`N min`)

#Partial Correlations

library(ppcor)

#C mineralization partial correlations

soil.4 <- soil.3[,c(6,7,8,10,13,15,19,9,11,26,39,40,49,32,33,31,34,35,36,37,38,21,24,30,44,45,46,47,48)]

pcor.test(soil.4[,10],soil.4[,14],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,10],soil.4[,15],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,10],soil.4[,18],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,10],soil.4[,22],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,10],soil.4[,23],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,10],soil.4[,24],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,10],soil.4[,25],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,10],soil.4[,26],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,10],soil.4[,27],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,10],soil.4[,28],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,10],soil.4[,29],soil.4[,c(1:7)],method="pearson")

pcor2_coef <- read.csv("C:/Users/ernie/Desktop/FD_isotopes/pcorr_coef2.csv")
pcor2_pval <- read.csv("C:/Users/ernie/Desktop/FD_isotopes/pcorr_pvalue2.csv")

names(pcor2_coef)[names(pcor2_coef) == "MBC.N"] <- "MBC:MBN"

names(pcor2_coef)[names(pcor2_coef) == "FB"] <- "ITS:16S"

names(pcor2_coef)[names(pcor2_coef) == "bac_alpha"] <- "16S (H')"

names(pcor2_coef)[names(pcor2_coef) == "fung_alpha"] <- "ITS (H')"

names(pcor2_coef)[names(pcor2_coef) == "bac_beta"] <- "16S (PCoA 1)"

names(pcor2_coef)[names(pcor2_coef) == "fung_beta"] <- "ITS (PCoA 1)"

names(pcor2_coef)[names(pcor2_coef) == "r.K"] <- "16S (r:K)"

matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m
}

pcor2_c <- matrix.please(pcor2_coef)
pcor2_p <- matrix.please(pcor2_pval)

jpeg(filename="CO2correlogram2.jpeg", bg="transparent", res=500, units = "in", height=1.5, width=6) 

corrplot(pcor2_c[1,1:12,drop=F], method="color", col=scalebluered, diag=FALSE, addCoef.col = "black",  type="upper",  tl.col="black",tl.srt=45, sig.level=0.05,p.mat = pcor2_p, insig = "blank", cl.pos='n')

dev.off()

#N mineralization partial correlations

pcor.test(soil.4[,12],soil.4[,16],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,12],soil.4[,17],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,12],soil.4[,22],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,12],soil.4[,23],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,12],soil.4[,24],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,12],soil.4[,25],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,12],soil.4[,26],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,12],soil.4[,27],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,12],soil.4[,28],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,12],soil.4[,29],soil.4[,c(1:7)],method="pearson")

pcor3_coef <- read.csv("C:/Users/ernie/Desktop/FD_isotopes/pcorr_coef3.csv")
pcor3_pval <- read.csv("C:/Users/ernie/Desktop/FD_isotopes/pcorr_pvalue3.csv")

names(pcor3_coef)[names(pcor3_coef) == "MBC.N"] <- "MBC:MBN"

names(pcor3_coef)[names(pcor3_coef) == "FB"] <- "ITS:16S"

names(pcor3_coef)[names(pcor3_coef) == "bac_alpha"] <- "16S (H')"

names(pcor3_coef)[names(pcor3_coef) == "fung_alpha"] <- "ITS (H')"

names(pcor3_coef)[names(pcor3_coef) == "bac_beta"] <- "16S (PCoA 1)"

names(pcor3_coef)[names(pcor3_coef) == "fung_beta"] <- "ITS (PCoA 1)"

names(pcor3_coef)[names(pcor3_coef) == "r.K"] <- "16S (r:K)"

pcor3_c <- matrix.please(pcor3_coef)
pcor3_p <- matrix.please(pcor3_pval)

jpeg(filename="NH4correlogram2.jpeg", bg="transparent", res=500, units = "in", height=1.5, width=6) 

corrplot(pcor3_c[1,1:11,drop=F], method="color", col=scalebluered, diag=FALSE, addCoef.col = "black",  type="upper",  tl.col="black",tl.srt=45, sig.level=0.05,p.mat = pcor3_p, insig = "blank", cl.pos='n')

dev.off()

#Nitrification partial correlations

pcor.test(soil.4[,11],soil.4[,19],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,11],soil.4[,20],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,11],soil.4[,21],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,11],soil.4[,22],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,11],soil.4[,23],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,11],soil.4[,24],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,11],soil.4[,25],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,11],soil.4[,26],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,11],soil.4[,27],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,11],soil.4[,28],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,11],soil.4[,29],soil.4[,c(1:7)],method="pearson")

pcor4_coef <- read.csv("C:/Users/ernie/Desktop/FD_isotopes/pcorr_coef4.csv")
pcor4_pval <- read.csv("C:/Users/ernie/Desktop/FD_isotopes/pcorr_pvalue4.csv")

names(pcor4_coef)[names(pcor4_coef) == "MBC.N"] <- "MBC:MBN"

names(pcor4_coef)[names(pcor4_coef) == "FB"] <- "ITS:16S"

names(pcor4_coef)[names(pcor4_coef) == "bac_alpha"] <- "16S (H')"

names(pcor4_coef)[names(pcor4_coef) == "fung_alpha"] <- "ITS (H')"

names(pcor4_coef)[names(pcor4_coef) == "bac_beta"] <- "16S (PCoA 1)"

names(pcor4_coef)[names(pcor4_coef) == "fung_beta"] <- "ITS (PCoA 1)"

names(pcor4_coef)[names(pcor4_coef) == "r.K"] <- "16S (r:K)"

pcor4_c <- matrix.please(pcor4_coef)
pcor4_p <- matrix.please(pcor4_pval)

jpeg(filename="NO3correlogram2.jpeg", bg="transparent", res=500, units = "in", height=1.5, width=6.7) 

corrplot(pcor4_c[1,1:12,drop=F], method="color", col=scalebluered, diag=FALSE, addCoef.col = "black",  type="upper",  tl.col="black",tl.srt=45, sig.level=0.05,p.mat = pcor4_p, insig = "blank", cl.pos='n')

dev.off()

#Multifunctionality partial correlations

pcor.test(soil.4[,13],soil.4[,22],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,13],soil.4[,23],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,13],soil.4[,24],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,13],soil.4[,25],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,13],soil.4[,26],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,13],soil.4[,27],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,13],soil.4[,28],soil.4[,c(1:7)],method="pearson")
pcor.test(soil.4[,13],soil.4[,29],soil.4[,c(1:7)],method="pearson")

pcor5_coef <- read.csv("C:/Users/ernie/Desktop/FD_isotopes/pcorr_coef5.csv")
pcor5_pval <- read.csv("C:/Users/ernie/Desktop/FD_isotopes/pcorr_pvalue5.csv")

names(pcor5_coef)[names(pcor5_coef) == "MBC.N"] <- "MBC:MBN"

names(pcor5_coef)[names(pcor5_coef) == "FB"] <- "ITS:16S"

names(pcor5_coef)[names(pcor5_coef) == "bac_alpha"] <- "16S (H')"

names(pcor5_coef)[names(pcor5_coef) == "fung_alpha"] <- "ITS (H')"

names(pcor5_coef)[names(pcor5_coef) == "bac_beta"] <- "16S (PCoA 1)"

names(pcor5_coef)[names(pcor5_coef) == "fung_beta"] <- "ITS (PCoA 1)"

names(pcor5_coef)[names(pcor5_coef) == "r.K"] <- "16S (r:K)"

pcor5_c <- matrix.please(pcor5_coef)
pcor5_p <- matrix.please(pcor5_pval)

jpeg(filename="multicorrelogram2.jpeg", bg="transparent", res=500, units = "in", height=1.5, width=6.7) 

corrplot(pcor5_c[1,1:9,drop=F], method="color", col=scalebluered, diag=FALSE, addCoef.col = "black",  type="upper",  tl.col="black",tl.srt=45, sig.level=0.05,p.mat = pcor5_p, insig = "blank", cl.pos='n')

dev.off()

#Multifunctionality regression models

library(MASS)
library(leaps)

soil.2$Multifunctionality <- soil_multi6$Multifunctionality

multireg1 <- data.frame(soil.2$Moisture, soil.2$pH, soil.2$NH4, soil.2$NO3, soil.2$DOC, soil.2$DOC.TDN, soil.2$C.N,soil.2$Multifunctionality) 
library(bestglm)
multi.bestglm1 <-
  bestglm(Xy = multireg1,
          family = gaussian,
          IC = "AIC",                 # Information criteria 
          method = "exhaustive")
summary(multi.bestglm1$BestModel)

multireg2 <- data.frame(soil.2$MBC, soil.2$MBC.N, soil.2$fung_alpha, soil.2$fung_beta, soil.2$bac_alpha,soil.2$bac_beta,soil.2$FB,soil.2$rK, soil.2$Multifunctionality) 
multi.bestglm2 <-
  bestglm(Xy = multireg2,
          family = gaussian,
          IC = "AIC",                 # Information criteria 
          method = "exhaustive")
summary(multi.bestglm2$BestModel)

multireg3 <- data.frame(soil.2$Moisture, soil.2$pH, soil.2$NH4, soil.2$NO3, soil.2$DOC, soil.2$DOC.TDN, soil.2$C.N,soil.2$MBC, soil.2$MBC.N, soil.2$fung_alpha, soil.2$fung_beta, soil.2$bac_alpha,soil.2$bac_beta,soil.2$FB,soil.2$rK, soil.2$Multifunctionality) 
multi.bestglm3 <-
  bestglm(Xy = multireg3,
          family = gaussian,
          IC = "AIC",                 # Information criteria 
          method = "exhaustive")
summary(multi.bestglm3$BestModel)

library(AICcmodavg)
library(arm)

aictab(cand.set=list(multi.bestglm1$BestModel,multi.bestglm2$BestModel,multi.bestglm3$BestModel),modnames=c("soil","microbes","full"))

multifullreg <- lm(Multifunctionality ~ Moisture+pH+NH4+NO3+DOC+DOC.TDN+C.N+MBC+MBC.N+FB+bac_alpha+fung_alpha+bac_beta+fung_beta+rK, data=soil.2,na.action=na.fail)
stdz.model4 <- standardize(multifullreg,standardize.y = TRUE)
res4<- dredge(stdz.model4, trace=2)
avg4 <- model.avg(subset(res4, delta <= 4))
#Turn off AICcmodavg package
importance(res4)
avg4

#Variable importance (dredge results) plot. Results in the .csv file were extracted from the dredge models above

var_imp <- read.csv("C:/Users/ernie/Desktop/FD_isotopes/FD_varimp.csv")

var_impCO2 <- var_imp[c(1:6),]

var_impCO2$var <- factor(var_impCO2$var, levels = var_impCO2$var[order(1-var_impCO2$imp)])

bar_plotCO2varimp <- ggplot(var_impCO2, aes(x=as.factor(var_impCO2$var), y=var_impCO2$imp,fill=as.factor(var_impCO2$vartype))) +
  scale_fill_manual(values=c("white", "grey50")) +
  geom_bar(stat="identity", color="black", size=.5, position="dodge") +
  labs(list(x ="")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title=element_text(size=16)) +
  theme(text=element_text(size=18)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(bquote('Variable Importance')) +
  theme(legend.position="none") +
  coord_cartesian(ylim=c(0,1))+
  theme(plot.margin = unit(c(1,1,0,0), "cm"))
bar_plotCO2varimp

var_impNH4 <- var_imp[c(7:12),]

var_impNH4$var <- factor(var_impNH4$var, levels = var_impNH4$var[order(1-var_impNH4$imp)])

bar_plotNH4varimp <- ggplot(var_impNH4, aes(x=as.factor(var_impNH4$var), y=var_impNH4$imp,fill=as.factor(var_impNH4$vartype))) +
  scale_fill_manual(values=c("white", "grey50")) +
  geom_bar(stat="identity", color="black", size=.5, position="dodge") +
  labs(list(x ="")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=10)) +
  theme(axis.title=element_text(size=16)) +
  theme(text=element_text(size=18)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  theme(axis.title.x=element_blank()) +
  theme(legend.position="none") +
  scale_y_continuous(bquote('Variable Importance')) +
  coord_cartesian(ylim=c(0,1))+
  theme(plot.margin = unit(c(1,1,0,0), "cm"))
bar_plotNH4varimp

var_impNO3 <- var_imp[c(13:18),]

var_impNO3$var <- factor(var_impNO3$var, levels = var_impNO3$var[order(1-var_impNO3$imp)])

bar_plotNO3varimp <- ggplot(var_impNO3, aes(x=as.factor(var_impNO3$var), y=var_impNO3$imp,fill=as.factor(var_impNO3$vartype))) +
  scale_fill_manual(values=c("white", "grey50")) +
  geom_bar(stat="identity", color="black", size=.5, position="dodge") +
  labs(list(x ="")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title=element_text(size=16)) +
  theme(text=element_text(size=18)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  theme(axis.title.x=element_blank()) +
  theme(legend.position="none") +
  scale_y_continuous(bquote('Variable Importance')) +
  coord_cartesian(ylim=c(0,1))
bar_plotNO3varimp

var_impmulti <- var_imp[c(19:24),]

var_impmulti$var <- factor(var_impmulti$var, levels = var_impmulti$var[order(1-var_impmulti$imp)])

bar_plotmultivarimp <- ggplot(var_impmulti, aes(x=as.factor(var_impmulti$var), y=var_impmulti$imp,fill=as.factor(var_impmulti$vartype))) +
  scale_fill_manual(values=c("white", "grey50")) +
  geom_bar(stat="identity", color="black", size=.5, position="dodge") +
  labs(list(x ="")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title=element_text(size=16)) +
  theme(text=element_text(size=18)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(bquote('Variable Importance')) +
  theme(legend.position="none") +
  coord_cartesian(ylim=c(0,1))
bar_plotmultivarimp

jpeg(filename="Figure3.jpeg", bg="transparent", res=500, units = "in", height=8, width=8) 

f4 <- plot_grid(bar_plotCO2varimp, bar_plotNH4varimp, bar_plotNO3varimp,bar_plotmultivarimp, ncol = 2, align="hv", label_size=20,label_y=1.05)

f4

dev.off()

#structural equation models

library(lavaan)

sem <-      
'gNmin ~  LAP + rK + MBC
CO2 ~rK +XYL +MBC
gNit ~ logCAOB
LAP ~  rK + MBC
XYL ~ rK + MBC 
rK ~  pH  + DOC.TDN
MBC ~ DOC.TDN + pH
logCAOB ~ rK + pH + DOC.TDN
pH ~ Disturbance
DOC.TDN ~ Disturbance'

semFit<-sem(sem, data=soil.2, estimator = "ML") 

summary(semFit,rsquare=T)

T.boot <- bootstrapLavaan(semFit, R=1000, type="bollen.stine",
                          FUN=fitMeasures, fit.measures=c("chisq","rmsea","cfi","srmr"))
summary(T.boot)

semstd <- standardizedSolution(semFit)

semstd <- subset(semstd, pvalue<=0.05)

semstd$est.std <- round(semstd$est.std, 3)
semstd$se <- round(semstd$se, 3)
semstd$pvalue <- round(semstd$pvalue, 3)
semstd$z <- round(semstd$z, 3)
semstd$ci.lower <- round(semstd$ci.lower, 3)
semstd$ci.upper <- round(semstd$ci.upper, 3)

write.csv(semstd,"C:\\Users\\ernie\\Desktop\\FD_isotopes\\SEM.csv", row.names = FALSE)

#Modification Indices
modI<-modificationIndices(semFit, standardized=F)
modI[which(modI$mi>3),]

sem2 <-      
'Multifunctionality ~  rK + MBC  
rK ~ DOC.TDN + pH 
MBC ~ pH+DOC.TDN 
DOC.TDN ~ Disturbance
pH~Disturbance'

semFit2<-sem(sem2, data=soil.2, estimator = "ML")

summary(semFit2,rsquare=T)

T.boot2 <- bootstrapLavaan(semFit2, R=1000, type="bollen.stine",
                          FUN=fitMeasures, fit.measures=c("chisq","rmsea","cfi","srmr"))
summary(T.boot2)

modI2<-modificationIndices(semFit2, standardized=F)
modI2[which(modI2$mi>3),]

semstd2 <- standardizedSolution(semFit2)

semstd2 <- subset(semstd2, pvalue<=0.05)

semstd2$est.std <- round(semstd2$est.std, 3)
semstd2$se <- round(semstd2$se, 3)
semstd2$pvalue <- round(semstd2$pvalue, 3)
semstd2$z <- round(semstd2$z, 3)
semstd2$ci.lower <- round(semstd2$ci.lower, 3)
semstd2$ci.upper <- round(semstd2$ci.upper, 3)

write.csv(semstd2,"C:\\Users\\ernie\\Desktop\\FD_isotopes\\SEM2.csv", row.names = FALSE)

