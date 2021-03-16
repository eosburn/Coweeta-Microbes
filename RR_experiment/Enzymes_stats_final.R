#Read csv

enz <- read.csv("C:/Users/ernie/Desktop/Data/Rhodo Data/Enzyme Data/Enzymes Master Files/RhodoEnzymesMaster.csv")

scare <- read.csv("C:/Users/ernie/Desktop/Data/Rhodo Data/SCARE.perma.csv")

#Log transform data

info <- enz[,c(1:5)]

enz2 <- log10(enz[,c(6:19,24:32)])

enz2 <- cbind(info, enz2)

enz4 <- enz[-17,]

enz5 <- enz[,c(25:32)]

enz6 <- enz[-17,]

enz5 <- enz5[-17,]

enz7 <- log10(enz4[,25:32])

info2 <- info[-17,]

enz7 <- cbind(info2, enz7)

#Convert milligrams to micrograms

enz$DOC <- enz$DOC * 1000

enz$TN <- enz$TN * 1000

enz$MBC <- enz$MBC * 1000

enz$MBN <- enz$MBN * 1000

#Calculate DON

enz$DON <- enz$TN - (enz$NH4 + enz$NO3)

enz2$DON <- log10(enz$DON)

#Enzyme ratios

enz$LAPNAG <- enz$LAP + enz$NAG

enz$CN <- (enz$BG/enz$LAPNAG)

enz$CP <- (enz$BG/enz$AP)

enz2$CP <- log10(enz$CP)

#mixed models

library(lme4)
library(car)
library(lsmeans)
library(AICcmodavg)
library(ggplot2)

c1 <- lmer(DOC ~ Treatment + (1| Plot), data = enz, REML=FALSE)
c1.1 <- lmer(DOC ~ Treatment*Date + (1| Plot), data = enz, REML=FALSE)
aictab(cand.set=list(c1, c1.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(c1))
Anova(c1)
lsmeans(c1, list(pairwise ~ Treatment), adjust = "tukey")

c2 <- lmer(MBC ~ Treatment + (1| Plot), data = enz2, REML=FALSE)
c2.1 <- lmer(MBC ~ Treatment*Date + (1| Plot), data = enz2, REML=FALSE)
aictab(cand.set=list(c2, c2.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(c2))
Anova(c2)
lsmeans(c2, list(pairwise ~ Treatment), adjust = "tukey")

c3 <- lmer(TN ~ Treatment + (1| Plot), data = enz, REML=FALSE)
c3.1 <- lmer(TN ~ Treatment*Date + (1| Plot), data = enz, REML=FALSE)
aictab(cand.set=list(c3, c3.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(c3.1))
Anova(c3.1)
lsmeans(c3.1, list(pairwise ~ Treatment*Date), adjust = "tukey")

c4 <- lmer(MBN ~ Treatment + (1| Plot), data = enz2, REML=FALSE)
c4.1 <- lmer(MBN ~ Treatment*Date + (1| Plot), data = enz2, REML=FALSE)
aictab(cand.set=list(c4, c4.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(c4))
Anova(c4)
lsmeans(c4, list(pairwise ~ Treatment), adjust = "tukey")

c5 <- lmer(NH4 ~ Treatment + (1| Plot), data = enz, REML=FALSE)
c5.1 <- lmer(NH4 ~ Treatment*Date + (1| Plot), data = enz, REML=FALSE)
aictab(cand.set=list(c5, c5.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(c5))
Anova(c5)
lsmeans(c5, list(pairwise ~ Treatment), adjust = "tukey")

c6 <- lmer(NO3 ~ Treatment + (1| Plot), data = enz, REML=FALSE)
c6.1 <- lmer(NO3 ~ Treatment*Date + (1| Plot), data = enz, REML=FALSE)
aictab(cand.set=list(c6, c6.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(c6))
Anova(c6)
lsmeans(c6, list(pairwise ~ Treatment), adjust = "tukey")

c7 <- lmer(DON ~ Treatment + (1| Plot), data = enz, REML=FALSE)
c7.1 <- lmer(DON ~ Treatment*Date + (1| Plot), data = enz, REML=FALSE)
aictab(cand.set=list(c7, c7.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(c7.1))
Anova(c7.1)
lsmeans(c7.1, list(pairwise ~ Treatment*Date), adjust = "tukey")

c8 <- lmer(ph ~ Treatment + (1| Plot), data = enz, REML=FALSE)
c8.1 <- lmer(ph ~ Treatment*Date + (1| Plot), data = enz, REML=FALSE)
aictab(cand.set=list(c8, c8.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(c8))
Anova(c8)
lsmeans(c8, list(pairwise ~ Treatment), adjust = "tukey")

e1 <- lmer(BG ~ Treatment + (1| Plot), data = enz, REML=FALSE)
e1.1 <- lmer(BG ~ Treatment*Date + (1| Plot), data = enz, REML=FALSE)
aictab(cand.set=list(e1, e1.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e1.1))
Anova(e1)
lsmeans(e1, list(pairwise ~ Treatment), adjust = "tukey")
MeanBG <- aggregate(enz$BG, by = list(Treatment = enz$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanBG

e2 <- lmer(CHB ~ Treatment + (1| Plot), data = enz2, REML=FALSE)
e2.1 <- lmer(CHB ~ Treatment*Date + (1| Plot), data = enz2, REML=FALSE)
aictab(cand.set=list(e2, e2.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e2))
Anova(e2)
lsmeans(e2, list(pairwise ~ Treatment), adjust = "tukey")
MeanCHB <- aggregate(enz$CHB, by = list(Treatment = enz$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanCHB

e3 <- lmer(XYL ~ Treatment + (1| Plot), data = enz, REML=FALSE)
e3.1 <- lmer(XYL ~ Treatment*Date + (1| Plot), data = enz, REML=FALSE)
aictab(cand.set=list(e3, e3.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e3))
MeanXYL <- aggregate(enz$XYL, by = list(Treatment = enz$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanXYL
Anova(e3)
lsmeans(e3, list(pairwise ~ Treatment), adjust = "tukey")

e4 <- lmer(AP ~ Treatment + (1| Plot), data = enz, REML=FALSE)
e4.1 <- lmer(AP ~ Treatment*Date + (1| Plot), data = enz, REML=FALSE)
aictab(cand.set=list(e4, e4.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e4))
Anova(e4)
lsmeans(e4, list(pairwise ~ Treatment), adjust = "tukey")
MeanAP <- aggregate(enz$AP, by = list(Treatment = enz$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanAP


e5 <- lmer(LAP ~ Treatment + (1| Plot), data = enz, REML=FALSE)
e5.1 <- lmer(LAP ~ Treatment*Date + (1| Plot), data = enz, REML=FALSE)
aictab(cand.set=list(e5, e5.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e5))
Anova(e5)
lsmeans(e5, list(pairwise ~ Treatment), adjust = "tukey")
MeanLAP <- aggregate(enz$LAP, by = list(Treatment = enz$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanLAP

e6 <- lmer(NAG ~ Treatment + (1| Plot), data = enz2, REML=FALSE)
e6.1 <- lmer(NAG ~ Treatment*Date + (1| Plot), data = enz2, REML=FALSE)
aictab(cand.set=list(e6, e6.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e6))
Anova(e6)
lsmeans(e6, list(pairwise ~ Treatment), adjust = "tukey")

e7 <- lmer(POX ~ Treatment + (1| Plot), data = enz2, REML=FALSE)
e7.1 <- lmer(POX ~ Treatment*Date + (1| Plot), data = enz2, REML=FALSE)
aictab(cand.set=list(e7, e7.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e7))
Anova(e7)
lsmeans(e7, list(pairwise ~ Treatment), adjust = "tukey")

e8 <- lmer(PER ~ Treatment + (1| Plot), data = enz2, REML=FALSE)
e8.1 <- lmer(PER ~ Treatment*Date + (1| Plot), data = enz2, REML=FALSE)
aictab(cand.set=list(e8, e8.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e8))
Anova(e8)
lsmeans(e8, list(pairwise ~ Treatment), adjust = "tukey")

e9 <- lmer(BG2 ~ Treatment + (1| Plot), data = enz7, REML=FALSE)
e9.1 <- lmer(BG2 ~ Treatment*Date + (1| Plot), data = enz7, REML=FALSE)
aictab(cand.set=list(e9, e9.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e9))
Anova(e9)
lsmeans(e9, list(pairwise ~ Treatment), adjust = "tukey")

e10 <- lmer(CHB2 ~ Treatment + (1| Plot), data = enz7, REML=FALSE)
e10.1 <- lmer(CHB2 ~ Treatment*Date + (1| Plot), data = enz7, REML=FALSE)
aictab(cand.set=list(e10, e10.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e10))
Anova(e10)
lsmeans(e10, list(pairwise ~ Treatment), adjust = "tukey")

e11 <- lmer(XYL2 ~ Treatment + (1| Plot), data = enz7, REML=FALSE)
e11.1 <- lmer(XYL2 ~ Treatment*Date + (1| Plot), data = enz7, REML=FALSE)
aictab(cand.set=list(e11, e11.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e11))
Anova(e11)
lsmeans(e11, list(pairwise ~ Treatment), adjust = "tukey")

e12 <- lmer(AP2 ~ Treatment + (1| Plot), data = enz7, REML=FALSE)
e12.1 <- lmer(AP2 ~ Treatment*Date + (1| Plot), data = enz7, REML=FALSE)
aictab(cand.set=list(e12, e12.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e12))
Anova(e12)
lsmeans(e12, list(pairwise ~ Treatment), adjust = "tukey")

e13 <- lmer(LAP2 ~ Treatment + (1| Plot), data = enz7, REML=FALSE)
e13.1 <- lmer(LAP2 ~ Treatment*Date + (1| Plot), data = enz7, REML=FALSE)
aictab(cand.set=list(e13, e13.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e13))
Anova(e13)
lsmeans(e13, list(pairwise ~ Treatment), adjust = "tukey")

e14 <- lmer(NAG2 ~ Treatment + (1| Plot), data = enz7, REML=FALSE)
e14.1 <- lmer(NAG2 ~ Treatment*Date + (1| Plot), data = enz7, REML=FALSE)
aictab(cand.set=list(e14, e14.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e14))
Anova(e14)
lsmeans(e14, list(pairwise ~ Treatment), adjust = "tukey")

e15 <- lmer(POX2 ~ Treatment + (1| Plot), data = enz7, REML=FALSE)
e15.1 <- lmer(POX2 ~ Treatment*Date + (1| Plot), data = enz7, REML=FALSE)
aictab(cand.set=list(e15, e15.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e15))
Anova(e15)
lsmeans(e15, list(pairwise ~ Treatment), adjust = "tukey")

e16 <- lmer(PER2 ~ Treatment + (1| Plot), data = enz7, REML=FALSE)
e16.1 <- lmer(PER2 ~ Treatment*Date + (1| Plot), data = enz7, REML=FALSE)
aictab(cand.set=list(e16, e16.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(e16))
Anova(e16)
lsmeans(e16, list(pairwise ~ Treatment), adjust = "tukey")

m1 <- lmer(LogBac ~ Treatment + (1| Plot), data = enz, REML=FALSE)
m1.1 <- lmer(LogBac ~ Treatment*Date + (1| Plot), data = enz, REML=FALSE)
aictab(cand.set=list(m1, m1.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(m1))
Anova(m1.1)
lsmeans(m1.1, list(pairwise ~ Treatment*Date), adjust = "tukey")
MeanBac <- aggregate(enz$Bac, by = list(Treatment = enz$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanBac

MeanBac2 <- aggregate(enz$Bac, by = list(Date = enz$Date), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanBac2

m2 <- lmer(LogITS ~ Treatment + (1| Plot), data = enz, REML=FALSE)
m2.1 <- lmer(LogITS ~ Treatment*Date + (1| Plot), data = enz, REML=FALSE)
aictab(cand.set=list(m2, m2.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(m2))
Anova(m2.1)
lsmeans(m2.1, list(pairwise ~ Date), adjust = "tukey")
MeanITS2 <- aggregate(enz$ITS, by = list(Date = enz$Date), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanITS2

m3 <- lmer(F.B ~ Treatment + (1| Plot), data = enz, REML=FALSE)
m3.1 <- lmer(F.B ~ Treatment*Date + (1| Plot), data = enz, REML=FALSE)
aictab(cand.set=list(m3, m3.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(m3.1))
Anova(m3.1)
lsmeans(m3.1, list(pairwise ~ Treatment*Date), adjust = "tukey")
MeanFB2 <- aggregate(enz$F.B, by = list(Date = enz$Date), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanFB2

r1 <- lmer(CN ~ Treatment + (1| Plot), data = enz, REML=FALSE)
r1.1 <- lmer(CN ~ Treatment*Date + (1| Plot), data = enz, REML=FALSE)
aictab(cand.set=list(r1, r1.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(r1.1))
Anova(r1.1)
lsmeans(r1.1, list(pairwise ~ Treatment*Date), adjust = "tukey")

r2 <- lmer(CP ~ Treatment + (1| Plot), data = enz2, REML=FALSE)
r2.1 <- lmer(CP ~ Treatment*Date + (1| Plot), data = enz2, REML=FALSE)
aictab(cand.set=list(r2, r2.1),modnames=c("Treatment","Treatment*Date"))
shapiro.test(resid(r2.1))
Anova(r2.1)
lsmeans(r2.1, list(pairwise ~ Treatment*Date), adjust = "tukey")

#Regressions and Correlations

reg1 <- lm(enz$CN ~ enz$F.B, data=enz)
shapiro.test(resid(reg1))
summary(reg1)

reg2 <- lm(enz2$CP ~ enz2$F.B, data=enz2)
shapiro.test(resid(reg2))
summary(reg2)

enz8 <- enz2[,c(6:13)]

enz8$LogBac <- enz$LogBac 
enz8$LogITS <- enz$LogITS

library("Hmisc")

cor1 <- rcorr(as.matrix(enz8))

cor1

round(cor1$r, 3)

round(cor1$P, 3)

#Multivariate analysis of Enzyme data

library(vegan)
library(BiodiversityR)

#Remove non-enzyme activity info

enz3 <- enz[,c(6:13)]

#colmax function for relativising enzyme data

colmax <- function(enz3) sapply(enz3, max, na.rm = TRUE)

ENZMAX <- colmax(enz3)

ENZMAX

ENZMAX <- as.data.frame(ENZMAX)

#Relativize data based on maximum activity for each enzyme

enz3[,1] <- enz3[,1]/ENZMAX[1,1]
enz3[,2] <- enz3[,2]/ENZMAX[2,1]
enz3[,3] <- enz3[,3]/ENZMAX[3,1]
enz3[,4] <- enz3[,4]/ENZMAX[4,1]
enz3[,5] <- enz3[,5]/ENZMAX[5,1]
enz3[,6] <- enz3[,6]/ENZMAX[6,1]
enz3[,7] <- enz3[,7]/ENZMAX[7,1]
enz3[,8] <- enz3[,8]/ENZMAX[8,1]

enz3 <- cbind(enz3, enz[,c(1:5)])

#PCA on relativized enzyme data

P <- princomp(enz3[,-c(9:13)], cor=FALSE, scores=TRUE)

summary(P)

P$loadings

#Testing for dispersion effects

Y <- vegdist(enz3[,-c(9:13)], method="euclidean") 

treat <- enz$Treatment

disp <- betadisper(Y, treat)

anova(disp)

#PERMANOVA 

adonis2(enz3[,-c(9:13)] ~ Treatment*Date, data = enz3, permutations = 999, method = "euclidean", strata=Plot)

library(tcltk2)

a <- enz3[,c(1:8)]

enz3$Plot <- as.factor(enz3$Plot)

b <- nested.npmanova(a ~ Treatment + Plot, data = enz3, method="euc", permutations=999) 
b

#Multivariate analysis of biomass-corrected enzyme activities

enz5 <- enz[,c(25:32)]

enz6 <- enz[-17,]

enz5 <- enz5[-17,]

#colmax function for relativising enzyme data

colmax <- function(enz5) sapply(enz5, max, na.rm = TRUE)

ENZMAX <- colmax(enz5)

ENZMAX

ENZMAX <- as.data.frame(ENZMAX)

#Relativize data based on maximum activity for each enzyme

enz5[,1] <- enz5[,1]/ENZMAX[1,1]
enz5[,2] <- enz5[,2]/ENZMAX[2,1]
enz5[,3] <- enz5[,3]/ENZMAX[3,1]
enz5[,4] <- enz5[,4]/ENZMAX[4,1]
enz5[,5] <- enz5[,5]/ENZMAX[5,1]
enz5[,6] <- enz5[,6]/ENZMAX[6,1]
enz5[,7] <- enz5[,7]/ENZMAX[7,1]
enz5[,8] <- enz5[,8]/ENZMAX[8,1]

enz5 <- cbind(enz5, enz6[,c(1:5)])

#PCA on biomass-corrected enzymes

P2 <- princomp(enz5[,-c(9:16)], cor=FALSE, scores=TRUE)

summary(P2)

P2$loadings

par(mfrow=c(2,1)) 

plot(P$scores[,1],P$scores[,2], col="white", cex.main=1.5, cex.lab=1.5, cex.axis=1.5, xlab='PCA Axis 1 (51.7%)', ylab='PCA Axis 2 (23.6%)')
points(P$scores[enz3$Treatment=='CR',1], P$scores[enz3$Treatment=='CR',2], pch=21,bg='#56B4E9', cex=1.5)
points(P$scores[enz3$Treatment=='FF',1], P$scores[enz3$Treatment=='FF',2], pch=22,bg='#F0E442', cex=1.5)
points(P$scores[enz3$Treatment=='CFFR',1], P$scores[enz3$Treatment=='CFFR',2], pch=25,bg='#D55E00', cex=1.5)
points(P$scores[enz3$Treatment=='REF',1], P$scores[enz3$Treatment=='REF',2], pch=24,bg='#999999', cex=1.5)
ordiellipse(cbind(P$scores[,1], P$scores[,2]), enz3$Treatment, kind="se", conf=0.95,
            lwd=2, alpha=90, draw = "polygon", col=c('#999999', '#F0E442','#56B4E9','#D55E00'))
par(xpd = TRUE)
legend(-.8,1.1, bty = 'n', c('REF','FF', 'CR', 'CFFR'),pch=c(24,22,21,25), ncol=4, pt.bg=c('#999999', '#F0E442', '#56B4E9', '#D55E00'), cex=1.5)
legend(-1.5,1.3, bty = 'n', c('A'), text.font=2, cex=1.5)

plot(P2$scores[,1],P2$scores[,2], col="white", cex.main=1.5, cex.lab=1.5, cex.axis=1.5, xlab='PCA Axis 1 (58.7%)', ylab='PCA Axis 2 (22.3%)')
points(P2$scores[enz5$Treatment=='CR',1], P2$scores[enz5$Treatment=='CR',2], pch=21,bg='#56B4E9', cex=1.5)
points(P2$scores[enz5$Treatment=='FF',1], P2$scores[enz5$Treatment=='FF',2], pch=22,bg='#F0E442', cex=1.5)
points(P2$scores[enz5$Treatment=='CFFR',1], P2$scores[enz5$Treatment=='CFFR',2], pch=25,bg='#D55E00', cex=1.5)
points(P2$scores[enz5$Treatment=='REF',1], P2$scores[enz5$Treatment=='REF',2], pch=24,bg='#999999', cex=1.5)
ordiellipse(cbind(P2$scores[,1], P2$scores[,2]), enz5$Treatment, kind="se", conf=0.95,
            lwd=2, draw = "polygon", col=c('#D55E00', '#56B4E9', '#F0E442', '#999999'))
par(xpd=TRUE)
#legend(-1.4,1, bty = 'n', c('REF','FF', 'CR', 'CFFR'),pch=c(24,22,21,25), ncol=4, pt.bg=c('#999999', '#F0E442', '#56B4E9', '#D55E00'), cex=1.5)

#Testing for dispersion effects

Y2 <- vegdist(enz5[,-c(9:16)], method="euclidean") 

treat2 <- enz5$Treatment

disp2 <- betadisper(Y2, treat2)

anova(disp2)

#PERMANOVA 

adonis2(enz5[,-c(9:16)] ~ Treatment * Date, data = enz5, permutations = 999, method = "euclidean", strata=Plot)

a2 <- enz5[,c(1:8)]

enz5$Plot <- as.factor(enz5$Plot)

b2 <- nested.npmanova(a2 ~ Treatment + Plot, data = enz5, method="euc", permutations=999) 
b2

#TN Plot

library(ggplot2)

MeanTN <- aggregate(enz$TN, by = list(Treatment = enz$Treatment, enz$Date), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))

MeanTN <- do.call(data.frame, MeanTN)

MeanTN$se <- MeanTN$x.sd / sqrt(MeanTN$x.n)

MeanTN

colnames(MeanTN) <- c("Treatment", "Date", "mean", "sd", "n", "se")

MeanTN$Date <- as.character(MeanTN$Date)

MeanTN[MeanTN=="4/1/2017"] <- "April 2017"
MeanTN[MeanTN=="7/1/2017"] <- "July 2017"

MeanTN$Treatment <- factor(MeanTN$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanTN$Treatment

setwd("C:/Users/ernie/Desktop/Rhodo Paper")

jpeg(filename="Figure1.jpeg", bg="transparent", res=1500, units = "in", height=3, width=5) 

bar_plotTN <- ggplot(MeanTN, aes(x=as.factor(MeanTN$Date), y=MeanTN$mean, fill=as.factor(MeanTN$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=MeanTN$mean-MeanTN$se, ymax=MeanTN$mean+MeanTN$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate("text", label="a", x=.66, y=82, size=5) +
  annotate("text", label="a", x=.881, y=95, size=5) +
  annotate("text", label="a", x=1.115, y=83, size=5) +
  annotate("text", label="a", x=1.33, y=118, size=5) +
  annotate("text", label="a", x=1.66, y=80, size=5) +
  annotate("text", label="a", x=1.881, y=82, size=5) +
  annotate("text", label="ab", x=2.115, y=113, size=5) +
  annotate("text", label="b", x=2.34, y=143, size=5) +
  annotate('text', .83, 142, label="Treatment: ~italic(P) == 0.016", size=3.5, parse=TRUE)+
  annotate('text', .83, 134, label="Date: italic(P) < 0.001", size = 3.5, parse=TRUE)+
  annotate('text', .83, 126, label="Treatment %*% Date: italic(P) < 0.001", size=3.5, parse=TRUE)+
  scale_fill_manual(values=c("white", "gray70","gray28", "gray10")) +
  labs(list(x ="")) +
  theme_classic() +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=20)) +
  theme(legend.position="none") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_text(colour="white")) +
  scale_y_continuous(bquote(~'Total Extractable N'~~ (mu*'g N'~~gdw^-1))) 
plot(bar_plotTN)

dev.off()

#Aggregated Soil Chem data

#Remove observation #17 due to screwed up fumigation extraction

enz6 <- enz2[-17,]

MeanMBC <- aggregate(enz6$MBC, by = list(Treatment = enz6$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanMBC <- do.call(data.frame, MeanMBC)
MeanMBC$se <- MeanMBC$x.sd / sqrt(MeanMBC$x.n)
MeanMBC

MeanMBN <- aggregate(enz6$MBN, by = list(Treatment = enz6$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanMBN <- do.call(data.frame, MeanMBN)
MeanMBN$se <- MeanMBN$x.sd / sqrt(MeanMBN$x.n)
MeanMBN

MeanNH4 <- aggregate(enz$NH4, by = list(Treatment = enz$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanNH4 <- do.call(data.frame, MeanNH4)
MeanNH4$se <- MeanNH4$x.sd / sqrt(MeanNH4$x.n)
MeanNH4

MeanNO3 <- aggregate(enz$NO3, by = list(Treatment = enz$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanNO3 <- do.call(data.frame, MeanNO3)
MeanNO3$se <- MeanNO3$x.sd / sqrt(MeanNO3$x.n)
MeanNO3

MeanDON <- aggregate(enz$DON, by = list(Treatment = enz$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanDON <- do.call(data.frame, MeanDON)
MeanDON$se <- MeanDON$x.sd / sqrt(MeanDON$x.n)
MeanDON

MeanDOC <- aggregate(enz$DOC, by = list(Treatment = enz$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanDOC <- do.call(data.frame, MeanDOC)
MeanDOC$se <- MeanDOC$x.sd / sqrt(MeanDOC$x.n)
MeanDOC

Meanph <- aggregate(enz$ph, by = list(Treatment = enz$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
Meanph <- do.call(data.frame, Meanph)
Meanph$se <- Meanph$x.sd / sqrt(Meanph$x.n)
Meanph

#Enzyme Plots


#AP

MeanAP <- aggregate(enz$AP, by = list(Treatment = enz$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanAP

MeanAP <- do.call(data.frame, MeanAP)

MeanAP

MeanAP$se <- MeanAP$x.sd / sqrt(MeanAP$x.n)

MeanAP

colnames(MeanAP) <- c("Treatment", "mean", "sd", "n", "se")

MeanAP$Treatment <- factor(MeanAP$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanAP$Treatment

MeanAP

bar_plotAP <- ggplot(MeanAP, aes(x=as.factor(MeanAP$Treatment), y=MeanAP$mean, fill=as.factor(MeanAP$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1) +
  geom_errorbar(aes(ymin=MeanAP$mean-MeanAP$se, ymax=MeanAP$mean + MeanAP$se), width=0.2, size=1, position = position_dodge(width=0.9)) +
  annotate("text", label="ab", x=1, y=(10^4), size=5) +
  annotate("text", label="ab", x=2, y=10^3.99, size=5) +
  annotate("text", label="a", x=3, y=10^3.85, size=5) +
  annotate("text", label="b", x=4, y=10^4.16, size=5) +
  annotate("text", 2, 10^4.16, label="Treatment: ~italic(P) == 0.003", parse=TRUE, size=3.5) +
  labs(list(x ="", y =expression(""), title="AP")) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999"))+
  theme_classic() +
  theme(axis.text=element_text(size=12)) +
  theme(text=element_text(size=12)) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.position="none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) 
plot(bar_plotAP)

MeanAP2 <- aggregate(enz7$AP2, by = list(Treatment = enz7$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanAP2

MeanAP2 <- do.call(data.frame, MeanAP2)

MeanAP2

MeanAP2$se <- MeanAP2$x.sd / sqrt(MeanAP2$x.n)

MeanAP2$se1 <- 10^(MeanAP2$x.mean + MeanAP2$se)

MeanAP2$se2 <- 10^(MeanAP2$x.mean - MeanAP2$se)

MeanAP2$x.mean <- 10^(MeanAP2$x.mean)

MeanAP2

colnames(MeanAP2) <- c("Treatment", "mean", "sd", "n", "se", "se1", "se2")


MeanAP2$Treatment <- factor(MeanAP2$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanAP2$Treatment

bar_plotAP2 <- ggplot(MeanAP2, aes(x=as.factor(MeanAP2$Treatment), y=MeanAP2$mean, fill=as.factor(MeanAP2$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1, position = "dodge") +
  geom_errorbar(aes(ymin=MeanAP2$se1, ymax=MeanAP2$se2), width=0.2, size=1,  position = position_dodge(0.9)) +
  annotate("text", label="a", x=1, y=31000, size=5) +
  annotate("text", label="a", x=2, y=30000, size=5) +
  annotate("text", label="a", x=3, y=24000, size=5) +
  annotate("text", label="a", x=4, y=30000, size=5) +
  annotate("text", label="Treatment: P = 0.353", x=3, y=36000, size=3.5) + 
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999")) +
  labs(list(x ="", y =expression(""), title="AP")) +
  theme_classic() +
  theme(axis.text=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position="none")
plot(bar_plotAP2)

#NAG

MeanNAG <- aggregate(enz2$NAG, by = list(Treatment = enz2$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanNAG

MeanNAG <- do.call(data.frame, MeanNAG)

MeanNAG$se <- MeanNAG$x.sd / sqrt(MeanNAG$x.n)

MeanNAG$se1 <- 10^(MeanNAG$x.mean + MeanNAG$se)

MeanNAG$se2 <- 10^(MeanNAG$x.mean - MeanNAG$se)

MeanNAG$x.mean <- 10^(MeanNAG$x.mean)

MeanNAG

colnames(MeanNAG) <- c("Treatment", "mean", "sd", "n", "se", "se1", "se2")

MeanNAG

MeanNAG$Treatment <- factor(MeanNAG$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanNAG$Treatment

bar_plotNAG <- ggplot(MeanNAG, aes(x=as.factor(MeanNAG$Treatment), y=MeanNAG$mean, fill=as.factor(MeanNAG$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1) +
  geom_errorbar(aes(ymin=MeanNAG$se1, ymax=MeanNAG$se2), width=0.2, size=1) +
  annotate("text", label="a", x=1, y=1100, size=5) +
  annotate("text", label="a", x=2, y=1000, size=5) +
  annotate("text", label="a", x=3, y=970, size=5) +
  annotate("text", label="a", x=4, y=1500, size=5) +
  annotate("text", 2, 1500, label="Treatment: ~italic(P) == 0.269", parse=TRUE, size=3.5) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999")) +
  labs(list(x ="", y =expression(""), title="NAG")) +
  theme_classic() +
  theme(axis.text=element_text(size=12)) +
  theme(text=element_text(size=12)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position="none") 
plot(bar_plotNAG)

MeanNAG2 <- aggregate(enz7$NAG2, by = list(Treatment = enz7$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanNAG2

MeanNAG2 <- do.call(data.frame, MeanNAG2)

MeanNAG2

MeanNAG2$se <- MeanNAG2$x.sd / sqrt(MeanNAG2$x.n)

MeanNAG2$se1 <- 10^(MeanNAG2$x.mean + MeanNAG2$se)

MeanNAG2$se2 <- 10^(MeanNAG2$x.mean - MeanNAG2$se)

MeanNAG2$x.mean <- 10^(MeanNAG2$x.mean)

MeanNAG2

colnames(MeanNAG2) <- c("Treatment", "mean", "sd", "n", "se", "se1", "se2")

MeanNAG2$Treatment <- factor(MeanNAG2$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanNAG2$Treatment

bar_plotNAG2 <- ggplot(MeanNAG2, aes(x=as.factor(MeanNAG2$Treatment), y=MeanNAG2$mean, fill=as.factor(MeanNAG2$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1, position = "dodge") +
  geom_errorbar(aes(ymin=MeanNAG2$se1, ymax=MeanNAG2$se2), width=0.2, size=1,  position = position_dodge(0.9)) +
  annotate("text", label="a", x=1, y=4100, size=5) +
  annotate("text", label="a", x=2, y=3200, size=5) +
  annotate("text", label="a", x=3, y=3500, size=5) +
  annotate("text", label="a", x=4, y=3400, size=5) +
  annotate("text", label="Treatment: P = 0.750", x=3, y = 4300, size = 3.5) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999")) +
  labs(list(x ="", y = "", title="NAG")) +
  theme_classic() +
  theme(axis.text=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position="none") 
plot(bar_plotNAG2)

#LAP

MeanLAP <- aggregate(enz$LAP, by = list(Treatment = enz$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))

MeanLAP <- do.call(data.frame, MeanLAP)

MeanLAP$se <- MeanLAP$x.sd / sqrt(MeanLAP$x.n)

MeanLAP

colnames(MeanLAP) <- c("Treatment", "mean", "sd", "n", "se")

MeanLAP$Treatment <- factor(MeanLAP$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanLAP$Treatment

bar_plotLAP <- ggplot(MeanLAP, aes(x=as.factor(MeanLAP$Treatment), y=MeanLAP$mean, fill=as.factor(MeanLAP$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1) +
  geom_errorbar(aes(ymin=MeanLAP$mean-MeanLAP$se, ymax=MeanLAP$mean+MeanLAP$se), width=0.2, size=1) +
  annotate("text", label="ab", x=1, y=10^1.69, size=5) +
  annotate("text", label="ab", x=2, y=10^1.68, size=5) +
  annotate("text", label="a*", x=3, y=10^1.62, size=5) +
  annotate("text", label="b", x=4, y=10^1.75, size=5) +
  annotate("text", 2, 10^1.8, label="Treatment: ~italic(P) == 0.092", parse=TRUE, size=3.5) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999")) +
  labs(list(x ="", y ="", title="LAP")) +
  theme_classic() +
  theme(axis.text=element_text(size=12)) +
  theme(text=element_text(size=12)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position="none") 
plot(bar_plotLAP)

MeanLAP2 <- aggregate(enz7$LAP2, by = list(Treatment = enz7$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))

MeanLAP2 <- do.call(data.frame, MeanLAP2)

MeanLAP2$se <- MeanLAP2$x.sd / sqrt(MeanLAP2$x.n)

MeanLAP2$e1 <- 10^(MeanLAP2$x.mean + MeanLAP2$se)

MeanLAP2$se2 <- 10^(MeanLAP2$x.mean - MeanLAP2$se)

MeanLAP2$x.mean <- 10^(MeanLAP2$x.mean)

MeanLAP2

colnames(MeanLAP2) <- c("Treatment", "mean", "sd", "n", "se", "se1", "se2")

MeanLAP2$Treatment <- factor(MeanLAP2$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanLAP2$Treatment

bar_plotLAP2 <- ggplot(MeanLAP2, aes(x=as.factor(MeanLAP2$Treatment), y=MeanLAP2$mean, fill=as.factor(MeanLAP2$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=MeanLAP2$se1, ymax=MeanLAP2$se2), width=0.2, size=1, position = position_dodge(0.9)) +
  annotate("text", label="a", x=1, y=185, size=5) +
  annotate("text", label="a", x=2, y=160, size=5) +
  annotate("text", label="a", x=3, y=155, size=5) +
  annotate("text", label="a", x=4, y=155, size=5) +
  annotate("text", label="Treatment: P = 0.312", x = 3, y = 200, size = 3.5) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999")) +
  labs(list(x ="", y ="", title="LAP")) +
  theme_classic() +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
  theme(axis.text=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position="none") 
plot(bar_plotLAP2)


#POX

MeanPOX <- aggregate(enz2$POX, by = list(Treatment = enz2$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanPOX

MeanPOX <- do.call(data.frame, MeanPOX)

MeanPOX$se <- MeanPOX$x.sd / sqrt(MeanPOX$x.n)

MeanPOX$se1 <- 10^(MeanPOX$x.mean + MeanPOX$se)

MeanPOX$se2 <- 10^(MeanPOX$x.mean - MeanPOX$se)

MeanPOX$x.mean <- 10^(MeanPOX$x.mean)

MeanPOX

colnames(MeanPOX) <- c("Treatment", "mean", "sd", "n", "se", "se1", "se2")

MeanPOX

MeanPOX$Treatment <- factor(MeanPOX$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanPOX$Treatment

bar_plotPOX <- ggplot(MeanPOX, aes(x=as.factor(MeanPOX$Treatment), y=MeanPOX$mean, fill=as.factor(MeanPOX$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1) +
  geom_errorbar(aes(ymin=MeanPOX$se1, ymax=MeanPOX$se2), width=0.2, size=1) +
  annotate("text", label="a", x=1, y=10^3.96, size=5) +
  annotate("text", label="a", x=2, y=10^3.91, size=5) +
  annotate("text", label="a", x=3, y=10^3.99, size=5) +
  annotate("text", label="a", x=4, y=10^3.89, size=5) +
  annotate("text", 2, 10^4.07, label="Treatment: ~italic(P) == 0.961", parse=TRUE, size=3.5) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999")) +
  labs(list(x ="", y ="", title="POX")) +
  theme_classic() +
  theme(axis.text=element_text(size=12)) +
  theme(text=element_text(size=12)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position="none") 
plot(bar_plotPOX)

MeanPOX2 <- aggregate(enz7$POX2, by = list(Treatment = enz7$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanPOX2

MeanPOX2 <- do.call(data.frame, MeanPOX2)

MeanPOX2$se <- MeanPOX2$x.sd / sqrt(MeanPOX2$x.n)

MeanPOX2$se1 <- 10^(MeanPOX2$x.mean + MeanPOX2$se)

MeanPOX2$se2 <- 10^(MeanPOX2$x.mean - MeanPOX2$se)

MeanPOX2$x.mean <- 10^(MeanPOX2$x.mean)

MeanPOX2

colnames(MeanPOX2) <- c("Treatment", "mean", "sd", "n", "se", "se1", "se2")


MeanPOX2$Treatment <- factor(MeanPOX2$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanPOX2$Treatment

bar_plotPOX2 <- ggplot(MeanPOX2, aes(x=as.factor(MeanPOX2$Treatment), y=MeanPOX2$mean, fill=as.factor(MeanPOX2$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1, position = "dodge") +
  geom_errorbar(aes(ymin=MeanPOX2$se1, ymax=MeanPOX2$se2), width=0.2, size=1, position = position_dodge(0.9)) +
  annotate("text", label="a", x=1, y=39000, size=5) +
  annotate("text", label="a", x=2, y=30000, size=5) +
  annotate("text", label="a", x=3, y=40000, size=5) +
  annotate("text", label="a", x=4, y=25000, size=5) +
  annotate("text", label="Treatment: P = 0.214", x=2, y=48000, size=3.5) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999")) +
  labs(list(x ="", y ="", title="POX")) +
  theme_classic() +
  theme(axis.text=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position="none") 
plot(bar_plotPOX2)

#PER

MeanPER <- aggregate(enz2$PER, by = list(Treatment = enz2$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))

MeanPER <- do.call(data.frame, MeanPER)

MeanPER$se <- MeanPER$x.sd / sqrt(MeanPER$x.n)

MeanPER$se1 <- 10^(MeanPER$x.mean + MeanPER$se)

MeanPER$se2 <- 10^(MeanPER$x.mean - MeanPER$se)

MeanPER$x.mean <- 10^(MeanPER$x.mean)

MeanPER

colnames(MeanPER) <- c("Treatment", "mean", "sd", "n", "se", "se1", "se2")

MeanPER$Treatment <- factor(MeanPER$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanPER$Treatment

bar_plotPER <- ggplot(MeanPER, aes(x=as.factor(MeanPER$Treatment), y=MeanPER$mean, fill=as.factor(MeanPER$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1) +
  geom_errorbar(aes(ymin=MeanPER$se1, ymax=MeanPER$se2), width=0.2, size=1) +
  annotate("text", label="a", x=1, y=10^4.49, size=5) +
  annotate("text", label="a", x=2, y=10^4.37, size=5) +
  annotate("text", label="a", x=3, y=10^4.45, size=5) +
  annotate("text", label="a", x=4, y=10^4.47, size=5) +
  annotate("text", 2, 10^4.55, label="Treatment: ~italic(P) == 0.465", parse=TRUE, size=3.5) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999")) +
  labs(list(x ="", y ="", title="PER")) +
  theme_classic() +
  theme(axis.text=element_text(size=12)) +
  theme(text=element_text(size=12)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position="none") 
plot(bar_plotPER)

MeanPER2 <- aggregate(enz7$PER2, by = list(Treatment = enz7$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))

MeanPER2 <- do.call(data.frame, MeanPER2)

MeanPER2$se <- MeanPER2$x.sd / sqrt(MeanPER2$x.n)

MeanPER2$se1 <- 10^(MeanPER2$x.mean + MeanPER2$se)

MeanPER2$se2 <- 10^(MeanPER2$x.mean - MeanPER2$se)

MeanPER2$x.mean <- 10^(MeanPER2$x.mean)

MeanPER2

colnames(MeanPER2) <- c("Treatment", "mean", "sd", "n", "se", "se1", "se2")

MeanPER2$Treatment <- factor(MeanPER2$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanPER2$Treatment

bar_plotPER2 <- ggplot(MeanPER2, aes(x=as.factor(MeanPER2$Treatment), y=MeanPER2$mean, fill=as.factor(MeanPER2$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1, position = "dodge") +
  geom_errorbar(aes(ymin=MeanPER2$se1, ymax=MeanPER2$se2), width=0.2, size=1, position = position_dodge(0.9)) +
  annotate("text", label="a", x=1, y=130000, size=5) +
  annotate("text", label="a", x=2, y=85000, size=5) +
  annotate("text", label="a", x=3, y=110000, size=5) +
  annotate("text", label="a", x=4, y=80000, size=5) +
  annotate("text", label="Treatment: P = 0.161", x=3, y=130000, size=3.5) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999")) +
  labs(list(x ="", y ="", title="PER")) +
  theme_classic() +
  theme(axis.text=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position="none") 
plot(bar_plotPER2)

#BG

MeanBG <- aggregate(enz$BG, by = list(Treatment = enz$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanBG

MeanBG <- do.call(data.frame, MeanBG)

MeanBG$se <- MeanBG$x.sd / sqrt(MeanBG$x.n)

MeanBG

colnames(MeanBG) <- c("Treatment", "mean", "sd", "n", "se")

MeanBG$Treatment <- factor(MeanBG$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanBG$Treatment

bar_plotBG <- ggplot(MeanBG, aes(x=as.factor(MeanBG$Treatment), y=MeanBG$mean, fill=MeanBG$Treatment)) + 
  geom_bar(stat="identity", color="black", size=1) +
  geom_errorbar(aes(ymin=MeanBG$mean-MeanBG$se, ymax=MeanBG$mean+MeanBG$se), width=0.2, size=1) +
  annotate("text", label="a*", x=1, y=10^3.1, size=5) +
  annotate("text", label="a*", x=2, y=10^3.1, size=5) +
  annotate("text", label="a", x=3, y=10^3.03, size=5) +
  annotate("text", label="b", x=4, y=10^3.26, size=5) +
  annotate("text", 2, 10^3.23, label="Treatment: ~italic(P) == 0.003", parse=TRUE, size=3.5) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999")) +
  labs(list(x ="", y ="", title="BG")) +
  theme_classic() +
  theme(axis.text=element_text(size=12)) +
  theme(text=element_text(size=12)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position="none") 
plot(bar_plotBG)

MeanBG2 <- aggregate(enz7$BG2, by = list(Treatment = enz7$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanBG2

MeanBG2 <- do.call(data.frame, MeanBG2)

MeanBG2$se <- MeanBG2$x.sd / sqrt(MeanBG2$x.n)

MeanBG2$se1 <- 10^(MeanBG2$x.mean + MeanBG2$se)

MeanBG2$se2 <- 10^(MeanBG2$x.mean - MeanBG2$se)

MeanBG2$x.mean <- 10^(MeanBG2$x.mean)

MeanBG2

colnames(MeanBG2) <- c("Treatment", "mean", "sd", "n", "se", "se1", "se2")

MeanBG2$Treatment <- factor(MeanBG2$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanBG2$Treatment

bar_plotBG2 <- ggplot(MeanBG2, aes(x=as.factor(MeanBG2$Treatment), y=MeanBG2$mean, fill=as.factor(MeanBG2$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=MeanBG2$se1, ymax=MeanBG2$se2), width=0.2, size=1, position = position_dodge(0.9)) +
  annotate("text", label="a", x=1, y=4000, size=5) +
  annotate("text", label="a", x=2, y=3300, size=5) +
  annotate("text", label="a", x=3, y=3500, size=5) +
  annotate("text", label="a", x=4, y=4000, size=5) +
  annotate("text", label = "Treatment: P = 0.709", x = 3, y = 4500, size = 3.5) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999")) +
  labs(list(x ="", y ="", title="BG")) +
  theme_classic() +
  theme(axis.text=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position="none") 
plot(bar_plotBG2)

#XYL

MeanXYL <- aggregate(enz$XYL, by = list(Treatment = enz$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanXYL

MeanXYL <- do.call(data.frame, MeanXYL)

MeanXYL$se <- MeanXYL$x.sd / sqrt(MeanXYL$x.n)

colnames(MeanXYL) <- c("Treatment", "mean", "sd", "n", "se")

MeanXYL$Treatment <- factor(MeanXYL$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanXYL$Treatment

bar_plotXYL <- ggplot(MeanXYL, aes(x=as.factor(MeanXYL$Treatment), y=MeanXYL$mean, fill=as.factor(MeanXYL$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1) +
  geom_errorbar(aes(ymin=MeanXYL$mean-MeanXYL$se, ymax=MeanXYL$mean+MeanXYL$se), width=0.2, size=1) +
  annotate("text", label="a", x=1, y=435, size=5) +
  annotate("text", label="a", x=2, y=385, size=5) +
  annotate("text", label="a", x=3, y=340, size=5) +
  annotate("text", label="b", x=4, y=695, size=5) +
  annotate("text", 2, 600, label="Treatment: ~italic(P) < 0.001", parse=TRUE, size=3.5) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999")) +
  labs(list(x ="", y = "", title="XYL")) +
  theme_classic() +
  theme(axis.text=element_text(size=12)) +
  theme(text=element_text(size=12)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position="none")
plot(bar_plotXYL)

MeanXYL2 <- aggregate(enz7$XYL2, by = list(Treatment = enz7$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanXYL2

MeanXYL2 <- do.call(data.frame, MeanXYL2)

MeanXYL2$se <- MeanXYL2$x.sd / sqrt(MeanXYL2$x.n)

MeanXYL2$se1 <- 10^(MeanXYL2$x.mean + MeanXYL2$se)

MeanXYL2$se2 <- 10^(MeanXYL2$x.mean - MeanXYL2$se)

MeanXYL2$x.mean <- 10^(MeanXYL2$x.mean)

MeanXYL2

colnames(MeanXYL2) <- c("Treatment", "mean", "sd", "n", "se", "se1", "se2")

MeanXYL2$Treatment <- factor(MeanXYL2$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanXYL2$Treatment

bar_plotXYL2 <- ggplot(MeanXYL2, aes(x=as.factor(MeanXYL2$Treatment), y=MeanXYL2$mean, fill=as.factor(MeanXYL2$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=MeanXYL2$se1, ymax=MeanXYL2$se2), width=0.2, size=1, position = position_dodge(0.9)) +
  annotate("text", label="a", x=1, y=1450, size=5) +
  annotate("text", label="a", x=2, y=1100, size=5) +
  annotate("text", label="a", x=3, y=1050, size=5) +
  annotate("text", label="a", x=4, y=1500, size=5) +
  annotate("text", label="Treatment: P = 0.082", x=2, y=1700, size=3.5) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999")) +
  labs(list(x ="", y ="", title="XYL")) +
  theme_classic() +
  theme(axis.text=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position="none") 
plot(bar_plotXYL2)

#CHB

MeanCHB <- aggregate(enz2$CHB, by = list(Treatment = enz2$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanCHB

MeanCHB <- do.call(data.frame, MeanCHB)

MeanCHB$se <- MeanCHB$x.sd / sqrt(MeanCHB$x.n)

MeanCHB$se1 <- 10^(MeanCHB$x.mean + MeanCHB$se)

MeanCHB$se2 <- 10^(MeanCHB$x.mean - MeanCHB$se)

MeanCHB$x.mean <- 10^(MeanCHB$x.mean)

MeanCHB

colnames(MeanCHB) <- c("Treatment", "mean", "sd", "n", "se", "se1", "se2")

MeanCHB

MeanCHB$Treatment <- factor(MeanCHB$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanCHB$Treatment

bar_plotCHB <- ggplot(MeanCHB, aes(x=as.factor(MeanCHB$Treatment), y=MeanCHB$mean, fill=as.factor(MeanCHB$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1) +
  geom_errorbar(aes(ymin=MeanCHB$se1 , ymax=MeanCHB$se2), width=0.2, size=1) +
  annotate("text", label="ab", x=1, y=210, size=5) +
  annotate("text", label="ab", x=2, y=200, size=5) +
  annotate("text", label="a", x=3, y=160, size=5) +
  annotate("text", label="b", x=4, y=315, size=5) +
  annotate("text", 2, 300, label="Treatment: ~italic(P) == 0.025", parse=TRUE, size=3.5) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999")) +
  labs(list(x ="", y ="", title="CHB")) +
  theme_classic() +
  theme(axis.text=element_text(size=12)) +
  theme(text=element_text(size=12)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position="none") +
  theme(legend.text=element_text(size=20)) +
  theme(legend.title=element_text(colour="white")) 
plot(bar_plotCHB)

MeanCHB2 <- aggregate(enz7$CHB2, by = list(Treatment = enz7$Treatment), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanCHB2

MeanCHB2 <- do.call(data.frame, MeanCHB2)

MeanCHB2$se <- MeanCHB2$x.sd / sqrt(MeanCHB2$x.n)

MeanCHB2$se1 <- 10^(MeanCHB2$x.mean + MeanCHB2$se)

MeanCHB2$se2 <- 10^(MeanCHB2$x.mean - MeanCHB2$se)

MeanCHB2$x.mean <- 10^(MeanCHB2$x.mean)

MeanCHB2

colnames(MeanCHB2) <- c("Treatment", "mean", "sd", "n", "se", "se1", "se2")

MeanCHB2

MeanCHB2$Treatment <- factor(MeanCHB2$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanCHB2$Treatment

bar_plotCHB2 <- ggplot(MeanCHB2, aes(x=as.factor(MeanCHB2$Treatment), y=MeanCHB2$mean, fill=as.factor(MeanCHB2$Treatment))) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=MeanCHB2$se1, ymax=MeanCHB2$se2), width=0.2, size=1, position = position_dodge(0.9)) +
  annotate("text", label="a", x=1, y=750, size=5) +
  annotate("text", label="a", x=2, y=570, size=5) +
  annotate("text", label="a", x=3, y=520, size=5) +
  annotate("text", label="a", x=4, y=680, size=5) +
  annotate("text", label="Treatment: P = 0.301", x=2.5, y=700, size=3.5) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999")) +
  labs(list(x ="", y ="", title="CHB")) +
  theme_classic() +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
  theme(axis.text=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position="none") 
plot(bar_plotCHB2)

library(gridExtra)

library(grid)

library(ggpubr)

library(gtable)

library(cowplot)

jpeg(filename="Figure3.jpeg", bg="transparent", res=1500, units = "in", height=8, width=6.5)

p <- plot_grid(bar_plotBG, bar_plotCHB, bar_plotXYL, bar_plotAP, bar_plotLAP, bar_plotNAG, bar_plotPER, bar_plotPOX, ncol = 2, labels=c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'), align="hv")

q <- grid.arrange(p, left = textGrob(expression("Enzyme activity" ~~("nmol" ~"gdw"^-1* ~"h"^-1)), rot = 90, vjust = .5, gp=gpar(fontsize=20,font=8)))

dev.off() 

r <- plot_grid(bar_plotBG2, bar_plotCHB2, bar_plotXYL2, bar_plotAP2, bar_plotLAP2, bar_plotNAG2, bar_plotPER2, bar_plotPOX2, ncol = 2, labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'), align="hv")

s <- grid.arrange(r, left = textGrob(expression("Enzyme activity" ~~("nmol" ~"mg"^-1* ~"h"^-1)), rot = 90, vjust = .5, gp=gpar(fontsize=20,font=8)))

#Enzyme Ratio Plots

a2 <- aggregate(enz$CN, by = list(Treatment = enz$Treatment, Date = enz$Date), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
a2

a2 <- do.call(data.frame, a2)

a2$se <- a2$x.sd / sqrt(a2$x.n)


a2$Date <- as.character(a2$Date)
a2$Treatment <- as.character(a2$Treatment)


a2[a2=="4/1/2017"] <- "April 2017"
a2[a2=="7/1/2017"] <- "July 2017"

a2[a2=="REF"] <- "Reference"
a2[a2=="FF"] <- "Burn"
a2[a2=="CR"] <- "Cut"
a2[a2=="CFFR"] <- "Cut+Burn"

a2

a2$Treatment <- factor(a2$Treatment, levels = c("Reference", "Burn", "Cut", "Cut+Burn"))
a2$Treatment

plota2 <- ggplot(a2, aes(x=as.factor(a2$Date), y=a2$x.mean, fill=a2$Treatment)) + 
  geom_bar(stat = "identity", color="black", size=1, position = "dodge") +
  geom_errorbar(aes(ymin=a2$x.mean-a2$se, ymax=a2$x.mean + a2$se), width=0.2, size=1, position = position_dodge(0.9)) +
  scale_fill_manual(values=c("#999999", "#F0E442", '#56B4E9', '#D55E00')) +
  labs(list(x ="", y = "BG:(NAG + LAP)", title="")) +
  theme_classic() +
  annotate("text", label="a", x=.66, y=1.4, size=10) +
  annotate("text", label="ab", x=.881, y=1.9, size=10) +
  annotate("text", label="a", x=1.115, y=1.4, size=10) +
  annotate("text", label="b", x=1.33, y=2, size=10) +
  annotate("text", label="a", x=1.66, y=1.4, size=10) +
  annotate("text", label="a", x=1.881, y=1.2, size=10) +
  annotate("text", label="a", x=2.115, y=1.2, size=10) +
  annotate("text", label="a", x=2.34, y=1.1, size=10) +
  theme(text=element_text(size=30)) +
  theme(axis.text=element_text(size=30)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position=c(.75,.9)) +
  theme(legend.title=element_text(colour="white")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  coord_cartesian(ylim=c(0,2.2)) 
plot(plota2)

a3 <- aggregate(enz$CP, by = list(Treatment = enz$Treatment, Date = enz$Date), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
a3

a3 <- do.call(data.frame, a3)

a3$se <- a3$x.sd / sqrt(a3$x.n)

colnames(a3) <- c("Treatment", "date", "mean", "sd", "n", "se")

a3$date <- as.character(a3$date)

a3[a3=="4/1/2017"] <- "April 2017"
a3[a3=="7/1/2017"] <- "July 2017"

a3

a3$Treatment <- factor(a3$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
a3$Treatment

plota3 <- ggplot(a3, aes(x=as.factor(a3$date), y=a3$mean, fill=a3$Treatment)) + 
  geom_bar(stat = "identity", color="black", size=.25, position = "dodge") +
  geom_errorbar(aes(ymin=a3$mean-a3$se, ymax=a3$mean + a3$se), width=0.2, size=.25, position = position_dodge(0.9)) +
  scale_fill_manual(values=c("white", "gray70","gray28", "gray10")) +
  labs(list(x ="", y = "BG:AP", title="")) +
  theme_classic() +
  annotate("text", label="a", x=.66, y=.23, size=4) +
  annotate("text", label="a", x=.881, y=.17, size=4) +
  annotate("text", label="a", x=1.115, y=.19, size=4) +
  annotate("text", label="a", x=1.33, y=.28, size=4) +
  annotate("text", label="a", x=1.66, y=.16, size=4) +
  annotate("text", label="a", x=1.881, y=.155, size=4) +
  annotate("text", label="a", x=2.115, y=.18, size=4) +
  annotate("text", label="a", x=2.34, y=.195, size=4) +
  annotate("text", 2, .25, label="Treatment: ~italic(P) == 0.725", parse=TRUE, size=2.5, fontface="bold")+  
  annotate("text", 2, .235, label="Date: ~italic(P) < 0.001", parse=TRUE, size=2.5, fontface="bold")+           
  annotate("text", 2, .22, label = "Treatment %*% Date: ~italic(P) == 0.051", parse=TRUE, size=2.5, fontface="bold")+         
  theme(text=element_text(size=14)) +
  theme(axis.text=element_text(size=14)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  theme(legend.position="none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.title=element_text(colour="white")) 
plot(plota3)

FBreg <- ggplot(enz, aes(x=enz$F.B, y=enz$CN, fill=enz$CN)) + 
  geom_point(size = 1) +    # Use hollow circles
  geom_smooth(method=lm, color="black",size = .5) +
  labs(list(y ="BG:(NAG + LAP)", x ="ITS:16s")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text=element_text(size=14)) +
  theme(axis.text=element_text(size=11.5)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  theme(legend.position="none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  annotate("text", label="R^2==0.31~~~~italic(P) < 0.001", x=.14, y=2, size=4, parse=TRUE)
FBreg

FBreg2 <- ggplot(enz, aes(x=enz$F.B, y=enz$CP, fill=enz$CP)) + 
  geom_point(size = 1) +    
  geom_smooth(method=lm, color="black",size = .5) +
  labs(list(y ="BG:AP", x ="ITS:16s")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text=element_text(size=14)) +
  theme(axis.text=element_text(size=11.5)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  theme(legend.position="none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  annotate("text", label="R^2==0.07~~~~~italic(P)== 0.07", x=.14, y=.35, size=4, parse=TRUE)
FBreg2

#qPCR plots

MeanITS <- aggregate(enz$LogITS, by = list(Treatment = enz$Treatment, Date = enz$Date), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanITS

MeanITS <- do.call(data.frame, MeanITS)

MeanITS$se <- MeanITS$x.sd / sqrt(MeanITS$x.n)

MeanITS

colnames(MeanITS) <- c("Treatment","Date", "mean", "sd", "n", "se")

MeanITS$Date <- as.character(MeanITS$Date)

MeanITS[MeanITS=="4/1/2017"] <- "April 2017"
MeanITS[MeanITS=="7/1/2017"] <- "July 2017"

MeanITS$Treatment <- factor(MeanITS$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanITS$Treatment
MeanITS

bar_plotITS2 <- ggplot(MeanITS, aes(x=as.factor(MeanITS$Date), y=MeanITS$mean, fill=MeanITS$Treatment)) + 
  geom_bar(stat="identity", color="black", size=.25, position="dodge") +
  geom_errorbar(aes(ymin=MeanITS$mean-MeanITS$se, ymax=MeanITS$mean + MeanITS$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate("text", label="Treatment: P = 0.179 \n Date: P < 0.001 \n Treatment x Date: P = 0.934", x=.83, y=9.4, size=4) +
  scale_fill_manual(values=c("white", "gray70","gray28", "gray10")) +
  labs(list(x ="", y ="Log(Gene copies g soil-1)", title="ITS Gene Copies")) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(text=element_text(size=15)) +
  theme(axis.text=element_text(size=15)) +
  theme(axis.line.y=element_line(colour="white", size=.5)) + 
  theme(axis.line.y=element_line(colour="white", size=.5)) + 
  theme(legend.position=c(.5,.5)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  coord_cartesian(ylim=c(8.5,9.5)) 
plot(bar_plotITS2)


MeanBac <- aggregate(enz$LogBac, by = list(Treatment = enz$Treatment, Date = enz$Date), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanBac

MeanBac <- do.call(data.frame, MeanBac)

MeanBac$se <- MeanBac$x.sd / sqrt(MeanBac$x.n)

MeanBac$Date <- as.character(MeanBac$Date)

MeanBac[MeanBac=="4/1/2017"] <- "April 2017"
MeanBac[MeanBac=="7/1/2017"] <- "July 2017"

MeanBac

colnames(MeanBac) <- c("Treatment", "Date", "mean", "sd", "n", "se")

MeanBac$Treatment <- factor(MeanBac$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanBac$Treatment
MeanBac

bar_plotBac <- ggplot(MeanBac, aes(x=as.factor(MeanBac$Date), y=MeanBac$mean, fill=MeanBac$Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=MeanBac$mean-MeanBac$se, ymax=MeanBac$mean + MeanBac$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate("text", label="Treatment: P = 0.014 \n Date: P < 0.001 \n Treatment x Date: P = 0.619", x=2, y=10.6, size=4) +
    scale_fill_manual(values=c("white", "gray70","gray28", "gray10")) +
  labs(list(x ="", y ="Log(Gene copies g soil-1)", title="16s Gene Copies")) +
  theme_classic() +
  theme(text=element_text(size=16)) +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.position=c(1,1)) +
  coord_cartesian(ylim=c(10,10.7))
plot(bar_plotBac)

legend <- get_legend(bar_plotITS2)

MeanFB <- aggregate(enz$F.B, by = list(Treatment = enz$Treatment, Date = enz$Date), FUN =function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
MeanFB

MeanFB <- do.call(data.frame, MeanFB)

MeanFB$se <- MeanFB$x.sd / sqrt(MeanFB$x.n)

MeanFB$Date <- as.character(MeanFB$Date)

MeanFB[MeanFB=="4/1/2017"] <- "April 2017"
MeanFB[MeanFB=="7/1/2017"] <- "July 2017"

colnames(MeanFB) <- c("Treatment", "Date", "mean", "sd", "n", "se")

MeanFB$Treatment <- factor(MeanFB$Treatment, levels = c("REF", "FF", "CR", "CFFR"))
MeanFB

bar_plotFB <- ggplot(MeanFB, aes(x=as.factor(MeanFB$Date), y=MeanFB$mean, fill=MeanFB$Treatment)) + 
  geom_bar(stat="identity", color="black", size=.25, position="dodge") +
  geom_errorbar(aes(ymin=MeanFB$mean-MeanFB$se, ymax=MeanFB$mean + MeanFB$se), width=0.2, size=.25, position=position_dodge(0.9)) +
  annotate("text", .94, .14, label="Treatment: ~italic(P) == 0.804", parse=TRUE, size=2.5, fontface="bold")+ 
  annotate("text", .94, 0.13, label="Date: ~italic(P) < 0.001", parse=TRUE, size=2.5, fontface="bold")+         
  annotate("text", .94, 0.12, label="Treatment %*% Date: ~italic(P) == 0.82", parse=TRUE, size=2.5, fontface="bold")+        
  scale_fill_manual(values=c("white", "gray70","gray28", "gray10")) +
  annotate("text", label="a", x=.66, y=.087, size=4) +
  annotate("text", label="a", x=.881, y=.07, size=4) +
  annotate("text", label="a", x=1.115, y=.065, size=4) +
  annotate("text", label="a", x=1.33, y=.089, size=4) +
  annotate("text", label="a", x=1.66, y=.165, size=4) +
  annotate("text", label="a", x=1.881, y=.185, size=4) +
  annotate("text", label="a", x=2.115, y=.155, size=4) +
  annotate("text", label="a", x=2.34, y=.175, size=4) +
  labs(list(x ="", y ="ITS:16s", title="")) +
  theme_classic() +
  theme(text=element_text(size=14)) +
  theme(axis.text=element_text(size=14)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position="none") 
plot(bar_plotFB)

#Create panel figure 4

p <- plot_grid(plota2, plota3, bar_plotFB, NULL, FBreg, FBreg2, ncol = 2, labels=c('A', 'B', 'C', '', 'D', 'E'), align="hv")
p + draw_grob(legend, x = .25, y = 0)
