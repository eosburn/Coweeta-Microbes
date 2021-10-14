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

din <- read.csv("C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/drought_chem.csv")

din$Day <- as.factor(din$Day)

plfa1 <- read.csv("C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/PLFA_groups.csv")

plfa1$Day <- as.factor(plfa1$Day)

glm1 <- lmer(Total_FAME ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm1))
Anova(glm1)
em_glm1 <- emmeans(glm1, ~ LandUse * Drought | Day)
contrast(em_glm1, method = "pairwise")

fame <- aggregate(Total_FAME~Treatment+LandUse+Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(fame)

fame2 <- do.call(data.frame, fame)
fame2

fame2$se <- fame2$Total_FAME.sd / sqrt(fame2$Total_FAME.n)
head(fame2)

fame2$Treatment <- factor(fame2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
fame2$Treatment
head(fame2)

line_plot_FAME <- ggplot(fame2, aes(x=fame2$Day, y=fame2$Total_FAME.mean, group=fame2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=fame2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=fame2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=fame2$Total_FAME.mean-fame2$se, ymax=fame2$Total_FAME.mean+fame2$se), width=.25, size=.75) +
  annotate('text', 1.5, 60000, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 1.5, 50000, label="Land~Use: italic(P) < 0.001", size = 5, parse=TRUE)+
  annotate('text', 1.5, 40000, label="Land~Use %*% Drought: italic(P) == 0.559", size=5, parse=TRUE)+
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
  scale_y_continuous(expression("Total FAME")) 
plot(line_plot_FAME)

glm12 = lmer(MBC~ LandUse*Drought*Day+(1|Plot), data=din)
qqPlot(resid(glm12))
Anova(glm12)
as.data.frame(Anova(glm12))
write.csv(as.data.frame(Anova(glm12)),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\MBC_glm.csv", row.names = TRUE)


em_glm12 <- emmeans(glm12, ~ LandUse * Drought | Day)
contrast(em_glm12, method = "pairwise")
write.csv(as.data.frame(contrast(em_glm12, method = "pairwise")),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\MBC_emmeans.csv", row.names = TRUE)

MBC.0 <- aggregate(MBC~Drought+LandUse+Day, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
MBC.0

MBC.0 <- aggregate(MBC~Drought+Day, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
MBC.0

MBC.1 <- aggregate(MBC~Treatment+Day, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
MBC.1

MBC <- aggregate(MBC~Treatment+Drought+Day+LandUse, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(MBC)

MBC <- do.call(data.frame, MBC)
MBC

MBC$se <- MBC$MBC.sd / sqrt(MBC$MBC.n)
head(MBC)

MBC$Treatment <- factor(MBC$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
MBC$Treatment

line_plot_MBC <- ggplot(MBC, aes(x=MBC$Day, y=MBC$MBC.mean, group=MBC$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=MBC$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=MBC$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=MBC$MBC.mean-MBC$se, ymax=MBC$MBC.mean+MBC$se), width=.25, size=.75) +
  annotate('text', 3.5, 250, label="Drought: ~italic(P) == 0.007", size=5, parse=TRUE)+
  annotate('text', 3.5, 225, label="Land~Use: italic(P) < 0.001", size = 5, parse=TRUE)+
  annotate('text', 3.5, 200, label="Land~Use %*% Drought: italic(P) == 0.008", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.88,.81)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("Microbial Biomass C")) 
plot(line_plot_MBC)


plfa1$MBC <- din$MBC

plfa3 <- plfa1[-29,]

cor.test(plfa3$Total_FAME, plfa3$MBC)


MBreg <- ggplot(plfa3, aes(x=plfa3$MBC, y=plfa3$Total_FAME)) + 
  geom_smooth(method='lm', colour="black") +
  geom_point(size=3, aes(color=plfa3$Treatment)) +   
  scale_color_brewer(palette="Set1") +
  xlab('Microbial Biomass C') +
  theme_classic() +
  theme(text=element_text(size=16)) +
  theme(axis.text=element_text(size=15)) +
  theme(legend.title=element_blank()) +
  theme(legend.position=c(.8,.22)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) +
  ylab('Total FAME')+
  annotate("text", label="r == 0.553~~~~~~~P < 0.001", x=100, y=220000, size=5.5, parse=TRUE) 
MBreg


glm2 <- lmer(AM_Fungi ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm2))
Anova(glm2)
em_glm2 <- emmeans(glm2, ~ LandUse * Drought | Day)
contrast(em_glm2, method = "pairwise")

am <- aggregate(AM_Fungi~Treatment+LandUse+Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(am)

am2 <- do.call(data.frame, am)
am2

am2$se <- am2$AM_Fungi.sd / sqrt(am2$AM_Fungi.n)
head(am2)

am2$Treatment <- factor(am2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
am2$Treatment
head(am2)

line_plot_am <- ggplot(am2, aes(x=am2$Day, y=am2$AM_Fungi.mean, group=am2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=am2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=am2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=am2$AM_Fungi.mean-am2$se, ymax=am2$AM_Fungi.mean+am2$se), width=.25, size=.75) +
  annotate('text', 5, 5500, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 5, 5250, label="Land~Use: italic(P) < 0.001", size = 5, parse=TRUE)+
  annotate('text', 5, 5000, label="Land~Use %*% Drought: italic(P) == 0.209", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.14,.1)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("Arbuscular Mycorrhizae")) 
plot(line_plot_am)

plfa1$Total_Bac <- plfa1$Gram_neg + plfa1$Gram_pos

glm3 <- lmer(Gram_neg ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm3))
Anova(glm3)
em_glm3 <- emmeans(glm3, ~ LandUse * Drought | Day)
contrast(em_glm3, method = "pairwise")

neg <- aggregate(Gram_neg~Treatment+LandUse+Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(neg)

neg2 <- do.call(data.frame, neg)
neg2

neg2$se <- neg2$Gram_neg.sd / sqrt(neg2$Gram_neg.n)
head(neg2)

neg2$Treatment <- factor(neg2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
neg2$Treatment
head(neg2)

line_plot_neg <- ggplot(neg2, aes(x=neg2$Day, y=neg2$Gram_neg.mean, group=neg2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=neg2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=neg2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=neg2$Gram_neg.mean-neg2$se, ymax=neg2$Gram_neg.mean+neg2$se), width=.25, size=.75) +
  annotate('text', 5.5, 70000, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 5.5, 67000, label="Land~Use: italic(P) < 0.001", size = 5, parse=TRUE)+
  annotate('text', 5.5, 64000, label="Land~Use %*% Drought: italic(P) == 0.392", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.5,.95)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("Gram negative")) 
plot(line_plot_neg)

glm4 <- lmer(Gram_pos ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm4))
Anova(glm4)
em_glm4 <- emmeans(glm4, ~ LandUse * Drought | Day)
contrast(em_glm4, method = "pairwise")

pos <- aggregate(Gram_pos~Treatment+LandUse+Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(pos)

pos2 <- do.call(data.frame, pos)
pos2

pos2$se <- pos2$Gram_pos.sd / sqrt(pos2$Gram_pos.n)
head(pos2)

pos2$Treatment <- factor(pos2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
pos2$Treatment
head(pos2)

line_plot_pos <- ggplot(pos2, aes(x=pos2$Day, y=pos2$Gram_pos.mean, group=pos2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=pos2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=pos2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=pos2$Gram_pos.mean-pos2$se, ymax=pos2$Gram_pos.mean+pos2$se), width=.25, size=.75) +
  annotate('text', 5, 34000, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 5, 32000, label="Land~Use: italic(P) < 0.001", size = 5, parse=TRUE)+
  annotate('text', 5, 30000, label="Land~Use %*% Drought: italic(P) == 0.899", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.3,.95)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("Gram positive")) 
plot(line_plot_pos)

plfa1$Total_Bac <- plfa1$Gram_neg + plfa1$Gram_pos

glm5 <- lmer(Total_Bac ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm5))
Anova(glm5)
em_glm5 <- emmeans(glm5, ~ LandUse * Drought | Day)
contrast(em_glm5, method = "pairwise")

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
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.4,.95)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("Total bacteria")) 
plot(line_plot_bac)

plfa1$bac <- din$X16S

glm5 <- lmer(bac ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm5))
Anova(glm5)
em_glm5 <- emmeans(glm5, ~ LandUse * Drought | Day)
contrast(em_glm5, method = "pairwise")

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
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.2,.2)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("16S gene copies")) 
plot(line_plot_bac)

cor.test(plfa1$Total_Bac, plfa1$bac)

cor.test(plfa1$Gram_neg, plfa1$bac)

cor.test(plfa1$Gram_pos, plfa1$bac)

plfa_d1 <- plfa1[plfa1$Day == 1,]

cor.test(plfa_d1$bac,plfa_d1$Total_Bac)

cor.test(plfa_d1$bac,plfa_d1$Gram_pos)

cor.test(plfa_d1$bac,plfa_d1$Gram_neg)


Bacreg <- ggplot(plfa1, aes(x=plfa1$bac, y=plfa1$Total_Bac)) + 
  geom_smooth(method='lm', colour="black") +
  geom_point(size=3, aes(color=plfa1$Treatment)) +   
  scale_color_brewer(palette="Set1") +
  xlab('16S gene copies') +
  theme_classic() +
  theme(text=element_text(size=16)) +
  theme(axis.text=element_text(size=15)) +
  theme(legend.title=element_blank()) +
  theme(legend.position=c(.8,.15)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) +
  ylab('Total Bacteria (FAME)')+
  annotate("text", label="r == 0.13~~~~~~~P == 0.121", x=9, y=126000, size=5.5, parse=TRUE) 
Bacreg

plfa1$posneg <- plfa1$Gram_neg/plfa1$Gram_pos

glm5 <- lmer(posneg ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm5))
Anova(glm5)
em_glm5 <- emmeans(glm5, ~ LandUse * Drought | Day)
contrast(em_glm5, method = "pairwise")

posneg <- aggregate(posneg~Treatment+LandUse+Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(posneg)

posneg2 <- do.call(data.frame, posneg)
posneg2

posneg2$se <- posneg2$posneg.sd / sqrt(posneg2$posneg.n)
head(posneg2)

posneg2$Treatment <- factor(posneg2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
posneg2$Treatment
head(posneg2)

line_plot_posneg <- ggplot(posneg2, aes(x=posneg2$Day, y=posneg2$posneg.mean, group=posneg2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=posneg2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=posneg2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=posneg2$posneg.mean-posneg2$se, ymax=posneg2$posneg.mean+posneg2$se), width=.25, size=.75) +
  annotate('text', 5, 2.3, label="Drought: ~italic(P) < 0.477", size=5, parse=TRUE)+
  annotate('text', 5, 2.25, label="Land~Use: italic(P) == 0.015", size = 5, parse=TRUE)+
  annotate('text', 5, 2.2, label="Land~Use %*% Drought: italic(P) == 0.064", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.5,.8)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("Gram- : Gram+")) 
plot(line_plot_posneg)



plfa1$ITS <- din$ITS

plfa2 <- plfa1

plfa2 <- plfa2[-30,]


glm6 <- lmer(Total_Fungi ~ Drought*LandUse*Day + (1|Plot), data=plfa2)
qqPlot(resid(glm6))
Anova(glm6)
em_glm6 <- emmeans(glm6, ~ LandUse * Drought | Day)
contrast(em_glm6, method = "pairwise")

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
  annotate('text', 5, 6000, label="Drought: ~italic(P) == 0.740", size=5, parse=TRUE)+
  annotate('text', 5, 5800, label="Land~Use: italic(P) < 0.001", size = 5, parse=TRUE)+
  annotate('text', 5, 5600, label="Land~Use %*% Drought: italic(P) == 0.807", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.5,.65)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("Total Fungi")) 
plot(line_plot_fung)

glm7 <- lmer(ITS ~ Drought*LandUse*Day + (1|Plot), data=din)
qqPlot(resid(glm7))
Anova(glm7)
em_glm7 <- emmeans(glm7, ~ LandUse * Drought | Day)
contrast(em_glm7, method = "pairwise")

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
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.2,.15)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("ITS gene copies")) 
plot(line_plot_ITS)


cor.test(plfa2$ITS, plfa2$Total_Fungi)


Fungreg <- ggplot(plfa2, aes(x=plfa2$ITS, y=plfa2$Total_Fungi)) + 
  geom_smooth(method='lm', colour="black") +
  geom_point(size=3, aes(color=plfa2$Treatment)) +   
  scale_color_brewer(palette="Set1") +
  xlab('ITS gene copies') +
  theme_classic() +
  theme(text=element_text(size=16)) +
  theme(axis.text=element_text(size=15)) +
  theme(legend.title=element_blank()) +
  theme(legend.position=c(.2,.6)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) +
  ylab('Total Fungi (FAME)')+
  annotate("text", label="r == 0.498~~~~~~~P < 0.001", x=7.8, y=7000, size=5.5, parse=TRUE) 
Fungreg

plfa1$FB <- plfa1$Total_Fungi/plfa1$Total_Bac

plfa1$ITS.16S <- din$ITS.16S

glm7 <- lmer(FB ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm7))
Anova(glm7)
em_glm7 <- emmeans(glm7, ~ LandUse * Drought | Day)
contrast(em_glm7, method = "pairwise")

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
  annotate('text', 5, .07, label="Drought: ~italic(P) == 0.003", size=5, parse=TRUE)+
  annotate('text', 5, .065, label="Land~Use: italic(P) < 0.219", size = 5, parse=TRUE)+
  annotate('text', 5, .06, label="Land~Use %*% Drought: italic(P) == 0.855", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.6,.5)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("Fungal:Bacterial (FAME)")) 
plot(line_plot_FB)

glm7 <- lmer(ITS.16S ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm7))
Anova(glm7)
em_glm7 <- emmeans(glm7, ~ LandUse * Drought | Day)
contrast(em_glm7, method = "pairwise")

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
  annotate('text', 5, .17, label="Drought: ~italic(P) == 0.07", size=5, parse=TRUE)+
  annotate('text', 5, .16, label="Land~Use: italic(P) == .008", size = 5, parse=TRUE)+
  annotate('text', 5, .15, label="Land~Use %*% Drought: italic(P) == 0.821", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.2,.86)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("Fungal:Bacterial (qPCR)")) 
plot(line_plot_ITS.16S)


cor.test(plfa1$ITS.16S, plfa1$FB)


Fungreg <- ggplot(plfa1, aes(x=plfa1$ITS.16S, y=plfa1$FB)) + 
  geom_smooth(method='lm', colour="black") +
  geom_point(size=3, aes(color=plfa1$Treatment)) +   
  scale_color_brewer(palette="Set1") +
  xlab('Fungal:Bacterial (qPCR)') +
  theme_classic() +
  theme(text=element_text(size=16)) +
  theme(axis.text=element_text(size=15)) +
  theme(legend.title=element_blank()) +
  theme(legend.position=c(.7,.6)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) +
  ylab('Fungal:Bacterial (FAME)')+
  annotate("text", label="r == 0.179~~~~~~~P == 0.03", x=.1, y=.1, size=5.5, parse=TRUE) 
Fungreg

glm7 <- lmer(Actino ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm7))
Anova(glm7)
em_glm7 <- emmeans(glm7, ~ LandUse * Drought | Day)
contrast(em_glm7, method = "pairwise")

Actino <- aggregate(Actino~Treatment+LandUse+Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(Actino)

Actino2 <- do.call(data.frame, Actino)
Actino2

Actino2$se <- Actino2$Actino.sd / sqrt(Actino2$Actino.n)
head(Actino2)

Actino2$Treatment <- factor(Actino2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
Actino2$Treatment
head(Actino2)

line_plot_Actino <- ggplot(Actino2, aes(x=Actino2$Day, y=Actino2$Actino.mean, group=Actino2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=Actino2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=Actino2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=Actino2$Actino.mean-Actino2$se, ymax=Actino2$Actino.mean+Actino2$se), width=.25, size=.75) +
  annotate('text', 5, 18000, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 5, 17500, label="Land~Use: italic(P) < .001", size = 5, parse=TRUE)+
  annotate('text', 5, 17000, label="Land~Use %*% Drought: italic(P) == 0.268", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.3,.95)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("Actinobacteria")) 
plot(line_plot_Actino)

glm7 <- lmer(Protozoa ~ Drought*LandUse*Day + (1|Plot), data=plfa1)
qqPlot(resid(glm7))
Anova(glm7)
em_glm7 <- emmeans(glm7, ~ LandUse * Drought | Day)
contrast(em_glm7, method = "pairwise")

Protozoa <- aggregate(Protozoa~Treatment+LandUse+Drought+Day, data=plfa1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(Protozoa)

Protozoa2 <- do.call(data.frame, Protozoa)
Protozoa2

Protozoa2$se <- Protozoa2$Protozoa.sd / sqrt(Protozoa2$Protozoa.n)
head(Protozoa2)

Protozoa2$Treatment <- factor(Protozoa2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
Protozoa2$Treatment
head(Protozoa2)

line_plot_Protozoa <- ggplot(Protozoa2, aes(x=Protozoa2$Day, y=Protozoa2$Protozoa.mean, group=Protozoa2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=Protozoa2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=Protozoa2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=Protozoa2$Protozoa.mean-Protozoa2$se, ymax=Protozoa2$Protozoa.mean+Protozoa2$se), width=.25, size=.75) +
  annotate('text', 4, 900, label="Drought: ~italic(P) == 0.05", size=5, parse=TRUE)+
  annotate('text', 4, 825, label="Land~Use: italic(P) == 0.881", size = 5, parse=TRUE)+
  annotate('text', 4, 750, label="Land~Use %*% Drought: italic(P) == 0.465", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.3,.95)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("Protozoa")) 
plot(line_plot_Protozoa)

#Multivariate stuff

plfa_multi1 <- read.csv("C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/PLFA1.csv")
plfa_multi2 <- read.csv("C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/PLFA2.csv")

#Merging the two datasets (too large to download all of it in one file from MG-RAST)

plfa_multi <- merge(plfa_multi1, plfa_multi2,all.x=TRUE, all.y=TRUE)

plfa_multi[is.na(plfa_multi)] <- 0

plfa_multi <- plfa_multi[order(plfa_multi$ID_Only),]

plfa_multi <- plfa_multi[-c(1,14,27,40,53,66,79,92,105,118,131,144),]

plfa1 <- plfa1[order(plfa1$ID_Only),]

meta <- plfa1[,c(5:13)]

plfa_multi <- cbind(plfa_multi, meta)

plfa_matrix <- plfa_multi[,c(7:110)]

plfa_matrix2 <- decostand(plfa_matrix, "total")

rowSums(plfa_matrix2)

plfa_multi2 <- decostand(plfa_multi[,c(7:110)], "total")

a <- vegdist(plfa_multi2, method="bray")

perm <- how(nperm=999)
setBlocks(perm) <- with(plfa_multi, Plot)
adonis2(a~LandUse * Drought, data=plfa_multi, permutations=perm)

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

ord_2 <- ggplot(multi_ord6, aes(x=MDS1.mean, y=MDS2.mean,color=Treatment,shape=Day)) + 
  geom_errorbar(aes(ymin=multi_ord6$MDS2.mean-multi_ord6$MDS2se, ymax=multi_ord6$MDS2.mean+multi_ord6$MDS2se), width=0.01, size=.75, position=position_dodge(0.9)) +
  geom_errorbarh(aes(xmin=multi_ord6$MDS1.mean-multi_ord6$MDS1se, xmax=multi_ord6$MDS1.mean+multi_ord6$MDS1se), height=0.01, size=.75) +
  geom_point(data=multi_ord6, aes(x=MDS1.mean,y=MDS2.mean,color=Treatment, shape=Day),size=5) + # add the point markers
  annotate("text", -.07, -.03, label="Land~Use: ~italic(P) == 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.07, -.06, label="'2D Stress' == 0.13", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.07, -.04, label="Drought: ~italic(P) == 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.07, -.05, label="`Land Use x Drought`: ~italic(P) == 0.324", parse=TRUE, size=4.5, fontface="bold")+
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
  ggtitle("All PLFAs") +
  theme(legend.text=element_text(size=15)) 
ord_2

plfa_multi1 <- read.csv("C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/PLFA1_bac.csv")
plfa_multi2 <- read.csv("C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/PLFA2_bac.csv")

#Merging the two datasets (too large to download all of it in one file from MG-RAST)

plfa_multi <- merge(plfa_multi1, plfa_multi2,all.x=TRUE, all.y=TRUE)

plfa_multi[is.na(plfa_multi)] <- 0

plfa_multi <- plfa_multi[order(plfa_multi$ID_Only),]

plfa_multi <- plfa_multi[-c(1,14,27,40,53,66,79,92,105,118,131,144),]

plfa1 <- plfa1[order(plfa1$ID_Only),]

meta <- plfa1[,c(5:13)]

plfa_multi <- cbind(plfa_multi, meta)

plfa_matrix <- plfa_multi[,c(7:110)]

plfa_matrix2 <- decostand(plfa_matrix, "total")

rowSums(plfa_matrix2)

plfa_multi2 <- decostand(plfa_multi[,c(7:110)], "total")

a <- vegdist(plfa_multi2, method="bray")

perm <- how(nperm=999)
setBlocks(perm) <- with(plfa_multi, Plot)
adonis2(a~LandUse * Drought, data=plfa_multi, permutations=perm)

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

ord_2 <- ggplot(multi_ord6, aes(x=MDS1.mean, y=MDS2.mean,color=Treatment,shape=Day)) + 
  geom_errorbar(aes(ymin=multi_ord6$MDS2.mean-multi_ord6$MDS2se, ymax=multi_ord6$MDS2.mean+multi_ord6$MDS2se), width=0.01, size=.75, position=position_dodge(0.9)) +
  geom_errorbarh(aes(xmin=multi_ord6$MDS1.mean-multi_ord6$MDS1se, xmax=multi_ord6$MDS1.mean+multi_ord6$MDS1se), height=0.01, size=.75) +
  geom_point(data=multi_ord6, aes(x=MDS1.mean,y=MDS2.mean,color=Treatment, shape=Day),size=5) + # add the point markers
  annotate("text", -.07, -.03, label="Land~Use: ~italic(P) == 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.07, -.06, label="'2D Stress' == 0.13", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.07, -.04, label="Drought: ~italic(P) == 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.07, -.05, label="`Land Use x Drought`: ~italic(P) == 0.319", parse=TRUE, size=4.5, fontface="bold")+
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
  ggtitle("Bacterial PLFAs") +
  theme(legend.text=element_text(size=15)) 
ord_2

#Stress

plfa_stress1<- read.csv("C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/PLFA1_stress.csv")
plfa_stress2 <- read.csv("C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/PLFA2_stress.csv")

#Merging the two datasets (too large to download all of it in one file from MG-RAST)

plfa_stress <- merge(plfa_stress1, plfa_stress2,all.x=TRUE, all.y=TRUE)

plfa_stress[is.na(plfa_stress)] <- 0

plfa_stress <- plfa_stress[order(plfa_stress$ID_Only),]

plfa_stress <- plfa_stress[-c(1,14,27,40,53,66,79,92,105,118,131,144),]

plfa1 <- plfa1[order(plfa1$ID_Only),]

meta <- plfa1[,c(5:13)]

plfa_stress <- cbind(plfa_stress, meta)

glm7 <- glmer(stress ~ Drought*LandUse*Day + (1|Plot), data=plfa_stress, family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(glm7))
Anova(glm7)

stress <- aggregate(stress~Treatment+LandUse+Drought+Day, data=plfa_stress, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(stress)

stress2 <- do.call(data.frame, stress)
stress2

stress2$se <- stress2$stress.sd / sqrt(stress2$stress.n)
head(stress2)

stress2$Treatment <- factor(stress2$Treatment, levels = c("Reference-Control", "Reference-Drought", "Disturbed-Control", "Disturbed-Drought"))
stress2$Treatment
head(stress2)

line_plot_stress <- ggplot(stress2, aes(x=stress2$Day, y=stress2$stress.mean, group=stress2$Treatment)) + 
  geom_point(position=position_dodge(width=0.25),aes(color=stress2$Treatment), size=4) +
  geom_line(position=position_dodge(width=0.25),aes(color=stress2$Treatment),size=3) +
  geom_errorbar(position=position_dodge(width=0.25),aes(ymin=stress2$stress.mean-stress2$se, ymax=stress2$stress.mean+stress2$se), width=.25, size=.75) +
  #annotate('text', 4, 900, label="Drought: ~italic(P) == 0.05", size=5, parse=TRUE)+
  #annotate('text', 4, 825, label="Land~Use: italic(P) == 0.881", size = 5, parse=TRUE)+
  #annotate('text', 4, 750, label="Land~Use %*% Drought: italic(P) == 0.465", size=5, parse=TRUE)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.7,.85)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(expression("PLFA Stress Index")) 
plot(line_plot_stress)
