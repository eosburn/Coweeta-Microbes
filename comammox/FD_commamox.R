#Import and prepare data 

setwd("C:/Users/ernie/Desktop/")

soil <- read.csv("C:/Users/ernie/Desktop/Data/Chapter 2 Data/FD_commamox.csv")

names(soil)[names(soil) == "C.N"] <- "C:N"

names(soil)[names(soil) == "DOCC.TDN"] <- "DOC:TDN"

library(lme4)
library(car)

#Statistical tests for disturbance effects on amoA groups and linear regressions with soil NO3

m1 <- lmer(logAOB~ Use+(1|Pair), data=soil)
shapiro.test(resid(m1))
Anova(m1)
m1.1 <- lm(log10(NO3) ~ logAOB, data=soil)
summary(m1.1)
plot(log10(soil$NO3), soil$logAOB)
a1.1 <- aggregate(soil$AOB ~ soil$Use, FUN = mean)
a1.1
#AOB is 5-fold higher in disturbed

m2 <- lmer(logAOA~ Use+(1|Pair), data=soil)
shapiro.test(resid(m2))
Anova(m2)
m2.1 <- lm(log10(NO3) ~ logAOA, data=soil)
summary(m2.1)
plot(log10(soil$NO3), soil$logAOA)
a2.1 <- aggregate(soil$AOA ~ soil$Use, FUN = mean)
a2.1
#AOA is 2.6-fold higher in disturbed

m3 <- lmer(logCAOB~ Use+(1|Pair), data=soil)
shapiro.test(resid(m3))
summary(m3)
m3.1 <- lm(log10(NO3) ~ logCAOB, data=soil)
summary(m3.1)
plot(log10(soil$NO3), soil$logCAOB)
a3.1 <- aggregate(soil$CAOB ~ soil$Use, FUN=mean)
a3.1
#CAOB is 2.9-fold higher in disturbed

m4.1 <- lm(log10(NO3) ~ logAOM, data=soil)
summary(m4.1)
plot(log10(soil$NO3), soil$logAOM)

m4.2 <- lm(log10(NO3) ~ logTAOB, data=soil)
summary(m4.2)
plot(log10(soil$NO3), soil$logTAOB)

m4.3 <- lm(log10(NO3) ~ logAOAAOB, data=soil)
summary(m4.3)
plot(log10(soil$NO3), soil$logAOAAOB)

#Model selection to determine best combination of amoA gropus to predict soil NO3

library(AICcmodavg)

mreg1 <- lm(log10(NO3) ~ logCAOB, data = soil)
summary(mreg1)
mreg2 <- lm(log10(NO3) ~ logAOB, data = soil)
summary(mreg2)
mreg3 <- lm(log10(NO3) ~ logAOA, data = soil)
summary(mreg3)
mreg4 <- lm(log10(NO3) ~ logAOA+logAOB, data = soil)
summary(mreg4)
mreg5 <- lm(log10(NO3) ~ logAOA+logCAOB, data = soil)
summary(mreg5)
mreg6 <- lm(log10(NO3) ~ logCAOB+logAOB, data = soil)
summary(mreg6)
mreg7 <- lm(log10(NO3) ~ logAOB+logCAOB+logAOA, data = soil)
summary(mreg7)

aictab(cand.set=list(mreg1,mreg2,mreg3,mreg4,mreg5,mreg6,mreg7))

#Correlations for amoA groups and soil variables

library(Hmisc)

soil2 <- soil[,c(5:21)]

cor1 <- rcorr(as.matrix(soil2), type="spearman")

cor1

round(cor1$r, 3)

round(cor1$P, 3)

#amoA box plots

CAOBbox <- ggplot(soil, aes(x=Use, y=logCAOB)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.size=10) +
  geom_point(position=position_jitter(.1),size=2) +
  labs(list(x ="")) +
  theme_classic() +
  theme(axis.text=element_text(size=18)) +
  theme(text=element_text(size=14)) +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  theme(axis.title.x=element_blank()) +
  annotate("text", .7, 7.7, label="CAOB", parse=TRUE, size=7, fontface="bold")+
  annotate("text", 2, 7.5, label="Disturbance: ~italic(P) < 0.001", parse=TRUE, size=5, fontface="bold")+  
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  scale_y_continuous(bquote(~'Log'[10]~'('~italic(amoA)~~ 'gene copies'~~gdw^-1*')')) +
  coord_cartesian(ylim=c(6,7.7)) 
CAOBbox

AOBbox <- ggplot(soil, aes(x=Use, y=logAOB)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.size=.25) +
  geom_point(position=position_jitter(.1),size=2) +
  labs(list(x ="")) +
  theme_classic() +
  theme(axis.text=element_text(size=18)) +
  theme(text=element_text(size=14)) +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  theme(axis.title.x=element_blank()) +
  annotate("text", .65, 5.5, label="AOB", parse=TRUE, size=7, fontface="bold")+
  annotate("text", 2, 5, label="Disturbance: ~italic(P) == 0.007", parse=TRUE, size=5, fontface="bold")+  
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  scale_y_continuous(bquote(~'Log'[10]~'('~italic(amoA)~~ 'gene copies'~~gdw^-1*')')) +
  coord_cartesian(ylim=c(1.75,5.5)) 
AOBbox

AOAbox <- ggplot(soil, aes(x=Use, y=logAOA)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.size=.25) +
  geom_point(position=position_jitter(.1),size=2) +
  labs(list(x ="")) +
  theme_classic() +
  theme(axis.text=element_text(size=18)) +
  theme(text=element_text(size=14)) +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  theme(axis.title.x=element_blank()) +
  annotate("text", .65, 7.75, label="AOA", parse=TRUE, size=7, fontface="bold")+
  annotate("text", 2, 7.25, label="Disturbance: ~italic(P) == 0.013", parse=TRUE, size=5, fontface="bold")+  
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) + 
  scale_y_continuous(bquote(~'Log'[10]~'('~italic(amoA)~~ 'gene copies'~~gdw^-1*')')) +
  coord_cartesian(ylim=c(3.5,7.75)) 
AOAbox

#Create new dataframe to make a barplot with all 3 amoA genes

library(anchors)

b1 <- as.data.frame(soil[,c(16,22)])

names(b1)[names(b1) == "logAOA"] <- "logamoA"

names(b1)[names(b1) == "AOA"] <- "amoA"

b1$gene <- "gene"

b1 <- replace.value(b1, "gene", from="gene", to="AOA")

b2 <- as.data.frame(soil[,c(17,23)])

names(b2)[names(b2) == "logAOB"] <- "logamoA"

names(b2)[names(b2) == "AOB"] <- "amoA"

b2$gene <- "gene"

b2 <- replace.value(b2, "gene", from="gene", to="AOB")

b3 <- as.data.frame(soil[,c(18,24)])

names(b3)[names(b3) == "logCAOB"] <- "logamoA"

names(b3)[names(b3) == "CAOB"] <- "amoA"

b3$gene <- "gene"

b3 <- replace.value(b3, "gene", from="gene", to="CAOB")

b4 <- as.data.frame(rbind(b1,b2,b3))

a <- aggregate(b4$logamoA, by = list(gene=b4$gene), FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x))) 
a

a <- do.call(data.frame, a)

a$se <- a$x.sd / sqrt(a$x.n)

a

a$gene <- factor(a$gene, levels = c("CAOB", "AOB", "AOA"))
a$gene

#ANOVA and pairwise comparison for amoA abundance 

library(emmeans)

amoa <- aov(logamoA ~ gene, data=b4)
shapiro.test(resid(amoa))
Anova(amoa)
em_amoa <- emmeans(amoa, ~ gene)
contrast(em_amoa, method = "pairwise")

a.1 <- aggregate(b4$amoA ~ b4$gene, FUN = mean)
a.1
#CAOB is 2.2-fold higher than AOA over all and 573-fold higher than AOB 
#AOA is 261-fold higher than AOB

#Create figure 1

jpeg(filename="Figure1.jpeg", bg="transparent", res=300, units = "in", height=4, width=5) 

bar_plotamoa <- ggplot(a, aes(x=a$gene, y=a$x.mean, fill=a$gene)) + 
  geom_bar(stat="identity", color="black", size=.25, position="dodge") +
  geom_errorbar(aes(ymin=a$x.mean-a$se, ymax=a$x.mean+a$se), width=0.2, size=.25, position=position_dodge(0.9)) +
  theme_classic() +
  scale_fill_manual(values=c("grey25","white","gray80")) +
  theme(legend.position="none") +
  annotate("text", label="a", x=1, y=7.1, size=5) +
  annotate("text", label="b", x=2, y=4, size=5) +
  annotate("text", label="c", x=3, y=6.4, size=5) +
  theme(legend.title = element_blank()) +
  theme(axis.title=element_text(size=14)) +
  theme(text=element_text(size=17)) +
  theme(axis.line.x=element_line(colour="black", size=.25)) + 
  theme(axis.line.y=element_line(colour="black", size=.25)) +
  labs(list(x ="")) +
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(bquote(~'Log'[10]~'('~italic(amoA)~~ 'gene copies'~~gdw^-1*')')) +
  coord_cartesian(ylim=c(2,7)) 
plot(bar_plotamoa)

dev.off()

#Linear Regression plots

CAOBreg <- ggplot(soil, aes(x=soil$logCAOB, y=log10(soil$NO3))) + 
  geom_point(size = 2) +    
  geom_smooth(method='lm', colour="black") +
  scale_x_continuous(bquote(~'Log'[10]~'(CAOB'~italic(amoA)~~ 'gene copies'~~gdw^-1*')')) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text=element_text(size=16)) +
  theme(axis.text=element_text(size=15)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) +
  theme(legend.position="none") +
  ylab(bquote(~'Log'[10]~(mu*'g'~'NO'[3]~'-N'~~gdw^-1))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  annotate("text", label="R^2 == 0.65~~~~~~~P < 0.001", x=6.55, y=1, size=5.5, parse=TRUE) 
CAOBreg

AOBreg <- ggplot(soil, aes(x=soil$logAOB, y=log10(soil$NO3))) + 
  geom_point(size = 2) +    
  geom_smooth(method='lm', colour="black") +
  scale_x_continuous(bquote(~'Log'[10]~'(AOB'~italic(amoA)~~ 'gene copies'~~gdw^-1*')')) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text=element_text(size=16)) +
  theme(axis.text=element_text(size=15)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) +
  theme(legend.position="none") +
  ylab(bquote(~'Log'[10]~(mu*'g'~'NO'[3]~'-N'~~gdw^-1))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  annotate("text", label="R^2 == 0.26~~~~~~~P < 0.001", x=2.9, y=1, size=5.5, parse=TRUE) 
AOBreg

AOAreg <- ggplot(soil, aes(x=soil$logAOA, y=log10(soil$NO3))) + 
  geom_point(size = 2) +    
  geom_smooth(method='lm', colour="black") +
  scale_x_continuous(bquote(~'Log'[10]~'(AOA'~italic(amoA)~~ 'gene copies'~~gdw^-1*')')) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text=element_text(size=16)) +
  theme(axis.text=element_text(size=15)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) +
  theme(legend.position="none") +
  ylab(bquote(~'Log'[10]~(mu*'g'~'NO'[3]~'-N'~~gdw^-1))) +
  annotate("text", label="R^2 == 0.32~~~~~~~P < 0.001", x=4.6, y=1, size=5.5, parse=TRUE) 
AOAreg

#Make panel figure 2

library(cowplot)

jpeg(filename="Figure2.jpeg", bg="transparent", res=300, units = "in", height=9, width=15) 

f1 <- plot_grid(CAOBbox, AOBbox, AOAbox, CAOBreg, AOBreg, AOAreg, ncol = 3, labels=c('A', 'B','C', 'D', 'E', 'F'), label_y = 1.02, align="hv", label_size=20)

f1

dev.off()
