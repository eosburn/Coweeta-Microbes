#Analyze soil properties in pre-treatment soils

setwd("C:/Users/ernie/OneDrive/Desktop/Coweeta_drought")

t0 <- read.csv("C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/Drought_t0.csv")

library(car)

#Analyze Pre-Treatment soil properties

Anova(aov(t0$pH ~ t0$LandUse))
sc1 <- aggregate(pH~LandUse, data=t0, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc1

Anova(aov(t0$NH4 ~ t0$LandUse))
sc2 <- aggregate(NH4~LandUse, data=t0, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc2

Anova(aov(t0$NO3 ~ t0$LandUse))
sc3 <- aggregate(NO3~LandUse, data=t0, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc3

Anova(aov(t0$TDN ~ t0$LandUse))
sc4 <- aggregate(TDN~LandUse, data=t0, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc4

Anova(aov(t0$DON ~ t0$LandUse))
sc4 <- aggregate(DON~LandUse, data=t0, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc4

Anova(aov(t0$DOC ~ t0$LandUse))
sc5 <- aggregate(DOC~LandUse, data=t0, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc5

Anova(aov(t0$DOC.TDN ~ t0$LandUse))
sc6 <- aggregate(DOC.TDN~LandUse, data=t0, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc6

Anova(aov(t0$MBC ~ t0$LandUse))
sc7 <- aggregate(MBC~LandUse, data=t0, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc7

Anova(aov(t0$TC ~ t0$LandUse))
sc8 <- aggregate(TC~LandUse, data=t0, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc8

Anova(aov(t0$TN ~ t0$LandUse))
sc9 <- aggregate(TN~LandUse, data=t0, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc9

Anova(aov(t0$CN ~ t0$LandUse))
sc10 <- aggregate(CN~LandUse, data=t0, FUN =function(x) c(mean=mean(x), stderr=sd(x)/sqrt(length(x))))
sc10

#Water content analyses and visualization

water <- read.csv("C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/Drought_water.csv")

water2 <- aggregate(GVM~Treatment+Day, data=water, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(water2)

water2 <- do.call(data.frame, water2)
water2

water2$se <- water2$GVM.sd / sqrt(water2$GVM.n)
head(water2)

water2$Treatment <- as.character(water2$Treatment)

water2[water2=="Drought"] <- "Drought-Rewet"

library(ggplot2)

jpeg(filename="water.jpeg", bg="transparent", res=1500, units = "in", height=4, width=8) 

H2Otime1 <- ggplot(water2, aes(x=Day, y=GVM.mean,  shape=Treatment)) + 
  geom_point(size=3) +
  geom_line(size=.5) +
  geom_hline(yintercept=.126, linetype="dashed", color = "black",size=1.5) +
  theme_classic() +
  theme(axis.title=element_text(size=16)) +
  annotate("text", 70,.145 , label="Wilting Point", size=5)+
  theme(text=element_text(size=18)) + 
  theme(legend.position=c(.75,.7)) +
  scale_y_continuous(bquote(~'Water Content'~~ ('g H'[2]~'O'~~gdw^-1))) 
H2Otime1

dev.off()

library(lme4)
library(car)
library(emmeans)

#Cumulative Respiration responses

resp2 <- read.csv("C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/Drought_resp2.csv")

glm2 <- lmer(Resp~Drought*LandUse+(1|Plot), data=resp2)
qqPlot(resid(glm2))
shapiro.test(resid(glm2))

#use glm instead, residuals deviate from normality

glm2 <- glmer(Resp~LandUse*Drought+(1|Plot), data=resp2, family=Gamma(link=log))

#Note: The process above was used throughout the code, though hereafter the original lmer models are deleted when residuals were found to deviate from normality to reduce confusion about which final model was reported

Anova(glm2)
em_glm2 <- emmeans(glm2, ~ LandUse * Drought)
contrast(em_glm2, method = "pairwise")

resp2.0 <- aggregate(Resp~Drought, data=resp2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(resp2.0)

resp2.0 <- aggregate(Resp~Treatment, data=resp2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(resp2.0)

resp2$Treatment <- factor(resp2$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
resp2$Treatment
head(resp2)

box_plotResp2 <- ggplot(resp2, aes(x=resp2$Treatment, y=resp2$Resp, fill=resp2$Treatment)) + 
  stat_boxplot(geom = "errorbar", size=.5,width = 0.5, position=position_dodge(.75)) +
  geom_boxplot(inherit.aes = TRUE,size=.5, outlier.size=.25) +
  stat_summary(fun.y=mean,geom="point",shape=23,size=3,color="black",fill="yellow")+
  #annotate('text', x=8, y=70, label="Drought: ~italic(P) == 0.020", size=4, parse=TRUE)+
  #annotate('text', 8, 68, label="Disturbance: italic(P) == 0.782", size = 4, parse=TRUE)+
  #annotate('text', 8, 66, label="Disturbance %*% Drought: italic(P) == 0.143", size=4, parse=TRUE)+
  annotate("text", label="ab", x=1, y=70, size=5) +
  annotate("text", label="ab", x=2, y=70, size=5) +
  annotate("text", label="a", x=3, y=70, size=5) +
  annotate("text", label="b", x=4, y=70, size=5) +
  labs(list(x ="")) +
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank()) +
  ylab(bquote(~mu*'g CO'[2]~'-C'~~gdw^-1)) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_blank())+
  coord_cartesian(ylim=c(43,70),clip='off') 
plot(box_plotResp2)

#Respiration time series analysis

resp <- read.csv("C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/Drought_resp.csv")

resp$Day <- as.factor(resp$Day)
glm3 <- glmer(Cmin~Drought*LandUse*Day+(1|Plot), data=resp, family=Gamma(link=log),nAGQ=0,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1000000)))
Anova(glm3,type="III")
glm3.1 <- aov(Cmin~Drought*LandUse*Day+Error(Plot), data=resp)
summary(glm3.1)
em_glm3 <- emmeans(glm3, ~ LandUse * Drought | Day)
as.data.frame(contrast(em_glm3, method = "pairwise"))

resp2.2 <- aggregate(Cmin~Treatment+LandUse+Drought+Day, data=resp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(resp2.2)

resp2.3 <- aggregate(Cmin~Drought+Day, data=resp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
resp2.3

resp2.3 <- aggregate(Cmin~Treatment+Day, data=resp, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
resp2.3

resp2.2 <- do.call(data.frame, resp2.2)
resp2.2

resp2.2$se <- resp2.2$Cmin.sd / sqrt(resp2.2$Cmin.n)
head(resp2.2)

resp2.2$Treatment <- factor(resp2.2$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
resp2.2$Treatment
head(resp2.2)

line_plotCmin <- ggplot(resp2.2, aes(x=resp2.2$Day, y=resp2.2$Cmin.mean, group=resp2.2$Treatment)) + 
  geom_point(position=position_dodge(width=0.9),aes(color=resp2.2$Treatment), size=3) +
  geom_line(position=position_dodge(width=0.9),aes(color=resp2.2$Treatment),size=1) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=resp2.2$Cmin.mean-resp2.2$se, ymax=resp2.2$Cmin.mean+resp2.2$se), width=3, size=.75) +
  annotate('text', 70, 3.75, label="Drought: ~italic(P) == 0.018", size=4, parse=TRUE)+
  annotate('text', 70, 3.5, label="Disturbance: italic(P) == 0.77", size = 4, parse=TRUE)+
  annotate('text', 70, 3.25, label="Disturbance %*% Drought: italic(P) == 0.51", size=4, parse=TRUE)+
  annotate('text', 44, 3.9, label="*", size=9)+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.25,.85)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  ylab(bquote(~mu*'g CO'[2]~'-C'~~gdw^-1~~d^-1)) 
plot(line_plotCmin)

library(cowplot)

jpeg(filename="Figure2.jpeg", bg="transparent", res=500, units = "in", height=5, width=13) 

f2 <- plot_grid(line_plotCmin, box_plotResp2,  ncol = 2, labels=c('A', 'B'),align="hv", label_size=25)

f2

dev.off()

#Analyze and visualize responses of soil properties and qPCR data

din <- read.csv("C:/Users/ernie/OneDrive/Desktop/Coweeta_drought/drought_chem.csv")

din$Day <- as.factor(din$Day)

#Relationships of AOB with NO3 on day 84 (fig. S4)

d84 <- din[din$Day=='84',]

d84$Treatment <- factor(d84$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
d84$Treatment

summary(lm(log10(NO3+.01)~AOB, data=d84))

jpeg(filename="figs5.jpeg", bg="transparent", res=500, units = "in", height=4, width=6) 

AOBreg <- ggplot(d84, aes(x=d84$AOB, y=log10(d84$NO3+.01))) + 
  geom_smooth(method='lm', colour="black") +
  geom_point(size=3, aes(color=d84$Treatment)) +   
  scale_color_brewer(palette="Set1") +
  scale_x_continuous(bquote(~'Log'[10]~'(AOB'~italic(amoA)~~ 'gene copies'~~gdw^-1*')')) +
  theme_classic() +
  theme(text=element_text(size=16)) +
  theme(axis.text=element_text(size=15)) +
  theme(legend.title=element_blank()) +
  theme(legend.position=c(.8,.22)) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) +
  ylab(bquote(~'Log'[10]~(mu*'g'~'NO'[3]^{"-"}~'-N'~~gdw^-1))) +
  annotate("text", label="R^2 == 0.771~~~~~~~P < 0.001", x=2.9, y=1.2, size=5.5, parse=TRUE) 
AOBreg

dev.off()

#Multiple regression and 3D plot of AOA, C:N, and NO3 

NO3reg <- lm(log10(NO3+.01) ~ DOC.TDN + AOA, data=din )
summary(NO3reg)
Anova(NO3reg)
write.csv(as.data.frame(Anova(NO3reg)),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\no3_lm.csv", row.names = TRUE)


din2<-din[,c(8,17,21,12)]

din2 <- na.omit(din2)

din2$Treatment <- factor(din2$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
din2$Treatment

library(scatterplot3d)
library(RColorBrewer)

brewer.pal(n = 4, name = "Set1")

colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
colors <- colors[as.numeric(din2$Treatment)]

din2$NO3 <- log10(din2$NO3+.01)

jpeg(filename="3d.jpeg", bg="transparent", res=1500, units = "in", height=5, width=5.5) 

s3d <- scatterplot3d(din2$NO3 ~ din2$DOC.TDN + din2$AOA, angle=50, color=colors, pch=16,grid=TRUE, box=TRUE,lab.z=3,type="h",scale.y=.75,                
                     xlab = "",
                     ylab = "",
                     zlab = "")

legend("bottom", legend = levels(din2$Treatment),
       col =  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"), pch = 16, bty="n",inset = -0.35, xpd = TRUE, ncol=2)

dev.off()

#Spearman rank correlations (supplemental figure)

library(corrplot)
library(RColorBrewer)
scalebluered <- colorRampPalette(brewer.pal(8, "RdBu"))(8)

Corr <- cor(din[,c(11:26)], use="complete.obs", method="spearman")

res1 <- cor.mtest(din[,c(11:26)], conf.level = .95)

jpeg(filename="correlogram.jpeg", bg="transparent", res=500, units = "in", height=7, width=7) 

corrplot(Corr, method="color", col=scalebluered,number.cex=.75, addCoef.col = "black",  type="upper",cl.cex=1, tl.cex=.75, tl.col="black", tl.srt=45,sig.level=0.05,p.mat = res1$p, insig = "blank")

dev.off()

#fungal:bacterial ratios

glmfb <- glmer(ITS.16S ~ Drought*LandUse*Day+ (1|Plot), data=din, family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(glmfb))
Anova(glmfb)
write.csv(as.data.frame(Anova(glmfb)),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\fb_glm.csv", row.names = TRUE)

em_glmfb <- emmeans(glmfb, ~ LandUse * Drought | Day)
contrast(em_glmfb, method = "pairwise")
write.csv(as.data.frame(contrast(em_glmfb, method = "pairwise")),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\fb_emmeans.csv", row.names = TRUE)


FB1 <- aggregate(ITS.16S~Drought+Day, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
FB1

FB1.1 <- aggregate(ITS.16S~LandUse, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
FB1.1

FB2 <- aggregate(ITS.16S~Treatment+Drought+Day+LandUse, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
FB2

FB2 <- do.call(data.frame, FB2)
FB2

FB2$se <- FB2$ITS.16S.sd / sqrt(FB2$ITS.16S.n)
head(FB2)

FB2$Treatment <- factor(FB2$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
FB2$Treatment

library(grid)

text_1 <- textGrob("Day 1:\nPre-drought", gp=gpar(fontsize=15))
text_42 <- textGrob("Day 42:\nPost-drought", gp=gpar(fontsize=15))
text_43 <- textGrob("Day 43: 1 Day\nafter re-wet", gp=gpar(fontsize=15))
text_45 <- textGrob("Day 45: 3 Days\nafter re-wet", gp=gpar(fontsize=15))
text_56 <- textGrob("Day 56: 14 Days\nafter re-wet", gp=gpar(fontsize=15))
text_84 <- textGrob("Day 84: 42 Days\nafter re-wet", gp=gpar(fontsize=15))

library(ggplot2)
bar_plotFB<- ggplot(FB2, aes(x=FB2$Day, y=FB2$ITS.16S.mean, fill=FB2$Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=FB2$ITS.16S.mean-FB2$se, ymax=FB2$ITS.16S.mean+FB2$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate('text', 3.5, .18, label="Drought: ~italic(P) == 0.067", size=5, parse=TRUE)+
  annotate('text', 3.5, .17, label="Disturbance: italic(P) == 0.021", size = 5, parse=TRUE)+
  annotate('text', 3.5, .16, label="Disturbance %*% Drought: italic(P) == 0.703", size=5, parse=TRUE)+
  #1
  annotate("text", label="a", x=.67, y=.14, size=5) +
  annotate("text", label="a", x=.89, y=.14, size=5) +
  annotate("text", label="a", x=1.12, y=.085, size=5) +
  annotate("text", label="a", x=1.34, y=.095, size=5) +
  #42
  annotate("text", label="ab", x=1.67, y=.095, size=5) +
  annotate("text", label="a", x=1.89, y=.15, size=5) +
  annotate("text", label="b", x=2.12, y=.065, size=5) +
  annotate("text", label="ab", x=2.34, y=.09, size=5) +
  #43
  annotate("text", label="a", x=2.67, y=.145, size=5) +
  annotate("text", label="ab", x=2.89, y=.13, size=5) +
  annotate("text", label="b", x=3.12, y=.085, size=5) +
  annotate("text", label="b", x=3.34, y=.08, size=5) +
  #45
  annotate("text", label="a", x=3.67, y=.142, size=5) +
  annotate("text", label="ab", x=3.89, y=.1, size=5) +
  annotate("text", label="a", x=4.12, y=.126, size=5) +
  annotate("text", label="b", x=4.34, y=.055, size=5) +
  #56
  annotate("text", label="a", x=4.67, y=.12, size=5) +
  annotate("text", label="a", x=4.89, y=.1, size=5) +
  annotate("text", label="a", x=5.12, y=.07, size=5) +
  annotate("text", label="a", x=5.34, y=.08, size=5) +
  #84
  annotate("text", label="a", x=5.67, y=.185, size=5) +
  annotate("text", label="ab", x=5.89, y=.13, size=5) +
  annotate("text", label="ab", x=6.12, y=.14, size=5) +
  annotate("text", label="b", x=6.34, y=.075, size=5) +
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position="none") +
  theme(legend.text = element_text(color = "black", size = 15)) +
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  ylab('ITS:16S') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  annotation_custom(text_1,xmin=1,xmax=1,ymin=-.02,ymax=-.03) + 
  annotation_custom(text_42,xmin=2,xmax=2,ymin=-.02,ymax=-.03) +
  annotation_custom(text_43,xmin=3,xmax=3,ymin=-.02,ymax=-.03) + 
  annotation_custom(text_45,xmin=4,xmax=4,ymin=-.02,ymax=-.03) +
  annotation_custom(text_56,xmin=5,xmax=5,ymin=-.02,ymax=-.03) + 
  annotation_custom(text_84,xmin=6,xmax=6,ymin=-.02,ymax=-.03) +
  coord_cartesian(clip="off") 
plot(bar_plotFB)

#Extractable DOC

glm8 = glmer(DOC~ LandUse*Drought*Day+(1|Plot), data=din, family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(glm8))
Anova(glm8)
write.csv(as.data.frame(Anova(glm8)),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\doc_glm.csv", row.names = TRUE)

em_glm8 <- emmeans(glm8, ~ LandUse * Drought | Day)
contrast(em_glm8, method = "pairwise")
write.csv((contrast(em_glm8, method = "pairwise")),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\doc_emmeans.csv", row.names = TRUE)


DOC <- aggregate(DOC~Treatment+Drought+Day+LandUse, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(DOC)

DOC <- do.call(data.frame, DOC)
DOC

DOC$se <- DOC$DOC.sd / sqrt(DOC$DOC.n)
head(DOC)

DOC$Treatment <- factor(DOC$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
DOC$Treatment

bar_plotDOC <- ggplot(DOC, aes(x=DOC$Day, y=DOC$DOC.mean, fill=DOC$Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=DOC$DOC.mean-DOC$se, ymax=DOC$DOC.mean+DOC$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate('text', 3.5, 800, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 3.5, 750, label="Disturbance: italic(P)< 0.001", size = 5, parse=TRUE)+
  annotate('text', 3.5, 700, label="Disturbance %*% Drought: italic(P) == 0.219", size=5, parse=TRUE)+
  #Day 1
  annotate("text", label="a", x=.67, y=725, size=5) +
  annotate("text", label="ab", x=.89, y=635, size=5) +
  annotate("text", label="b", x=1.12, y=460, size=5) +
  annotate("text", label="b", x=1.34, y=450, size=5) +
  # Day 42
  annotate("text", label="a", x=1.67, y=400, size=5) +
  annotate("text", label="b", x=1.89, y=825, size=5) +
  annotate("text", label="c", x=2.12, y=290, size=5) +
  annotate("text", label="a", x=2.34, y=460, size=5) +
  #Day 43
  annotate("text", label="a", x=2.67, y=400, size=5) +
  annotate("text", label="b", x=2.89, y=650, size=5) +
  annotate("text", label="c", x=3.12, y=270, size=5) +
  annotate("text", label="a", x=3.34, y=340, size=5) +
  #Day 45
  annotate("text", label="ab", x=3.67, y=380, size=5) +
  annotate("text", label="d", x=3.89, y=530, size=5) +
  annotate("text", label="c", x=4.12, y=250, size=5) +
  annotate("text", label="bc", x=4.34, y=330, size=5) +
  #Day 56
  annotate("text", label="a", x=4.67, y=400, size=5) +
  annotate("text", label="a", x=4.89, y=415, size=5) +
  annotate("text", label="b", x=5.12, y=250, size=5) +
  annotate("text", label="b", x=5.34, y=250, size=5) +
  #Day 84
  annotate("text", label="a", x=5.67, y=415, size=5) +
  annotate("text", label="a", x=5.89, y=480, size=5) +
  annotate("text", label="b", x=6.12, y=280, size=5) +
  annotate("text", label="b", x=6.34, y=290, size=5) +
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.85,.9)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  ylab(bquote(~'Extractable DOC'~~'('~mu*'g'~~'C'~~gdw^-1~')')) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  annotation_custom(text_1,xmin=1,xmax=1,ymin=-200,ymax=-2) + 
  annotation_custom(text_42,xmin=2,xmax=2,ymin=-200,ymax=-2) +
  annotation_custom(text_43,xmin=3,xmax=3,ymin=-200,ymax=-2) + 
  annotation_custom(text_45,xmin=4,xmax=4,ymin=-200,ymax=-2) +
  annotation_custom(text_56,xmin=5,xmax=5,ymin=-200,ymax=-2) + 
  annotation_custom(text_84,xmin=6,xmax=6,ymin=-200,ymax=-2) +
  coord_cartesian(clip="off") 
plot(bar_plotDOC)

#Microbial Biomass C

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

library(ggplot2)

bar_plotMBC <- ggplot(MBC, aes(x=MBC$Day, y=MBC$MBC.mean, fill=MBC$Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=MBC$MBC.mean-MBC$se, ymax=MBC$MBC.mean+MBC$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate('text', 3.5, 250, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 3.5, 230, label="Disturbance: italic(P) == 0.007", size = 5, parse=TRUE)+
  annotate('text', 3.5, 210, label="Disturbance %*% Drought: italic(P) == 0.008", size=5, parse=TRUE)+
  #Day 1
  annotate("text", label="a", x=.67, y=175, size=5) +
  annotate("text", label="a", x=.89, y=180, size=5) +
  annotate("text", label="a", x=1.12, y=130, size=5) +
  annotate("text", label="a", x=1.34, y=150, size=5) +
  #Day 42
  annotate("text", label="a", x=1.67, y=260, size=5) +
  annotate("text", label="b", x=1.89, y=145, size=5) +
  annotate("text", label="b", x=2.12, y=130, size=5) +
  annotate("text", label="b", x=2.34, y=90, size=5) +
  #Day 43
  annotate("text", label="a", x=2.67, y=195, size=5) +
  annotate("text", label="ab", x=2.89, y=150, size=5) +
  annotate("text", label="b", x=3.12, y=105, size=5) +
  annotate("text", label="ab", x=3.34, y=130, size=5) +
  #Day 45
  annotate("text", label="a", x=3.67, y=165, size=5) +
  annotate("text", label="a", x=3.89, y=125, size=5) +
  annotate("text", label="a", x=4.12, y=120, size=5) +
  annotate("text", label="a", x=4.34, y=90, size=5) +
  #Day 56
  annotate("text", label="a", x=4.67, y=185, size=5) +
  annotate("text", label="b", x=4.89, y=110, size=5) +
  annotate("text", label="b", x=5.12, y=105, size=5) +
  annotate("text", label="b", x=5.34, y=95, size=5) +
  #Day 84
  annotate("text", label="a", x=5.67, y=180, size=5) +
  annotate("text", label="a", x=5.89, y=135, size=5) +
  annotate("text", label="a", x=6.12, y=125, size=5) +
  annotate("text", label="a", x=6.34, y=105, size=5) +
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  theme(legend.position=c(.85,.9)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  ylab(bquote(~'Microbial Biomass C'~~'('~mu*'g'~~'C'~~gdw^-1~')')) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  annotation_custom(text_1,xmin=1,xmax=1,ymin=-75,ymax=-2) + 
  annotation_custom(text_42,xmin=2,xmax=2,ymin=-75,ymax=-2) +
  annotation_custom(text_43,xmin=3,xmax=3,ymin=-75,ymax=-2) + 
  annotation_custom(text_45,xmin=4,xmax=4,ymin=-75,ymax=-2) +
  annotation_custom(text_56,xmin=5,xmax=5,ymin=-75,ymax=-2) + 
  annotation_custom(text_84,xmin=6,xmax=6,ymin=-75,ymax=-2) +
  coord_cartesian(clip="off") 
plot(bar_plotMBC)

#Microbial Biomass N

glm13 = lmer(MBN~ LandUse*Drought*Day+(1|Plot), data=din)
qqPlot(resid(glm13))
Anova(glm13)
em_glm13 <- emmeans(glm13, ~ LandUse * Drought | Day)
contrast(em_glm13, method = "pairwise")

MBN <- aggregate(MBN~Treatment+Drought+Day+LandUse, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(MBN)

MBN <- do.call(data.frame, MBN)
MBN

MBN$se <- MBN$MBN.sd / sqrt(MBN$MBN.n)
head(MBN)

MBN$Treatment <- factor(MBN$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
MBN$Treatment

library(ggplot2)

bar_plotMBN <- ggplot(MBN, aes(x=MBN$Day, y=MBN$MBN.mean, fill=MBN$Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=MBN$MBN.mean-MBN$se, ymax=MBN$MBN.mean+MBN$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate('text', 4.5, 50, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 4.5, 47, label="Disturbance: italic(P) == 0.012", size = 5, parse=TRUE)+
  annotate('text', 4.5, 44, label="Disturbance %*% Drought: italic(P) == 0.004", size=5, parse=TRUE)+
  
  annotate("text", label="a", x=.67, y=43, size=5) +
  annotate("text", label="a", x=.89, y=38, size=5) +
  annotate("text", label="a", x=1.12, y=25, size=5) +
  annotate("text", label="a", x=1.34, y=26, size=5) +
  
  annotate("text", label="a", x=1.67, y=39, size=5) +
  annotate("text", label="b", x=1.89, y=12, size=5) +
  annotate("text", label="ab", x=2.12, y=25, size=5) +
  annotate("text", label="b", x=2.34, y=18, size=5) +
  
  annotate("text", label="a", x=2.67, y=43, size=5) +
  annotate("text", label="b", x=2.89, y=19, size=5) +
  annotate("text", label="ab", x=3.12, y=32, size=5) +
  annotate("text", label="b", x=3.34, y=23, size=5) +
  
  annotate("text", label="a", x=3.67, y=31, size=5) +
  annotate("text", label="a", x=3.89, y=27, size=5) +
  annotate("text", label="a", x=4.12, y=23, size=5) +
  annotate("text", label="a", x=4.34, y=20, size=5) +
  
  annotate("text", label="a", x=4.67, y=38, size=5) +
  annotate("text", label="b", x=4.89, y=19, size=5) +
  annotate("text", label="ab", x=5.12, y=26, size=5) +
  annotate("text", label="b", x=5.34, y=22, size=5) +
  
  annotate("text", label="a", x=5.67, y=45, size=5) +
  annotate("text", label="b", x=5.89, y=29, size=5) +
  annotate("text", label="ab", x=6.12, y=33, size=5) +
  annotate("text", label="ab", x=6.34, y=28, size=5) +
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.16,.9)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  ylab(bquote(~'Microbial Biomass N'~~'('~mu*'g'~~'N'~~gdw^-1~')')) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  annotation_custom(text_1,xmin=1,xmax=1,ymin=-9,ymax=-2) + 
  annotation_custom(text_42,xmin=2,xmax=2,ymin=-9,ymax=-2) +
  annotation_custom(text_43,xmin=3,xmax=3,ymin=-9,ymax=-2) + 
  annotation_custom(text_45,xmin=4,xmax=4,ymin=-9,ymax=-2) +
  annotation_custom(text_56,xmin=5,xmax=5,ymin=-9,ymax=-2) + 
  annotation_custom(text_84,xmin=6,xmax=6,ymin=-9,ymax=-2) +
  coord_cartesian(clip="off") 
plot(bar_plotMBN)

#Microbial Biomass C:N

glm14 = lmer(MBCN~ LandUse*Drought*Day+(1|Plot), data=din)
qqPlot(resid(glm14))
Anova(glm14)
write.csv(as.data.frame(Anova(glm14)),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\MBCN_glm.csv", row.names = TRUE)

em_glm14 <- emmeans(glm14, ~ LandUse * Drought | Day)
contrast(em_glm14, method = "pairwise")
write.csv(as.data.frame(contrast(em_glm14, method = "pairwise")),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\MBCN_emmeans.csv", row.names = TRUE)

MBCN <- aggregate(MBCN~Drought+Day, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(MBCN)

MBCN <- aggregate(MBCN~Treatment+Drought+Day+LandUse, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(MBCN)

MBCN <- do.call(data.frame, MBCN)
MBCN

MBCN$se <- MBCN$MBCN.sd / sqrt(MBCN$MBCN.n)
head(MBCN)

MBCN$Treatment <- factor(MBCN$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
MBCN$Treatment

library(ggplot2)

bar_plotMBCN <- ggplot(MBCN, aes(x=MBCN$Day, y=MBCN$MBCN.mean, fill=MBCN$Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=MBCN$MBCN.mean-MBCN$se, ymax=MBCN$MBCN.mean+MBCN$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate('text', 3.5, 15, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 3.5, 13.75, label="Disturbance: italic(P) == 0.32", size = 5, parse=TRUE)+
  annotate('text', 3.5, 12.5, label="Disturbance %*% Drought: italic(P) == 0.842", size=5, parse=TRUE)+
  #Day 1
  annotate("text", label="a", x=.67, y=5, size=5) +
  annotate("text", label="a", x=.89, y=6, size=5) +
  annotate("text", label="a", x=1.12, y=5.5, size=5) +
  annotate("text", label="a", x=1.34, y=8, size=5) +
  #Day 42
  annotate("text", label="ac", x=1.67, y=8, size=5) +
  annotate("text", label="b", x=1.89, y=15, size=5) +
  annotate("text", label="a", x=2.12, y=6.5, size=5) +
  annotate("text", label="bc", x=2.34, y=16, size=5) +
  #Day 43
  annotate("text", label="a", x=2.67, y=6, size=5) +
  annotate("text", label="b", x=2.89, y=10.5, size=5) +
  annotate("text", label="a", x=3.12, y=5, size=5) +
  annotate("text", label="ab", x=3.34, y=7.5, size=5) +
  #Day 45
  annotate("text", label="a", x=3.67, y=6.5, size=5) +
  annotate("text", label="a", x=3.89, y=7.25, size=5) +
  annotate("text", label="a", x=4.12, y=6.5, size=5) +
  annotate("text", label="a", x=4.34, y=8, size=5) +
  #Day 56
  annotate("text", label="a", x=4.67, y=7.25, size=5) +
  annotate("text", label="a", x=4.89, y=7.25, size=5) +
  annotate("text", label="a", x=5.12, y=5.75, size=5) +
  annotate("text", label="a", x=5.34, y=6.75, size=5) +
  #Day 84
  annotate("text", label="a", x=5.67, y=5, size=5) +
  annotate("text", label="a", x=5.89, y=5.5, size=5) +
  annotate("text", label="a", x=6.12, y=5, size=5) +
  annotate("text", label="a", x=6.34, y=5, size=5) +
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position="none") +
  theme(legend.text = element_text(color = "black", size = 15)) +
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  ylab('Microbial Biomass C:N') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  annotation_custom(text_1,xmin=1,xmax=1,ymin=-2.5,ymax=-2) + 
  annotation_custom(text_42,xmin=2,xmax=2,ymin=-2.5,ymax=-2) +
  annotation_custom(text_43,xmin=3,xmax=3,ymin=-2.5,ymax=-2) + 
  annotation_custom(text_45,xmin=4,xmax=4,ymin=-2.5,ymax=-2) +
  annotation_custom(text_56,xmin=5,xmax=5,ymin=-2.5,ymax=-2) + 
  annotation_custom(text_84,xmin=6,xmax=6,ymin=-2.5,ymax=-2) +
  coord_cartesian(clip="off") 
plot(bar_plotMBCN)

#Total Extractable N

glm9 = glmer(TDN~ LandUse*Drought*Day+(1|Plot), data=din, family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(glm9))
Anova(glm9)
write.csv(as.data.frame(Anova(glm9)),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\tdn_glm.csv", row.names = TRUE)

em_glm9 <- emmeans(glm9, ~ LandUse * Drought | Day)
contrast(em_glm9, method = "pairwise")
write.csv(as.data.frame(contrast(em_glm9, method = "pairwise")),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\tdn_emmeans.csv", row.names = TRUE)

TDN <- aggregate(TDN~Treatment+Drought+Day+LandUse, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(TDN)

TDN.1 <- aggregate(TDN~Drought+Day, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
TDN.1

TDN.2 <- aggregate(TDN~Treatment+Day, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
TDN.2

TDN <- do.call(data.frame, TDN)
TDN

TDN$se <- TDN$TDN.sd / sqrt(TDN$TDN.n)
head(TDN)

TDN$Treatment <- factor(TDN$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
TDN$Treatment

library(ggplot2)

bar_plotTDN <- ggplot(TDN, aes(x=TDN$Day, y=TDN$TDN.mean, fill=TDN$Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=TDN$TDN.mean-TDN$se, ymax=TDN$TDN.mean+TDN$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate('text', 1.5, 150, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 1.5, 140, label="Disturbance: italic(P) == 0.054", size = 5, parse=TRUE)+
  annotate('text', 1.5, 130, label="Disturbance %*% Drought: italic(P) < 0.001", size=5, parse=TRUE)+
  #1
  annotate("text", label="a", x=.67, y=50, size=5) +
  annotate("text", label="a", x=.89, y=50, size=5) +
  annotate("text", label="a", x=1.12, y=45, size=5) +
  annotate("text", label="a", x=1.34, y=45, size=5) +
  #42
  annotate("text", label="a", x=1.67, y=50, size=5) +
  annotate("text", label="b", x=1.89, y=105, size=5) +
  annotate("text", label="a", x=2.12, y=55, size=5) +
  annotate("text", label="a", x=2.34, y=65, size=5) +
  #43
  annotate("text", label="a", x=2.67, y=57, size=5) +
  annotate("text", label="c", x=2.89, y=145, size=5) +
  annotate("text", label="ab", x=3.12, y=95, size=5) +
  annotate("text", label="bc", x=3.34, y=100, size=5) +
  #45
  annotate("text", label="a", x=3.67, y=62, size=5) +
  annotate("text", label="b", x=3.89, y=135, size=5) +
  annotate("text", label="a", x=4.12, y=57, size=5) +
  annotate("text", label="b", x=4.34, y=105, size=5) +
  #56
  annotate("text", label="ac", x=4.67, y=60, size=5) +
  annotate("text", label="b", x=4.89, y=146, size=5) +
  annotate("text", label="a", x=5.12, y=52, size=5) +
  annotate("text", label="c", x=5.34, y=80, size=5) +
  #84
  annotate("text", label="ab", x=5.67, y=78, size=5) +
  annotate("text", label="c", x=5.89, y=136, size=5) +
  annotate("text", label="a", x=6.12, y=70, size=5) +
  annotate("text", label="bc", x=6.34, y=112, size=5) +
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=20)) +
  theme(legend.position="none") +
  theme(legend.text = element_text(color = "black", size = 15)) +
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  ylab(bquote(~'Total Extractable N'~~'('~mu*'g'~~'N'~~gdw^-1~')')) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  annotation_custom(text_1,xmin=1,xmax=1,ymin=-40,ymax=-2) + 
  annotation_custom(text_42,xmin=2,xmax=2,ymin=-40,ymax=-2) +
  annotation_custom(text_43,xmin=3,xmax=3,ymin=-40,ymax=-2) + 
  annotation_custom(text_45,xmin=4,xmax=4,ymin=-40,ymax=-2) +
  annotation_custom(text_56,xmin=5,xmax=5,ymin=-40,ymax=-2) + 
  annotation_custom(text_84,xmin=6,xmax=6,ymin=-40,ymax=-2) +
  coord_cartesian(clip="off") 
plot(bar_plotTDN)

#Extractable dissolved organic N

glm10 = glmer(DON~ LandUse*Drought*Day+(1|Plot), data=din, family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(glm10))
Anova(glm10)
em_glm10 <- emmeans(glm10, ~ LandUse * Drought | Day)
contrast(em_glm10, method = "pairwise")

DON <- aggregate(DON~Treatment+Drought+Day+LandUse, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(DON)

DON <- do.call(data.frame, DON)
DON

DON$se <- DON$DON.sd / sqrt(DON$DON.n)
head(DON)

DON$Treatment <- factor(DON$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
DON$Treatment

library(ggplot2)

bar_plotDON <- ggplot(DON, aes(x=DON$Day, y=DON$DON.mean, fill=DON$Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=DON$DON.mean-DON$se, ymax=DON$DON.mean+DON$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate('text', 1.5, 150, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 1.5, 140, label="Disturbance: italic(P) == 0.02", size = 5, parse=TRUE)+
  annotate('text', 1.5, 130, label="Disturbance %*% Drought: italic(P) == 0.002", size=5, parse=TRUE)+
  #1
  annotate("text", label="a", x=.67, y=48, size=5) +
  annotate("text", label="a", x=.89, y=48, size=5) +
  annotate("text", label="a", x=1.12, y=43, size=5) +
  annotate("text", label="a", x=1.34, y=43, size=5) +
  #42
  annotate("text", label="a", x=1.67, y=48, size=5) +
  annotate("text", label="b", x=1.89, y=105, size=5) +
  annotate("text", label="a", x=2.12, y=52, size=5) +
  annotate("text", label="a", x=2.34, y=65, size=5) +
  #43
  annotate("text", label="a", x=2.67, y=53, size=5) +
  annotate("text", label="c", x=2.89, y=134, size=5) +
  annotate("text", label="ab", x=3.12, y=87, size=5) +
  annotate("text", label="bc", x=3.34, y=95, size=5) +
  #45
  annotate("text", label="a", x=3.67, y=60, size=5) +
  annotate("text", label="b", x=3.89, y=128, size=5) +
  annotate("text", label="a", x=4.12, y=53, size=5) +
  annotate("text", label="b", x=4.34, y=98, size=5) +
  #56
  annotate("text", label="ac", x=4.67, y=57, size=5) +
  annotate("text", label="b", x=4.89, y=133, size=5) +
  annotate("text", label="a", x=5.12, y=48, size=5) +
  annotate("text", label="c", x=5.34, y=73, size=5) +
  #84
  annotate("text", label="ab", x=5.67, y=72, size=5) +
  annotate("text", label="c", x=5.89, y=127, size=5) +
  annotate("text", label="a", x=6.12, y=63, size=5) +
  annotate("text", label="bc", x=6.34, y=104, size=5) +
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.11,.65)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  ylab(bquote(~'Extractable Organic N'~~'('~mu*'g'~~'N'~~gdw^-1~')')) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  annotation_custom(text_1,xmin=1,xmax=1,ymin=-30,ymax=-2) + 
  annotation_custom(text_42,xmin=2,xmax=2,ymin=-30,ymax=-2) +
  annotation_custom(text_43,xmin=3,xmax=3,ymin=-30,ymax=-2) + 
  annotation_custom(text_45,xmin=4,xmax=4,ymin=-30,ymax=-2) +
  annotation_custom(text_56,xmin=5,xmax=5,ymin=-30,ymax=-2) + 
  annotation_custom(text_84,xmin=6,xmax=6,ymin=-30,ymax=-2) +
  coord_cartesian(clip="off") 
plot(bar_plotDON)

#Extractable C:N

glm11 = glmer(DOC.TDN~ LandUse*Drought*Day+(1|Plot), data=din, family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(glm11))
Anova(glm11)
em_glm11 <- emmeans(glm11, ~ LandUse * Drought | Day)
contrast(em_glm11, method = "pairwise")

DOC.TDN <- aggregate(DOC.TDN~Treatment+Drought+Day+LandUse, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(DOC.TDN)

DOC.TDN <- do.call(data.frame, DOC.TDN)
DOC.TDN

DOC.TDN$se <- DOC.TDN$DOC.TDN.sd / sqrt(DOC.TDN$DOC.TDN.n)
head(DOC.TDN)

DOC.TDN$Treatment <- factor(DOC.TDN$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
DOC.TDN$Treatment

library(ggplot2)

bar_plotDOC.TDN <- ggplot(DOC.TDN, aes(x=DOC.TDN$Day, y=DOC.TDN$DOC.TDN.mean, fill=DOC.TDN$Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=DOC.TDN$DOC.TDN.mean-DOC.TDN$se, ymax=DOC.TDN$DOC.TDN.mean+DOC.TDN$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate('text', 3, 20, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 3, 18, label="Disturbance: italic(P) == 0.02", size = 5, parse=TRUE)+
  annotate('text', 3, 16, label="Disturbance %*% Drought: italic(P) == 0.003", size=5, parse=TRUE)+
  #1
  annotate("text", label="a", x=.67, y=20, size=5) +
  annotate("text", label="a", x=.89, y=19, size=5) +
  annotate("text", label="a", x=1.12, y=16.5, size=5) +
  annotate("text", label="a", x=1.34, y=15, size=5) +
  #42
  annotate("text", label="a", x=1.67, y=10, size=5) +
  annotate("text", label="a", x=1.89, y=9, size=5) +
  annotate("text", label="a", x=2.12, y=6.5, size=5) +
  annotate("text", label="a", x=2.34, y=8.5, size=5) +
  #43
  annotate("text", label="a", x=2.67, y=8.5, size=5) +
  annotate("text", label="b", x=2.89, y=5.5, size=5) +
  annotate("text", label="b", x=3.12, y=5, size=5) +
  annotate("text", label="b", x=3.34, y=4.5, size=5) +
  #45
  annotate("text", label="a", x=3.67, y=10, size=5) +
  annotate("text", label="b", x=3.89, y=6, size=5) +
  annotate("text", label="b", x=4.12, y=6.5, size=5) +
  annotate("text", label="b", x=4.34, y=5, size=5) +
  #56
  annotate("text", label="a", x=4.67, y=8.75, size=5) +
  annotate("text", label="b", x=4.89, y=4.5, size=5) +
  annotate("text", label="a", x=5.12, y=6, size=5) +
  annotate("text", label="b", x=5.34, y=4, size=5) +
  #84
  annotate("text", label="a", x=5.67, y=8, size=5) +
  annotate("text", label="bc", x=5.89, y=5, size=5) +
  annotate("text", label="ac", x=6.12, y=5, size=5) +
  annotate("text", label="b", x=6.34, y=3.5, size=5) +
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.85,.85)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  ylab('Extractable C:N') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  annotation_custom(text_1,xmin=1,xmax=1,ymin=-3.5,ymax=-2) + 
  annotation_custom(text_42,xmin=2,xmax=2,ymin=-3.5,ymax=-2) +
  annotation_custom(text_43,xmin=3,xmax=3,ymin=-3.5,ymax=-2) + 
  annotation_custom(text_45,xmin=4,xmax=4,ymin=-3.5,ymax=-2) +
  annotation_custom(text_56,xmin=5,xmax=5,ymin=-3.5,ymax=-2) + 
  annotation_custom(text_84,xmin=6,xmax=6,ymin=-3.5,ymax=-2) +
  coord_cartesian(clip="off") 
plot(bar_plotDOC.TDN)

#NH4-N

glm2 <- glmer(NH4~ LandUse*Drought*Day+(1|Plot), data=din, family=Gamma(link=log), glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(glm2))
Anova(glm2)
write.csv(as.data.frame(Anova(glm2)),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\nh4_glm.csv", row.names = TRUE)

em_glm2 <- emmeans(glm2, ~  LandUse * Drought | Day)
contrast(em_glm2, method = "pairwise")
write.csv(as.data.frame(contrast(em_glm2, method = "pairwise")),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\nh4_emmeans.csv", row.names = TRUE)

NH4 <- aggregate(NH4~Treatment+Drought+Day+LandUse, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(NH4)

NH4.1 <- aggregate(NH4~Drought+Day, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
NH4.1

NH4.2 <- aggregate(NH4~Treatment+Day, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
NH4.2

NH4 <- do.call(data.frame, NH4)
NH4

NH4$se <- NH4$NH4.sd / sqrt(NH4$NH4.n)
head(NH4)

NH4$Treatment <- factor(NH4$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
NH4$Treatment

bar_plotNH4 <- ggplot(NH4, aes(x=NH4$Day, y=NH4$NH4.mean, fill=NH4$Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=NH4$NH4.mean-NH4$se, ymax=NH4$NH4.mean+NH4$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate('text', 1.5, 16, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 1.5, 15, label="Disturbance: italic(P) == 0.491", size = 5, parse=TRUE)+
  annotate('text', 1.5, 14, label="Disturbance %*% Drought: italic(P) == 0.003", size=5, parse=TRUE)+
  #Day 1
  annotate("text", label="a", x=.67, y=2.25, size=5) +
  annotate("text", label="a", x=.89, y=2.25, size=5) +
  annotate("text", label="a", x=1.12, y=2.25, size=5) +
  annotate("text", label="a", x=1.34, y=2.25, size=5) +
  #Day 42
  annotate("text", label="a", x=1.67, y=4, size=5) +
  annotate("text", label="a", x=1.89, y=4, size=5) +
  annotate("text", label="a", x=2.12, y=4, size=5) +
  annotate("text", label="a", x=2.34, y=4, size=5) +
  #Day 43
  annotate("text", label="a", x=2.67, y=3, size=5) +
  annotate("text", label="c", x=2.89, y=11.5, size=5) +
  annotate("text", label="b", x=3.12, y=7.25, size=5) +
  annotate("text", label="bc", x=3.34, y=8, size=5) +
  #Day 45
  annotate("text", label="a", x=3.67, y=4.5, size=5) +
  annotate("text", label="b", x=3.89, y=14, size=5) +
  annotate("text", label="a", x=4.12, y=4, size=5) +
  annotate("text", label="b", x=4.34, y=13, size=5) +
  #Day 56
  annotate("text", label="a", x=4.67, y=4, size=5) +
  annotate("text", label="b", x=4.89, y=17.25, size=5) +
  annotate("text", label="a", x=5.12, y=3.25, size=5) +
  annotate("text", label="b", x=5.34, y=9, size=5) +
  #Day 84
  annotate("text", label="ab", x=5.67, y=5.25, size=5) +
  annotate("text", label="c", x=5.89, y=13.5, size=5) +
  annotate("text", label="a", x=6.12, y=3.5, size=5) +
  annotate("text", label="bc", x=6.34, y=9, size=5) +
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=18)) +
  theme(text=element_text(size=20)) +
  theme(legend.position="none") +
  theme(legend.text = element_text(color = "black", size = 15)) +
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  ylab(bquote(~mu*'g'~'NH'[4]^{"+"}~~'-N'~~gdw^-1)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  annotation_custom(text_1,xmin=1,xmax=1,ymin=-3,ymax=-2) + 
  annotation_custom(text_42,xmin=2,xmax=2,ymin=-3,ymax=-2) +
  annotation_custom(text_43,xmin=3,xmax=3,ymin=-3,ymax=-2) + 
  annotation_custom(text_45,xmin=4,xmax=4,ymin=-3,ymax=-2) +
  annotation_custom(text_56,xmin=5,xmax=5,ymin=-3,ymax=-2) + 
  annotation_custom(text_84,xmin=6,xmax=6,ymin=-3,ymax=-2) +
coord_cartesian(clip="off") 
plot(bar_plotNH4)

#NO3 -N

glm3 = glmer(NO3+.01~ LandUse*Drought*Day+(1|Plot), data=din, family=Gamma(link=log),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
qqPlot(resid(glm3))
Anova(glm3)
write.csv(as.data.frame(Anova(glm3)),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\no3_glm.csv", row.names = TRUE)

em_glm3 <- emmeans(glm3, ~ LandUse * Drought | Day)
contrast(em_glm3, method = "pairwise")
write.csv(as.data.frame(contrast(em_glm3, method = "pairwise")),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\no3_emmeans.csv", row.names = TRUE)

NO3 <- aggregate(NO3~Treatment+Drought+Day+LandUse, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(NO3)

NO3.1 <- aggregate(NO3~LandUse, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(NO3.1)

NO3.2 <- aggregate(NO3~Drought+Day, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(NO3.2)


NO3 <- do.call(data.frame, NO3)
NO3

NO3$se <- NO3$NO3.sd / sqrt(NO3$NO3.n)
head(NO3)

NO3$Treatment <- factor(NO3$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
NO3$Treatment

bar_plotNO3 <- ggplot(NO3, aes(x=NO3$Day, y=NO3$NO3.mean, fill=NO3$Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=NO3$NO3.mean-NO3$se, ymax=NO3$NO3.mean+NO3$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate('text', 1.5, 5, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 1.5, 4.6, label="Disturbance: italic(P) < 0.001", size = 5, parse=TRUE)+
  annotate('text', 1.5, 4.2, label="Disturbance %*% Drought: italic(P) == 0.237", size=5, parse=TRUE)+
  #Day 1
  annotate("text", label="ab", x=.67, y=.5, size=5) +
  annotate("text", label="a", x=.89, y=.5, size=5) +
  annotate("text", label="b", x=1.12, y=1, size=5) +
  annotate("text", label="b", x=1.34, y=1, size=5) +
  #Day 42
  annotate("text", label="a", x=1.67, y=.5, size=5) +
  annotate("text", label="b", x=1.89, y=.5, size=5) +
  annotate("text", label="c", x=2.12, y=2.75, size=5) +
  annotate("text", label="a", x=2.34, y=.75, size=5) +
  #Day 43
  annotate("text", label="a", x=2.67, y=.5, size=5) +
  annotate("text", label="a", x=2.89, y=.5, size=5) +
  annotate("text", label="b", x=3.12, y=5.6, size=5) +
  annotate("text", label="b", x=3.34, y=1.75, size=5) +
  #Day 45
  annotate("text", label="a", x=3.67, y=.5, size=5) +
  annotate("text", label="a", x=3.89, y=.5, size=5) +
  annotate("text", label="b", x=4.12, y=3.25, size=5) +
  annotate("text", label="b", x=4.34, y=2, size=5) +
  #Day 56
  annotate("text", label="a", x=4.67, y=.5, size=5) +
  annotate("text", label="a", x=4.89, y=.5, size=5) +
  annotate("text", label="b", x=5.12, y=4.1, size=5) +
  annotate("text", label="b", x=5.34, y=3, size=5) +
  #Day 84
  annotate("text", label="b", x=5.67, y=2.5, size=5) +
  annotate("text", label="a", x=5.89, y=.5, size=5) +
  annotate("text", label="b", x=6.12, y=4.75, size=5) +
  annotate("text", label="b", x=6.34, y=7.75, size=5) +
  labs(list(x ="")) +
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.14,.9)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  ylab(bquote(~mu*'g'~'NO'[3]^{"-"}~~'-N'~~gdw^-1)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  annotation_custom(text_1,xmin=1,xmax=1,ymin=-1.25,ymax=-.7) + 
  annotation_custom(text_42,xmin=2,xmax=2,ymin=-1.25,ymax=-.7) +
  annotation_custom(text_43,xmin=3,xmax=3,ymin=-1.25,ymax=-.7) + 
  annotation_custom(text_45,xmin=4,xmax=4,ymin=-1.25,ymax=-.7) +
  annotation_custom(text_56,xmin=5,xmax=5,ymin=-1.25,ymax=-.7) + 
  annotation_custom(text_84,xmin=6,xmax=6,ymin=-1.25,ymax=-.7) +
  coord_cartesian(clip="off") 
plot(bar_plotNO3)

#AOA amoA

glm5 = lmer(AOA~ LandUse*Drought*Day+(1|Plot), data=din)
qqPlot(resid(glm5))
Anova(glm5)
write.csv(as.data.frame(Anova(glm5)),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\aoa_glm.csv", row.names = TRUE)

em_glm5 <- emmeans(glm5, ~ LandUse * Drought | Day)
contrast(em_glm5, method = "pairwise")
write.csv(as.data.frame(contrast(em_glm5, method = "pairwise")),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\aoa_emmeans.csv", row.names = TRUE)

AOA <- aggregate(AOA~Treatment+Drought+Day+LandUse, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(AOA)

AOA <- do.call(data.frame, AOA)
AOA

AOA$se <- AOA$AOA.sd / sqrt(AOA$AOA.n)
head(AOA)

AOA$Treatment <- factor(AOA$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
AOA$Treatment

library(ggplot2)
bar_plotAOA <- ggplot(AOA, aes(x=AOA$Day, y=AOA$AOA.mean, fill=AOA$Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=AOA$AOA.mean-AOA$se, ymax=AOA$AOA.mean+AOA$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate('text', 3.5, 10, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 3.5, 9.5, label="Disturbance: italic(P) < 0.001", size = 5, parse=TRUE)+
  annotate('text', 3.5, 9, label="Disturbance %*% Drought: italic(P) == 0.57", size=5, parse=TRUE)+
  #Day 1
  annotate("text", label="a", x=.67, y=4.5, size=5) +
  annotate("text", label="a", x=.89, y=4, size=5) +
  annotate("text", label="b", x=1.12, y=6.8, size=5) +
  annotate("text", label="b", x=1.34, y=7.6, size=5) +
  #Day 42
  annotate("text", label="a", x=1.67, y=5.1, size=5) +
  annotate("text", label="a", x=1.89, y=4.8, size=5) +
  annotate("text", label="a", x=2.12, y=7.4, size=5) +
  annotate("text", label="a", x=2.34, y=7, size=5) +
  #Day 43
  annotate("text", label="ac", x=2.67, y=5.6, size=5) +
  annotate("text", label="a", x=2.89, y=4.8, size=5) +
  annotate("text", label="b", x=3.12, y=7.6, size=5) +
  annotate("text", label="bc", x=3.34, y=7, size=5) +
  #Day 45
  annotate("text", label="a", x=3.67, y=4.9, size=5) +
  annotate("text", label="a", x=3.89, y=4.9, size=5) +
  annotate("text", label="b", x=4.12, y=8.1, size=5) +
  annotate("text", label="b", x=4.34, y=7.6, size=5) +
  #Day 56
  annotate("text", label="a", x=4.67, y=4.8, size=5) +
  annotate("text", label="a", x=4.89, y=4.8, size=5) +
  annotate("text", label="b", x=5.12, y=8.4, size=5) +
  annotate("text", label="a", x=5.34, y=7, size=5) +
  #Day 84
  annotate("text", label="ac", x=5.67, y=6.6, size=5) +
  annotate("text", label="b", x=5.89, y=4.4, size=5) +
  annotate("text", label="a", x=6.12, y=7.4, size=5) +
  annotate("text", label="bc", x=6.34, y=6.5, size=5) +
  
  labs(list(x ="")) +
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position="none") +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(bquote(~'Log'[10]~'('~AOA~italic(amoA)~~gdw^-1*')'))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  annotation_custom(text_1,xmin=1,xmax=1,ymin=-1.5,ymax=-1) + 
  annotation_custom(text_42,xmin=2,xmax=2,ymin=-1.5,ymax=-1) +
  annotation_custom(text_43,xmin=3,xmax=3,ymin=-1.5,ymax=-1) + 
  annotation_custom(text_45,xmin=4,xmax=4,ymin=-1.5,ymax=-1) +
  annotation_custom(text_56,xmin=5,xmax=5,ymin=-1.5,ymax=-1) + 
  annotation_custom(text_84,xmin=6,xmax=6,ymin=-1.5,ymax=-1) +
  coord_cartesian(clip="off") 
plot(bar_plotAOA)

#AOB amoA

glm6 = lmer(AOB~ LandUse*Drought*Day+(1|Plot), data=din)
qqPlot(resid(glm6))
Anova(glm6)
write.csv(as.data.frame(Anova(glm6)),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\aob_glm.csv", row.names = TRUE)

em_glm6 <- emmeans(glm6, ~ LandUse * Drought | Day)
contrast(em_glm6, method = "pairwise")
write.csv(as.data.frame(contrast(em_glm6, method = "pairwise")),"C:\\Users\\ernie\\OneDrive\\Desktop\\Coweeta_drought\\aob_emmeans.csv", row.names = TRUE)

AOB <- aggregate(AOB~Treatment+Drought+Day+LandUse, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(AOB)

AOB <- do.call(data.frame, AOB)
AOB

AOB$se <- AOB$AOB.sd / sqrt(AOB$AOB.n)
head(AOB)

AOB$Treatment <- factor(AOB$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
AOB$Treatment

bar_plotAOB <- ggplot(AOB, aes(x=AOB$Day, y=AOB$AOB.mean, fill=AOB$Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=AOB$AOB.mean-AOB$se, ymax=AOB$AOB.mean+AOB$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate('text', 3.5, 7.5, label="Drought: ~italic(P) == 0.765", size=5, parse=TRUE)+
  annotate('text', 3.5, 7, label="Disturbance: italic(P) == 0.013", size = 5, parse=TRUE)+
  annotate('text', 3.5, 6.5, label="Disturbance %*% Drought: italic(P) == 0.99", size=5, parse=TRUE)+
  
  annotate("text", label="a", x=.67, y=5, size=5) +
  annotate("text", label="a", x=.89, y=4.75, size=5) +
  annotate("text", label="a", x=1.12, y=5.4, size=5) +
  annotate("text", label="a", x=1.34, y=5.15, size=5) +
  
  annotate("text", label="a", x=1.67, y=4, size=5) +
  annotate("text", label="a", x=1.89, y=4.5, size=5) +
  annotate("text", label="a", x=2.12, y=6.2, size=5) +
  annotate("text", label="a", x=2.34, y=5.75, size=5) +
  
  annotate("text", label="a", x=2.67, y=5.3, size=5) +
  annotate("text", label="a", x=2.89, y=4, size=5) +
  annotate("text", label="a", x=3.12, y=5.8, size=5) +
  annotate("text", label="a", x=3.34, y=5.3, size=5) +
  
  annotate("text", label="a", x=3.67, y=4.25, size=5) +
  annotate("text", label="a", x=3.89, y=4.5, size=5) +
  annotate("text", label="a", x=4.12, y=5.75, size=5) +
  annotate("text", label="a", x=4.34, y=5.5, size=5) +
  
  annotate("text", label="a", x=4.67, y=3.5, size=5) +
  annotate("text", label="ab", x=4.89, y=5.2, size=5) +
  annotate("text", label="ab", x=5.12, y=5.6, size=5) +
  annotate("text", label="b", x=5.34, y=6.2, size=5) +
  
  annotate("text", label="ab", x=5.67, y=4.75, size=5) +
  annotate("text", label="a", x=5.89, y=3.75, size=5) +
  annotate("text", label="b", x=6.12, y=5.75, size=5) +
  annotate("text", label="b", x=6.34, y=6.5, size=5) +
  
  labs(list(x ="")) +
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position="none") +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(bquote(~'Log'[10]~'('~AOB~italic(amoA)~~gdw^-1*')'))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  annotation_custom(text_1,xmin=1,xmax=1,ymin=-1.5,ymax=-.5) + 
  annotation_custom(text_42,xmin=2,xmax=2,ymin=-1.5,ymax=-.5) +
  annotation_custom(text_43,xmin=3,xmax=3,ymin=-1.5,ymax=-.5) + 
  annotation_custom(text_45,xmin=4,xmax=4,ymin=-1.5,ymax=-.5) +
  annotation_custom(text_56,xmin=5,xmax=5,ymin=-1.5,ymax=-.5) + 
  annotation_custom(text_84,xmin=6,xmax=6,ymin=-1.5,ymax=-.5) +
  coord_cartesian(clip="off") 
plot(bar_plotAOB)

#CAOB amoA

glm7<- lmer(CAOB~Drought*LandUse*Day+(1|Plot), data=din)
qqPlot(resid(glm7))
Anova(glm7)
em_glm7 <- emmeans(glm7, ~ LandUse * Drought | Day)
contrast(em_glm7, method = "pairwise")

CAOB <- aggregate(CAOB~Treatment+Drought+Day+LandUse, data=din, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))
head(CAOB)

CAOB <- do.call(data.frame, CAOB)
CAOB

CAOB$se <- CAOB$CAOB.sd / sqrt(CAOB$CAOB.n)
head(CAOB)

CAOB$Treatment <- factor(CAOB$Treatment, levels = c("Reference - Control", "Reference - Drought", "Disturbed - Control", "Disturbed - Drought"))
CAOB$Treatment

bar_plotCAOB <- ggplot(CAOB, aes(x=CAOB$Day, y=CAOB$CAOB.mean, fill=CAOB$Treatment)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin=CAOB$CAOB.mean-CAOB$se, ymax=CAOB$CAOB.mean+CAOB$se), width=0.2, size=1, position=position_dodge(0.9)) +
  annotate('text', 3.5, 8.5, label="Drought: ~italic(P) < 0.001", size=5, parse=TRUE)+
  annotate('text', 3.5, 8, label="Disturbance: italic(P) == 0.18", size = 5, parse=TRUE)+
  annotate('text', 3.5, 7.5, label="Disturbance %*% Drought: italic(P) == 0.36", size=5, parse=TRUE)+
  
  annotate("text", label="a", x=.67, y=5.7, size=5) +
  annotate("text", label="a", x=.89, y=5.7, size=5) +
  annotate("text", label="a", x=1.12, y=6.5, size=5) +
  annotate("text", label="a", x=1.34, y=6.75, size=5) +
  
  annotate("text", label="a", x=1.67, y=6, size=5) +
  annotate("text", label="a", x=1.89, y=5.5, size=5) +
  annotate("text", label="a", x=2.12, y=7, size=5) +
  annotate("text", label="a", x=2.34, y=6.5, size=5) +
  
  annotate("text", label="a", x=2.67, y=6, size=5) +
  annotate("text", label="a", x=2.89, y=5.5, size=5) +
  annotate("text", label="a", x=3.12, y=6.8, size=5) +
  annotate("text", label="a", x=3.34, y=6.5, size=5) +
  
  annotate("text", label="a", x=3.67, y=5.8, size=5) +
  annotate("text", label="a", x=3.89, y=5.3, size=5) +
  annotate("text", label="a", x=4.12, y=6.6, size=5) +
  annotate("text", label="a", x=4.34, y=6.1, size=5) +
  
  annotate("text", label="ab", x=4.67, y=6, size=5) +
  annotate("text", label="b", x=4.89, y=4.8, size=5) +
  annotate("text", label="a", x=5.12, y=6.7, size=5) +
  annotate("text", label="b", x=5.34, y=5.5, size=5) +
  
  annotate("text", label="ab", x=5.67, y=6, size=5) +
  annotate("text", label="b", x=5.89, y=5.2, size=5) +
  annotate("text", label="a", x=6.12, y=6.9, size=5) +
  annotate("text", label="b", x=6.34, y=5.7, size=5) +
  
  labs(list(x ="")) +
  theme_classic() +
  scale_fill_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position=c(.16,.9)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Day") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank(),
        legend.margin=margin(c(1,5,5,5))) +
  scale_y_continuous(bquote(~'Log'[10]~'('~CAOB~italic(amoA)~~ 'gene copies'~~gdw^-1*')'))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(plot.margin = unit(c(1,1,4,1), "lines")) +
  annotation_custom(text_1,xmin=1,xmax=1,ymin=-.75,ymax=-1) + 
  annotation_custom(text_42,xmin=2,xmax=2,ymin=-.75,ymax=-1) +
  annotation_custom(text_43,xmin=3,xmax=3,ymin=-.75,ymax=-1) + 
  annotation_custom(text_45,xmin=4,xmax=4,ymin=-.75,ymax=-1) +
  annotation_custom(text_56,xmin=5,xmax=5,ymin=-.75,ymax=-1) + 
  annotation_custom(text_84,xmin=6,xmax=6,ymin=-.75,ymax=-1) +
  coord_cartesian(clip="off") 
plot(bar_plotCAOB)

library(cowplot)

jpeg(filename="fig3.jpeg", bg="transparent", res=500, units = "in", height=16, width=12) 

f3 <- plot_grid(bar_plotMBC, bar_plotMBCN, bar_plotFB, ncol = 1, labels=c('A', 'B', 'C'),align="hv", label_size=30)

f3

dev.off()

jpeg(filename="fig4.jpeg", bg="transparent", res=500, units = "in", height=16, width=12) 

f4 <- plot_grid(bar_plotDOC, bar_plotTDN, bar_plotNH4, ncol = 1, labels=c('A', 'B', 'C'),align="hv", label_size=30)

f4

dev.off()

jpeg(filename="fig5.jpeg", bg="transparent", res=500, units = "in", height=16, width=12) 

f5 <- plot_grid(bar_plotNO3, bar_plotAOA,bar_plotAOB, ncol = 1, labels=c('A', 'B','C'),align="hv", label_size=30)

f5

dev.off()

