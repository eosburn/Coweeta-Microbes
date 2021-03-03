###This script runs all analyses using the fungal phylogenetic tree from PASTA-aligned ITS sequences, i.e., 'fungtree.tre' 
###These analyses are the ones presented in the main manuscript body. 

library(picante)
library(iCAMP)
library(biomformat)

Bac_otus <- "C:/Users/ernie/OneDrive/Desktop/phylogenetic analyses/FD-16S-table-97-rarefied.biom"

x1 <- read_biom(Bac_otus)

meta <- read.csv("C:/Users/ernie/OneDrive/Desktop/Data/Chapter 2 Data/FD_16s_metadata.csv")

bac <- as(biom_data(x1), "matrix")

bac <- t(bac)

library(ape)

tree <- read.tree("C:/Users/ernie/OneDrive/Desktop/phylogenetic analyses/16S_fasttree.tree")

samporder <- bac[, tree$tip.label]
samporder

meta<- meta[order(meta$id),]

fung_otus <- "C:/Users/ernie/OneDrive/Desktop/phylogenetic analyses/FD-ITS-Table-rarefied-97.biom"

x2 <- read_biom(fung_otus)

fung <- as(biom_data(x2), "matrix")

fung <- t(fung)

fungtree <-  read.tree("C:/Users/ernie/OneDrive/Desktop/phylogenetic analyses/fungtree.tre")

samporder2 <- fung[, fungtree$tip.label]
samporder2

phydist <- cophenetic(tree)

phydist2 <- cophenetic(fungtree)

#Testing for phylogenetic signal

library(funrar)

soil <- read.csv("C:/Users/ernie/OneDrive/Desktop/Data/Chapter 2 Data/FD_soil.csv")

veg <- read.csv("C:/Users/ernie/OneDrive/Desktop/Data/Chapter 2 Data/FD_veg.csv")

soil <- soil[order(soil$id),]

soil_vp <- soil[,c(6:19)]

library(caret)

soil.impute<- preProcess(soil_vp, method = "bagImpute")

soil_vp2 <- predict(soil.impute, soil_vp)

bac2 <- make_relative(bac)

soil2 <- soil_vp2[,c(1:7,12:14)]

bac3 <- cbind(bac2,soil2)

ph <- t(data.frame(lapply(bac3[ ,1:3551], weighted.mean, x = bac3$pH),check.names=FALSE))

colnames(ph)[1] <- "ph"

NH4 <- t(data.frame(lapply(bac3[ ,1:3551], weighted.mean, x = bac3$NH4),check.names=FALSE))

colnames(NH4)[1] <- "NH4"

NO3 <- t(data.frame(lapply(bac3[ ,1:3551], weighted.mean, x = bac3$NO3),check.names=FALSE))

colnames(NO3)[1] <- "NO3"

TDN <- t(data.frame(lapply(bac3[ ,1:3551], weighted.mean, x = bac3$TDN),check.names=FALSE))

colnames(TDN)[1] <- "TDN"

DOC <- t(data.frame(lapply(bac3[ ,1:3551], weighted.mean, x = bac3$DOC),check.names=FALSE))

colnames(DOC)[1] <- "DOC"

DOC.TDN <- t(data.frame(lapply(bac3[ ,1:3551], weighted.mean, x = bac3$DOC.TDN),check.names=FALSE))

colnames(DOC.TDN)[1] <- "DOC.TDN"

C.N <- t(data.frame(lapply(bac3[ ,1:3551], weighted.mean, x = bac3$C.N),check.names=FALSE))

colnames(C.N)[1] <- "C.N"

soil_otus <- cbind(ph,NH4,TDN,NO3,DOC,C.N,DOC.TDN)

soil_otus <- scale(soil_otus)

soil_dist <- vegdist(soil_otus, method="euclidean")

soil_dist.2 <- as.matrix(soil_dist)

library(reshape2)
phydist.2 <- setNames(melt(phydist), c('taxa1', 'taxa2', 'dist'))

phydist.3 <- phydist.2[order(phydist.2[,1], phydist.2[,2]),]

library(tidyr)
phydist.4 <- spread(phydist.3, taxa1, dist)

phydist.4 <- data.frame(phydist.4, row.names = 1,check.names=FALSE)

phydist.4 <- as.matrix(phydist.4)

soil_dist.3<- setNames(melt(soil_dist.2), c('taxa1', 'taxa2', 'dist'))

soil_dist.4 <- soil_dist.3[order(soil_dist.3[,1], soil_dist.3[,2]), ]

soil_dist.5 <- spread(soil_dist.4, taxa1, dist)

soil_dist.5 <- data.frame(soil_dist.5, row.names = 1,check.names=FALSE)

soil_dist.5 <- as.matrix(soil_dist.5)

phydist.5 <- decostand(phydist.4, method="range")

cor1 <- mantel.correlog(soil_dist.5, phydist.5, n.class=100, mult="bonferroni",nperm=999)

plot(cor1)

#Phylogenetic signal - fungi

fung2 <- make_relative(fung)

fung3 <- cbind(fung2,soil2)

ph2 <- t(data.frame(lapply(fung3[ ,1:1197], weighted.mean, x = fung3$pH),check.names=FALSE))

colnames(ph2)[1] <- "ph"

NH42 <- t(data.frame(lapply(fung3[ ,1:1197], weighted.mean, x = fung3$NH4),check.names=FALSE))

colnames(NH42)[1] <- "NH4"

NO32 <- t(data.frame(lapply(fung3[ ,1:1197], weighted.mean, x = fung3$NO3),check.names=FALSE))

colnames(NO32)[1] <- "NO3"

TDN2 <- t(data.frame(lapply(fung3[ ,1:1197], weighted.mean, x = fung3$TDN),check.names=FALSE))

colnames(TDN2)[1] <- "TDN"

DOC2 <- t(data.frame(lapply(fung3[ ,1:1197], weighted.mean, x = fung3$DOC),check.names=FALSE))

colnames(DOC2)[1] <- "DOC"

DOC.TDN2 <- t(data.frame(lapply(fung3[ ,1:1197], weighted.mean, x = fung3$DOC.TDN),check.names=FALSE))

colnames(DOC.TDN2)[1] <- "DOC.TDN"

C.N2 <- t(data.frame(lapply(fung3[ ,1:1197], weighted.mean, x = fung3$C.N),check.names=FALSE))

colnames(C.N2)[1] <- "C.N"

soil_otus2 <- cbind(ph2,NH42,TDN2,NO32,DOC2,C.N2,DOC.TDN2)

soil_otus2 <- scale(soil_otus2)

soil_dist2 <- vegdist(soil_otus2, method="euclidean")

soil_dist2.2 <- as.matrix(soil_dist2)

library(reshape2)
phydist2.2 <- setNames(melt(phydist2), c('taxa1', 'taxa2', 'dist'))

phydist2.2$taxa1 <- as.character(phydist2.2$taxa1)

phydist2.2$taxa2 <- as.character(phydist2.2$taxa2)

phydist2.3 <- phydist2.2[order(phydist2.2$taxa1, phydist2.2$taxa2),]

library(tidyr)
phydist2.4 <- spread(phydist2.3, taxa1, dist)

phydist2.4 <- data.frame(phydist2.4, row.names = 1,check.names=FALSE)

phydist2.4 <- as.matrix(phydist2.4)

soil_dist2.3<- setNames(melt(soil_dist2.2), c('taxa1', 'taxa2', 'dist'))

soil_dist2.3$taxa1 <- as.character(soil_dist2.3$taxa1)

soil_dist2.3$taxa2 <- as.character(soil_dist2.3$taxa2)

soil_dist2.4 <- soil_dist2.3[order(soil_dist2.3[,1], soil_dist2.3[,2]), ]

soil_dist2.5 <- spread(soil_dist2.4, taxa1, dist)

soil_dist2.5 <- data.frame(soil_dist2.5, row.names = 1,check.names=FALSE)

soil_dist2.5 <- as.matrix(soil_dist2.5)

phydist2.5 <- decostand(phydist2.4, method="range")

cor2 <- mantel.correlog(soil_dist2.5, phydist2.5,mult="bonferroni", n.class=100,nperm=999)

plot(cor2)

#Community Assembly

#All bacteria

assem_bac <- qpen(comm = bac, pd = phydist, pd.big.wd = NULL,
                  pd.big.spname = NULL, tree = NULL,
                  bNTI = NULL, RC = NULL, ab.weight = TRUE,
                  meta.ab = NULL, exclude.conspecifics = FALSE,
                  rand.time = 1000, sig.bNTI = 1.96, sig.rc = 0.95,
                  nworker = 4, memory.G = 50, project = NA,
                  wd = getwd(), output.detail = FALSE, save.bNTIRC = FALSE)

bac_assem_ratio <- data.frame(assem_bac[["ratio"]])

bac_assem_result <- data.frame(assem_bac[["result"]])

bac_assem_result['Taxon']= 'Bacteria'


#All fungi

assem_fung <- qpen(comm = fung, pd = phydist2, pd.big.wd = NULL,
                   pd.big.spname = NULL, tree = NULL,
                   bNTI = NULL, RC = NULL, ab.weight = TRUE,
                   meta.ab = NULL, exclude.conspecifics = FALSE,
                   rand.time = 1000, sig.bNTI = 1.96, sig.rc = 0.95,
                   nworker = 4, memory.G = 50, project = NA,
                   wd = getwd(), output.detail = FALSE, save.bNTIRC = FALSE)

fung_assem_ratio <- data.frame(assem_fung[["ratio"]])

fung_assem_result <- data.frame(assem_fung[["result"]])

fung_assem_result['Taxon']= 'Fungi'


bnti <- rbind(bac_assem_result,fung_assem_result)



library(car)
m5 <- lm(bNTI~Taxon,data=bnti)
shapiro.test(resid(m5))
Anova(m5)
kruskal.test(bNTI~Taxon,data=bnti)

bnti2<-bnti[!(bnti$bNTI<6),]

bnti3<-bnti[!(bnti$bNTI>6),]

library(ggplot2)
viol <-ggplot(bnti3, aes(x=Taxon, y=bNTI)) + 
  geom_violin(trim=TRUE,fill="gray")+
  labs(x="Taxon", y ="ß-NTI")+
  geom_boxplot(width=0.1,outlier.shape=NA)+
  coord_cartesian(ylim = c(-6,6)) +
  theme_classic()+
  geom_hline(yintercept=2, linetype="dashed", color = "black",size=1.5) +
  geom_hline(yintercept=-2, linetype="dashed", color = "black",size=1.5) +
  theme(axis.title=element_text(size=20)) +
  annotate("text", x=1.85,y=-4.5,label="Taxon: ~italic(P) < 0.001", parse=TRUE, size=5, fontface="bold")+
  annotate("rect", xmin=-Inf, xmax=Inf, ymin = 2, ymax = Inf, alpha = .3, fill = "#0571B0") +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin = -Inf, ymax = -2, alpha = .3, fill = "#0571B0") +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin = -2, ymax = 2, alpha = .3, fill = "#CA0020") +
  theme(text=element_text(size=20)) 
viol

#Reference bacteria

bac2 <- cbind(bac, meta)

bac_ref <- subset(bac2, Treatment=="Reference")

bac_ref2 <- data.frame(bac_ref[,c(1:3551)], check.names=FALSE)

bac_ref3 <- bac_ref2[, colSums(bac_ref2 != 0) > 0]

bac_ref4 <- as.matrix(bac_ref3)

reftree <- prune.sample(bac_ref4, tree)

distref <- cophenetic(reftree)

samporderref <- bac_ref4[, reftree$tip.label]
samporderref

assem_ref <- qpen(comm = bac_ref4, pd = distref, pd.big.wd = NULL,
                  pd.big.spname = NULL, tree = NULL,
                  bNTI = NULL, RC = NULL, ab.weight = TRUE,
                  meta.ab = NULL, exclude.conspecifics = FALSE,
                  rand.time = 1000, sig.bNTI = 1.96, sig.rc = 0.95,
                  nworker = 4, memory.G = 50, project = NA,
                  wd = getwd(), output.detail = FALSE, save.bNTIRC = FALSE)

ref_assem_ratio <- data.frame(assem_ref[["ratio"]])

ref_assem_result <- data.frame(assem_ref[["result"]])

ref_assem_ratio['Treatment']= 'Reference'

#Disturbed bacteria

bac_dis <- subset(bac2, Treatment=="Disturbed")

bac_dis2 <- data.frame(bac_dis[,c(1:3551)], check.names=FALSE)

bac_dis3 <- bac_dis2[, colSums(bac_dis2 != 0) > 0]

bac_dis4 <- as.matrix(bac_dis3)

distree <- prune.sample(bac_dis4, tree)

distdis <- cophenetic(distree)

samporderdis <- bac_dis4[, distree$tip.label]
samporderdis

assem_dis <- qpen(comm = bac_dis4, pd = distdis, pd.big.wd = NULL,
                  pd.big.spname = NULL, tree = NULL,
                  bNTI = NULL, RC = NULL, ab.weight = TRUE,
                  meta.ab = NULL, exclude.conspecifics = FALSE,
                  rand.time = 1000, sig.bNTI = 1.96, sig.rc = 0.95,
                  nworker = 4, memory.G = 50, project = NA,
                  wd = getwd(), output.detail = FALSE, save.bNTIRC = FALSE)


dis_assem_ratio <- data.frame(assem_dis[["ratio"]])

dis_assem_result <- data.frame(assem_dis[["result"]])

dis_assem_ratio['Treatment']= 'Disturbed'

#z tests for assembly processes in bacteria

#Selection overall

n_s_r=sum((ref_assem_result$process == "Homogeneous.Selection")+(ref_assem_result$process == "Heterogeneous.Selection")) #total number of negative edges for the disturbed
t_s_r=nrow(ref_assem_result)#total edges for the disturbed
n_s_d=sum((dis_assem_result$process == "Homogeneous.Selection")+(dis_assem_result$process == "Heterogeneous.Selection"))  #total number of positive edges for the disturbed
t_s_d=nrow(dis_assem_result)#total edges for the reference

prop.test(c(n_s_r,n_s_d),c(t_s_r,t_s_d))

#Homo selection

n_hs_r=sum(ref_assem_result$process == "Homogeneous.Selection") #total number of negative edges for the disturbed
t_hs_r=nrow(ref_assem_result)#total edges for the disturbed
n_hs_d=sum(dis_assem_result$process == "Homogeneous.Selection") #total number of positive edges for the disturbed
t_hs_d=nrow(dis_assem_result)#total edges for the reference

prop.test(c(n_hs_r,n_hs_d),c(t_hs_r,t_hs_d))

#Homo dispersal

n_hd_r=sum(ref_assem_result$process == "Homogenizing.Dispersal") #total number of negative edges for the disturbed
t_hd_r=nrow(ref_assem_result)#total edges for the disturbed
n_hd_d=sum(dis_assem_result$process == "Homogenizing.Dispersal") #total number of positive edges for the disturbed
t_hd_d=nrow(dis_assem_result)#total edges for the reference

prop.test(c(n_hd_r,n_hd_d),c(t_hd_r,t_hd_d))

#Dispersal Limitation

n_dl_r=sum(ref_assem_result$process == "Dispersal.Limitation") #total number of negative edges for the disturbed
t_dl_r=nrow(ref_assem_result)#total edges for the disturbed
n_dl_d=sum(dis_assem_result$process == "Dispersal.Limitation") #total number of positive edges for the disturbed
t_dl_d=nrow(dis_assem_result)#total edges for the reference

prop.test(c(n_dl_r,n_dl_d),c(t_dl_r,t_dl_d))

#Drift

n_d_r=sum(ref_assem_result$process == "Undominated") #total number of negative edges for the disturbed
t_d_r=nrow(ref_assem_result)#total edges for the disturbed
n_d_d=sum(dis_assem_result$process == "Undominated") #total number of positive edges for the disturbed
t_d_d=nrow(dis_assem_result)#total edges for the reference

prop.test(c(n_d_r,n_d_d),c(t_d_r,t_d_d))

bac_assem <- rbind(ref_assem_ratio,dis_assem_ratio)


library(reshape2)
bac_assem2 <- melt(bac_assem, id.vars = c("Treatment"))

bac_assem2 <- bac_assem2[c(1:10),]

bac_assem2$variable <- as.character(bac_assem2$variable)

bac_assem2[bac_assem2 =="Undominated"] <- "Drift"
bac_assem2[bac_assem2 =="Dispersal.Limitation"] <- "Dispersal Limitation"
bac_assem2[bac_assem2 =="Homogenizing.Dispersal"] <- "Homogenizing Dispersal"
bac_assem2[bac_assem2 =="Homogeneous.Selection"] <- "Homogeneous Selection"
bac_assem2[bac_assem2 =="Heterogeneous.Selection"] <- "Heterogeneous Selection"

bac_assem2$variable <- factor(bac_assem2$variable, levels = c("Homogenizing Dispersal", "Dispersal Limitation", "Drift", "Heterogeneous Selection","Homogeneous Selection"))
bac_assem2$variable

bac_bar <-ggplot(bac_assem2, aes(fill=variable, y=value, x=Treatment)) + 
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous("Proportion") +
  scale_fill_manual(values = c("#CA0020", "#F4A582", "#FDDBC7","#92C5DE", "#0571B0")) +
  theme_classic() +
  #theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill="Assembly Process") +
  ggtitle("Bacteria") +
  theme(legend.position = "none")+
  theme(axis.title=element_text(size=20)) +
  annotate("text", x=2,y=.88,label='*', size=8, fontface="bold")+
  #annotate("text", x=1,y=.82,label='*', size=8, fontface="bold")+
  annotate("text", x=1,y=.23,label='***', size=8, fontface="bold")+
  theme(text=element_text(size=20)) 
bac_bar

#Reference Fungi

fung2 <- cbind(fung, meta)

fung_ref <- subset(fung2, Treatment=="Reference")

fung_ref2 <- data.frame(fung_ref[,c(1:1197)],check.names=FALSE)

fung_ref3 <- fung_ref2[, colSums(fung_ref2 != 0) > 0]

fung_ref4 <- as.matrix(fung_ref3)

reftree2 <- prune.sample(fung_ref4, fungtree)

distref2 <- cophenetic(reftree2)

samporderref2 <- fung_ref4[, reftree2$tip.label]
samporderref2

assem_ref2 <- qpen(comm = fung_ref4, pd = distref2, pd.big.wd = NULL,
                   pd.big.spname = NULL, tree = NULL,
                   bNTI = NULL, RC = NULL, ab.weight = TRUE,
                   meta.ab = NULL, exclude.conspecifics = FALSE,
                   rand.time = 1000, sig.bNTI = 2, sig.rc = 0.95,
                   nworker = 4, memory.G = 50, project = NA,
                   wd = getwd(), output.detail = FALSE, save.bNTIRC = FALSE)

ref_assem_ratio2 <- data.frame(assem_ref2[["ratio"]])

ref_assem_result2 <- data.frame(assem_ref2[["result"]])

#Disturbed Fungi

fung_dis <- subset(fung2, Treatment=="Disturbed")

fung_dis2 <- data.frame(fung_dis[,c(1:1197)], check.names=FALSE)

fung_dis3 <- fung_dis2[, colSums(fung_dis2 != 0) > 0]

fung_dis4 <- as.matrix(fung_dis3)

distree2 <- prune.sample(fung_dis4, fungtree)

distdis2 <- cophenetic(distree2)

samporderdis2 <- fung_dis4[, distree2$tip.label]
samporderdis2

assem_dis2 <- qpen(comm = fung_dis4, pd = distdis2, pd.big.wd = NULL,
                   pd.big.spname = NULL, tree = NULL,
                   bNTI = NULL, RC = NULL, ab.weight = TRUE,
                   meta.ab = NULL, exclude.conspecifics = FALSE,
                   rand.time = 1000, sig.bNTI = 1.96, sig.rc = 0.95,
                   nworker = 2, memory.G = 50, project = NA,
                   wd = getwd(), output.detail = FALSE, save.bNTIRC = FALSE)


dis_assem_ratio2 <- data.frame(assem_dis2[["ratio"]])

dis_assem_result2 <- data.frame(assem_dis2[["result"]])

#z tests for assembly processes in fungi

#Heterogeneous Selection

n_hd_r2=sum(ref_assem_result2$process == "Heterogeneous.Selection") #total number of negative edges for the disturbed
t_hd_r2=nrow(ref_assem_result2)#total edges for the disturbed
n_hd_d2=sum(dis_assem_result2$process == "Heterogeneous.Selection") #total number of positive edges for the disturbed
t_hd_d2=nrow(dis_assem_result2)#total edges for the reference

prop.test(c(n_hd_r2,n_hd_d2),c(t_hd_r2,t_hd_d2))

#Dispersal Limitation

n_dl_r2=sum(ref_assem_result2$process == "Dispersal.Limitation") #total number of negative edges for the disturbed
t_dl_r2=nrow(ref_assem_result2)#total edges for the disturbed
n_dl_d2=sum(dis_assem_result2$process == "Dispersal.Limitation") #total number of positive edges for the disturbed
t_dl_d2=nrow(dis_assem_result2)#total edges for the reference

prop.test(c(n_dl_r2,n_dl_d2),c(t_dl_r2,t_dl_d2))

#Drift

n_d_r2=sum(ref_assem_result2$process == "Undominated") #total number of negative edges for the disturbed
t_d_r2=nrow(ref_assem_result2)#total edges for the disturbed
n_d_d2=sum(dis_assem_result2$process == "Undominated") #total number of positive edges for the disturbed
t_d_d2=nrow(dis_assem_result2)#total edges for the reference

prop.test(c(n_d_r2,n_d_d2),c(t_d_r2,t_d_d2))

dis_assem_ratio2['Treatment']= 'Disturbed'

ref_assem_ratio2['Treatment']= 'Reference'

fung_assem <- rbind(ref_assem_ratio2,dis_assem_ratio2)


library(reshape2)
fung_assem2 <- melt(fung_assem, id.vars = c("Treatment"))

fung_assem2 <- fung_assem2[c(1:10),]

fung_assem2$variable <- as.character(fung_assem2$variable)

fung_assem2[fung_assem2 =="Undominated"] <- "Drift"
fung_assem2[fung_assem2=="Dispersal.Limitation"] <- "Dispersal Limitation"
fung_assem2[fung_assem2 =="Homogenizing.Dispersal"] <- "Homogenizing Dispersal"
fung_assem2[fung_assem2 =="Homogeneous.Selection"] <- "Homogeneous Selection"
fung_assem2[fung_assem2 =="Heterogeneous.Selection"] <- "Heterogeneous Selection"

fung_assem2$variable <- factor(fung_assem2$variable, levels = c("Homogenizing Dispersal", "Dispersal Limitation", "Drift", "Heterogeneous Selection","Homogeneous Selection"))
fung_assem2$variable

fung_bar <-ggplot(fung_assem2, aes(fill=variable, y=value, x=Treatment)) + 
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous("Proportion") +
  scale_fill_manual(values = c("#CA0020", "#F4A582", "#FDDBC7","#92C5DE", "#0571B0")) +
  theme_classic() +
  #theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill="Assembly Process") +
  ggtitle("Fungi") +
  annotate("text", x=2,y=.6,label='*', size=8, fontface="bold")+
  annotate("text", x=1,y=.04,label='***', size=8, fontface="bold")+
  theme(legend.position = "none")+
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) 
fung_bar

fung_bar2 <-ggplot(fung_assem2, aes(fill=variable, y=value, x=Treatment)) + 
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous("Proportion") +
  scale_fill_manual(values = c("#CA0020", "#F4A582", "#FDDBC7","#92C5DE", "#0571B0")) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill="Assembly Process") +
  ggtitle("Fungi") +
  annotate("text", x=2,y=.6,label='*', size=8, fontface="bold")+
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=15)) 
fung_bar2

library(cowplot)
legend <- get_legend(fung_bar2+ theme(plot.margin = margin(-100, -100, -100, -100)))

plot(legend)

setwd("C:/Users/ernie/OneDrive/Desktop/phylogenetic analyses")

jpeg(filename="fig3.jpeg", bg="transparent", res=500, units = "in", height=6, width=15) 

f3 <- plot_grid(viol, bac_bar, fung_bar, legend,ncol = 4, labels=c('A', 'B', 'C',''),align="hv", label_size=25)
f3

dev.off()

#Variance partitioning - Bacteria 

#spatial stuff

library(Imap)

ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}

coords <- read.csv("C:/Users/ernie/OneDrive/Desktop/phylogenetic analyses/coords.csv")

coords$name <- as.character(coords$name)

df.cities <- coords

coords2 <- GeoDistanceInMetresMatrix(df.cities)

space <- pcnm(coords2)

space2 <- scores(space)

a <- plot(x=coords$lat, y=coords$lon, data=coords)

ordisurf(coords[,c(2:3)],scores(space,choi=1),bubble=4,
         main='PCNM 1')

ordisurf(coords[,c(2:3)],scores(space,choi=2),bubble=4,
         main='PCNM 2')

ordisurf(coords[,c(2:3)],scores(space,choi=3),bubble=4,
         main='PCNM 3')

ordisurf(coords[,c(2:3)],scores(space,choi=4),bubble=4,
         main='PCNM 4')


veg <- veg[order(veg$id),]

veg_vp <- veg[,c(7:22)]

veg_vp2 <- decostand(veg_vp, method="hellinger")

bac_bnti <- bac_assem_result[,c(1,2,5)]

bac_bnti <- bac_bnti[order(bac_bnti[,1], bac_bnti[,2]),]

library(tidyr)
bac_bnti  <- spread(bac_bnti , sample1, bNTI)

bac_bnti <- data.frame(bac_bnti, row.names = 1,check.names=FALSE)

bac_bnti <- as.matrix(bac_bnti)

bac_bnti <- cbind(FD1 = NA, bac_bnti)

bac_bnti <- rbind(bac_bnti,FD9 = NA )

diag(bac_bnti) <- 0

library(Matrix)

bac_bnti<- as.matrix(forceSymmetric(bac_bnti))

as.dist(bac_bnti)

bac_bnti1 <- decostand(bac_bnti, method="range")

soil_vp3 <- data.frame(scale(soil_vp2))

cap <- capscale(as.dist(bac_bnti1) ~ .,data=soil_vp3)

ord <- ordistep(cap)

anova(ord)

anova.cca(ord,by="term")

soil_vp4 <- soil_vp3[,c(5,6,7,10,13)]

cap <- capscale(as.dist(bac_bnti1) ~ .,data=veg_vp2)

ord <- ordistep(cap)

anova(ord)

anova.cca(ord, by="term")

veg_vp3 <- veg_vp2[,c(2,3,8,12)]

cap <- capscale(as.dist(bac_bnti1) ~ .,data=as.data.frame(space2))

ord <- ordistep(cap)

anova(ord)

anova.cca(ord, by="term")

space3 <- space2[,c(1,3,4)]

varpart_16s <- varpart(as.dist(bac_bnti1),soil_vp4,veg_vp3)

varpart_16s

plot(varpart_16s)

ord1 <- capscale(as.dist(bac_bnti1)~as.matrix(soil_vp4))

anova(ord1)

ord2 <- capscale(as.dist(bac_bnti1)~as.matrix(soil_vp4)+Condition(as.matrix(veg_vp3)))

anova(ord2)

ord3 <- capscale(as.dist(bac_bnti1)~as.matrix(veg_vp3))

anova(ord3)

ord4 <- capscale(as.dist(bac_bnti1)~as.matrix(veg_vp3)+Condition(as.matrix(soil_vp4)))

anova(ord4)

ord5 <- capscale(as.dist(bac_bnti1)~as.matrix(veg_vp3)+as.matrix(soil_vp4))

anova(ord5)

vars <- cbind(soil_vp4, veg_vp3)

ord6 <- capscale(as.dist(bac_bnti1)~vars$MBN+vars$DOC+vars$TDN+vars$DOC.TDN+vars$TN+Condition(vars$qrub+vars$blen+vars$ltul+vars$pstr))

anova(ord6)

anova.cca(ord6, by="term")

ord7 <- capscale(as.dist(bac_bnti1)~vars$qrub+vars$blen+vars$ltul+vars$pstr+Condition(vars$MBN+vars$DOC+vars$TDN+vars$DOC.TDN+vars$TN))

anova(ord7)

anova.cca(ord7, by="term")

#Variance partitioning - Fungal 

fung_bnti <-fung_assem_result[,c(1,2,5)]

fung_bnti <- fung_bnti[order(fung_bnti[,1], fung_bnti[,2]),]

library(tidyr)
fung_bnti  <- spread(fung_bnti , sample1, bNTI)

fung_bnti <- data.frame(fung_bnti, row.names = 1,check.names=FALSE)

fung_bnti <- as.matrix(fung_bnti)

fung_bnti <- cbind(FD1 = NA, fung_bnti)

fung_bnti <- rbind(fung_bnti,FD9 = NA )

diag(fung_bnti) <- 0

library(Matrix)

fung_bnti<- as.matrix(forceSymmetric(fung_bnti))

fung_bnti1 <- decostand(fung_bnti, method="range")

cap <- capscale(as.dist(fung_bnti1) ~ .,data=soil_vp3)

ord <- ordistep(cap)

anova(ord)

anova.cca(ord, by="term")

soil_vp5 <- soil_vp3[,c(6,7,9,10,13)]

cap <- capscale(as.dist(fung_bnti1) ~ .,data=veg_vp2)

ord <- ordistep(cap)

anova(ord)

anova.cca(ord, by="term")

veg_vp4 <- veg_vp2[,c(1,7,8,12)]

cap <- capscale(as.dist(fung_bnti1) ~ .,data=as.data.frame(space2))

ord <- ordistep(cap)

anova(ord)

anova.cca(ord, by="term")

space3 <- space2[,c(1,3)]

varpart_its <- varpart(as.dist(fung_bnti1),soil_vp5,veg_vp4)

varpart_its

plot(varpart_its)

ord1 <- capscale(as.dist(fung_bnti1)~as.matrix(soil_vp5))

anova(ord1)

ord2 <- capscale(as.dist(fung_bnti1)~as.matrix(soil_vp5)+Condition(as.matrix(veg_vp4)))

anova(ord2)

ord3 <- capscale(as.dist(fung_bnti1)~as.matrix(veg_vp4))

anova(ord3)

ord4 <- capscale(as.dist(fung_bnti1)~as.matrix(veg_vp4)+Condition(as.matrix(soil_vp5)))

anova(ord4)

ord5 <- capscale(as.dist(fung_bnti1)~as.matrix(veg_vp4)+as.matrix(soil_vp5))

anova(ord5)

vars <- cbind(soil_vp5, veg_vp4)

ord6 <- capscale(as.dist(fung_bnti1)~vars$MBN+vars$MBC+vars$TDN+vars$DOC.TDN+vars$TN+Condition(vars$rmax+vars$carya+vars$ltul+vars$pstr))

anova(ord6)

anova.cca(ord6, by="term")

ord7 <- capscale(as.dist(fung_bnti1)~vars$rmax+vars$carya+vars$ltul+vars$pstr+Condition(vars$MBN+vars$MBC+vars$TDN+vars$DOC.TDN+vars$TN))

anova(ord7)

anova.cca(ord7, by="term")

