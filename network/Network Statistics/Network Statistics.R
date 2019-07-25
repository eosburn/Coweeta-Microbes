###Network Analysis###
###Network Statistics for Fungi and Bacteria###
###Coweeta Microbe project###
###Modified from https://github.com/ryanjw/co-occurrence###
###Modified by Steve McBride; mcbridsg@vt.edu
###last Modified 07/25/2019###
###generates network graphs##
###Data from Ernest Osburn Forest  Disturbance study at Coweeta###
###Data collected in Summer 2018###
###Files Needed to begin analysis: created in Sequence Processing and Network Correlation Matrix Generation.R
###"Reference_OTU_correlation_matrix_min10_16S.csv", "Disturbed_OTU_correlation_matrix_min10_16S.csv"
###"Reference_OTU_correlation_matrix_min10_ITS.csv", "Disturbed_OTU_correlation_matrix_min10_ITS.csv"
###Files Generated and imported 
###"dis16snetstats.csv","ref16snetstats.csv","disITSnetstats.csv","refITSnetstats.csv"
###All files generated and used in the following code are available https://github.com/eosburn/Coweeta-Microbes



#####Create environment and load necessary packages#####

### Set your working directory. This should be a folder that you have access to that
### all of your files will be stored in.
setwd('your/workingdirectory/')

library(igraph)

######Import Correlation matrix data####

##16S
refcormat_16S=read.csv(file="Reference_OTU_correlation_matrix_min10_16S.csv", header=TRUE)
discormat_16S=read.csv(file="Disturbed_OTU_correlation_matrix_min10_16S.csv", header=TRUE)

#selects all values with p.value column above 0.01
sigdis16s=subset(discormat_16S,p.value <=0.01)
sigref16s=subset(refcormat_16S,p.value<=0.01)

#Remove the unncessary first columen
sigdis16s=sigdis16s[,-1]
sigref16s=sigref16s[,-1]

##ITS
refcormat_ITS=read.csv(file="Reference_OTU_correlation_matrix_min10_ITS.csv", header=TRUE)
discormat_ITS=read.csv(file="Disturbed_OTU_correlation_matrix_min10_ITS.csv", header=TRUE)

#selects all values with p.value column above 0.01
sigdisITS=subset(discormat_ITS,p.value <=0.01)
sigrefITS=subset(refcormat_ITS,p.value<=0.01)

#Remove the unncessary first columen
sigdisITS=sigdisITS[,-1]
sigrefITS=sigrefITS[,-1]


#####Create Network Statistics 16S####
#Disturbed 16S
dis16S<-simplify(graph_from_edgelist(as.matrix(sigdis16s[,c(2,3)]),directed=FALSE))# The columns selected here '[,c(2,3)]' should match the columns with the different taxa
dis16S_degree=(as.matrix(degree(dis16S)))
dis16S_betweenness=(as.matrix(betweenness(dis16S)))
dis16S_closeness=(as.matrix(closeness(dis16S)))
dis16S_transitivity=as.matrix(transitivity(dis16S,type="local"))
dis16S_norm_degree=(as.matrix(degree(dis16S, normalized=TRUE)))
dis16S_norm_between=(as.matrix(betweenness(dis16S, normalized=TRUE)))
dis16S_norm_close=(as.matrix(closeness(dis16S, normalized=TRUE)))
dis16S_all=cbind(dis16S_degree,dis16S_betweenness,dis16S_closeness,dis16S_transitivity,dis16S_norm_degree, 
                 dis16S_norm_between, dis16S_norm_close)#binds the columns for each of the different network stats
#write.csv(dis16S_all, "dis16snetstats.csv") ##I had trouble using ggplot on the data without first converting to .csv then importing that .csv
dis16S_all=read.csv(file="dis16snetstats.csv", header=TRUE)

#Reference 16S
ref16S<-simplify(graph_from_edgelist(as.matrix(sigref16s[,c(2,3)]),directed=FALSE))# The columns selected here '[,c(2,3)]' should match the columns with the different taxa
ref16S_degree=(as.matrix(degree(ref16S)))
ref16S_betweenness=(as.matrix(betweenness(ref16S)))
ref16S_closeness=(as.matrix(closeness(ref16S)))
ref16S_transitivity=as.matrix(transitivity(ref16S,type="local"))
ref16S_norm_degree=(as.matrix(degree(ref16S, normalized=TRUE)))
ref16S_norm_between=(as.matrix(betweenness(ref16S, normalized=TRUE)))
ref16S_norm_close=(as.matrix(closeness(ref16S, normalized=TRUE)))
ref16S_all=cbind(ref16S_degree,ref16S_betweenness,ref16S_closeness,ref16S_transitivity,ref16S_norm_degree, 
                 ref16S_norm_between, ref16S_norm_close)#binds the columns for each of the different network stats

#write.csv(ref16S_all, "ref16snetstats.csv") ##I had trouble using ggplot on the data without first converting to .csv then importing that .csv
ref16S_all=read.csv(file="ref16snetstats.csv", header=TRUE)

###Add treatment to individual dataframes
dis16S_all$trt="Disturbed"
ref16S_all$trt="Reference"

#combine dataframes
netstats_16S=rbind(ref16S_all,dis16S_all)

#name columns and convert trt to a factor
colnames(netstats_16S)=c("X","degree", "betweenness","closeness","clustering_coeff","norm_degree","norm_between","norm_close","trt")
netstats_16S$trt=as.factor(netstats_16S$trt)

#Stats for 16S
kruskal.test(data=netstats_16S,norm_degree~trt)
kruskal.test(data=netstats_16S,norm_between~trt)
kruskal.test(data=netstats_16S,norm_close~trt)
kruskal.test(data=netstats_16S,clustering_coeff~trt)

#z test for negative edges
nd16s=sum(sigdis16s$rho < 0) #total number of negative edges for the disturbed
td16s=nrow(sigdis16s)#total edges for the disturbed
nr16s=sum(sigref16s$rho < 0) #total number of positive edges for the disturbed
tr16s=nrow(sigref16s)#total edges for the reference

prop.test(c(nd16s,nr16s),c(td16s,tr16s), correct=FALSE)


#####Create Network Statistics ITS#####
#Disturbed ITS
disITS<-simplify(graph_from_edgelist(as.matrix(sigdisITS[,c(2,3)]),directed=FALSE))# The columns selected here '[,c(2,3)]' should match the columns with the different taxa
disITS_degree=(as.matrix(degree(disITS)))
disITS_betweenness=(as.matrix(betweenness(disITS)))
disITS_closeness=(as.matrix(closeness(disITS)))
disITS_transitivity=as.matrix(transitivity(disITS,type="local"))
disITS_norm_degree=(as.matrix(degree(disITS, normalized=TRUE)))
disITS_norm_between=(as.matrix(betweenness(disITS, normalized=TRUE)))
disITS_norm_close=(as.matrix(closeness(disITS, normalized=TRUE)))
disITS_all=cbind(disITS_degree,disITS_betweenness,disITS_closeness,disITS_transitivity,disITS_norm_degree, 
                 disITS_norm_between, disITS_norm_close)#binds the columns for each of the different network stats

#write.csv(disITS_all, "disITSnetstats.csv") ##I had trouble using ggplot on the data without first converting to .csv then importing that .csv
disITS_all=read.csv(file="disITSnetstats.csv", header=TRUE)

#Reference ITS
refITS<-simplify(graph_from_edgelist(as.matrix(sigrefITS[,c(2,3)]),directed=FALSE))# The columns selected here '[,c(2,3)]' should match the columns with the different taxa
refITS_degree=(as.matrix(degree(refITS)))
refITS_betweenness=(as.matrix(betweenness(refITS)))
refITS_closeness=(as.matrix(closeness(refITS)))
refITS_transitivity=as.matrix(transitivity(refITS,type="local"))
refITS_norm_degree=(as.matrix(degree(refITS, normalized=TRUE)))
refITS_norm_between=(as.matrix(betweenness(refITS, normalized=TRUE)))
refITS_norm_close=(as.matrix(closeness(refITS, normalized=TRUE)))
refITS_all=cbind(refITS_degree,refITS_betweenness,refITS_closeness,refITS_transitivity,refITS_norm_degree, 
                 refITS_norm_between, refITS_norm_close)#binds the columns for each of the different network stats


refITS_all$trt=reference
#write.csv(refITS_all, "refITSnetstats.csv") ##I had trouble manipulating the data without first converting to .csv then importing that .csv
refITS_all=read.csv(file="refITSnetstats.csv", header=TRUE)

###Add treatment to individual dataframes
disITS_all$trt="Disturbed"
refITS_all$trt="Reference"

#combine dataframes
ITSnetstats=rbind(refITS_all,disITS_all)

#name columns and convert trt to a factor
colnames(ITSnetstats)=c("X","degree", "betweenness","closeness","clustering_coeff","norm_degree","norm_between","norm_close","trt")
ITSnetstats$trt=as.factor(ITSnetstats$trt)

#Stats for ITS
kruskal.test(data=ITSnetstats,norm_degree~trt)
kruskal.test(data=ITSnetstats,norm_between~trt)
kruskal.test(data=ITSnetstats,norm_close~trt)
kruskal.test(data=ITSnetstats,clustering_coeff~trt)


#z test for negative edges
ndITS=sum(sigdisITS$rho < 0) #total number of negative edges for the disturbed
tdITS=nrow(sigdisITS)#total edges for the disturbed
nrITS=sum(sigrefITS$rho < 0) #total number of positive edges for the disturbed
trITS=nrow(sigrefITS)#total edges for the reference

prop.test(c(ndITS,nrITS),c(tdITS,trITS), correct=FALSE)
