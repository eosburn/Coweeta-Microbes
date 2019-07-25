###Network Analysis###
###Network Generation for Fungi and Bacteria###
###Coweeta Microbe project###
###Modified from https://github.com/ryanjw/co-occurrence###
###Modified by Steve McBride; mcbridsg@vt.edu
###last Modified 07/23/2019###
###generates network graphs##
###Data from Ernest Osburn Forest  Disturbance study at Coweeta###
###Data collected in Summer 2018###
###Files Needed to begin analysis: 
###Bacteria correlation data "Reference_OTU_correlation_matrix_min10_16Sp01rho7.csv", "Disturbed_OTU_correlation_matrix_min10_16Sp01rho7.csv"
###Fungi correlation data "Reference_OTU_correlation_matrix_min10_ITSp01.csv", "Disturbed_OTU_correlation_matrix_min10_ITSp01.csv"
###Files Generated and imported after processing outside of R
###Files for adding colors to vertecies: 16sfeatures.csv, ITSdata.csv
###Files for setting size of verticies refotu_ITS_min10samples_size.csv, disotu_ITS_min10samples_size.csv,
###refotu_16S_min10samples_size.csv, refotu_16S_min10samples_size.csv
###Files for adding color to edges posneg.csv
###All files generated and used in the following code are available https://github.com/eosburn/Coweeta-Microbes



#####Create environment and load necessary packages#####

library(igraph)


### Set your working directory. This should be a folder that you have access to that
### all of your files will be stored in.
setwd('your/workingdirectory/')

######Import feature data####

##16S features
features_16S=read.csv(file="16sfeatures.csv",header=TRUE)
features_16S$OTU2=paste(features_16S$XXX,features_16S$OTU, sep="")

refsizefeatures_16S=read.csv(file="refotu_16S_min10samples_size.csv",header=TRUE)
refsizefeatures_16S$Sample2=paste(refsizefeatures_16S$XXX,refsizefeatures_16S$Sample, sep="")

dissizefeatures_16S=read.csv(file="disotu_16S_min10samples_size.csv",header=TRUE)
dissizefeatures_16S$Sample2=paste(dissizefeatures_16S$XXX,dissizefeatures_16S$Sample, sep="")

##ITS features
features_ITS=read.csv(file="ITSfeatures.csv",header=TRUE)
features_ITS$OTU2=paste(features_ITS$XXX,features_ITS$OTU, sep="")


refsizefeatures_ITS=read.csv(file="refotu_ITS_min10samples_size.csv",header=TRUE)
refsizefeatures_ITS$Sample2=paste(refsizefeatures_ITS$XXX,refsizefeatures_ITS$Sample, sep="")

dissizefeatures_ITS=read.csv(file="disotu_ITS_min10samples_size.csv",header=TRUE)
dissizefeatures_ITS$Sample2=paste(dissizefeatures_ITS$XXX,dissizefeatures_ITS$Sample, sep="")

edgefeatures=read.csv(file="posneg.csv", header=TRUE)#Used for both 16s and ITS

#####Assigned colors to 16Sfeatures.csv and ITSfeatures.csv based on taxonomy.
#####Assigned size to each to the appropriate files by setting the minimum size to 4 then adding that
#####to the relative (abundance *100)
#####edge colors in posneg.csv are just red and blue. This is paired with the $cor data below


#####Import Correlation Data####
###Import data that will be used for creating network graphs. These will be your correlation matricies
###I modified the Disturbed_OTU_correlation_matrix_min10_16S.csv, and Reference_OTU_correlation_matrix_min10_16S.csv
###in excel to only include p.value < 0.01, and |rho| >0.7
refdata_16S=read.csv(file="Reference_OTU_correlation_matrix_min10_16Sp01rho7.csv",header=TRUE)
disdata_16S=read.csv(file="Disturbed_OTU_correlation_matrix_min10_16Sp01rho7.csv",header=TRUE)

###I modified the Disturbed_OTU_correlation_matrix_min10_ITS.csv, and Reference_OTU_correlation_matrix_min10_ITS.csv
###in excel to only include p.value < 0.01, and |rho| > 0.5
refdata_ITS=read.csv(file="Reference_OTU_correlation_matrix_min10_ITSp01.csv",header=TRUE)
disdata_ITS=read.csv(file="Disturbed_OTU_correlation_matrix_min10_ITSp01.csv",header=TRUE)

#####Subset data so you can create networks for the negative and positive correlations separately#####

###16S first  
refnegdata_16S=subset(refdata_16S,rho<0)
refnegdata_16S$cor="negative"
refposdata_16S=subset(refdata_16S,rho>0)
refposdata_16S$cor="positive"
refalldata_16S=rbind(refnegdata_16S,refposdata_16S)

disnegdata_16S=subset(disdata_16S,rho<0)
disnegdata_16S$cor="negative"
disposdata_16S=subset(disdata_16S,rho>0)
disposdata_16S$cor="positive"
disalldata_16S=rbind(disnegdata_16S,disposdata_16S)

###SubsetITS
refnegdata_ITS=subset(refdata_ITS,rho<0)
refnegdata_ITS$cor="negative"
refposdata_ITS=subset(refdata_ITS,rho>0)
refposdata_ITS$cor="positive"
refalldata_ITS=rbind(refnegdata_ITS,refposdata_ITS)

disnegdata_ITS=subset(disdata_ITS,rho<0)
disnegdata_ITS$cor="negative"
disposdata_ITS=subset(disdata_ITS,rho>0)
disposdata_ITS$cor="positive"
disalldata_ITS=rbind(disnegdata_ITS,disposdata_ITS)

######Reference 16S Graphs#####
refgraph1_16S<-simplify(graph.edgelist(as.matrix(refnegdata_16S[,c(3,4)]),directed=FALSE)) #The columns being called here "[,c(3,4)]" should be the columns with your different taxa
refgraph2_16S<-simplify(graph.edgelist(as.matrix(refposdata_16S[,c(3,4)]),directed=FALSE)) #same as above
E(refgraph2_16S)$posneg="blue" #will use for coloring edges
refsoils_graph_16S<-graph.union(refgraph1_16S,refgraph2_16S) #combines the two graphs

V(refsoils_graph_16S)$color = as.character(features_16S$color[match(V(refsoils_graph_16S)$name, features_16S$OTU2)])#Colors verticies based on 
V(refsoils_graph_16S)$size = as.character(refsizefeatures_16S$size[match(V(refsoils_graph_16S)$name, refsizefeatures_16S$Sample2)])#Adds colors to the different vertcies based on taxa
#V(refsoils_graph_16S)$Phylum = as.character(features_16S$Phylum[match(V(refsoils_graph_16S)$name, features_16S$OTU2)])#Adds taxanomic information to the different verticies
refvsize_16S=as.matrix(V(refsoils_graph_16S)$size)
refvsize_16S=as.numeric(refvsize_16S)

E(refsoils_graph_16S)$color=as.character(edgefeatures$color[match(E(refsoils_graph_16S)$posneg, edgefeatures$V1)])

#png("Reference_16S_Network.png",units="in",width=5, height=5,res=1500) ##This can be used to write to PNG if you want to create a file on your computer
set.seed(145)
par(mar=c(0,0,0,0)+.1)
plot.igraph(refsoils_graph_16S, vertex.size=refvsize_16S, vertex.label=NA,
            edge.width=c(1), margin=-.1, vertex.frame.color=NA)
#dev.off()##If you remove the # in front of png above you will need to do the same here.

######Disturbed 16S Graphs#####
disgraph1_16S<-simplify(graph.edgelist(as.matrix(disnegdata_16S[,c(3,4)]),directed=FALSE)) #The columns being called here "[,c(3,4)]" should be the columns with your different taxa
disgraph2_16S<-simplify(graph.edgelist(as.matrix(disposdata_16S[,c(3,4)]),directed=FALSE)) #same as above
E(disgraph2_16S)$posneg="blue" #will use for coloring edges
dissoils_graph_16S<-graph.union(disgraph1_16S,disgraph2_16S) #combines the two graphs

V(dissoils_graph_16S)$color = as.character(features_16S$color[match(V(dissoils_graph_16S)$name, features_16S$OTU2)])#Colors verticies based on 
V(dissoils_graph_16S)$size = as.character(dissizefeatures_16S$size[match(V(dissoils_graph_16S)$name, dissizefeatures_16S$Sample2)])#Adds colors to the different vertcies based on taxa
#V(dissoils_graph_16S)$Phylum = as.character(features_16S$Phylum[match(V(dissoils_graph_16S)$name, features_16S$OTU2)])#Adds taxanomic information to the different verticies
disvsize_16S=as.matrix(V(dissoils_graph_16S)$size)
disvsize_16S=as.numeric(disvsize_16S)

E(dissoils_graph_16S)$color=as.character(edgefeatures$color[match(E(dissoils_graph_16S)$posneg, edgefeatures$V1)])

#png("Disturbed_16S_Network.png",units="in",width=5, height=5,res=1500) ##This can be used to write to PNG if you want to create a file on your computer
set.seed(1)
par(mar=c(0,0,0,0)+.1)
plot.igraph(dissoils_graph_16S, vertex.size=disvsize_16S, vertex.label=NA,
            edge.width=c(1), margin=-.1, vertex.frame.color=NA)

#dev.off()##If you remove the # in front of png above you will need to do the same here.




#####Reference ITS Graphs#####
refgraph1_ITS<-simplify(graph.edgelist(as.matrix(refnegdata_ITS[,c(3,4)]),directed=FALSE)) #The columns being called here "[,c(3,4)]" should be the columns with your different taxa
refgraph2_ITS<-simplify(graph.edgelist(as.matrix(refposdata_ITS[,c(3,4)]),directed=FALSE)) #same as above
E(refgraph2_ITS)$posneg="blue" #will use for coloring edges
refsoils_graph_ITS<-graph.union(refgraph1_ITS,refgraph2_ITS) #combines the two graphs

V(refsoils_graph_ITS)$color = as.character(features_ITS$color[match(V(refsoils_graph_ITS)$name, features_ITS$OTU2)])#Colors verticies based on 
V(refsoils_graph_ITS)$size = as.character(refsizefeatures_ITS$size[match(V(refsoils_graph_ITS)$name, refsizefeatures_ITS$Sample2)])#Adds colors to the different vertcies based on taxa
#V(refsoils_graph_ITS)$Phylum = as.character(features_ITS$Phylum[match(V(refsoils_graph_ITS)$name, features_ITS$OTU2)])#Adds taxanomic information to the different verticies
refvsize_ITS=as.matrix(V(refsoils_graph_ITS)$size)
refvsize_ITS=as.numeric(refvsize_ITS)

E(refsoils_graph_ITS)$color=as.character(edgefeatures$color[match(E(refsoils_graph_ITS)$posneg, edgefeatures$V1)])

#png("Reference_ITS_Network.png",units="in",width=5, height=5,res=1500) ##This can be used to write to PNG if you want to create a file on your computer
set.seed(145)
par(mar=c(0,0,0,0)+.1)
plot.igraph(refsoils_graph_ITS, vertex.size=refvsize_ITS, vertex.label=NA,
            edge.width=c(1), margin=-.1, vertex.frame.color=NA)

#dev.off()##If you remove the # in front of png above you will need to do the same here.

#####Disturbed ITS Graphs######
disgraph1_ITS<-simplify(graph.edgelist(as.matrix(disnegdata_ITS[,c(3,4)]),directed=FALSE)) #The columns being called here "[,c(3,4)]" should be the columns with your different taxa
disgraph2_ITS<-simplify(graph.edgelist(as.matrix(disposdata_ITS[,c(3,4)]),directed=FALSE)) #same as above
E(disgraph2_ITS)$posneg="blue" #will use for coloring edges
dissoils_graph_ITS<-graph.union(disgraph1_ITS,disgraph2_ITS) #combines the two graphs

V(dissoils_graph_ITS)$color = as.character(features_ITS$color[match(V(dissoils_graph_ITS)$name, features_ITS$OTU2)])#Colors verticies based on 
V(dissoils_graph_ITS)$size = as.character(dissizefeatures_ITS$size[match(V(dissoils_graph_ITS)$name, dissizefeatures_ITS$Sample2)])#Adds colors to the different vertcies based on taxa
#V(dissoils_graph_ITS)$Phylum = as.character(features_ITS$Phylum[match(V(dissoils_graph_ITS)$name, features_ITS$OTU2)])#Adds taxanomic information to the different verticies
disvsize_ITS=as.matrix(V(dissoils_graph_ITS)$size)
disvsize_ITS=as.numeric(disvsize_ITS)

E(dissoils_graph_ITS)$color=as.character(edgefeatures$color[match(E(dissoils_graph_ITS)$posneg, edgefeatures$V1)])


#png("Disturbed_ITS_Network.png",units="in",width=5, height=5,res=1500) ##This can be used to write to PNG if you want to create a file on your computer
set.seed(1)
par(mar=c(0,0,0,0)+.1)
plot.igraph(dissoils_graph_ITS, vertex.size=disvsize_ITS, vertex.label=NA,
            edge.width=c(1), margin=-.1, vertex.frame.color=NA)

#dev.off()##If you remove the # in front of png above you will need to do the same here.