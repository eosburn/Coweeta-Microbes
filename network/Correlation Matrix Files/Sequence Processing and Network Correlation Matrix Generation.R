###Network Analysis###
###Sequence Processesing and Correlation Matrix Generation###
###Coweeta Microbe project###
###Modified from https://github.com/ryanjw/co-occurrence###
###Modified by Steve McBride; mcbridsg@vt.edu
###last Modified 07/22/2019###
###generates Correlation Matrix for data###
###Data from Ernest Osburn Forest  Disturbance study at Coweeta###
###Data collected in Summer 2018###
###Files Needed to begin analysis: metadata = 16s_metadata.txt, sequence data - 
###ITS=FD-ITs-table-with-taxonomy-final.biom, 16S=FD-16s-table-with-taxonomy-final.biom###
###Files Generated and imported after processingoutside of R
###refotu_ITS.csv, disotu_ITS.csv, refotu_16S.csv, disotu_16S.csv
###refotu_ITS_min10samples.csv, disotu_ITS_min10samples.csv, refotu_16S_min10samples.csv, disotu_16S_min10samples.csv
###Reference_OTU_correlation_matrix_min10_ITS.csv, Disturbed_OTU_correlation_matrix_min10_ITS.csv
###Reference_OTU_correlation_matrix_min10_16S.csv, Disturbed_OTU_correlation_matrix_min10_16S.csv###
###All files generated and used in the following code are available https://github.com/eosburn/Coweeta-Microbes



#####Create environment and load necessary packages#####

library(dplyr)
library(phyloseq)
library(tidyr)

### Set your working directory. This should be a folder that you have access to that
### all of your files will be stored in.
setwd('your/workingdirectory/')

#####Import Data and format for analysis#####
##these files need to be in your working directory for this code to work.
meta2 <- import_qiime_sample_data("16s_metadata.txt") ###imports metadata, The ITSmetadat file in this GitHub repository is only necessary for demultiplexing raw sequences, either 16s_metadata.txt or ITS_metadata.txt can be used for network analysis
ernITS=import_biom('FD-ITs-table-with-taxonomy-final.biom') ###imports ITS sequence data, if using your own data be sure to change file name
ern16s=import_biom('FD-16s-table-with-taxonomy-final.biom') ###imports 16S sequence data, if using your own data be sure to change file name

##create phyloseq objects
OTU2_ITS <- merge_phyloseq(ernITS, meta2) 
OTU2_16S <- merge_phyloseq(ernITS, meta2)

#####subset data, transform for cooccurence analysis #####
set.seed(TRUE)
OTU3_ITS=rarefy_even_depth(OTU2_ITS, sample.size = min(sample_sums(OTU2_ITS)),
                       rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
OTU3_16S=rarefy_even_depth(OTU2_16S, sample.size = min(sample_sums(OTU2_16S)),
                           rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

##adds metadata to OTUs
OTU4_ITS=psmelt(OTU3_ITS)
OTU4_16S=psmelt(OTU3_16S)

#Select from the dplyr package allows me to select the columns Abundance, treatment, sample, and OTU from the OTU4 object
OTUselect_ITS=select(OTU4_ITS,Abundance,Treatment,Sample,OTU) 
OTUselect_16S=select(OTU4_16S,Abundance,Treatment,Sample,OTU)

#Spread from the tidyr package allows me to take my long data (taxa as rows) and convert it to wide (taxa as columns)
OTU5_ITS=spread(data = OTUselect_ITS,key=OTU,value=Abundance)
OTU5_16S=spread(data = OTUselect_16S,key=OTU,value=Abundance)

##subset by treatment ###Datasets for processing - networks will be created for each treatment
refOTU_ITS=subset(OTU5_ITS, Treatment=="Reference")
disOTU_ITS=subset(OTU5_ITS, Treatment=="Disturbed")

refOTU_16S=subset(OTU5_16S, Treatment=="Reference")
disOTU_16S=subset(OTU5_16S, Treatment=="Disturbed")

###Although this could probably be automated in R, I created .CSV files to remove taxa that appeared
###in a minimum number of samples. I keep the code disable until I create these files so they aren't 
###accidentally overwritte. remove the "#" from in front of the code below to create these .csv files

#write.csv(refOTU_ITS,file="refotu_ITS.csv")
#write.csv(disOTU_ITS,file="disotu_ITS.csv")

#write.csv(refOTU_16S,file="refotu_16S.csv")
#write.csv(disOTU_16S,file="disotu_16S.csv")

###These are the new dataframes that only contain taxa that appeared in a minimum of 10 samples.
###This allows us to ensure the robustness of our correlations per Shi et a. 2016 doi:10.1111/ele.12630.
###I did this by transposing the rows to columns. Then, I used the "countif" command in excel to sum
###the number of samples each taxa appeared in.Then, I removed all taxa that did not appear in at least 10 samples
refotumin10_ITS=read.csv(file="refotu_ITS_min10samples.csv", header=TRUE,row.names = 1)
disotumin10_ITS=read.csv(file="disotu_ITS_min10samples.csv", header=TRUE,row.names = 1)
refotumin10_16S=read.csv(file="refotu_16S_min10samples.csv", header=TRUE,row.names = 1)
disotumin10_16S=read.csv(file="disotu_16S_min10samples.csv", header=TRUE,row.names = 1)


#####GENERATE CORRELATION MATRICIES#####
####set comm.data to the dataset that you want to calculate these metrics for
comm.data<- refotumin10_ITS # will need to repeat with each relevant dataframe disotumin10_ITS, refotumin10_16S, disotumin10_16S

trts<-as.vector((unique(comm.data$Treatment)))
results<-matrix(nrow=0,ncol=7)
options(warnings=-1)

###The for loop will need to be run for each new comm.data entry.
for(a in 1:length(trts)){
  #pull the first element from the vector of treatments
  trt.temp<-trts[a]
  #subset the dataset for those treatments
  temp<-subset(comm.data, Treatment==trt.temp)
  #in this case the community data started at column 6, so the loop for co-occurrence has to start at that point
  for(b in 3:(dim(temp)[2]-1)){
    #every species will be compared to every other species, so there has to be another loop that iterates down the rest of the columns
    for(c in (b+1):(dim(temp)[2])){
      
      #summing the abundances of species of the columns that will be compared
      species1.ab<-sum(temp[,b])
      species2.ab<-sum(temp[,c])
      #if the column is all 0's no co-occurrence will be performed
      if(species1.ab >1 & species2.ab >1){
        test<-cor.test(temp[,b],temp[,c],method="spearman",na.action=na.rm)
        rho<-test$estimate
        p.value<-test$p.value
      }
      
      if(species1.ab <=1 | species2.ab <= 1){
        rho<-0
        p.value<-1
      }	
      
      new.row<-c(trts[a],names(temp)[b],names(temp)[c],rho,p.value,species1.ab,species2.ab)
      results<-rbind(results,new.row)			
      
    }      }
  
  
  print(a/length(trts))
  
}

head(results)
results<-data.frame(data.matrix(results))
names(results)<-c("trt","taxa1","taxa2","rho","p.value","ab1","ab2")

###Write data to .csv
#write.csv(results,file="Reference_OTU_correlation_matrix_min10_ITS.csv")
#write.csv(results,file="Disturbed_OTU_correlation_matrix_min10_ITS.csv")
#write.csv(results,file="Reference_OTU_correlation_matrix_min10_16S.csv")
#write.csv(results,file="Disturbed_OTU_correlation_matrix_min10_16S.csv")

###Once these .csv files were created I deleted the column of row names e.g. "new.row.1"