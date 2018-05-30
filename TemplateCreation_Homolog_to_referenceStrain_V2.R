########################################################################################
###############
############### R script 
############## create template to match blast query and result. 
############## Modify locus Tag with homolog to N16961 using blast .xml (output6)
########## 12 oct 2017

########################################################################################

# Libraries
#install.packages('genoPlotR')
library('genoPlotR')
library(RColorBrewer)
library(pheatmap)
#######################################
# Set working directory

setwd('~/Documents/Noemie/NoemiePaper_VF/Strains_paper/db_Prokka_protein/')

###############################
# Preléiminary 
#IN BASH

# 1. Create PROKKA ANNOTATION OF NEW GENOMES

#for i in $(ls Vibrio_*.fasta); do echo $i ; ~/software/prokka-1.12/prokka/bin/prokka --outdir Annotation_$(echo $i | cut -d'_' -f3) --genus Vibrio --species cholerae --strain $(echo $i | cut -d'_' -f3) --locustag VC-$(echo $i | cut -d'_' -f3) --prefix $(echo $i | cut -d'_' -f3)_Prokka --rfam --usegenus $i ; done

# 2. Make protein db with the new anotations
#a=0;for i in $(ls Annotation*/*.faa); do echo $(echo $i | cut -d'/' -f2|cut -d'_' -f1) ;makeblastdb -dbtype prot -in $i -parse_seqids -out db_Prokka_protein/$(echo $i | cut -d'/' -f2|cut -d'_' -f1)_db ; done

# 3. Blastp REFERENCE proteins TO new assemblies proteins.
#a=0;for i in $(ls *.faa); do echo $(echo $i | cut -d'_' -f1,2,3) ;blastp -db db_Prokka_protein/$(echo $i | cut -d'_' -f1,2,3)_db -outfmt 6 -evalue 1e-8 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out db_Prokka_protein/blast_$(echo $i | cut -d'_' -f1,2,3).xml -query REF_*/*.faa ; done

# Import xml files (honmology info)
# xml best blast hit per gene.
# 4 CREATE GENENAMES FILES PER STRAIN
# for i in $(ls *_Prokka.faa); do cat $i | grep '>' | sed 's/ /\t/' | sed 's/>//' | cut -f1 > db_Prokka_protein/$(echo $i | cut -d'_' -f3).GeneNames ; done

###############################################################3333
filesToProcess <- dir(pattern = "*\\.xml$")  #files to pro# if event 3 merged
Homolog_info <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                            error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                 V4="NA",V5="NA",V6="NA",
                                                                                                 V7=0,V8=0,V9="NA",
                                                                                                 V10="NA",V11="NA",V12="NA")))

names(Homolog_info)<-gsub(".xml","",gsub("blast_","",filesToProcess))

# Create column header

colnam<-strsplit(names(Homolog_info), "_" )
colnam<-lapply(colnam, function (x) x[c(2,1)])
colnam<-lapply(colnam, function (x) gsub(".xml","",x))

for (i in 1:length(filesToProcess)) {
  NAME<-c(unlist(colnam[[i]]),"Identity","Length","Mnn","Humm","Start_strain1","End_strain1","Start_strain2","End_strain2","Evalue","Bitscore")
  
  colnames(Homolog_info[[i]]) <- NAME }

Homolog_info[[1]]
names(Homolog_info)

###############################################################################################3
########### FILTERING THE DATA
# filter 100% identity
Homolog_info<-lapply(Homolog_info, function(x) x[x$Identity>79.9,] )

# length no higher than 20% 

Homolog_info<-lapply(Homolog_info, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) >0.8 ,] )

Homolog_info<-lapply(Homolog_info, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) <1.2 ,] )

# select only names of homology
Homolog_info<-lapply(Homolog_info, function(x) x[,c(1,2)])

Homolog_info<-lapply(Homolog_info, function(x) x[!duplicated(x),])

##############################################################################################

###############################

# Import Gene names for file
# xml best blast hit per gene.

filesToProcess <- dir(pattern = "*\\.GeneNames$")  #files to pro# if event 3 merged
Gene_Names <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                          error= function (e) cbind.data.frame(V1="NA")))

names(Gene_Names)<-gsub(".xml","",gsub("blast","",filesToProcess))

# Create column header
head(Gene_Names[[1]])
colnam<-strsplit(names(Gene_Names), "\\." )
colnam<-lapply(colnam, function (x) x[c(1)])

for (i in 1:length(filesToProcess)) {
  NAME<-c(unlist(colnam[[i]]))
  
  colnames(Gene_Names[[i]]) <- NAME }

head(Gene_Names[[1]])



##############################################################################################
names(Gene_Names)<-gsub("\\.GeneNames","",names(Gene_Names))
names(Homolog_info)<-gsub("_N16961VCheidelbergNames","", names(Homolog_info) )

head(Homolog_info[[i]])
head(Gene_Names[[i]])

Info_homo<-list()
Merged_name_homolog<- list()
for (i in names(Gene_Names) ){
  
  Info_homo[[i]]<-merge (Gene_Names[[i]],Homolog_info[[i]],by=gsub(".GeneNames","",i), all=T)
  paste("temp",seq(1,dim(Info_homo[[i]])[1]),sep="")
  Info_homo[[i]][,2]<-gsub("^","inASM1542v1_",Info_homo[[i]][,2])
  Info_homo[[i]][,3]<-paste("temp",seq(1,dim(Info_homo[[i]])[1]),sep="")
  c("file",paste("temp",seq(1,dim(Info_homo[[i]])[1]  -1),sep=""))
  Info_homo[[i]]<-cbind.data.frame(c("file",paste("temp",seq(1,dim(Info_homo[[i]])[1]  -1),sep="")),Info_homo[[i]])
  
  write.table(Info_homo[[i]],paste("~/Documents/Noemie/NoemiePaper_VF/Strains_paper/db_Prokka_protein/NewLocusTag/LocusHomology_",gsub(".GeneNames","",i),".txt",sep="" ),sep = "\t",quote = F,row.names = F, col.names = F)
  
  Merged_name_homolog[[i]]<-merge (Gene_Names[[i]],Homolog_info[[i]],by=gsub(".GeneNames","",i), all=T)
  
  Merged_name_homolog[[i]]<-gsub("^","inASM1542v1_",Merged_name_homolog[[i]][,2])
  Merged_name_homolog[[i]]<- as.data.frame(Merged_name_homolog[[i]])
  
  
  write.table(Merged_name_homolog[[i]],paste("~/Documents/Noemie/NoemiePaper_VF/Strains_paper/db_Prokka_protein/NewLocusTag/LocusTag_",gsub(".GeneNames","",i),".txt",sep="" ),sep = "\t",quote = F,row.names = F, col.names = F)
  
  
}

