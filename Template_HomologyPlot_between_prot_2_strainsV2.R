################################################################################################
########################################################################################
###############
############### Make homology plot between strains with pheatmap
############### Blast all proteins of a strain to the other.
########## 12.04 2018

########################################################################################
# Libraries
#install.packages('genoPlotR')
library('genoPlotR')
library(RColorBrewer)
library(pheatmap)

###############################
# USe protein .faa file 
# 2. Make protein db with the new anotations
#a=0;for i in $(ls *.faa); do echo $(echo $i | cut -d'_' -f1) ;makeblastdb -dbtype prot -in $i -parse_seqids -out db_protein/$(echo $i |cut -d'_' -f1)_db ; done

# 3. Blastp REFERENCE proteins TO new assemblies proteins. Change reference
#a=0;for i in $(ls *.faa); do echo $(echo $i | cut -d'_' -f1) ;blastp -db db_protein/$(echo $i | cut -d'_' -f1)_db -outfmt 6 -evalue 1e-8 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out db_protein/Blast_$(echo $i | cut -d'_' -f1)_IN_N16961Blokesch_.xml -query N16961Blokesch_Prokka.faa ; done

# 4. Extract Gene names 
# for i in $(ls *.faa); do cat $i | grep ">" | cut -d' ' -f1 | sed 's/>//' > $(echo $i | cut -d'_' -f1).GeneNames ; done



###########################################################################################
########### First command of Rscript after library import.


# 5. Extracting Gene Name of 1st gene on 2nd chromosome
# look at gbk file and extract the name
Chr2A1552_1stGene<-"VC-A1552Blokesch_02796"
Chr2N16961_1stGene<-"VC-N16961Blokesch_02736"
Chr2Sa5Y_1stGene<-"VC-Sa5YBlokech_02729"


#############################################
setwd('~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/4_BlastHomologyStrains/db_protein/')

# Import xml files
# xml best blast hit per gene.

filesToProcess <- dir(pattern = "*\\.xml$")  #files to pro# if event 3 merged



####################################################
############  N16961-REF A1552        ##############
####################################################

filesToProcess<-filesToProcess[grep("IN_N16961Blo",filesToProcess)]
filesToProcess<-filesToProcess[grep("Blast_A1552Blo",filesToProcess)]

listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                           error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                V4="NA",V5="NA",V6="NA",
                                                                                                V7=0,V8=0,V9="NA",
                                                                                                V10="NA",V11="NA",V12="NA")))
names(listOfFiles)<-gsub(".xml","",gsub("blast","",filesToProcess))


# Create column header

colnam<-strsplit(names(listOfFiles), "_" )
colnam<-lapply(colnam, function (x) x[c(4,2)])
colnam<-lapply(colnam, function (x) gsub(".xml","",x))

for (i in 1:length(filesToProcess)) {
  NAME<-c(unlist(colnam[[i]]),"Identity","Length","Mnn","Humm","Start_strain1","End_strain1","Start_strain2","End_strain2","Evalue","Bitscore")
  
  colnames(listOfFiles[[i]]) <- NAME }

listOfFiles[[1]]
names(listOfFiles)

###############################################################################################3
########### FILTERING THE DATA
# filter 100% identity
listOfFiles<-lapply(listOfFiles, function(x) x[x$Identity>79.9,] )

# length no higher than 20% 

listOfFiles<-lapply(listOfFiles, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) >0.8 ,] )

listOfFiles<-lapply(listOfFiles, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) <1.2 ,] )


###############################################################################################

# not reciprocal hits
Comp1<-listOfFiles

# Merge names and blast
labREFN16961<-read.table("../N16961Blokesch.GeneNames",h=F)
colnames(labREFN16961)<-"N16961Blokesch"

merged<-merge(Comp1[[1]],labREFN16961,by="N16961Blokesch",all=T)



Comp1<-cbind.data.frame(N16961=as.character(merged$N16961Blokesch),  A1552=as.character(merged$A1552Blokesch) ,stringsAsFactors = FALSE)

summary(Comp1)
legenda<-Comp1$N16961
Comp1[!is.na(Comp1)] <- 1
Comp1[is.na(Comp1)] <- 0


rownames(Comp1)<-legenda

# error if duplicate rownames 
# correct by
legenda[grep("VC-N16961Blokesch_01244",legenda)]<-c("VC-N16961Blokesch_01244","VC-N16961Blokesch_01244b")
rownames(Comp1)<-legenda

Comp1$N16961<-as.numeric(Comp1$N16961)
Comp1$A1552<-as.numeric(Comp1$A1552)

rownames(Comp1[Comp1$A1552==0,])

# Order Comp1 
Comp1  <- Comp1[order(row.names(Comp1 )),] 
# Order Comp1 vector

#

LargeChromosomeN16961<-Comp1[1:match(Chr2N16961_1stGene,rownames(Comp1)),]
SmallChromosomeN16961<-Comp1[match(Chr2N16961_1stGene,rownames(Comp1)):length(rownames(Comp1)),]

rownames(SmallChromosomeN16961)

pdf("~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/4_BlastHomologyStrains/N16961REF_A1552_LargeChromosome.pdf", height=9, width = 4)
pheatmap(LargeChromosomeN16961,cluster_rows = F,cluster_cols = F,color = c("white", "black"),cellwidth =90 ,cellheight = 0.2,fontsize_col = 20,legend = F, show_colnames = T, show_rownames = F)
dev.off()

pdf("~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/4_BlastHomologyStrains/N16961REF_A1552_SmallChromosome.pdf", height=9, width = 4)
pheatmap(SmallChromosomeN16961,cluster_rows = F,cluster_cols = F,color = c("white", "black"),cellwidth =90 ,cellheight = 0.2,fontsize_col = 20,legend = F, show_colnames = T, show_rownames = F)
dev.off()
#####################################################
