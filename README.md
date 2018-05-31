# Genebank_files
V1
IM

# Several scripts that allow the annotation of prokaryotic genomes, custom modification and finally submission on NCBI databases.

# 1. Annotation of prokaryotic genomes
at local:

    for i in $(ls *.fna); do echo $i;  ~/software/prokka-1.12/prokka/bin/prokka --outdir Annotation_$(echo $i | cut -d'_' -f3) --genus $(echo $i | cut -d'_' -f1) --species $(echo $i | cut -d'_' -f2) --strain $(echo $i | cut -d'_' -f3) --locustag Csp_$(echo $i | cut -d'_' -f3) --prefix $(echo $i | cut -d'_' -f1,2,3)_Prokka --rfam --usegenus $i  ;  done

on Vital-IT cluster:

    for i in $(ls *.fa); do echo $i; bsub -q normal -L /bin/bash -J $(echo $i | cut -d'_' -f3) -u ivan.mateusgonzalez@epfl.ch -n 8 -R "rusage[mem=2000]" -M 2000000  -N  " module add UHTS/Analysis/prokka/1.12; module add UHTS/Analysis/rnammer/1.2;  module add UHTS/Analysis/LMAT/1.2.6; module add SequenceAnalysis/HMM-Profile/hmmer/3.1b2; module add SequenceAnalysis/SequenceAlignment/tbl2asn/25.3; prokka --outdir Annotation_$(echo $i | cut -d'_' -f3) --genus $(echo $i | cut -d'_' -f1) --species $(echo $i | cut -d'_' -f2) --strain $(echo $i | cut -d'_' -f3)  --cpus 8 --locustag Sp_$(echo $i | cut -d'_' -f3) --prefix $(echo $i | cut -d'_' -f1,2,3)_Prokka --rfam --usegenus $i " ;  done


# 2. Modification of gbk files

Identification of the locus tag of the gene homologs on a related strain.

 - Make protein db with the new anotations
     
       a=0;for i in $(ls Annotation*/*.faa); do echo $(echo $i | cut -d'/' -f2|cut -d'_' -f1) ;makeblastdb -dbtype prot -in $i -parse_seqids -out db_Prokka_protein/$(echo $i | cut -d'/' -f2|cut -d'_' -f1)_db ; done
      
 - Blastp REFERENCE proteins TO new assemblies proteins.
     
       a=0;for i in $(ls *.faa); do echo $(echo $i | cut -d'_' -f1,2,3) ;blastp -db db_Prokka_protein/$(echo $i | cut -d'_' -f1,2,3)_db -outfmt 6 -evalue 1e-8 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out db_Prokka_protein/blast_$(echo $i | cut -d'_' -f1,2,3).xml -query REF_*/*.faa ; done

 - Create GENENAMES files per strains

       for i in $(ls *_Prokka.faa); do cat $i | grep '>' | sed 's/ /\t/' | sed 's/>//' | cut -f1 > db_Prokka_protein/$(echo $i | cut -d'_' -f3).GeneNames ; done


Inclusion of the homolog locus tag in the gbk file

 - In R use of script: 
      
     TemplateCreation_Homolog_to_referenceStrain_V2.R . This create a template of genome-wide homologs of each gene to a reference strain.
      
 - Modify original .gbk file in text editor (i.e. sublime). 
      
      - add a new space to replace info
                
                1. replace: /locus_tag="Ab_(\w+)"\n  for /locus_tag="Ab_\1"\n\t\t\t\t\t\t\t/note="Ab_\1"\n
                2. save with different name (i.e. ANNO IM)
                  
      - Final modification of .gbk file. Create modifications comands in a file to be applied in bash.
           
                1. Open in sublime the LOCUSomology file created by the Rscript. (i.e LocusHomology_VC-Sa5Y.txt )
                2. replace file by the respective annotation
                        replace :  
                                    #         ^             with    cat 
                                    #         \tVC          WITH     | sed 's/\\/note="VC
                                    #         \tin          with    "/\\/note="in
                                    #         \t            with    "/' > 
                                    #         ' > NA"/'     WITH    \\/note="NA"/'
                3. Add end of file : cat the last temp file and remove the temporary files temp*  
                                    # cat temp3610 > VC-A1552_Annotated.gbk 
                                    #  rm temp*

                4. Copy paste on bash
                    ex: cat Vibrio_cholerae_A1552-BlokeschLab_v4.gbk | sed 's/\/note="NC002505_VC0002"/\/note="NC002505_VC0002"\n\/gene_syn="mioC"/' > temp1
                    cat temp1 | sed 's/\/note="NC002505_VC0003"/\/note="NC002505_VC0003"\n\/gene_syn="trmE_mnmE"/' > temp2
                    cat temp2 > Vibrio_cholerae_A1552-BlokeschLab_v5.gbk


# 3. Create gene homology plot between two strains.

Blast each strain to the other

   - Make protein db with the new anotations
     
            a=0;for i in $(ls *.faa); do echo $(echo $i | cut -d'_' -f1) ;makeblastdb -dbtype prot -in $i -parse_seqids -out db_protein/$(echo $i |cut -d'_' -f1)_db ; done

   - Blastp REFERENCE proteins TO new assemblies proteins. Change reference

            a=0;for i in $(ls *.faa); do echo $(echo $i | cut -d'_' -f1) ;blastp -db db_protein/$(echo $i | cut -d'_' -f1)_db -outfmt 6 -evalue 1e-8 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out db_protein/Blast_$(echo $i | cut -d'_' -f1)_IN_N16961Blokesch_.xml -query N16961Blokesch_Prokka.faa ; done

   - Extract Gene names 

            for i in $(ls *.faa); do cat $i | grep ">" | cut -d' ' -f1 | sed 's/>//' > $(echo $i | cut -d'_' -f1).GeneNames ; done
           
   - Extracting Gene Name of 1st gene on 2nd chromosome. Look at gbk file and extract the name of the first gene of chromosome 1 and chromosome 2. 
    This will be the first command in r:
    
            Chr2A1552_1stGene<-"VC-A1552Blokesch_02796"
            Chr2N16961_1stGene<-"VC-N16961Blokesch_02736"
            Chr2Sa5Y_1stGene<-"VC-Sa5YBlokech_02729"

In R use of Script:
        
  Template_HomologyPlot_between_prot_2_strainsV2.R This script will : filter the blast matchs, and produce a heatmap of gene presence/absence between two strains, for two different chromosomes.
    


# 3. Submission to NCBI 

[Guidelines from Genome NCBI submission]
    - gene names: no repeated names ex: "tfoX1_1 (tfoY)"  ->  gene_name ="tfoX_1" gene_syn= "TfoY"  
    - locus tag: no "-"  ex: Sa5Y-chr1_0001 -> Sa5Y_0001      Sa5y-chr2_0001 -> Sa5Y_A0001
    - misc_RNA feature do not exist. change to one category from http://www.insdc.org/rna_vocab.html. or misc_feature
    - note: change NC_0024.._001 to NC0023.._001

-  use gbf2tbl.pl to create annotated file to put in ncbi 
                
                perl gbf2tbl.pl Vibrio_cholerae_N16961-BlokeschLab.gbk

- Modify template.sbt with info about authors of genome submission

- use tbl2asn to create .sqn file ready to upload and submit

        tbl2asn -a s -V v -c b -Z discrep -i Vibrio_cholerae_N16961-BlokeschLab.fsa -f Vibrio_cholerae_N16961-BlokeschLab.tbl -t template.sbt -X C

-Upload .sqn file with custom annotation on NCBI 
           
    



