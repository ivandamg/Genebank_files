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

- Inclusion of the locus tag of the gene homologs on a related strain.

      -Make protein db with the new anotations
a=0;for i in $(ls Annotation*/*.faa); do echo $(echo $i | cut -d'/' -f2|cut -d'_' -f1) ;makeblastdb -dbtype prot -in $i -parse_seqids -out db_Prokka_protein/$(echo $i | cut -d'/' -f2|cut -d'_' -f1)_db ; done
      
      -Blastp REFERENCE proteins TO new assemblies proteins.
a=0;for i in $(ls *.faa); do echo $(echo $i | cut -d'_' -f1,2,3) ;blastp -db db_Prokka_protein/$(echo $i | cut -d'_' -f1,2,3)_db -outfmt 6 -evalue 1e-8 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out db_Prokka_protein/blast_$(echo $i | cut -d'_' -f1,2,3).xml -query REF_*/*.faa ; done

      -Create GENENAMES files per strains
for i in $(ls *_Prokka.faa); do cat $i | grep '>' | sed 's/ /\t/' | sed 's/>//' | cut -f1 > db_Prokka_protein/$(echo $i | cut -d'_' -f3).GeneNames ; done


# 3. Submission to NCBI 



