#!/bin/bash

#ZMT_June_2023
#This script combines reads from 2 lanes of a novaseq S2 flowcell, performs pre-trimming fastqc, trims NEBNext adapters with Cutadapt, and runs post-trimming fastqc.
#cat lane 1 and 2 reads 

#mkdir cat_trim

R1L1=($(ls 06092023_Zach/*L001_R1*))
R1L2=($(ls 06092023_Zach/*L002_R1*))
R2L1=($(ls 06092023_Zach/*L001_R2*))
R2L2=($(ls 06092023_Zach/*L002_R2*))

#outnamesR1=($(ls 06092023_Zach/*L001_R1* | awk '{gsub("06092023_Zach/","") ; gsub ("L001_R1_001","R1") ; print $0}'))

#outnamesR2=($(ls 06092023_Zach/*L001_R2* | awk '{gsub("06092023_Zach/","") ; gsub ("L001_R2_001","R2") ; print $0}'))

#for ((a=0; a<${#R1L1[@]}; a++)) 
#	do
#	cat ${R1L1[$a]} ${R2L2[$a]} >> cat_trim/${outnamesR1[$a]}
#	cat ${R2L1[$a]} ${R2L2[$a]} >> cat_trim/${outnamesR2[$a]}
#done

#trim contaminating adapter sequences with cutadapt
        #usage cutadapt -a <adapterFwd> -A<adapterRev> -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq

        #create I/O filename arrays
#mkdir trim

R1=($(ls cat_trim/*R1*))
R2=($(ls cat_trim/*R2*))

R1out=($(echo ${R1[@]} | awk '{gsub("R1","R1_trim") ; print $0 }'))

R2out=($(echo ${R2[@]} | awk '{gsub("R2","R2_trim") ; print $0 }'))

#AdaptorRead 1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
#AdaptorRead 2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

for (( a=0 ; a<${#R1[@]} ; a++))
        do
        cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ${R1out[$a]} -p ${R2out[$a]} ${R1[$a]} ${R2[$a]}
        done

#run fastqc on trimmed reads

trim=($(for a in $(ls cat_trim/*trim*) ; do printf $a" " ; done))

mkdir ClipFASTQC
fastqc -t 6 -o ClipFASTQC ${trim[@]}
