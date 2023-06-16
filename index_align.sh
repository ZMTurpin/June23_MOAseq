#!/bin/bash

#ZMT_June_2023
#This Script prepares bwa indices and aligns MOA-seq reads. Input files are concatenated, adapter-trimmed gzipped-fastq format.

#build bwa indices for each genome if they don't already exist (B73, Oh7b, Ky21, Sorghum Bicolor)
        #B73 exists already
                B73index=/psi/basslab/zturpin/genomes/Zea_mays/B73_v5/bwa_index/Zm-B73-REFERENCE-NAM-5.0.fa

        #Oh7b
if [ -e /psi/basslab/zturpin/genomes/Zea_mays/Oh7b/Zm-Oh7B-REFERENCE-NAM-1.0.bwt ]
then
	echo "index exists"
else
	cd /psi/basslab/zturpin/genomes/Zea_mays/Oh7b
	bwa index -p Zm-Oh7B-REFERENCE-NAM-1.0 Zm-Oh7B-REFERENCE-NAM-1.0.fa
fi

oh7bIndex=/psi/basslab/zturpin/genomes/Zea_mays/Oh7b/Zm-Oh7B-REFERENCE-NAM-1.0

	#Ky21
if [ -e /psi/basslab/zturpin/genomes/Zea_mays/Ky21/Zm-Ky21-REFERENCE-NAM-1.0.bwt ]
then
        echo "index exists"
else
	cd /psi/basslab/zturpin/genomes/Zea_mays/Ky21
	bwa index -p Zm-Ky21-REFERENCE-NAM-1.0 Zm-Ky21-REFERENCE-NAM-1.0.fa
fi

Ky21Index=/psi/basslab/zturpin/genomes/Zea_mays/Ky21/Zm-Ky21-REFERENCE-NAM-1.0

	#SorBi
if [ -e /psi/basslab/zturpin/genomes/Zea_mays/SorBi_NCBIv3.55/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.bwt ]
then
        echo "index exists"
else
        cd /psi/basslab/zturpin/genomes/Zea_mays/SorBi_NCBIv3.55
	bwa index -p Sorghum_bicolor.Sorghum_bicolor_NCBIv3 Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa
fi


SbIndex=/psi/basslab/zturpin/genomes/Zea_mays/SorBi_NCBIv3.55/Sorghum_bicolor.Sorghum_bicolor_NCBIv3
cd /psi/bassdata/zturpin/Flooding_MOAseq

#align reads to respective genomes with bwa-mem

B73R1=($(ls cat_trim/*B73*R1*trim*))
B73R2=($(ls cat_trim/*B73*R2*trim*))

B73out=($(echo ${B73R1[@]} | awk '{gsub("cat_trim","sams") ; gsub("_R1_trim","") ; gsub(".fastq.gz",".sam") ; print $0}'))

Oh7bR1=($(ls cat_trim/*Oh7b*R1*trim*))
Oh7bR2=($(ls cat_trim/*Oh7b*R2*trim*))

Oh7bout=($(echo ${Oh7bR1[@]} | awk '{gsub("cat_trim","sams") ; gsub("_R1_trim","") ; gsub(".fastq.gz",".sam") ; print $0}'))

Ky21R1=($(ls cat_trim/*Ky21*R1*trim*))
Ky21R2=($(ls cat_trim/*Ky21*R2*trim*))

Ky21out=($(echo ${Ky21R1[@]} | awk '{gsub("cat_trim","sams") ; gsub("_R1_trim","") ; gsub(".fastq.gz",".sam") ; print $0}'))

SbR1=($(ls cat_trim/*Sorghum*R1*trim*))
SbR2=($(ls cat_trim/*Sorghum*R2*trim*))

Sbout=($(echo ${SbR1[@]} | awk '{gsub("cat_trim","sams") ; gsub("_R1_trim","") ; gsub(".fastq.gz",".sam") ; print $0}'))

#mkdir sams if does not already exist
if [ -e sams ]
then
	echo "directory exists"
else
	mkdir sams
fi

#Align reads to respective genome indices
for ((a=0; a<${#B73R1[@]}; a++))
do
        bwa mem -t 8 $B73index ${B73R1[$a]} ${B73R2[$a]} >> ${B73out[$a]}
        bwa mem -t 8 $oh7bIndex ${Oh7bR1[$a]} ${Oh7bR2[$a]} >> ${Oh7bout[$a]}
        bwa mem -t 8 $Ky21Index ${Ky21R1[$a]} ${Ky21R2[$a]} >> ${Ky21out[$a]}
        bwa mem -t 8 $SbIndex ${SbR1[$a]} ${SbR2[$a]} >> ${Sbout[$a]}
done
