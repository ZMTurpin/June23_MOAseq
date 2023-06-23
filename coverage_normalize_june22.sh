#!/bin/bash


#Compute scaling factor for each replicate. Output to text file. (number of raw reads per replicate, in millions)

files=($(ls cat_trim/*R1*trim*gz))
sample=($(echo ${files[@]} | awk '{gsub("cat_trim/","") ; gsub("_R1_trim.fastq.gz","") ; print $0}'))

#if the Scalefactors file already exists, this block will not execute.
#Else, count lines/4 (number of raw reads) in each file and scale to millions (/1,000,000)
if [ -e Scalefactors.txt ]
then
	echo "Scalefactors.txt already exists"
else
	for ((a=0; a<${#files[@]}; a++))
	do
		echo "$(zcat ${files[$a]} | wc -l)/4000000" | bc -l >> Scalefactors.txt
	done
fi

#input bed files (per genotype grouped)

B73=($(ls B73*q30.bed))
Ky21=($(ls Ky21*q30.bed))
Oh7b=($(ls Oh7b*q30.bed))
Sb=($(ls Sorghum*q30.bed))


#genome directory /psi/basslab/zturpin/genomes/Zea_mays/

B73Genome=/psi/basslab/zturpin/genomes/Zea_mays/B73_v5/Zm-B73-REFERENCE-NAM-5.0_chr.fa

#prepare chr only references for each genome
#prepare chromsizes files
if [ -e /psi/basslab/zturpin/genomes/Zea_mays/Ky21/Zm-Ky21-REFERENCE-NAM-1.0_chromsizes_chr.txt ]
then
	echo "Ky21 chromsizes exists"
else
samtools faidx -o /psi/basslab/zturpin/genomes/Zea_mays/Ky21/Zm-Ky21-REFERENCE-NAM-1.0_chr.fa /psi/basslab/zturpin/genomes/Zea_mays/Ky21/Zm-Ky21-REFERENCE-NAM-1.0.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
awk '{OFS="\t" ; print $1,$2}' /psi/basslab/zturpin/genomes/Zea_mays/Ky21/Zm-Ky21-REFERENCE-NAM-1.0_chr.fai >> /psi/basslab/zturpin/genomes/Zea_mays/Ky21/Zm-Ky21-REFERENCE-NAM-1.0_chromsizes_chr.txt
fi

if [ -e /psi/basslab/zturpin/genomes/Zea_mays/Oh7b/Zm-Oh7B-REFERENCE-NAM-1.0_chromsizes_chr.tx ]
then
	echo "Oh7b chromsizes exists"
else
samtools faidx -o /psi/basslab/zturpin/genomes/Zea_mays/Oh7b/Zm-Oh7B-REFERENCE-NAM-1.0_chr.fa /psi/basslab/zturpin/genomes/Zea_mays/Oh7b/Zm-Oh7B-REFERENCE-NAM-1.0.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
samtools faidx /psi/basslab/zturpin/genomes/Zea_mays/Oh7b/Zm-Oh7B-REFERENCE-NAM-1.0_chr.fa
awk '{OFS="\t" ; print $1,$2}' /psi/basslab/zturpin/genomes/Zea_mays/Oh7b/Zm-Oh7B-REFERENCE-NAM-1.0_chr.fai >> /psi/basslab/zturpin/genomes/Zea_mays/Oh7b/Zm-Oh7B-REFERENCE-NAM-1.0_chromsizes_chr.txt
fi

if [ -e /psi/basslab/zturpin/genomes/Zea_mays/SorBi_NCBIv3.55/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel_chromsizes_chr.txt ]
then
	echo "Sorghum chromsizes exists"
else
samtools faidx -o /psi/basslab/zturpin/genomes/Zea_mays/SorBi_NCBIv3.55/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel_chr.fa /psi/basslab/zturpin/genomes/Zea_mays/SorBi_NCBIv3.55/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa 1 2 3 4 5 6 7 8 9 10
samtools faidx /psi/basslab/zturpin/genomes/Zea_mays/SorBi_NCBIv3.55/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel_chr.fa
awk '{OFS="\t" ; print $1,$2}' /psi/basslab/zturpin/genomes/Zea_mays/SorBi_NCBIv3.55/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel_chr.fai >> /psi/basslab/zturpin/genomes/Zea_mays/SorBi_NCBIv3.55/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel_chromsizes_chr.txt
fi

#variables to chromsizes paths
B73_chromsizes=/psi/basslab/zturpin/genomes/Zea_mays/B73_v5/B73v5_chromsizes_chr
Ky21_chromsizes=/psi/basslab/zturpin/genomes/Zea_mays/Ky21/Zm-Ky21-REFERENCE-NAM-1.0_chromsizes_chr.txt
Oh7b_chromsizes=/psi/basslab/zturpin/genomes/Zea_mays/Oh7b/Zm-Oh7B-REFERENCE-NAM-1.0_chromsizes_chr.txt
SbGenome_chromsizes=/psi/basslab/zturpin/genomes/Zea_mays/SorBi_NCBIv3.55/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel_chromsizes_chr.txt

#generate genome window files (20bp) if they don't already exist

if [ -e /psi/basslab/zturpin/genomes/Zea_mays/B73_v5/B73v5_chr_w20.bed ]
then
	echo "B73w20 exists"
else
bedtools makewindows -g $B73_chromsizes -w 20 >> /psi/basslab/zturpin/genomes/Zea_mays/B73_v5/B73v5_chr_w20.bed
fi

if [ -e /psi/basslab/zturpin/genomes/Zea_mays/Ky21/Zm-Ky21-REFERENCE-NAM-1.0_chr_w20.bed ]
then
        echo "Ky21w20 exists"
else
bedtools makewindows -g $Ky21_chromsizes -w 20 >> /psi/basslab/zturpin/genomes/Zea_mays/Ky21/Zm-Ky21-REFERENCE-NAM-1.0_chr_w20.bed
fi

if [ -e /psi/basslab/zturpin/genomes/Zea_mays/Oh7b/Zm-Oh7B-REFERENCE-NAM-1.0_chr_w20.bed ]
then
        echo "Oh7bw20 exists"
else
bedtools makewindows -g $Oh7b_chromsizes -w 20 >> /psi/basslab/zturpin/genomes/Zea_mays/Oh7b/Zm-Oh7B-REFERENCE-NAM-1.0_chr_w20.bed
fi

if [ -e /psi/basslab/zturpin/genomes/Zea_mays/SorBi_NCBIv3.55/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel_chr_w20.bed ]
then
        echo "Sorghumw20 exists"
else
bedtools makewindows -g $SbGenome_chromsizes -w 20 >> /psi/basslab/zturpin/genomes/Zea_mays/SorBi_NCBIv3.55/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel_chr_w20.bed
fi

#count reads on w20 intervals and normalize counts using Scalefactors.txt

#usage
#bedtools intersect -c -a reads.bed -b windows.bed >> reads_w20.bg
B73=($(ls B73*q30.bed))
b73Out=($(echo ${B73[@]} | awk '{gsub(".bed","") ; print $0}'))
Ky21=($(ls Ky21*q30.bed))
KyOut=($(echo ${Ky21[@]} | awk '{gsub(".bed","") ; print $0}'))
Oh7b=($(ls Oh7b*q30.bed))
OhOut=($(echo ${Oh7b[@]} | awk '{gsub(".bed","") ; print $0}'))
Sb=($(ls Sorghum*q30.bed))
SbOut=($(echo ${Sb[@]} | awk '{gsub(".bed","") ; print $0}'))


for ((a=0 ; a<${#B73[@]} ; a++))
do
	bedtools intersect -c -a /psi/basslab/zturpin/genomes/Zea_mays/B73_v5/B73v5_chr_w20.bed -b ${B73[$a]} >> ${b73Out[$a]}_w20.bg 
done

for ((b=0 ; b<${#Ky21[@]} ; b++))
do
	bedtools intersect -c -a /psi/basslab/zturpin/genomes/Zea_mays/Ky21/Zm-Ky21-REFERENCE-NAM-1.0_chr_w20.bed -b ${Ky21[$b]} >> ${KyOut[$b]}_w20.bg
done

for ((c=0 ; c<${#Oh7b[@]} ; c++))
do
	bedtools intersect -c -a /psi/basslab/zturpin/genomes/Zea_mays/Oh7b/Zm-Oh7B-REFERENCE-NAM-1.0_chr_w20.bed -b ${Oh7b[$c]} >> ${OhOut[$c]}_w20.bg
done

for ((d=0 ; d<${#Sb[@]} ; d++))
do
	bedtools intersect -c -a /psi/basslab/zturpin/genomes/Zea_mays/SorBi_NCBIv3.55/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel_chr_w20.bed -b ${Sb[$d]} >> ${SbOut[$d]}_w20.bg
done

#normalize read counts (divide each count by the number of raw reads in the sample (in millions) ex. 100M reads = /100)

bgs=($(ls *w20.bg))
rpm=($(echo ${bgs[@]} | awk '{gsub(".bg","_rpm.bg") ; print $0}'))
factors=($(awk '{ORS=" " ; print $0}' Scalefactors.txt))

for ((e=0 ; e<${#factors[@]} ; e++))
do
	awk -v f=${factors[$e]} '{$5=$4/f ; OFS="\t" ; print $1,$2,$3,$5}' ${bgs[$e]} >> ${rpm[$e]}
done


#make bigwig files for browser inspection
#rpm=($(ls *_rpm.bg))
#File conversion requires sorted bed input
sorted=($(echo ${rpm[@]} | awk '{gsub(".bg","_sorted.bg") ; print $0}'))
for ((h=0; h<${#rpm[@]}; h++))
do
	sort -k1,1 -k2,2n ${rpm[$h]} >> ${sorted[$h]}
done	

	#Convert human readable bedgraphs to binary bigwig files
	#Usage bedGraphToBigWig input.bg chromsizes.txt output.bw

#B73_chromsizes=/psi/basslab/zturpin/genomes/Zea_mays/B73_v5/B73v5_chromsizes_chr
#Ky21_chromsizes=/psi/basslab/zturpin/genomes/Zea_mays/Ky21/Zm-Ky21-REFERENCE-NAM-1.0_chromsizes_chr.txt
#Oh7b_chromsizes=/psi/basslab/zturpin/genomes/Zea_mays/Oh7b/Zm-Oh7B-REFERENCE-NAM-1.0_chromsizes_chr.txt
#SbGenome_chromsizes=/psi/basslab/zturpin/genomes/Zea_mays/SorBi_NCBIv3.55/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel_chromsizes_chr.txt

#prepare file I/O arrays for each genotype
B73bg=($(ls B73*sorted.bg))
bout=($(echo ${B73bg[@]} | awk '{gsub(".bg",".bw") ; print $0}'))
Ky21bg=($(ls Ky21*sorted.bg))
kout=($(echo ${Ky21bg[@]} | awk '{gsub(".bg",".bw") ; print $0}'))
Oh7bbg=($(ls Oh*sorted.bg))
Oout=($(echo ${Oh7bbg[@]} | awk '{gsub(".bg",".bw") ; print $0}'))
SBbg=($(ls Sorgh*sorted.bg))
Sout=($(echo ${SBbg[@]} | awk '{gsub(".bg",".bw") ; print $0}'))

for ((i=0; i<${#B73bg[@]}; i++))
do
	bedGraphToBigWig ${B73bg[$i]} $B73_chromsizes ${bout[$i]}
        bedGraphToBigWig ${Ky21bg[$i]} $Ky21_chromsizes ${kout[$i]}
        bedGraphToBigWig ${Oh7bbg[$i]} $Oh7b_chromsizes ${Oout[$i]}
        bedGraphToBigWig ${SBbg[$i]} $SbGenome_chromsizes ${Sout[$i]}
done
