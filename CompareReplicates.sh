#!/bin/bash

#This script used multiBamSummary to construct correlation matrices across all replicates and conditions for a given genotype. Correlation values are produced for genome coverage per kilobase of reference sequence. 

#Compare Replicates within genotype groups

#usage multiBamSummary -bins --binSize 1000 -p 6 -b [bamfiles.bam] -o [results.npz]

		#set bin size to 1000 instead of default 10,000 
		#set option ignoreDuplicates
		#set number of processors = 6

#index all markdup files
mkdup=($(ls *mkdup.bam))
#samtools index -M ${mkdup[@]}

B73=($(ls B73*_chr_mkdup.bam))
b73Labels=($(echo ${B73[@]} | awk '{gsub("_sort_chr_chr_mkdup.bam","") ; print $0}'))
multiBamSummary bins --ignoreDuplicates --minMappingQuality 30 --binSize 1000 -p 6 -b ${B73[@]} --labels ${b73Labels[@]} -o B73_multiBamSummary.npz

Ky21=($(ls Ky21*_chr_mkdup.bam))
Ky21Labels=($(echo ${Ky21[@]} | awk '{gsub("_sort_chr_chr_mkdup.bam","") ; print $0}'))
multiBamSummary bins --ignoreDuplicates --minMappingQuality 30 --binSize 1000 -p 6 -b ${Ky21[@]} --labels ${Ky21Labels[@]} -o Ky21_multiBamSummary.npz

Oh7b=($(ls Oh7b*_chr_mkdup.bam))
Oh7bLabels=($(echo ${Oh7b[@]} | awk '{gsub("_sort_chr_chr_mkdup.bam","") ; print $0}'))
multiBamSummary bins --ignoreDuplicates --minMappingQuality 30 --binSize 1000 -p 6 -b ${Oh7b[@]} --labels ${Oh7bLabels[@]} -o Oh7b_multiBamSummary.npz


SorBi=($(ls Sorghum*_chr_mkdup.bam))
SbLabels=($(echo ${SorBi[@]} | awk '{gsub("_sort_chr_chr_mkdup.bam","") ; print $0}'))
multiBamSummary bins --ignoreDuplicates --minMappingQuality 30 --binSize 1000 -p 6 -b ${SorBi[@]} --labels ${SbLabels[@]} -o Sorghum_multiBamSummary.npz


#visualize Correlation matrix as a coverage heatmap
	#usage
		#plotCorrelation --corData <file.npz> --corMethod spearman --whatToPlot heatmap --plotFile 1Kb_genomeCov_spearmanCor_heatmap.pdf --labels ${names[@]} --plotTitle "1Kb genomeCoverage spearman Correlation of RNA-seq" --plotNumbers

#B73
plotCorrelation --corData B73_multiBamSummary.npz --corMethod spearman --whatToPlot heatmap --plotFile B73_MOA_1kb_genomeCov_spearmanCor_heatmap.pdf --labels ${b73Labels[@]} --plotTitle "1Kb genomeCoverage spearman Correlation of B73 MOA-seq" --plotNumbers
#Oh7b

plotCorrelation --corData Oh7b_multiBamSummary.npz --corMethod spearman --whatToPlot heatmap --plotFile Oh7b_MOA_1kb_genomeCov_spearmanCor_heatmap.pdf --labels ${Oh7bLabels[@]} --plotTitle "1Kb genomeCoverage spearman Correlation of Oh7b MOA-seq" --plotNumbers


#Ky21

plotCorrelation --corData Ky21_multiBamSummary.npz --corMethod spearman --whatToPlot heatmap --plotFile Ky21_MOA_1kb_genomeCov_spearmanCor_heatmap.pdf --labels ${Ky21Labels[@]} --plotTitle "1Kb genomeCoverage spearman Correlation of Ky21 MOA-seq" --plotNumbers


#Sorghum

plotCorrelation --corData Sorghum_multiBamSummary.npz --corMethod spearman --whatToPlot heatmap --plotFile Sorghum_MOA_1kb_genomeCov_spearmanCor_heatmap.pdf --labels ${SbLabels[@]} --plotTitle "1Kb genomeCoverage spearman Correlation of Sorghum MOA-seq" --plotNumbers


#Plot PCA

colors=(blue blue blue goldenrod goldenrod goldenrod)
markers=('s' 'o' '<' 's' 'o' '<')

	#usage 
		#plotPCA --corData <file.npz> --labels ${names[@]} --plotTitle "1Kb genomeCoverage RNA-seq PCA" --plotFile 1Kb_genomeCov_PCA.pdf --colors blue yellow


#B73

plotPCA --corData B73_multiBamSummary.npz --labels ${b73Labels[@]} --plotTitle "1Kb genomeCoverage B73 MOA-seq PCA" --plotFile B73_MOAseq_1Kb_genomeCov_PCA.pdf --colors ${colors[@]} --markers ${markers[@]}


#Oh7b

plotPCA --corData Oh7b_multiBamSummary.npz --labels ${Oh7bLabels[@]} --plotTitle "1Kb genomeCoverage Oh7b MOA-seq PCA" --plotFile Oh7b_MOAseq_1Kb_genomeCov_PCA.pdf --colors ${colors[@]} --markers ${markers[@]}

#Ky21

plotPCA --corData Ky21_multiBamSummary.npz --labels ${Ky21Labels[@]} --plotTitle "1Kb genomeCoverage Ky21 MOA-seq PCA" --plotFile Ky21_MOAseq_1Kb_genomeCov_PCA.pdf --colors ${colors[@]} --markers ${markers[@]}


#Sorghum

plotPCA --corData Sorghum_multiBamSummary.npz --labels ${SbLabels[@]} --plotTitle "1Kb genomeCoverage Sorghum MOA-seq PCA" --plotFile Sorghum_MOAseq_1Kb_genomeCov_PCA.pdf --colors ${colors[@]} --markers ${markers[@]}
