# June23_MOAseq
ZMT-written scripts catenate, trim, align, and quantify PE50 Flooding MOA-seq data
INPUT: demultiplexed fastq files
OUTPUTs: 
  catenated, gzipped fastq files
  Pre-adapter trimmed Fastqc reports
  Adapter-trimmed, gzipped fastq files
  Post-adapter trimmed Fastqc reports
  bwa genome indices for genome assemblies of interest
  Primary sequence alignments (.bam)
  
  
Scripts:
  merge.sh - This shell script combines reads from 2 lands of a novaseq S2 flowcell, performs pre-trimming fastqc, trims NEBNext adapters with Cutadapt, and runs post-trimming fastqc.
  
  index_align.sh - This shell script prepares bwa indices and aligns MOA-seq reads.

compareReplicates.sh - This shell script prepares 1kb-binned Pearson correlation matrices of genomic coverages for all control and flood replicates within each genotype. Resulting replicate-replicate correlations are then visualized as sample-sample distance heatmaps AND PCA plots with attached Scree plots (as .pdfs) 
  
