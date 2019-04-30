#!/usr/bin/zsh

set -ueo pipefail

outdir="./sRNA-seq"
rawout="${outdir}/Raw_sequences"
trimmedout="${outdir}/Trimmed_sequences"
mappedout="${outdir}/Trimmed_mapped_sequences"
ShortStackout="${outdir}/ShortStack"
genomefile="./Auxiliary_files/TAIR10_chr_nuclear.fa"
rfambed="./Auxiliary_files/TAIR10_fornorm.bed"
bamfile="${ShortStackout}/merged_alignments_fornorm_w_unmapped.bam"

#Creating the directory structure

mkdir -p "${rawout}/FastQC_${rawout##*/}" "${trimmedout}/FastQC_${trimmedout##*/}" "${mappedout}/FastQC_${mappedout##*/}"

#Determining the number of cores on the computer to set the number of threads for adapter trimming.

p=$(egrep -c '^processor' '/proc/cpuinfo')

#Determining the optimal memory for sorting the bam file.

m=$(egrep 'MemTotal' '/proc/meminfo' | awk '{if ($1 > 1048576) {printf "%iK\n", $2/8} else {print "768M"} }')

while read line
  do
    Library_Name=$(echo $line | cut -f6)
    Run=$(echo $line | cut -f9)
    rawname="${rawout}/${Library_Name}_raw.fastq.gz"
    trimmedname="${trimmedout}/${Library_Name}_trimmed.fastq.gz"

#Downloading raw sequences from the SRA database, renaming them by the sample name and placing them into the appropriate directory.
    if [[ ! -f $rawname ]]; then
      echo "Downloading $Run (${Library_Name}) from SRA"
      fastqerq-dump $Run -p -e $p -m $m && pigz -p $p ${Run}.fastq > $rawname
      fastqc -t $p -o "${rawout}/FastQC_${rawout##*/}" $rawname
    fi

#Trimming the Illumina TruSeq Small RNA 3' adapter (RA3).
#A quality check is performed before and after read processing using FastQC.
    if [[ ! -f $trimmedname ]]; then
      cutadapt -j $p -a 'TGGAATTCTCGGGTGCCAAGG' -m 20 -M 25 -q 20 --max-n=0 --discard-untrimmed $rawname | pigz -p $p > $trimmedname &&
      fastqc -t $p -o "${trimmedout}/FastQC_${trimmedout##*/}" $trimmedname
    fi

  done < <(awk 'NR>1{print}' './Auxiliary_files/SraRunTable.txt')

if [[ ! -f "${ShortStackout}/merged_alignments.bam" ]]; then
  ShortStack --align_only --bowtie_cores $p --sort_mem $m --bowtie_m 1000 --readfile *_trimmed.fastq.gz --outdir $ShortStackout &&
  samtools view -b -F3840 -L $rfambed "${ShortStackout}/merged_alignments.bam" > $bamfile
  samtools view -H $bamfile | awk -F "\t" '/^@RG/{print substr($2, 4, length($2))}' > "${ShortStackout}/rg_list.txt"
fi

while read rg
  do
    mappedname="${mappedout}/${rg%_trimmed}_trimmed_mapped.fastq.gz"
    samtools view -bu -F4 -r $rg $bamfile | samtools fastq - > $mappedname
    fastqc -t $p -o "${mappedout}/FastQC_${mappedout##*/}" $mappedname
  done < "${ShortStackout}/rg_list.txt"

#Creating the mapping statistics (Table S1)

if [[ ! -f "${outdir}/Table_S1.txt" ]]; then
  ./Scripts/mapping_statistics.sh
fi

#Creating a sequence count table from the ShortStack alignment file

if [[ ! -f "${outdir}/norm_count_table.txt" ]]; then
  ./Scripts/norm_count_table.sh
fi

#Annotating the expression table

#Creating genome browser tracks for the 21 and 24-nt sRNAs from the ShortStack alignment file

if [[ ! -d "${outdir}/Genome_browser_tracks" ]]; then
  mkdir "${outdir}/Genome_browser_tracks"
  ./Scripts/Genome_browser_tracks.sh
fi

#Creating heatmaps for the different sRNA classes

if [[ -f './sRNA-seq/miRBase_mature_sequences_norm_count_table.txt' ]]; then
  Rscript './Scripts/heatmaps.R miRNA'
  Rscript './Scripts/heatmaps.R siRNA21'
  Rscript './Scripts/heatmaps.R siRNA24'
fi
