#!/usr/bin/zsh

set -ueo pipefail

genomeurl='https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas'
geneurl='https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_blastsets/TAIR10_cdna_20101214_updated'
igurl='https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_blastsets/TAIR10_intergenic_20101028'
genomefile="./Auxiliary_files/TAIR10_chr_all.fas"
outdir='./sRNA-seq'
rawout="${outdir}/Raw_sequences"
trimmedout="${outdir}/Trimmed_sequences"
mappedout="${outdir}/Trimmed_mapped_sequences"
ShortStackout="${outdir}/ShortStack_results"
bamfile="${ShortStackout}/merged_alignments_filtered_w_unmapped.bam"
filter="./Auxiliary_files/filter.bed"

# Downloading the TAIR10 sequences from the TAIR site

if [[ ! -f $genomefile ]]; then
  cd './Auxiliary_files'
  wget $genomeurl
  wget $geneurl $igurl -O 'TAIR10_genes_intergenic_merged.fasta'
  cd ..
fi

# Creating the directory structure

if [[ ! -d $outdir ]]; then
  mkdir -p "${rawout}/FastQC_${rawout##*/}" "${trimmedout}/FastQC_${trimmedout##*/}" "${mappedout}/FastQC_${mappedout##*/}"
fi

# Determining the number of cores on the computer.

p=$(egrep -c '^processor' '/proc/cpuinfo')

# Determining the optimal memory for sorting the bam file.

m=$(egrep 'MemTotal' '/proc/meminfo' | awk '{if ($1 > 1048576) {printf "%iK\n", $2/8} else {print "768M"} }')

# Downloading the SRA metadata file

esearch -db 'sra' -query 'PRJNA540255' | efetch -format 'runinfo' | tr ',' '\t' > './Auxiliary_files/runinfo.txt'

while read line
  do
    Library_Name=$(echo $line | cut -f12 | sed 's/ /_/g')
    Run=$(echo $line | cut -f1)
    rawname="${rawout}/${Library_Name}_raw.fastq.gz"
    trimmedname="${trimmedout}/${Library_Name}_trimmed.fastq.gz"

    # Downloading raw sequences from the SRA database, renaming them by the sample name and placing them into the appropriate directory.
    # A quality check is performed before read processing using FastQC.

    if [[ ! -f $rawname ]]; then
      echo "Downloading $Run (${Library_Name}) from SRA"
      fasterq-dump $Run -p -e $p -m $m && pigz -p $p ${Run}.fastq > $rawname
      fastqc -t $p -o "${rawout}/FastQC_${rawout##*/}" $rawname
    fi

    # Trimming the Illumina TruSeq Small RNA 3' adapter (RA3).
    # A quality check is performed after read processing using FastQC.

    if [[ ! -f $trimmedname ]]; then
      cutadapt -j $p -a 'TGGAATTCTCGGGTGCCAAGG' -m 20 -M 25 -q 20 --max-n=0 --discard-untrimmed $rawname | pigz -p $p > $trimmedname &&
      fastqc -t $p -o "${trimmedout}/FastQC_${trimmedout##*/}" $trimmedname
    fi

  done < <(awk 'NR>1{print}' './Auxiliary_files/runinfo.txt')

# Mapping processed sequences to the Arabidopsis thaliana TAIR10 genome using ShortStack.
# After mapping, the sequences mapped to the Arabidopsis rRNA and tRNA genes are removed.

if [[ ! -f "${ShortStackout}/merged_alignments.bam" ]]; then
  ShortStack --align_only --bowtie_cores $p --sort_mem $m --bowtie_m 1000 --readfile *_trimmed.fastq.gz --outdir $ShortStackout &&
  samtools view -b -L $filter "${ShortStackout}/merged_alignments.bam" > $bamfile
  samtools view -H $bamfile | awk -F "\t" '/^@RG/{print substr($2, 4, length($2))}' > "${ShortStackout}/rg_list.txt"
fi

# The filtered sequences are extracted and a quality check is perfomed using FastQC.
# Normalisation factors are calculated using the FastQC results.

while read rg
  do
    mappedname="${mappedout}/${rg%_trimmed}_trimmed_mapped.fastq.gz"
    samtools view -bu -F4 -r $rg $bamfile | samtools fastq - > $mappedname
    fastqc -t $p -o "${mappedout}/FastQC_${mappedout##*/}" $mappedname
    fastqcdata="${mappedout}/FastQC_${mappedout##*/}/${mappedname%.fastq.gz}_fastqc/fastqc_data.txt"
    awk -v rg="$rg" 'BEGIN{FS=OFS="\t"}/"Total Sequences"/{print rg, (1000000/$2)}' $fastqcdata >> "${ShortStackout}/norm_factors.txt"
  done < "${ShortStackout}/rg_list.txt"

# Creating the mapping statistics (Table S1)

if [[ ! -f "${outdir}/Table_S1.txt" ]]; then
  ./Scripts/mapping_statistics.sh
fi

# Creating a sequence count table from the ShortStack alignment file

if [[ ! -f "${outdir}/norm_count_table.txt.gz" ]]; then
  ./Scripts/norm_count_table.sh
fi

# Getting the abundance data for miRBase mature miRNA sequences

miRBase='./Auxiliary_files/miRBase_mature_sequences.fasta'

awk 'BEGIN{RS="^>";FS="\n";OFS="\t"}NR>1{print $1, $2}' $miRBase > "${miRBase%.fasta}.txt"

awk 'BEGIN{FS=OFS="\t"}
     NR==FNR{a[$2]=$1; next}{
       if (FNR == 1) {print "Name", $0}
       else if ($1 in a) {print a[$1], $0}
     }' "./Auxiliary_files/miRBase_mature_sequences.txt" <(zcat "${outdir}/norm_count_table_sorted.txt.gz") > "${outdir}/miRBase_mature_sequences_norm_count_table.txt"

# Getting the top 5000 most abundant sequences

zcat "${outdir}/norm_count_table_sorted.txt.gz" | head -5000 > "${outdir}/Top_5000_sequences.txt"

# Principal component analysis using the top 5000 most abundant sequences

if [[ ! -f "PCA_plot.png" ]]; then
  Rscript 'PCA.R'
fi

# Annotating the top 5000 sequences

if [[ ! -f "${outdir}/Top_5000_sequences_miRBase_tasiRNA_TAIR10.txt" ]]; then
  ./Scripts/annotation.sh
fi

# Creating genome browser tracks for the 21 and 24-nt sRNAs from the ShortStack alignment file

if [[ ! -d "${outdir}/Genome_browser_tracks" ]]; then
  mkdir "${outdir}/Genome_browser_tracks"
  ./Scripts/Genome_browser_tracks.sh
fi

# Creating heatmaps for the different sRNA classes

if [[ -f './sRNA-seq/miRBase_mature_sequences_norm_count_table.txt' ]]; then
  Rscript './Scripts/miRNA_heatmaps.R miRNAs'
  Rscript './Scripts/miRNA_heatmaps.R miRNAs_5p-U'
  Rscript './Scripts/miRNA_heatmaps.R miRNAs_by_abundance'
  Rscript './Scripts/siRNA_heatmaps.R 21nt_siRNAs'
  Rscript './Scripts/siRNA_heatmaps.R 24nt_siRNAs'
fi
