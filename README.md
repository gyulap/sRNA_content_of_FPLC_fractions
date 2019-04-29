# sRNA content of FPLC fractions
Here we provide the scripts that reproduce the results and figures of the paper Dalmadi et al. 2019. These scripts are not meant to be reusable. Most of the figures and tables underwent post-processing in Excel, Powerpoint, or other GUI softwares, therefore the results will not exactly be the same as they appear in the paper.

## Requirements
* Linux operation system with zsh shell
* sra-toolkit 2.9.2 [link](https://github.com/ncbi/sra-tools)
* FastQC 0.11.3 [link](https://github.com/s-andrews/FastQC)
* cutadapt 1.9.1 [link](https://github.com/marcelm/cutadapt)
* pigz 2.4 [link](https://github.com/madler/pigz)
* ShortStack 3.4 [link](https://github.com/MikeAxtell/ShortStack)
* SAMtools 1.3.1 [link](https://github.com/samtools/samtools)
* BEDtools 2.26.0 [link](https://github.com/arq5x/bedtools2)
* R 3.4
* R libraries 'pheatmap', 'RColorBrewer', and 'ggfortify'

## Usage
Download this directory and run the `sRNA_analysis.sh` script which downloads sRNA-seq reads from the SRA database and processes them. The results are placed in the 'sRNA-seq' directory.
