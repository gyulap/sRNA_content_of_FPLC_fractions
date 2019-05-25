#!/usr/bin/zsh

wd=$PWD
outdir="${wd}/sRNA-seq"

miRBase="${wd}/Auxiliary_files/miRBase_mature_sequences.fasta"
tasiRNA="${wd}/Auxiliary_files/tasiRNA_sequences.fasta"
TAIR10="${wd}/Auxiliary_files/TAIR10_genes_intergenic_merged.fasta"
pattern="${outdir}/Top_5000_sequences.fasta"

cd $outdir

patman -D $miRBase -P $pattern -e 0 -s > 'miRBase.patman'
patman -D $tasiRNA -P $pattern -e 0 -s > 'tasiRNA.patman'
patman -D $TAIR10 -P $pattern -e 1 -s > 'TAIR10.patman'

annotation(){
            d=$1
            p=$2
            awk -v name="$d" 'BEGIN{FS=OFS="\t"}
                 NR==FNR{a[$2]=a[$2]";"$1"("$3"-"$4","$5", mm: "$6")"; next
                   }
                   {if ($1 == "Sequence")
                       {print $0, name}
                    else if ($1 in a)
                       {print $0, a[$1]}
                    else {print $0, "No hit"}
                   }' "${d}.patman" $p |\
            sed 's/\t[;]/\t/' > "${p%.txt}_${d}.txt"
            }

annotation 'miRBase' 'Top_5000_sequences.txt'
annotation 'tasiRNA' 'Top_5000_sequences_miRBase.txt'
annotation 'TAIR10' 'Top_5000_sequences_miRBase_tasiRNA.txt'

rm -f 'Top_5000_sequences_miRBase.txt' 'Top_5000_sequences_miRBase_tasiRNA.txt' *.patman

cd $wd
