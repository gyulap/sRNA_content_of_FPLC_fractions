#!/usr/bin/zsh

outdir='./sRNA-seq'

cd $outdir

miRBase='./Auxiliary_files/miRBase_mature_sequences.fasta'
tasiRNA='./Auxiliary_files/tasiRNA_sequences.fasta'
TAIR10='./Auxiliary_files/TAIR10_genes_intergenic_merged.fasta'
pattern='Top_5000_sequences.fasta'

patman -D $miRbase -P $pattern -e 0 -s > 'miRBase.patman'
patman -D $tasiRNA -P $pattern -e 0 -s > 'tasiRNA.patman'
patman -D $TAIR10 -P $pattern -e 1 -s > 'TAIR10.patman'

annot(){
       awk -v name="$1" 'BEGIN{FS=OFS="\t"}
       NR==FNR{a[$2]=a[$2]";
               "$1"("$3"-"$4","$5", mm: "$6")";
               next
              }
              {if ($2 == "Sequence")
                  {print $0, name}
               else if ($2 in a)
                  {print $0, a[$2]}
               else {print $0, "No hit"}
              }' "${1}.patman" $2 |\
       sed 's/\t[;]/\t/' > "${$2%.txt}_${$1}.txt"
       }

annot 'miRBase' 'Top_5000_sequences.txt'
annot 'tasiRNA' 'Top_5000_sequences_miRBase.txt'
annot 'TAIR10' 'Top_5000_sequences_miRBase_tasiRNA.txt'

rm -f 'Top_5000_sequences_miRBase.txt' 'Top_5000_sequences_miRBase_tasiRNA.txt'

cd ..
