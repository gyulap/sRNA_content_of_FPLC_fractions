#! /usr/bin/zsh

outdir='./sRNA-seq'
ShortStackout="${outdir}/ShortStack_results"
bamfile="${ShortStackout}/merged_alignments_fornorm_w_unmapped.bam"
normfile="${outdir}/norm_factors.txt"

# Determining the number of cores on the computer.

p=$(egrep -c '^processor' '/proc/cpuinfo')

# Extracting all the unique sequences from the ShortStack alignment file considering only the mapped reads and sorting them by total abundance.

samtools fasta -F4 $bamfile |\
awk '
  NR%2==0{a[$0]++}
  END{print "Sequence";
      for (i in a)
        {print i}
     }' > 'NR_sequences.txt'

# Counting every unique sequence per sample.

while read line
  do
    rg=$(echo $line | cut -f 1)
    normfactor=$(echo $line | cut -f 2)
    samtools view -bu -F4 -r $rg $bamfile |\
    samtools fasta - |\
    awk -v nf="$normfactor" '
      NR==FNR{if (NR%2==0)
               {a[$0]++}
             }
      END{for (i in a)
           {print i"\t"(nf*a[i])}
         }' > "${rg%_sRNA_processed}_NR.txt"
  done < $normfile

# Making a new sequence count table for every sample that matches the sequence content and order in the total unique sequence table.
# In case of a missing sequence, the count is set to 0.

for i in *NR.txt
  do
    awk -v i="${i%_NR.txt}" '
      BEGIN{FS=OFS="\t"}
           NR==FNR{a[$1]=$2; next}
           {if ($1 == "Sequence")
              {print i}
            else if ($1 in a)
              {print a[$1]}
            else
              {print 0}
           }' $i 'NR_sequences.txt' > "${i%.txt}2.txt"
  done

# Pasting all the samples together.

paste 'NR_sequences.txt' Leaf_input*NR2.txt Leaf_HMW*NR2.txt Leaf_LMW*NR2.txt Leaf_unbound*NR2.txt Flower_input*NR2.txt Flower_HMW*NR2.txt Flower_LMW*NR2.txt Flower_unbound*NR2.txt | pigz -p $p -c > "${outdir}/norm_count_table.txt.gz"

rm -f *NR*

# Sorting the table by mean normalized abundance

zcat "${outdir}/norm_count_table.txt.gz" |\
awk 'BEGIN{FS=OFS="\t"}{sum=0; for (i=2;i<=NF;i++) {sum=sum+$i}; mean=(sum/16); print $0, mean}' |\
sort -k18,18nr | pigz -p $p -c > "${outdir}/norm_count_table_sorted.txt.gz" &&
rm -f "${outdir}/norm_count_table.txt.gz"
