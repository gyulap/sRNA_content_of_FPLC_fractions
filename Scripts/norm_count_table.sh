#! /usr/bin/zsh

mydir=$PWD
ShortStackout="./sRNA-seq/ShortStack"
bamfile="${ShortStackout}/merged_alignments_fornorm_w_unmapped.bam"

cd './sRNA-seq/ShortStack_results'

#Extracting all the unique sequences from the ShortStack alignment file considering only the mapped reads and sorting them by total abundance.

samtools fasta -F4 $bamfile |\
awk '
  NR%2==0{a[$0]++}
  END{print "Sequence";
      PROCINFO["sorted_in"] = "@val_num_desc";
      for (i in a)
        {print i}
     }' > 'NR_sequences.txt'

#Counting every unique sequences per samples.

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
  done < 'norm_factors.txt'

#Making a new sequence count table for every sample that matches the sequence content and order in the total unique sequence table.
#In case of a missing sequence, the count is set to 0.

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

#Pasting all the samples together.

paste 'NR_sequences.txt' Leaf_input*NR2.txt Leaf_HMW*NR2.txt Leaf_LMW*NR2.txt Leaf_unbound*NR2.txt Flower_input*NR2.txt Flower_HMW*NR2.txt Flower_LMW*NR2.txt Flower_unbound*NR2.txt | pigz -c > 'norm_count_table_sorted.txt.gz'

rm -f *NR*

cd $mydir
