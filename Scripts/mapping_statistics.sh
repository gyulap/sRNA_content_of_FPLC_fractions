#!/usr/bin/zsh

# Summary statistics

globalqc=("Unprocessed" "Passed" "Discarded" "Mapped" "Unmapped" "Uniquely mapped" "Multi-mapped (placed by density)" "Multi-mapped (placed randomly)")
auc=("Redundant" "%Redundant" "Non-Redundant" "%Non-Redundant" "Complexity")
outdir='./sRNA-seq'
ShortStackout="${outdir}/ShortStack_results"
bamfile="${ShortStackout}/merged_alignments_filtered_w_unmapped.bam"
rgfile="${ShortStackout}/rg_list.txt"
normfile="${ShortStackout}/norm_factors.txt"
rawout="${outdir}/Raw_sequences"

# Determining the number of cores on the computer.

p=$(egrep -c '^processor' '/proc/cpuinfo')

outputfile="${outdir}/Summary_statistics.txt"
    {
    printf "%s\n\n" "Summary statistics"
    printf "\t%s\t\t\t\t\t" $globalqc
    printf "\n\t" ""
    for t in $globalqc
      do
        printf "%s\t" $auc
        printf "\t" ""
      done
    ;} | cut -f1-2,8- | tee $outputfile

  while read rg
    do
        {
        rgname="${rg%_trimmed_mapped}"
        printf "%s " $(echo $rgname | sed 's/_/ /g'); printf "%s\t" ""
        total_unprocessed_R=$(( $(pigz -d -p $p -c "${rawout}/${rgname}_raw.fastq.gz" | sed -n '2~4p' | wc -l) ))
        printf "%.0f\t" $total_unprocessed_R
        awk -v unproc="$total_unprocessed_R" '
        BEGIN {FS=OFS="\t"}
        {
          passedarray[$2]++
          if ($1 == "0" || $1 == "16") {mappedcount++; mappedarray[$2]++};
          if ($1 == "4") {unmappedcount++ ; unmappedarray[$2]++};
          if ($3 ~ /U/) {uniquecount++; uniquearray[$2]++};
          if ($3 ~ /P/) {multicount++; multiarray[$2]++};
          if ($3 ~ /R/) {randomcount++; randomarray[$2]++};
        }
        END {
          {printf "%.0f\t", NR}
          {if (unproc != 0) {printf "%.4f\t", ((NR*100)/unproc)} else {printf "%s\t", "NA"}}
          {printf "%.0f\t", length(passedarray)}
          {if (unproc !=0) {printf "%.4f\t", ((length(passedarray)*100)/unproc)} else {printf "%s\t", "NA"}}
          {if (NR != 0) {printf "%.4f\t\t", (length(passedarray)/NR)} else {printf "%s\t\t", "NA"}}

          {printf "%.0f\t", unproc-NR}
          {if (unproc != 0) {printf "%.4f\t\t\t\t\t", ((unproc-NR)*100/unproc)} else {printf "%s\t\t\t\t\t", "NA"}}

          {printf "%.0f\t", mappedcount}
          {if (NR != 0) {printf "%.4f\t", ((mappedcount*100)/NR)} else {printf "%s\t", "NA"}}
          {printf "%.0f\t", length(mappedarray)}
          {if (NR != 0) {printf "%.4f\t", ((length(mappedarray)*100)/NR)} else {printf "%s\t", "NA"}}
          {if (mappedcount != 0) {printf "%.4f\t\t", (length(mappedarray)/mappedcount)} else {printf "%s\t\t", "NA"}}

          {printf "%.0f\t", unmappedcount}
          {if (NR != 0) {printf "%.4f\t", ((unmappedcount*100)/NR)} else {printf "%s\t", "NA"}}
          {printf "%.0f\t", length(unmappedarray)}
          {if (NR != 0) {printf "%.4f\t", ((length(unmappedarray)*100)/NR)} else {printf "%s\t", "NA"}}
          {if (unmappedcount != 0) {printf "%.4f\t\t", (length(unmappedarray)/unmappedcount)} else {printf "%s\t\t", "NA"}}

          {printf "%.0f\t", uniquecount}
          {if (NR != 0) {printf "%.4f\t", ((uniquecount*100)/NR)} else {printf "s\t", "NA"}}
          {printf "%.0f\t", length(uniquearray)}
          {if (NR != 0) {printf "%.4f\t", ((length(uniquearray)*100)/NR)} else {printf "%s\t", "NA"}}
          {if (uniquecount != 0) {printf "%.4f\t\t", (length(uniquearray)/uniquecount)} else {printf "%s\t\t", "NA"}}

          {printf "%.0f\t", multicount}
          {if (NR !=0) {printf "%.4f\t", ((multicount*100)/NR)} else {printf "%s\t", "NA"}}
          {printf "%.0f\t", length(multiarray)}
          {if (NR != 0) {printf "%.4f\t", ((length(multiarray)*100)/NR)} else {printf "%s\t", "NA"}}
          {if (multicount != 0) {printf "%.4f\t\t", (length(multiarray)/multicount)} else {printf "%s\t\t", "NA"}}

          {printf "%.0f\t", randomcount}
          {if (NR !=0) {printf "%.4f\t", ((randomcount*100)/NR)} else {printf "%s\t", "NA"}}
          {printf "%.0f\t", length(randomarray)}
          {if (NR != 0) {printf "%.4f\t", ((length(randomarray)*100)/NR)} else {printf "%s\t", "NA"}}
          {if (randomcount != 0) {printf "%.4f\n", (length(randomarray)/randomcount)} else {printf "%s\n", "NA"}}
        }' <(samtools view -@ $p -r $rg $bamfile | cut -f2,10,16)

        ;} | tee -a $outputfile
    done < $rgfile

{
 while read rg
   do
     rgname="${rg%_trimmed_mapped}"
     tempfile="./stat_${rgname}.txt"
     if [ ! -f $tempfile ]; then
       normfactor=$(awk -v rg="$rg" 'BEGIN{FS=OFS="\t"} {if ($1 == rg) print $2}' $normfile)
       temp_file=$(samtools view -@ $p -F4 -r $rg $bamfile | cut -f10)

       printf "%s" "Creating ${tempfile}"

       (printf "%s\t" "Length(nt)" "Redundant (raw)" "Non-Redundant (raw)" "Redundant (normalised)" "Non-Redundant (normalised)"
       printf "%s\n" "Complexity") > $tempfile

       for readlength in {20..25}
         do
           (cat $temp_file | awk -v nf="$normfactor" 'BEGIN{FS=OFS="\t"}NR>1{$4=$2*nf; $5=$3*nf; printf "%s\t%s\t%s\t%.4f\t%.4f\t%.4f\n", $1, $2, $3, $4, $5, $6}') >> $tempfile
         done
       printf "\t%s\n" "Done."
       unset temp_file
     fi
   done < $rgfile
;}

unset bamfile rgname temp_file
{
 for nuc in "Redundant" "Non-Redundant" "Complexity"
   do
     for rn in "raw" "normalised"
       do
         if [ $nuc = "Complexity" ] && [ $rn = "raw" ]; then continue; fi
         printf "\n\n\n%s\n" "$nuc Read Length Distribution ($rn)"
         printf "%s\t" "Length(nt)" $(for n in {20..25}; do echo $n; done) "" "Sum"

         while read rg
           do
             rgname="${rg%_trimmed_mapped}"
             tempfile="./stat_${rgname}.txt"
             temp_file=$(cat $tempfile)
             sumlength=$(printf "%.0f" 0)
             printf "%s\n" ""
             printf "%s " $(echo $rgname | sed 's/_/ /g')
             for readlength in {20..25}
               do
                 case $rn in
                   "raw")
                      case $nuc in
                        "Redundant")
                          R=$(echo $temp_file | awk -v rl="$readlength" '$1 ~ rl {print $2}')
                          printf "\t%.0f" $R
                          sumlength=$(( $sumlength + $R ))
                        ;;
                        "Non-Redundant")
                          NR=$(echo $temp_file | awk -v rl="$readlength" '$1 ~ rl {print $3}')
                          printf "\t%.0f" $NR
                          sumlength=$(( $sumlength + $NR ))
                        ;;
                      esac
                      ;;
                   "normalised")
                      case $nuc in
                        "Redundant")
                          Rnorm=$(echo $temp_file | awk -v rl="$readlength" '$1 ~ rl {print $4}')
                          printf "\t%.4f" $Rnorm
                          sumlength=$(( $sumlength + $Rnorm ))
                        ;;
                        "Non-Redundant")
                          NRnorm=$(echo $temp_file | awk -v rl="$readlength" '$1 ~ rl {print $5}')
                          printf "\t%.4f" $NRnorm
                          sumlength=$(( $sumlength + $NRnorm ))
                        ;;
                        "Complexity")
                          complexity=$(echo $temp_file | awk -v rl="$readlength" '$1 ~ rl {print $6}')
                          printf "\t%.4f" $complexity
                        ;;
                      esac
                      ;;
                  esac
               done
               case $rn in
                 "raw")
                   printf "\t\t%.0f" $sumlength
                 ;;
                 "normalised")
                   if [ $nuc = "Complexity" ]; then continue; else printf "\t\t%.4f" $sumlength ; fi
                 ;;
               esac
               unset temp_file
               rm -f $tempfile
             done < $rgfile
       done
   done
;} | tee -a $outputfile
