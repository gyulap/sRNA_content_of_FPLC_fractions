#!/usr/bin/zsh

genomefile='./Auxiliary_files/TAIR10_nuclear.txt'
outdir='./sRNA-seq/ShortStack_results'
bamfile="${outdir}/merged_alignments_fornorm_w_unmapped.bam"
rgfile="${outdir}/norm_factors.txt"

if [[ ! -d "${outdir}/Genome_browser_tracks" ]]; then
  mkdir "${outdir}/Genome_browser_tracks"
fi

while read line
  do
    rg=$(echo $line | cut -f 1)
    echo "Processing $rg ..."
    normfactor=$(echo $line | cut -f 2)
    fraction=$(echo $rg | awk 'BEGIN{FS="_"}{print $2}')

# Sets the track colors by fractions. Color codes are RGB (Red, Green, Blue).

    case $fraction in
      "input")
        color='51,153,255'
      ;;
      "HMW")
        color='51,255,51'
      ;;
      "LMW")
        color='255,153,51'
      ;;
      "unbound")
        color='255,51,153'
      ;;
    esac

# Creates the tracks for every readlength, positive and negative strands separately and then merges them into one track file.
    for rl in {21,24}
      do
        trackname="$(echo ${rg%_processed_20_25} | sed 's/_/ /g') ${rl}nt"
        trackline="track type=bedGraph name=\"${trackname}\" visibility=full color=${color} graphType=bar viewLimits=-200.0:200.0"
        plus=$(bedtools genomecov -bg -strand + -ibam <(samtools view -h -F4 $bamfile | awk -v rg="$rg" -v rl="$rl" 'BEGIN{FS=OFS="\t"}{if ($1 ~ /^@/ || ($0 ~ rg && length($10) == rl)) {print $0}}' | samtools view -bu;) -g $genomefile -scale $normfactor;)
        minus=$(bedtools genomecov -bg -strand - -ibam <(samtools view -h -F4 $bamfile | awk -v rg="$rg" -v rl="$rl" 'BEGIN{FS=OFS="\t"}{if ($1 ~ /^@/ || ($0 ~ rg && length($10) == rl)) {print $0}}' | samtools view -bu;) -g $genomefile -scale $normfactor | awk 'BEGIN{FS=OFS="\t"}{$4=-$4; print $0}';)
        cat <(echo $plus;) <(echo $minus;) | bedtools sort | sed "1i$trackline" > "${outdir}/Genome_browser_tracks/${rg%_processed_20_25}_${rl}nt_norm.bedgraph" &&
        echo "$rg done."
      done
  done < $rgfile
