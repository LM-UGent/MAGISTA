#!/bin/bash
<<<ENV>>>
contigs=$1
falist=$2
reference=$3
prefix=$4

tfile=$(realpath $falist)
tdir=${tfile%/*}/

path2GenDisCal=<<<GDC_PATH>>>
#AMBER must be installed and on the PATH
magista_path=<<<MAGISTA_PATH>>>

echo "working in $tdir"

if [ ! -d "${tdir}${prefix}_MQassignment" ]; then
 echo "running MQ"
 ${magista_path}/make_MetaQUAST_bins.sh "$contigs" "$reference" "${tdir}${prefix}_MQassignment/"
fi

binline=$(while read line; do printf "$line "; done < "$falist")
refline=$(while read line; do printf "$line "; done < "${tdir}${prefix}_MQassignment/fa.list")

echo "converting to BB format"
$magista_path/dependencies/convert_fasta_bins_to_biobox_format.py -o "${tdir}/${prefix}.bins.profile.tmp" $binline
$magista_path/dependencies/convert_fasta_bins_to_biobox_format.py -o "${tdir}/${prefix}.reference.profile.tmp" $refline

sed -i '/^$/d' "$contigs"

echo "adding length column"
$magista_path/dependencies/add_length_column.py -g "${tdir}/${prefix}.reference.profile.tmp" -f "$contigs" > "${tdir}/${prefix}.refs.profile"
$magista_path/dependencies/add_length_column.py -g "${tdir}/${prefix}.bins.profile.tmp" -f "$contigs" > "${tdir}/${prefix}.bins.profile"

echo "running AMBER"
amber.py -g "${tdir}/${prefix}.refs.profile" -o "$tdir/${prefix}_AMBER" "${tdir}/${prefix}.bins.profile"

echo "Parsing results"
mkdir "$tdir/${prefix}_AMBER/bin_contents"

echo "AMBER LOG" > "$tdir/${prefix}_AMBER/script-log.txt"
while read line
do
  echo "$line" >> "$tdir/${prefix}_AMBER/script-log.txt"
  bname=$(basename $line)
  tname=${bname%.*}
  cat $line | grep '>' | sed "s/^>//g" > "$tdir/${prefix}_AMBER/bin_contents/${tname}.tmp"
  while read line2
  do
    cat "${tdir}/${prefix}.refs.profile" | grep $line2 | uniq
  done < "$tdir/${prefix}_AMBER/bin_contents/${tname}.tmp" > "$tdir/${prefix}_AMBER/bin_contents/${tname}.tsv"
  rm "$tdir/${prefix}_AMBER/bin_contents/${tname}.tmp"
  cat "$tdir/${prefix}_AMBER/bin_contents/${tname}.tsv" | cut -f 2 | sort | uniq > "$tdir/${prefix}_AMBER/bin_contents/${tname}.genomes.list"
  while read line2
  do
    i=0
    while read fraglen
    do
      i=$(expr $i + $fraglen)
    done < <(cat "$tdir/${prefix}_AMBER/bin_contents/${tname}.tsv" | grep "$line2" | cut -f 3)
    bname2=$(basename $line2)
    tname2=${bname2%.*}
    refgenome=$(cat "$reference" | grep "$tname2")
    taxonomy=$($path2GenDisCal "$ref_signatures" -s "$refgenome" -w 0.33 -r -a -t "$ref_taxonomy" -q all | head -2 | tail -1 | cut -d , -f 2)
    echo "$tname2,$taxonomy,$i"
  done < "$tdir/${prefix}_AMBER/bin_contents/${tname}.genomes.list" > "$tdir/${prefix}_AMBER/bin_contents/${tname}.bps.csv"
  i=$(nucops fastasummary $line -t -q)
  descriptor=Pure
  goodbps=$(cat "$tdir/${prefix}_AMBER/bin_contents/${tname}.bps.csv" | cut -d , -f 3 | sort -n | tail -1)
  while read line2
  do
    curbps=$(echo "$line2" | cut -d , -f 3)
    relbps=$(expr $curbps \* 100 / $goodbps )
    if [ $relbps -gt 10 ]
    then
      echo "$line2" >> "$tdir/${prefix}_AMBER/bin_contents/${tname}.sigfracs.csv"
    fi
    if [ $curbps -ne $goodbps ]
    then
      if [ $descriptor == Pure  ]
      then
        descriptor=Contaminated
      fi
      if [ $relbps -gt 10 ]
      then
        descriptor=Mix
      fi
    fi
  done < "$tdir/${prefix}_AMBER/bin_contents/${tname}.bps.csv"
  Gcount=$(cat "$tdir/${prefix}_AMBER/bin_contents/${tname}.sigfracs.csv" | cut -d , -f 2 | cut -d " " -f 1 | sort | uniq | wc -l)
  Scount=$(cat "$tdir/${prefix}_AMBER/bin_contents/${tname}.sigfracs.csv" | cut -d , -f 2 | cut -d " " -f 2 | sort | uniq | wc -l)
  Tcount=$(cat "$tdir/${prefix}_AMBER/bin_contents/${tname}.sigfracs.csv" | cut -d , -f 2 | cut -d " " -f 3 | sort | uniq | wc -l)
  echo "total,G${Gcount}-S${Scount}-T${Tcount}-${descriptor},$i" >> "$tdir/${prefix}_AMBER/bin_contents/${tname}.bps.csv"
done < "$falist"

tmppwd=$(pwd)
cd "$tdir/${prefix}_AMBER/bin_contents/"
for f in *.*.sigfracs.csv
do
  while read line
  do
    echo "$f,$line" | cut -d , -f 1,2,4
  done < $f
done > bin_genomes.csv
cd "$tmppwd"

mpb_path=$(realpath "${tdir}/${prefix}_AMBER/genome/${prefix}.bins.profile/metrics_per_bin.tsv")

while read line 
do
  bestref=$(echo "$line" | cut -f 3)
  bestgenome=$(basename $bestref)
  genomefile=$(cat "$reference" | grep "${bestgenome%.*}")
  genomesize=$(nucops fastasummary -t "$genomefile" -q)
  binfile=$(echo "$line" | cut -f 2)
  binfile_bn=$(basename $binfile)
  prefix1=$(echo "$line" | cut -f 2,6,7 | sed "s/\\t/,/g")
  D2M=$($path2GenDisCal $genomefile $binfile -q all | tail -1 | cut -d , -f 4)
  binclass=$(cat "${tdir}/${prefix}_AMBER/bin_contents/${binfile_bn%.*}.bps.csv" | tail -1 | cut -d , -f 2)
done < <( tail -n +2 "$mpb_path") > "${tdir}/${prefix}_AMBER/best_stats.tmp"

for i in 1
do
 header="$(head -1 "$mpb_path" | cut -f 2,6,7 | sed "s/\\t/,/g"),Genome_size,Best_match,distance_to_match,bin_class"
 echo $header
 cat "${tdir}/${prefix}_AMBER/best_stats.tmp"
done > "${tdir}/${prefix}_AMBER/best_stats.csv"

python ${magista_path}/AMBER_best_to_summary.py "${tdir}/${prefix}_AMBER/best_stats.csv" > "${tdir}/${prefix}_amber_summary.csv"

rm  "${tdir}/${prefix}_AMBER/best_stats.tmp"


source deactivate


