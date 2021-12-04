#!/bin/bash

contigs=$1
ref_list=$2
outdir=$3

nucops_path=<<<NUCOPS_PATH>>>

mkdir "$outdir"
rm -rf "$outdir/reference"
mkdir "$outdir/reference"

while read line
do
 >&2 echo $line
 cp -t "$outdir/reference" $line
done < "$ref_list"

>&2 echo "Running: metaquast.py -o $outdir/mqres -R $outdir/reference -t 16 $contigs --fragmented --unique-mapping &>/dev/null"
metaquast.py -o "$outdir/mqres" -R "$outdir/reference" -t 16 "$contigs" --fragmented --unique-mapping &>/dev/null

bcontigs=$(basename $contigs)
contigs_noext=${bcontigs%.*}

cat "$outdir/mqres/combined_reference/contigs_reports/all_alignments_${contigs_noext}.tsv" | grep -vP "\tunaligned" | grep -vP "\tmis_unaligned" > "$outdir/mqres/alns.tsv"
cat "$outdir/mqres/alns.tsv" | sed -n '3~2p' > "$outdir/mqres/row_tails.tmp"
cat "$outdir/mqres/alns.tsv" | sed -n '2~2p' > "$outdir/mqres/row_heads.tmp"
paste "$outdir/mqres/row_heads.tmp" "$outdir/mqres/row_tails.tmp" > "$outdir/mqres/rows.tsv"
rm "$outdir/mqres/row_tails.tmp"
rm "$outdir/mqres/row_heads.tmp"
cat "$outdir/mqres/rows.tsv" | grep CONTIG | cut -f 5,6 > "$outdir/contig_assignment.tsv"

cat "$outdir/mqres/combined_reference/contigs_reports/all_alignments_${contigs_noext}.tsv" | grep -P "\tunaligned" | grep "CONTIG" | cut -f 2 >  "$outdir/mqres/unalns.tsv"
cat "$outdir/mqres/combined_reference/contigs_reports/all_alignments_${contigs_noext}.tsv" | grep -P "\tmis_unaligned" | grep "CONTIG" | cut -f 2 >> "$outdir/mqres/unalns.tsv"
cat "$outdir/mqres/unalns.tsv" | sed "s/^/unknown_genome\t/g" >> "$outdir/contig_assignment.tsv"


while read line
do
 name=$(basename $line)
 base=${name%.f*}
 cat "$outdir/contig_assignment.tsv" | grep "$base" | cut -f 2 | uniq > "$outdir/$base.list"
 ${nucops_path} select -i "$contigs" -e "$outdir/$base.list" -o "$outdir/$base.fa" -q progress 2> tmp.err
 if [ $(cat tmp.err | grep "\[WARN\]" | wc -l) -gt 0 ]
 then
  >&2 echo $line
  >&2 cat tmp.err
 fi
 rm tmp.err
 rm "$outdir/$base.list"
 realpath "$outdir/$base.fa"
done < "$ref_list" > "$outdir/fa.list"

rm -rf "$outdir/reference"
#rm "$outdir/contig_assignment"


