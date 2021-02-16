#!/usr/bin/bash
cat $1 | grep -ve "^\[" | grep -v "^---" | sed "s/[ ]+/\t/g" | sed "s/\t\(/\(/g" | sed "s/\#\t/#/g" | sed "s/Marker\tlineage/checkm.Marker.lineage/" | sed "s/#genomes/checkm.ngenomes/" |  sed "s/#markers/checkm.nmarkers/" | sed "s/#marker\tsets/checkm.nmarker.sets/" | sed "s/Bin\tId/Bin.Id/" | sed "s/Completeness/checkm.Completeness/" | sed "s/Contamination/checkm.Contamination/" | sed "s/marker\tsets/marker.sets/" | sed "s/Strain\theterogeneity/checkm.Strain.heterogeneity/" | sed "s/\t/,/g" | sed "s/^,//g" | awk 'NR<2{print;next}{print| "sort "}' | cut -d , -f 1,2,3,4,5,12,13,14
