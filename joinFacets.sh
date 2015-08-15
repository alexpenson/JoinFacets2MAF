#!/bin/bash

if [ "$#" -ne "2" ]; then
	echo USAGE: joinFacets.sh MAF FACETS_DIR
	exit
fi

FACETS=/home/socci/Code/Pipelines/FACETS/FACETS.app

MAF=$1
FDIR=$2

TDIR=_scratch

mkdir -p $TDIR

$FACETS/out2tbl.py $FDIR/*out >$TDIR/facets.out

cat $FDIR/*_cncf.txt | egrep -v "^ID" \
	| awk '{print $2,$3-1,$4" "$0}' \
	| tr '\t' '|' | tr ' ' '\t' >$TDIR/facets.bed

cat $MAF | egrep -v "^(#|Hugo_Symbol)" \
	| awk '{print $5,$6-1,$7" "$0}' \
	| tr '\t' '|' | tr ' ' '\t' >$TDIR/maf.bed

rm -f $TDIR/join_maf_facets.bed missingCNCFData.txt

samples=$(cat $MAF | egrep -v "^(#|Hugo_Symbol)" | cut -f16 | sed 's/.*_//' | sort | uniq)
for sample in $samples; do
	echo $sample
	fgrep $sample $TDIR/maf.bed >$TDIR/sMaf.bed
	fgrep $sample $TDIR/facets.bed >$TDIR/sFacets.bed
	if [ -s "$TDIR/sFacets.bed" ]; then
		bedtools intersect -loj \
			-a $TDIR/sMaf.bed \
			-b $TDIR/sFacets.bed \
			>>$TDIR/join_maf_facets
	else
		echo $sample >>missingCNCFData.txt
	fi
done

./joinFacetsInfo.py $TDIR/join_maf_facets $TDIR/facets.out $MAF \
	>$(basename $MAF | sed 's/.maf//')__JoinFacets.maf


