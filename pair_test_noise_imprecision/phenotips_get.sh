#! /bin/bash

indir=$1
outdir=$2

for f in $1/*hpo.txt; do
	hps=`cat $f | sed 's/,/%20/g' | sed 's/:/%3A/g'`
	postdata=`cat $f | sed 's/,/\&symptom=/g' | sed 's/:/%3A/g'`
	wget -O $2/"`basename $f`".results "https://playground.phenotips.org/bin/get/PhenoTips/OmimPredictService?format=html&q=${hps}&reqNo=`grep -o , $f | wc -l`" --post-data "symptom=${postdata}" --no-check-certificate
  sleep 5
done
