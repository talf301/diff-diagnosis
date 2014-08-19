for s in $1/*_hpo.txt.results; do
	name=`basename $s _hpo.txt.results`
	omim=`grep $name /dupa-filer/talf/diff-diagnosis/diagnosed_patients/diagnoses.txt | cut -f 3`
	#echo $name
	echo `grep -n $omim $1/"$name"_hpo.txt.results | cut -f 1 -d ':'`
	#echo `cat $1/$s.results | head -n 1 | cut -f 1`
done | sort -n | uniq -c
