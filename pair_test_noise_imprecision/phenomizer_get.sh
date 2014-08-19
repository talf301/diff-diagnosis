for f in $1/*_hpo.txt; do
    wget -O $2/`basename $f`.results "http://compbio.charite.de/phenomizer/phenomizer/PhenomizerServiceURI?mobilequery=true&terms=`cat $f`&numres=20"
    sleep 5
done
