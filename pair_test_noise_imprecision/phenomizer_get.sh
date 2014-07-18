for f in First* Second*; do
    wget -O results/"$f".results "http://compbio.charite.de/phenomizer/phenomizer/PhenomizerServiceURI?mobilequery=true&terms=`cat $f`&numres=20"
    sleep 5
done
