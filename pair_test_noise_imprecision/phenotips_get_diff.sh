for f in $1/*.new; do
    echo $f
    arg="https://playground.phenotips.org/bin/get/PhenoTips/DiffDiagnosisService?format=html&q=`sed 's/:/%3A/g;s/,/%20/g' $f` --data `sed 's/:/%3A/g;s/H/symptom=H/g;s/,/\&/g' $f`"
    echo $arg
    curl -k -o $2/"`basename $f`".results $arg
    sleep 1
done
