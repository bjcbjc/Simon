
infile=$1
line=$2

split -l $line -d $infile $infile".p"
header=`egrep ^#[^#] ${infile}`
for f in `ls ${infile}.p*`; do
    if [ $f != ${infile}".p00" ]; then
	sed -i "1i ${header}" $f
    fi
done
