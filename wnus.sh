nline=`wc -l nusse.out | tr -d '[:alpha:][= =][=.=]'`
nline=`echo "$nline/2" | bc`
nu=`awk -v nl=$nline 'NR<=nl {x+=$2}END{print x/nl}' nusse.out`
echo "Half wT: $nu"
nu=`awk '{x+=$2}END{print x/NR}' nusse.out`
echo "Full wT: $nu"
nu=`awk '{x+=$3}END{print x/NR}' nusse.out`
echo "Full T: $nu"
nu=`awk '{x+=$2}END{print x/NR}' nusse3.out`
echo "Full kindiss: $nu"
nu=`awk '{x+=$3}END{print x/NR}' nusse3.out`
echo "Full thediss: $nu"

