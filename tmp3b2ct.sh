#!/bin/sh
rm tmp3.ct; rm tmp3.tmp; rm tmp3.ok
touch tmp3.ct
i="1"
while [ $i -lt 1001 ]
do
j=$i
# sed -n 'j,j p' < tmp3.b | tr -d '\012' > tmp3.tmp
head -n$i tmp3.b | tail -n1 | tr -d '\012' > tmp3.tmp
cat tmp3.seq tmp3.tmp end > tmp3.ok
cat tmp3.ok | b2ct >> tmp3.ct
i=$[$i+1]
# i='expr $i+1'
done
exit


