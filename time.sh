#!/bin/sh

N=14
loop=$1
st=15	#115
en=`expr $loop + 14`	#214

### initialize ###

cat hoge > ymm_mul.h

echo "N,Optimized,Normal,GMP">measure.csv
i=1
while [ $i -lt $N ]
do
	echo "$i,=average(A$st:A$en),=average(B$st:B$en),=average(C$st:C$en)">>measure.csv
	i=`expr $i + 1`
	st=`expr $st + $loop`
	en=`expr $en + $loop`
done

### measure ###

make clean all

outer=1
while [ $outer -lt $N ]
do
	i=0
	while [ $i -lt $loop ]
	do
		echo "$outer : $i"
		./mul `perl -e "print '11111111'x(4*$outer-1)"` >> measure.csv
		sleep 0.1s
		i=`expr $i + 1`
	done
	outer=`expr $outer + 1`
	cat hoge | sed "s/128 \* 1/128 \* $outer/" > ymm_mul.h
	make clean all
done
mv measure.csv data/
make clean


