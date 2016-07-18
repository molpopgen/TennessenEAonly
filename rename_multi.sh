#!/bin/bash

for lambda in 0.1 0.25 0.5
do
	for h in 0 1 
	do
		for i in $(seq 1 1 100)
		do
		fn=multi.lambda"$lambda".VA.h"$h"."$i".h5
		n=`basename $fn "$i".h5`
		#echo $fn "$n"0."$i".h5
		fn2="$n"0."$i".h5
		echo $fn $fn2
		cp $fn $fn2
		done
	done
done
	
