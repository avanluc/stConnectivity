#!/bin/sh

N=`head -1 $1 | cut -d" " -f 1`
for i in `seq 0 $N`
	do
		#../build/stConnectivity $1 1 0 $i
		o1=`../build/stConnectivity $1 1 0 $i `
		o2=`../build/stConnectivity $1 3 0 $i `

		if [ "$o1" != "$o2" ]; then
			echo $o1
			echo $o2
		fi
	done