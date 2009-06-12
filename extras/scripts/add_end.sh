#!/bin/sh
for j in 16 32 64 128 256 500 512 1024 ;
do
	if test -d t$j
	then
		for i in 0 1 2 3 4 5 6 7 8 9 
		do
			if test -z `tail t${j}/rep${i}t${j}/nuc.GTRIG.boot.tre | grep 'end'` 
			then
				echo 'end;' >> t${j}/rep${i}t${j}/nuc.GTRIG.boot.tre 
			fi
		done
	else
		echo "no directory t$j"
	fi
done
