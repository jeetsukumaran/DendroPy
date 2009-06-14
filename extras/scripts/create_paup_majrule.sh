#!/bin/sh
for i in 16 32 64 128 256 500 512  ;
do
	if test -d t$i
	then
		echo '#NEXUS'> t$i/paup.nex ; 
		echo "set maxtrees=2000 ; Begin Taxa ; dimensions ntax = $i ; taxlabels " >> t$i/paup.nex ;
		cat t$i/rep0t$i/nuc.GTRIG.boot.tre | grep -e '^  [0-9]' | sed -E "s/.+ ([-a-zA-Z\'0-9]+)[,;]/\1/" >> t$i/paup.nex ;
		echo '; end ; begin paup ; log start replace ; ' >> t$i/paup.nex  ;
		for j in 0 1 2 3 4 5 6 7 8 9 ; 
		do
			echo "gettrees file = rep${j}t${i}/nuc.GTRIG.boot.tre mode = 7; " >> t$i/paup.nex ; 
		done ;
		echo 'contree / nostrict majrule le50 grpfreq treefile=bootmajrule.tre replace;' >> t$i/paup.nex ;
		echo 'log stop ;  quit ; end ;' >> t$i/paup.nex ;
	else
		echo "directory t${i} not found"
	fi
done 
