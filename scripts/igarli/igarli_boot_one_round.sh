#!/bin/sh

#!/bin/sh
set -x

prevRound=$1
if ! test -f "${prevRound}"
then
	echo "expecting the file ${prevRound} to exist!"
	exit 1
fi 

wtsFile=$2
if ! test -f "${wtsFile}"
then
	echo "expecting the file ${wtsFile} to exist!"
	exit 1
fi 

ntaxa=$3
if test -z $ntaxa
then
	echo "expecting 2 arguments the filepath to the trees from the previous round and the number of taxa."
	exit 1
fi 

if test -z $DENDROPY_SCRIPTS_PAR
then
	DENDROPY_SCRIPTS_PAR="${HOME}/Documents/projects/dendropy/dendropy/scripts"
fi


################################################################################
# add the taxa to the trees from the previous round
echo "${DENDROPY_SCRIPTS_PAR}/igarli_boot_add_tree.py" "${prevRound}" "${wtsFile}"
(set -x ; time "${DENDROPY_SCRIPTS_PAR}/igarli_boot_add_tree.py" "${prevRound}" "${wtsFile}" >> boot_add_taxon_commands.txt) 2>time_igarli_boot_add_tree.txt 
if ! test $? -eq 0
then
	cat time_igarli_add_tree.txt 
	exit 1
fi
echo 'quit' >> boot_add_taxon_commands.txt 

iGarli ../boot_add_garli.conf < boot_add_taxon_commands.txt >igarli_boot_add_err_out.txt 2>&1
if ! test $? -eq 0
then
	cat igarli_boot_add_err_out.txt 
	exit 1
fi

echo '#NEXUS' > incrgarli.tre
echo 'begin trees;' >> incrgarli.tre
grep '\[iGar' igarli_boot_add_err_out.txt >>incrgarli.tre  || exit
echo 'end;' >> incrgarli.tre

