geodispersal-analysis.py

Usage: geodispersal-analysis.py is a script to assist in analyses that use the
modified Brooks Parsimony Analysis of Lieberman and Eldredge (1996).

The input is a NEXUS file with a data matrix and a single tree. The matrix 
should be coded with in the standard datatype, and should countain a single 
character. The different states of the character represent the different 
geographical regions in the analysis. The {}-braces are used to indicate the 
occurrence of a taxon in mulitple regions.  The tree depicts the evolutionary
relationships for the taxa.

For instance:

################################################################################
#NEXUS
BEGIN DATA;
	DIMENSIONS ntax = 4 nchar = 1;
	FORMAT datatype = standard symbols = "01234";
MATRIX
H_sapiens  {01234}
P_bonobo   0
G_gorilla  0
P_pygmaeus 2
;
END;
BEGIN TREES;
    TREE one = [&R] (P_pygmaeus,(G_gorilla,(P_bonobo,H_sapiens)));
END;
################################################################################
is a valid (but uninteresting) input file.

If you have PAUP installed (and on your PATH as 'paup'), then running
################################################################################
python geodispersal-analysis.py biogeoexample.nex --labels=labels.txt --vicariance=vic --dispersal=dis -p
################################################################################
will:
    1. read the area assignments and tree from biogeoexample.nex
    2. read the labels for the areas from labels.txt (one label per line, and the
        labels are interpretted as the names for the SYMBOLS in the FORMAT 
        command of the NEXUS file's DATA block).
    3. Lieberman and Eldredge's modified BPA analysis to produce a vicariance
        matrix (written to the file 'vic') and a dispersal matrix (written
        to the file 'dis').
    4. add PAUP commands to these files to perform a branch and bound search 
        using ordered parsimony, and save up to 100 most parsimonious trees to 
        either 'Vicariance.tre' or 'Dispersal.tre'
    5. invoke PAUP on the vic and dis files to produce the trees.

Omitting the -p option will skip the last step (so no trees will be produced).

Simply running:
################################################################################
python geodispersal-analysis.py biogeoexample.nex --labels=labels.txt
################################################################################
will print the vicariance and dispersal matrices to standard out.

Including multiple unnamed arguments (arguments which are not prefixed with 
dashes and option names, such as --vicariance) will indicate that multiple data
files are to be read and used to create a concatenated vicariance matrix and a
concatenated dispersal matrix. The character state symbols in each of the files
must use the same encoding of areas to SYMBOLS.
For example running
################################################################################
python geodispersal-analysis.py LiebermanECrassiproetus.nex LiebermanEBasidechenella.nex
################################################################################
produces matrices that contain the characters that would be created from
running LiebermanECrassiproetus.nex and LiebermanEBasidechenella.nex separately.



################################################################################
Optional arguments:
  -h, --help            show this help message and exit
  --vicariance=VICARIANCE
                        The name of the output file for the vicariance matrix
  --dispersal=DISPERSAL
                        The name of the output file for the dispersal matrix
  --labels=LABELS       Name of an (optional) file with labels for the areas.
                        The format for this file is simply one label per line.
                        Note that the labels should correspond to that the
                        state codes are listed in the SYMBOLS option of the
                        FORMAT command input files.
  -p, --paup            If specified, then PAUP* will be invoked to produce
                        Vicariance.tre and Dispersal.tre Note that both the
                        --vicariance and --dispersal options must also be used
                        if you use the --paup option
  --absent-as-missing   If specified, and an area is absent for all taxa in
                        one of the input analyses then the area will be coded
                        as missing (?).  If this option is not used, then the
                        area will be assigend a 0 for all characters generated
                        from that analysis.


Refer to the DendroPy instructions (two directories above this) for notes on
    how to install DendroPy
