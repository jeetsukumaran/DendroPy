**********************
Examining Data Objects
**********************

High-level summaries of the contents of DendroPy phylogenetic data objects are given by the :meth:`description()` instance method of the |Tree|, |TreeList|, |CharacterMatrix|-derived, and |DataSet| classes.
This method optionally takes a numeric value as its first argument that determines the level of detail (or depth) of the summary::

    >>> import dendropy
    >>> d = dendropy.DataSet.get_from_path('pythonidae.nex', 'nexus')
    >>> print(d.description())
    DataSet object at 0x79dd0: 1 Taxon Sets, 0 Tree Lists, 1 Character Matrices
    >>> print(d.description(3))
    DataSet object at 0x79dd0: 1 Taxon Sets, 0 Tree Lists, 1 Character Matrices
        [Taxon Sets]
            [0] TaxonSet object at 0x5a4a20 (TaxonSet5917216): 29 Taxa
                [0] Taxon object at 0x22c0fd0 (Taxon36442064): 'Python regius'
                [1] Taxon object at 0x22c0f10 (Taxon36441872): 'Python sebae'
                [2] Taxon object at 0x22c0ed0 (Taxon36441808): 'Python brongersmai'
                [3] Taxon object at 0x22c0f70 (Taxon36441968): 'Antaresia maculosa'
                [4] Taxon object at 0x22c0f30 (Taxon36441904): 'Python timoriensis'
                [5] Taxon object at 0x22c0f50 (Taxon36441936): 'Python molurus'
                [6] Taxon object at 0x22c0ff0 (Taxon36442096): 'Morelia carinata'
                [7] Taxon object at 0x23ae050 (Taxon37412944): 'Morelia boeleni'
                [8] Taxon object at 0x23ae030 (Taxon37412912): 'Antaresia perthensis'
                [9] Taxon object at 0x23ae070 (Taxon37412976): 'Morelia viridis'
                [10] Taxon object at 0x23ae090 (Taxon37413008): 'Aspidites ramsayi'
                [11] Taxon object at 0x23ae0b0 (Taxon37413040): 'Aspidites melanocephalus'
                [12] Taxon object at 0x22c0fb0 (Taxon36442032): 'Morelia oenpelliensis'
                [13] Taxon object at 0x23ae0d0 (Taxon37413072): 'Bothrochilus boa'
                [14] Taxon object at 0x23ae130 (Taxon37413168): 'Morelia bredli'
                [15] Taxon object at 0x23ae110 (Taxon37413136): 'Morelia spilota'
                [16] Taxon object at 0x23ae150 (Taxon37413200): 'Antaresia stimsoni'
                [17] Taxon object at 0x23ae0f0 (Taxon37413104): 'Antaresia childreni'
                [18] Taxon object at 0x23ae1b0 (Taxon37413296): 'Leiopython albertisii'
                [19] Taxon object at 0x23ae170 (Taxon37413232): 'Python reticulatus'
                [20] Taxon object at 0x23ae190 (Taxon37413264): 'Morelia tracyae'
                [21] Taxon object at 0x23ae1d0 (Taxon37413328): 'Morelia amethistina'
                [22] Taxon object at 0x23ae230 (Taxon37413424): 'Morelia nauta'
                [23] Taxon object at 0x23ae250 (Taxon37413456): 'Morelia kinghorni'
                [24] Taxon object at 0x23ae210 (Taxon37413392): 'Morelia clastolepis'
                [25] Taxon object at 0x23ae290 (Taxon37413520): 'Liasis fuscus'
                [26] Taxon object at 0x23ae2b0 (Taxon37413552): 'Liasis mackloti'
                [27] Taxon object at 0x23ae270 (Taxon37413488): 'Liasis olivaceus'
                [28] Taxon object at 0x23ae2f0 (Taxon37413616): 'Apodora papuana'
        [Character Matrices]
            [0] DnaCharacterMatrix object at 0x22c0f90 (DnaCharacterMatrix36442000):  29 Sequences
                [Taxon Set]
                    TaxonSet object at 0x5a4a20 (TaxonSet5917216): 29 Taxa
                [Characters]
                    [0] Python regius : 1114 characters
                    [1] Python sebae : 1114 characters
                    [2] Python brongersmai : 1114 characters
                    [3] Antaresia maculosa : 1114 characters
                    [4] Python timoriensis : 1114 characters
                    [5] Python molurus : 1114 characters
                    [6] Morelia carinata : 1114 characters
                    [7] Morelia boeleni : 1114 characters
                    [8] Antaresia perthensis : 1114 characters
                    [9] Morelia viridis : 1114 characters
                    [10] Aspidites ramsayi : 1114 characters
                    [11] Aspidites melanocephalus : 1114 characters
                    [12] Morelia oenpelliensis : 1114 characters
                    [13] Bothrochilus boa : 1114 characters
                    [14] Morelia bredli : 1114 characters
                    [15] Morelia spilota : 1114 characters
                    [16] Antaresia stimsoni : 1114 characters
                    [17] Antaresia childreni : 1114 characters
                    [18] Leiopython albertisii : 1114 characters
                    [19] Python reticulatus : 1114 characters
                    [20] Morelia tracyae : 1114 characters
                    [21] Morelia amethistina : 1114 characters
                    [22] Morelia nauta : 1114 characters
                    [23] Morelia kinghorni : 1114 characters
                    [24] Morelia clastolepis : 1114 characters
                    [25] Liasis fuscus : 1114 characters
                    [26] Liasis mackloti : 1114 characters
                    [27] Liasis olivaceus : 1114 characters
                    [28] Apodora papuana : 1114 characters

If you want to see the data in a particular schema, you can call the :meth:`as_string()` method, passing it a schema-specification string ("nexus", "newick", "fasta", "phylip", etc.), as well as other optional arguments specific to varous formats::

    >>> import dendropy
    >>> d = dendropy.DataSet.get_from_path('pythonidae.nex', 'nexus')
    >>> print(d.as_string("nexus"))
    >>> print(d.as_string("fasta"))
    >>> print(d.as_string("phylip"))

It is also possible to get an ASCII text plot of a tree object::

    >>> import dendropy
    >>> t = dendropy.Tree.get_from_path('pythonidae.mle.nex', 'nexus')
    >>> print(t.as_ascii_plot())
                              /---------- Python regius
                              |
       /----------------------+      /--- Python sebae
       |                      |   /--+
       |                      \---+  \--- Python molurus
       |                          |
       |                          \------ Python curtus
       |
       |                             /--- Morelia bredli
       |            /----------------+
       |            |                \--- Morelia spilota
       |            |
       |            |      /------------- Morelia tracyae
       |         /--+      |
       |         |  |   /--+      /------ Morelia clastolepis
       |         |  |   |  |  /---+
       |         |  |   |  |  |   |  /--- Morelia kinghorni
       |         |  |   |  \--+   \--+
       |         |  \---+     |      \--- Morelia nauta
       |         |      |     |
       |         |      |     \---------- Morelia amethistina
       |      /--+      |
       |      |  |      \---------------- Morelia oenpelliensis
       |      |  |
    /--+      |  |            /---------- Antaresia maculosa
    |  |      |  |         /--+
    |  |      |  |         |  |   /------ Antaresia perthensis
    |  |      |  |         |  \---+
    |  |      |  |         |      |  /--- Antaresia stimsoni
    |  |      |  \---------+      \--+
    |  |      |            |         \--- Antaresia childreni
    |  |      |            |
    |  |      |            |      /------ Morelia carinata
    |  |      |            \------+
    |  |      |                   |  /--- Morelia viridisN
    |  |  /---+                   \--+
    |  |  |   |                      \--- Morelia viridisS
    |  |  |   |
    |  |  |   |                      /--- Apodora papuana
    |  |  |   |                   /--+
    |  |  |   |                   |  \--- Liasis olivaceus
    |  |  |   |               /---+
    |  |  |   |               |   |  /--- Liasis fuscus
    |  |  |   |               |   \--+
    +  |  |   |            /--+      \--- Liasis mackloti
    |  |  |   |            |  |
    |  \--+   |            |  |      /--- Antaresia melanocephalus
    |     |   |         /--+  \------+
    |     |   |         |  |         \--- Antaresia ramsayi
    |     |   |         |  |
    |     |   \---------+  |         /--- Liasis albertisii
    |     |             |  \---------+
    |     |             |            \--- Bothrochilus boa
    |     |             |
    |     |             \---------------- Morelia boeleni
    |     |
    |     |                          /--- Python timoriensis
    |     \--------------------------+
    |                                \--- Python reticulatus
    |
    |                                /--- Xenopeltis unicolor
    |--------------------------------+
    |                                \--- Candola aspera
    |
    \------------------------------------ Loxocemus bicolor

