PEP_CLUSTER

This program reads a set of lists of detected peptides and outputs a
table summarising the families of proteins probably present. 

It tries to map peptides to proteins intelligently: it outputs every protein
that a perptide could be evidence for, and it looks for relationships between
protein sequences which indicate shared function.

SAMPLE DATA

Some sample data is included. Run with:

$ ./combine2.rb row1.csv 
loading row1.csv ...
found 644 peptides
loading human_protein2010.fasta ...
loaded 29180 proteins from human_protein2010.fasta
ignored  9603 'PREDICTED:' proteins
building protein table ...
indexing protein database ...
found 560 proteins
selecting and combining regions for each protein ...
removing duplicate peptides ...
(1674 peptides before dedupe)
(1674 peptides after dedupe)
removing proteins with less than two peps ...
(91 proteins remain)
building protein groups ...
(38 protein groups found)
linking peptides to groups ...
finding cross-group peptides ...
(104 cross-group peptides found)
linking groups ...
(4 protein clusters found)
finding all unused groups ...
(27 unused groups)
writing output.csv ...
writing output2.csv ...

The file "output.csv" looks something like this:

-------
** cluster of 4 proteins
32567786,NP_787028.1,"keratin, type II cytoskeletal 79"
119703753,NP_005546.2,"keratin, type II cytoskeletal 6B"
155969697,NP_775109.2,"keratin, type II cytoskeletal 6C"
5031839,NP_005545.1,"keratin, type II cytoskeletal 6A"
[.. peptide evidence supporting this cluster]
** cluster of 2 proteins
24430190,NP_002266.2,"keratin, type I cytoskeletal 15"
24430192,NP_005548.2,"keratin, type I cytoskeletal 16"
195972866,NP_000412.3,"keratin, type I cytoskeletal 10"
4557701,NP_000413.1,"keratin, type I cytoskeletal 17"
24234699,NP_002267.2,"keratin, type I cytoskeletal 19"
15431310,NP_000517.2,"keratin, type I cytoskeletal 14"
114431246,NP_853513.2,"keratin, type I cytoskeletal 28"
[.. peptide evidence supporting this cluster]
** cluster of 3 proteins
41322908,NP_958781.1,plectin isoform 1e
41322910,NP_958783.1,plectin isoform 1d
41322912,NP_958780.1,plectin isoform 1f
41322914,NP_958785.1,plectin isoform 1g
41322916,NP_958782.1,plectin isoform 1
41322919,NP_958784.1,plectin isoform 1b
207452735,NP_112598.2,epiplakin
47607492,NP_000436.2,plectin isoform 1c
[.. peptide evidence supporting this cluster]
** cluster of 2 proteins
160420317,NP_001104026.1,filamin-A isoform 2
116805322,NP_001449.3,filamin-C isoform a
188595687,NP_001120959.1,filamin-C isoform b
116063573,NP_001447.2,filamin-A isoform 1
[.. peptide evidence supporting this cluster]
[.. more proteins follow]
-------

The file output2.csv simply lists all detected proteins.



