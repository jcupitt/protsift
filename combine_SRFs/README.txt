This program was written for the paper ELUCIDATION OF THE PATHWAYS OF 
XENOBIOTIC METABOLISM IN HUMAN SKIN AND HUMAN SKIN MODELS BY PROTEOMIC 
PROFILING.

You give it a set of SRFs of a large number of experiments and it produces a
list of the xenobiotic metabolising proteins detected. It knows which fraction
(cytoplasm or membrane) the protein was found in. It 
tries to name proteins sensibly. It only reports the most significant
occurence of each protein. Where a peptide can appear in many proteins, it
reports all possible protein hits.

The program uses John Prince's mspire library to read SRF files.

http://mspire.rubyforge.org/
