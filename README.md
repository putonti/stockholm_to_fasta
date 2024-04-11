# stockholm_to_fasta
This code was developed to translate a Stockholm-format MSA with UniProt identifiers to FASTA-format file with NCBI taxonomy identifier.

Written in Python3 by Catherine Putonti

Command-line script takes in 2 parameters,
* -i, --input: A PFAM stockholm file with canonical Uniprot identifiers for each sequence
* -o, --output: FIlename to write out a Stockholm MSA file and unaligned FASTA file with NCBI taxonomy ids for each sequence.
