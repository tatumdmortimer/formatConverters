formatConverters
================

Scripts to convert between different file formats.

###BratStandardToTabular.py
Converts the standard output of BRATNextGen (http://www.helsinki.fi/bsg/software/BRAT-NextGen/) to the tabular output.

Usage: BratStandardtoTabular.py [inputfile] [outputfile]

###FastaToNexus.py
Converts a fasta format alignment to a non-interleaved nexus format alignment.

Requirements: Biopython (http://biopython.org/)

Current Versions: Python 2.7.3, Biopython 1.63

Usage: FastaToNexus.py [inputfile] [outputfile]

###gappedVCF.py
Replaces reference allele with a '-' in the VCF output from snp-sites at sites where the is a gap in the alignment.

Requirements: Biopython (http://biopython.org/)

Current Versions: Python 2.7.3, Biopython 1.63

Usage: gappedVCF.py [vcf] [fasta]
