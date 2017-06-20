formatConverters
================

Scripts to convert between different file formats.

###BratStandardToTabular.py
Converts the standard output of BRATNextGen (http://www.helsinki.fi/bsg/software/BRAT-NextGen/) to the tabular output.

Usage: BratStandardtoTabular.py [inputfile] [outputfile]

### EMBLtoBED.py
Converts EMBL annotation format to BED format (note score is automatically 1000).

Usage: EMBLtoBED.py [emblfile]

### EMBLtoGFF.py
Converts EMBL annotation header to GFF format. Also requires the core genome alignment with the "reference"
genome as the first sequence".

Usage: EMBLtoGFF.py [emblfile] [alignment]

###FastaReverseComplement.py
Reverse complements entries in a fasta file (can contain multiple records)

Requirements: Python 2.7.3, Biopython 1.63

Usage: FastaReverseComplement.py [inputfile] [outputfile]

###FastaToNexus.py
Converts a fasta format alignment to a non-interleaved nexus format alignment.

Requirements: Biopython (http://biopython.org/)

Current Versions: Python 2.7.3, Biopython 1.63

Usage: FastaToNexus.py [inputfile] [outputfile]

###FastaToNexusInterleave.py
Converts a fasta format alignmnent to an interleaved nexus format alignment.

Requirements: Biopython

Usage: FastaToNexusInterleaved.py [inputfile] [outputfile]

###FastaToPhylip.py
Converts a fasta format alignment to a sequential phylip format alignment. 

Requirements: Biopython

Usage: FastaToPhylip.py [inputfile] [outputfile]

###KodonToLDhat.py
Converts SNP table from Kodon to input for LDhat. **This script has not been tested recently.** 

Usage: KodonToLDhat.py [inputfile] [outputfile prefix] [reference name] [N or ALL]

###KodonToNexus.py
Converts SNP table from Kodon to Nexus format alignment. **This script has not been tested recently.**

Usage: Usage:  KodonToNexus.py  [input file] [outputfile] [name of reference sequence] [N or ALL]

###KodonToPhylip.py
Converts SNP table from Kodon to Phylip format alignment. **This script has not been tested recently.**

Usage: KodonToPhylip.py [inputfile] [outputfile] [namefile] [number of sequences] [number of nucleotides in alignment]

###gappedVCF.py
Replaces reference allele with a '-' in the VCF output from snp-sites at sites where the is a gap in the alignment.

Requirements: Biopython (http://biopython.org/)

Current Versions: Python 2.7.3, Biopython 1.63

Usage: gappedVCF.py [vcf] [fasta]
