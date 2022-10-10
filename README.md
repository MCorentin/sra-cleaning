# sra-cleaning

Python script to automatically parse a "Contamination.txt" file from the Sequence Read Archive (SRA) and correct the assembly FASTA and annotation GFF files.

## Motivation

When an assembly is submitted to the SRA, some checks are performed, including a Contamination check. The SRA then sends a "Contamination.txt" file containing a list of sequences to remove / trim. This script will parse the "Contamination.txt" file and output ready-to-submit assembly (FASTA) and annotation (GFF) files to submit to SRA.

**The script is performing some sanity checks to assess if the proper changes were made to the assembly, but it is highly recommended to manually check the output files!**

**note: sra-cleaning.py does not yet support multiple spans per sequence for the trimming step.**

Note: I could not find a "standard" format for the Contamination.txt files. The script will parse the file to find the sequences below the "Exclude:" and "Trim:" lines. If your Contamination.txt file is different, feel free to contact me or create an issue on GitHub to improve the parsing within this script.

## Dependencies

sra-cleaning was developed and tested on Linux (Ubuntu 16.04.7 LTS) using Python 3.7.4.

You will need to have the following dependencies:
- **Python v3 or higher**
- **The following python packages:** [biopython](https://biopython.org/ "biopython Homepage"), getopt, sys, os.

Note: the following packages should be installed by default:
```
getopt
sys
os 
```

## Usage

````
python sra-cleaning.py -a [assembly.fasta] -g [annotation.gff] -c [Contamination.txt] -o [output_folder]

Input:
     -a/--assembly           The assembly to clean in fasta format.
     -g/--gff                Optional: a GFF file containing annotation for the assembly.
     -c/--contamination      The Contamination.txt file given by the SRA.
Output:
     -o/--output             The path to the output folder, './sra-cleaning/' by default

Parameters:
     -n/--nosplit            If the contamination is found in the middle of a sequence, do not split the sequence.
                             (By default, the sequence is split in two and '_1' and '_2' are added to the IDs).
Other:
     -h/--help               Print the usage and help and exit.
     -v/--version            Print the version and exit.
````

## Workflow:

For the assembly FASTA:

1. Remove the sequences listed in the "Exclude:" list.
2. Trim the regions from the "Trim:" list.
     * By default, if the region to trim falls in the middle of a sequence (ie: the start is not 1 or the end does not correspond to the sequence's length) the sequence will be split in two and the resulting sequences will have "_1" and "_2" added to the end of their IDs. Use "--nosplit" if you want to avoid this.
     * If the region to trim is at the start or end of the sequence, it will just be trimmed from the sequence.
3. Polish the assembly as requested by the "Contamination.txt" file:
     * Removing trailing Ns.
     * Removing sequences < 200 nt.

For the annotation GFF:

1. Remove the features with _seqname_ listed in the "Exclude:" sequences
2. Remove the features which overlap with the "Trim:" regions
     * If the region falls in the middle of a sequence, add "_1" and "_2" to the *seqname* of features situated before and after the region respectively (unless --nosplit was used).
     * If a feature falls after a trimmed region, update the coordinates (if "--nosplit" is used, substract by the length of the trimmed region, else substract by the 'end' coordinate of the trimmed region (since "sequence_2" start again from 1)).
3. Remove the features with _seqname_ that were removed in the polishing step (sequence length < 200 nt).

## Output files:

- _cleaned.fasta_: contains the filtered FASTA file.
- _cleaned.gff_: contains the filtered GFF annotations.
- _removed\_features.gff_: contains the removed features.

## ToDo list:

- Implement parsing for multiple spans on the same sequence in the "Contamination.txt" file.
- Add more checks (check if GFF features have >= 9 columns and if the GFF and FASTA have corresponding sequence IDs)
