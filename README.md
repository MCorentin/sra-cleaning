# sra-cleaning

Python script to automatically parse a "Contamination.txt" file from the Sequence Read Archive (SRA) and correct the assembly FASTA  and annotation GFF files.

## Motivation

When an assembly is submitted to the SRA, some checks are performed, including a Contamination check. The SRA then sends a "Contamination.txt" file containing a list of sequences to remove / trim. This script will parse the "Contamination.txt" file and output ready-to-submit assembly (FASTA) and annotation (GFF) files to submit to SRA after the contamination cleaning.

**The script is performing some checks to assess if the correct length was trimmed from the assembly, but it is highly recommended to manually check if the proper changes were made in the output files!**

Note: I could not find a "standard" format for the Contamination.txt files. The script will parse the file to find the sequences below the "Exclude:" and "Trim:" lines. If your Contamination.txt file is different, feel free to contact me or create and issue on GitHub to improve the parsing within this script.

## Dependencies

sra-cleaning was developed and tested on Linux (Ubuntu 16.04.7 LTS) but it is pure python so I expect it should work on other OS as well.

You will need to have the following dependencies:
- **Python v3 or higher (tested with Python 3.7.4)**
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
````

````
Input:
     -a/--assembly          The assembly to clean in fasta format.
     -g/--gff                Optional: a GFF file containing annotation for the assembly.
     -c/--contamination      The Contamination.txt file given by the SRA.

Output:
     -o/--output             The path to the output folder, './sra-cleaning/' by default
Parameters:
    -m/--mask                If this option is selected, the script will mask the sequences instead of trimming them (the 'exclude' sequences will still be removed).
   
Other:
     -h/--help               Print the usage and help and exit.
     -v/--version            Print the version and exit.
````

Note: the -g/--gff option is not supported yet.

## TO DO

 - Add the support for gff files.
 - Add the option to mask (with Ns) instead of trimming.