# This script will parse a Contamination.txt file from the Sequence Read Archive (SRA)
# and clean the corresponding FASTA (assembly) and GFF (gene annotation, optional).

# This means removing / trimming sequences and adjusting the corresponding GFF
# file (removing genes if they overlap with the removed sequences).

# This script needs the following python libraries "getopt" / "sys" / "os"
# This script needs the following tools "bedtools" / 

global VERSION
VERSION = "0.1.0"

import getopt
import sys
import os

# Default values for the arguments:
assembly_fasta = None
annotation_gff = None
contamination_file = "Contamination.txt"
output_folder = "./sra-cleaning/"

# Help and usage:
def print_version():
    """Print the version.
    """
    global VERSION
    print("sra-cleaning v" + VERSION)

def usage():
    """Print the usage.
    """
    print("\npython sra-cleaning.py -a [assembly.fasta] -g [annotation.gff] -c [Contamination.txt] -o [output_folder]")
    print("\nInput:")
    print("     -a/--assembly               The assembly to clean in fasta format.")
    print("     -g/--gff                    Optional: a GFF file containing annotation for the assembly.")
    print("     -c/--contamination          The Contamination.txt file given by the SRA.")
    print("\nOutput:")
    print("     -o/--output                 The path to the output folder, './sra-cleaning/' by default")
    print("\nOther:")
    print("     -h/--help                   Print the usage and help and exit.")
    print("     -v/--version                Print the version and exit.")


# Get the arguments:
try:
    opts, args = getopt.getopt(sys.argv[1:], "a:g:c:o:hv", ["assembly=", "gff=", "contamination=", "output=", "help", "version"])
except getopt.GetoptError as err:
    print(str(err))
    usage()
    sys.exit(2)

for o,a in opts:
    if o in ("-a", "--assembly"):
        assembly_fasta = str(a)
    elif o in ("-g", "--gff"):
        annotation_gff = str(a)
    elif o in ("-c", "--contamination"):
        contamination_file = str(a)
    elif o in ("-o", "--output"):
        output_folder = str(a)
    elif o in ("-h", "--help"):
        usage()
        sys.exit(1)
    elif o in ("-v", "--version"):
        print_version()
        sys.exit(1)
    else:
        assert False, "Unhandled option !"

# Exit if output folder already exists (to avoid overwriting an already existing project)
if(os.path.isdir(output_folder)):
    print("\nFolder '"+output_folder+"' already exists, stopping now...\n")
    sys.exit(2)

# Checking if the FASTA is provided and exists:
if(assembly_fasta == None):
    print("The -a/--assembly argument is mandatory ('None' found).")
    sys.exit(2)
else:
    if not os.path.isfile(assembly_fasta):
        error_mssg = "Error: '"+str(assembly_fasta)+"' is not a file."
        sys.exit(2)

# Checking if the Contamination.txt is provided and exists:
if(contamination_file == None):
    print("The -c/--contamination argument is mandatory ('None' found).")
    sys.exit(2)
else:
    if not os.path.isfile(contamination_file):
        error_mssg = "Error: '"+str(contamination_file)+"' is not a file."
        sys.exit(2)

# Checking if the GFF is provided and exists:
if(annotation_gff == None):
    print("No -g/--gff given, the script will only correct the assembly.")
else:
    if not os.path.isfile(annotation_gff):
        error_mssg = "Error: '"+str(annotation_gff)+"' is not a file."
        sys.exit(2)


def parse_contamination(contamination_file):
    """Parse the contamination txt file and returns a list and a dict, one for the sequences to exclude,
    the other for the regions to mask/trim.
    Returns:
        A list and a dict, 'to_exclude' and 'to_trim'
        to_exclude is a list of sequence IDs to remove from the assembly /gff
        to_trim is a dict of sequences to mask/trim from the assembly / gf, with sequence IDs as key and regions as values.
    """
    to_exclude = list()
    to_trim = dict()

    exclude_start = -1
    in_exclude = False
    exclude_stop = -1

    trim_start = -1
    in_trim = False
    trim_stop = -1

    with open(contamination_file, 'r') as contamination:
        contamination_contents = contamination.readlines() 
        contamination_contents = [line.rstrip() for line in contamination_contents]
        index = -1

        for line in contamination_contents:
            index = index + 1
            if(line == "Exclude:"):
                exclude_start = index
                in_exclude = True
            elif(line == "Trim:"):
                trim_start = index
                in_trim = True
            elif(line == ""):
                if(in_exclude == True):
                    exclude_stop = index
                    in_exclude = False
                if(in_trim == True):
                    trim_stop = index
                    in_trim = False

        # start + 2 because we have the "Exlude:" line and a header line
        # stop corresponds to an empty line bu "[]" is exclusive for the last element (ie: "[[")
        if(exclude_start != -1):
            to_exclude = contamination_contents[(exclude_start + 2):(exclude_stop)]
        if(trim_start != -1):
            to_trim = contamination_contents[(trim_start + 2):(trim_stop)]

        print(to_exclude)
        print(to_trim)


parse_contamination(contamination_file)