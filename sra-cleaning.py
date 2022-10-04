# This script will parse a Contamination.txt file from the Sequence Read Archive (SRA)
# and clean the corresponding FASTA (assembly) and GFF (gene annotation, optional).
#
# This means removing / trimming sequences in the FASTA and adjusting the 
# corresponding GFF file (removing genes if they overlap with the removed sequences).
#
# This script needs the following python libraries "getopt" / "sys" / "os" / "Biopython"
#
# Author: Molitor Corentin, 2022

import getopt, sys, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

global VERSION
VERSION = "0.1.0"

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
    print("     -a/--assembly           The assembly to clean in fasta format.")
    print("     -g/--gff                Optional: a GFF file containing annotation for the assembly.")
    print("     -c/--contamination      The Contamination.txt file given by the SRA.")
    print("\nOutput:")
    print("     -o/--output             The path to the output folder, './sra-cleaning/' by default")
    print("\nOther:")
    print("     -h/--help               Print the usage and help and exit.")
    print("     -v/--version            Print the version and exit.")

def get_arguments():
    """Get the arguments from the command line."""
    # Default values for the arguments:
    assembly_fasta = None
    annotation_gff = None
    contamination_file = "Contamination.txt"
    output_folder = "./sra-cleaning/"
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:g:c:o:hv", ["assembly=", "gff=", "contamination=", 
                                                                "output=", "help", "version"])
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

    return assembly_fasta, annotation_gff, contamination_file, output_folder


def check_arguments(assembly_fasta, annotation_gff, contamination_file, output_folder):
    """Check the different arguments (if folder exists, file exists) etc... and 
    exit if the conditions are not fulfilled."""
    # Exit if output folder already exists (to avoid overwriting an already existing folder)
    if(os.path.isdir(output_folder)):
        print("\nFolder '"+output_folder+"' already exists, stopping now...\n")
        sys.exit(2)

    # Checking if the FASTA is provided and exists:
    if(assembly_fasta == None):
        print("The -a/--assembly argument is mandatory ('None' found).")
        sys.exit(2)
    elif not os.path.isfile(assembly_fasta):
        print("Error for option -a/--assembly: '"+str(assembly_fasta)+"' is not a file.")
        sys.exit(2)

    # Checking if the Contamination.txt is provided and exists:
    if(contamination_file == None):
        print("The -c/--contamination argument is mandatory ('None' found).")
        sys.exit(2)
    elif not os.path.isfile(contamination_file):
        print("Error for option -c/--contamination: '"+str(contamination_file)+"' is not a file.")
        sys.exit(2)

    # Checking if the GFF is provided (optional, just send a message) and exists:
    if(annotation_gff == None):
        print("No -g/--gff given, the script will only correct the assembly.")
    elif not os.path.isfile(annotation_gff):
        print("Error for option -g/--gff: '"+str(annotation_gff)+"' is not a file.")
        sys.exit(2)


def parse_contamination(contamination_file):
    """Parse the contamination txt file and returns a list and a dict, one for the sequences to exclude,
    the other for the regions to mask/trim.
    Returns:
        A list and a dict, 'to_exclude' and 'to_trim'
        'to_exclude' is a list of sequence IDs to remove from the assembly / gff
        'to_trim' is a dict of sequences to mask/trim from the assembly / gff, with sequence IDs as key and regions as values.
    """
    to_exclude = list()
    exclude_start = -1
    exclude_stop = -1
    in_exclude = False

    to_trim = dict()
    trim_start = -1
    trim_stop = -1
    in_trim = False

    with open(contamination_file, 'r') as contamination:
        contamination_contents = [line.rstrip() for line in contamination.readlines() ]

        for index, line in enumerate(contamination_contents):
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

        # start + 2 because we have the "Exclude:" line and + a header line
        # stop corresponds to an empty line bu "[]" is exclusive for the last element (ie: "[[")
        if(exclude_start != -1):
            exclude_lines = contamination_contents[(exclude_start + 2):(exclude_stop)]
            for line in exclude_lines:
                to_exclude.append(line.split("\t")[0])
        
        # Creates a dictonary of sequences to trim, with sequence ID as key 
        # and list of two values [start, stop] as value, in the contamination file
        # the start and stop values are registered as follow: "start..stop"
        if(trim_start != -1):
            trim_lines = contamination_contents[(trim_start + 2):(trim_stop)]
            for line in trim_lines:
                id = line.split("\t")[0]
                start = line.split("\t")[2].split("..")[0]
                stop = line.split("\t")[2].split("..")[1]
                to_trim.update( {id : [start,stop]} )

        return to_exclude, to_trim


def clean_assembly(assembly_fasta, to_exclude, to_trim):
    """This script will call exclude_sequences() and trim_sequences()"""
    # First, records contains the assembly:
    assembly_records = list(SeqIO.parse(assembly_fasta, 'fasta'))
    # If we have to exclude sequences, they will be removed from 'records' by exclude_sequences()
    if(len(to_exclude) > 0):
        assembly_records = exclude_sequences(assembly_records, to_exclude)
    # If we have to mask/trim sequences, they will be trimmed from 'records' by trim_sequences()
    if(len(to_trim) > 0):
        assembly_records = trim_sequences(assembly_records, to_trim)
    # We return the "cleaned" records
    return assembly_records


def exclude_sequences(records, to_exclude):
    """Get the a list of SeqRecords and a list of sequences IDs to exclude.
    Then use Biopython to remove the corresponding sequences from the original 
    records list and returns a list of 'cleaned' SeqRecords"""
    # Creating a new fasta file without the sequences in "to_exclude":
    after_exclude_records = [SeqRecord(record.seq, id = record.id, description = '', name = '') 
                             for record in records if record.id not in to_exclude]
    
    print("\nOriginal FASTA had ", len(records), " sequences.")
    print("The filtered FASTA has ", len(after_exclude_records), " sequences")
    if((len(records) - len(after_exclude_records)) != len(to_exclude)):
        print("Warning! We were supposed to remove " + str(len(to_exclude)) + " sequences "+
        "but " + str((len(records) - len(after_exclude_records))) + " sequences were removed !")
    
    return after_exclude_records


def trim_sequences(records, to_trim):
    """This function takes a list of SeqRecords and a dictionary of sequences to 
    trim (with seqID as key and a list of [start, stop] as values.
    The function will mask or trim the sequences from the dict in the FASTA file."""
    trimmed_records = list()

    for record in records:
        if(record.id in to_trim.keys()):
            # +1 because python is 0-based, but [] is exclusive, so no '+1' to stop
            start = int(to_trim.get(record.id)[0]) - 1
            stop = int(to_trim.get(record.id)[1])
            trimmed_records.append(SeqRecord(record.seq[:start] + record.seq[stop:], id = record.id, description = '', name = ''))
        else:
            trimmed_records.append(SeqRecord(record.seq, id = record.id, description = '', name = ''))

    len_input = sum([len(record.seq) for record in records])
    len_trimmed = sum([len(record.seq) for record in trimmed_records])
    len_cont_file = sum([(int(to_trim.get(key)[1]) - int(to_trim.get(key)[0]) + 1) for key in to_trim.keys()])

    print("\nBefore Trimming: ", str(len_input), " bases.")
    print("After Trimming: ", str(len_trimmed), " bases")
    print("The script trimmed " + str(len_input - len_trimmed) + " bases")

    if((len_input - len_trimmed) != len_cont_file):
        print("Warning! We were supposed to remove " + str(len_cont_file) + " bases "+
        "but " + str(len_input - len_trimmed) + " bases were removed !")
    
    return trimmed_records


def main():
    assembly_fasta, annotation_gff, contamination_file, output_folder = get_arguments()
    check_arguments(assembly_fasta, annotation_gff, contamination_file, output_folder)

    os.makedirs(output_folder)
    output_fasta = str(output_folder+"/cleaned.fasta")

    to_exclude, to_trim = parse_contamination(contamination_file)
    print(to_exclude, to_trim)
    
    cleaned_records = clean_assembly(assembly_fasta, to_exclude, to_trim)
    print("\nWriting cleaned fasta file to:", output_fasta)
    SeqIO.write(cleaned_records, output_fasta, "fasta")

    if(annotation_gff != None):
        print("\nSorry the GFF option is not supported yet!")

if __name__ == "__main__":
    main()
