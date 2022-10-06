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
from Bio.SeqRecord import SeqRecord

global VERSION
VERSION = "0.2.0"

# Help and usage:
def print_version():
    """Print the version"""
    global VERSION
    print("sra-cleaning v" + VERSION)

def usage():
    """Print the usage"""
    print("\npython sra-cleaning.py -a [assembly.fasta] -g [annotation.gff] -c [Contamination.txt] -o [output_folder]")
    print("\nInput:")
    print("     -a/--assembly           The assembly to clean in fasta format.")
    print("     -g/--gff                Optional: a GFF file containing annotation for the assembly.")
    print("     -c/--contamination      The Contamination.txt file given by the SRA.")
    print("Output:")
    print("     -o/--output             The path to the output folder, './sra-cleaning/' by default")
    print("\nParameters:")
    print("     -n/--nosplit            If the contamination is found in the middle of a sequence, do not split the sequence.\n"+ 
          "                             (By default, the sequence is split in two and '_1' and '_2' are added to the IDs).")
    print("Other:")
    print("     -h/--help               Print the usage and help and exit.")
    print("     -v/--version            Print the version and exit.")

def get_arguments():
    """Get the arguments from the command line."""
    # Default values for the arguments:
    assembly_fasta = None
    annotation_gff = None
    contamination_file = "Contamination.txt"
    output_folder = "./sra-cleaning/"
    split = True
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:g:c:o:nhv", ["assembly=", "gff=", "contamination=", 
                                                                 "output=", "nosplit", "help", "version"])
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
        elif o in ("-n", "--nosplit"):
            split = False
        elif o in ("-h", "--help"):
            usage()
            sys.exit(1)
        elif o in ("-v", "--version"):
            print_version()
            sys.exit(1)
        else:
            assert False, "Unhandled option !"

    return assembly_fasta, annotation_gff, contamination_file, output_folder, split


def check_arguments(assembly_fasta, annotation_gff, contamination_file, output_folder):
    """Check the different arguments (if folder exists, file exists) etc... and 
    exit if the conditions are not fulfilled."""
    # Exit if output folder already exists (to avoid overwriting an already existing folder)
    if(os.path.isdir(output_folder)):
        print("\nFolder '" + output_folder + "' already exists, stopping now...\n")
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
        contamination_contents = [line.rstrip() for line in contamination.readlines()]
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


def clean_assembly(assembly_fasta, to_exclude, to_trim, split):
    """This script will call exclude_sequences() and trim_sequences()
    """
    assembly_records = list(SeqIO.parse(assembly_fasta, 'fasta'))

    if(len(to_exclude) > 0):
        assembly_records = exclude_sequences(assembly_records, to_exclude)

    if(len(to_trim) > 0):
        assembly_records = trim_sequences(assembly_records, to_trim, split)

    assembly_records, ids_removed = polish_assembly(assembly_records)

    return assembly_records, ids_removed


def exclude_sequences(records, to_exclude):
    """Get the a list of SeqRecords and a list of sequences IDs to exclude.
    Then use Biopython to remove the corresponding sequences from the original 
    records list and returns a list of 'cleaned' SeqRecords
    """
    # Creating a new fasta file without the sequences in "to_exclude":
    after_exclude_records = [SeqRecord(record.seq, id = record.id, description = '', name = '') 
                             for record in records if record.id not in to_exclude]
    
    print("\nThe Contamination file listed ", len(to_exclude), " sequences to remove.")
    print("We removed ", str((len(records) - len(after_exclude_records))), " sequences")
    if((len(records) - len(after_exclude_records)) != len(to_exclude)):
        print("Warning! We were supposed to remove " + str(len(to_exclude)) + " sequences "+
        "but " + str((len(records) - len(after_exclude_records))) + " sequences were removed !")
    
    return after_exclude_records


def trim_sequences(records, to_trim, split):
    """This function takes a list of SeqRecords and a dictionary of sequences to 
    trim (with seqID as key and a list of [start, stop] as values.
    The function will mask or trim the sequences from the dict in the FASTA file."""
    trimmed_records = list()
    n_split = 0

    for record in records:
        if(record.id in to_trim.keys()):
            # trim_start-1 because python is 0-based, but [] is exclusive, so no need to trim_stop-1
            trim_start = int(to_trim.get(record.id)[0]) - 1
            trim_stop = int(to_trim.get(record.id)[1])
            # If the sequence is not at the start or end of the sequence and -n was not selected, we split the sequence:
            if(trim_start != 0 and trim_stop != len(record.seq) and split):
                trimmed_records.append(SeqRecord(record.seq[:trim_start], id = record.id + str("_1"), description = '', name = ''))
                trimmed_records.append(SeqRecord(record.seq[trim_stop:], id = record.id + str("_2"), description = '', name = ''))
                n_split = n_split + 1
            else:
                trimmed_records.append(SeqRecord(record.seq[:trim_start] + record.seq[trim_stop:], id = record.id, description = '', name = ''))
        # If not in list of sequences to trim, we just write the sequence as it is:
        else:
            trimmed_records.append(SeqRecord(record.seq, id = record.id, description = '', name = ''))

    len_input = sum([len(record.seq) for record in records])
    len_trimmed = sum([len(record.seq) for record in trimmed_records])
    len_cont_file = sum([(int(to_trim.get(key)[1]) - int(to_trim.get(key)[0]) + 1) for key in to_trim.keys()])

    print("\nBefore Trimming: ", str(len_input), " bases.")
    print("After Trimming: ", str(len_trimmed), " bases")
    print("The script removed " + str(len_input - len_trimmed) + " bases")
    if((len_input - len_trimmed) != len_cont_file):
        print("Warning! We were supposed to remove " + str(len_cont_file) + " bases "+
        "but " + str(len_input - len_trimmed) + " bases were removed !")
    if(split):
        print("Note: " + str(n_split) + " sequences were split (due to internal contamination, use --nosplit if you want to avoid this.)")

    return trimmed_records


def polish_assembly(assembly_records):
    """In this function we follow the recommendations from the Contamination.txt file:
     'After you remove the contamination, trim any Ns at the ends of the sequence
      and remove any sequences that are shorter than 200 nt'
      Returns the cleaned records and a list of IDs that have been removed"""
    polished_records = list()
    ids_removed = list()
    n_N = 0

    print("\nPolishing the assembly:")
    for record in assembly_records:
        # Remove trailing Ns:
        if(record.seq[-1] == 'N' or record.seq[-1] == 'N'):
            i = -1
            while(record.seq[i] == "N" or record.seq[i] == "n"):
                i = i - 1
            # To cancel the "last" loop (which corresponds to a non-N base):
            i = i + 1
            n_N = n_N + abs(i)
            record.seq = record.seq[:i]
        
        # Check sequence length (remove if < 200 nt):
        if(len(record.seq) < 200):
            print("Removing " + record.id + " (< 200 nt)")
            ids_removed.append(record.id)
        else:
            polished_records.append(SeqRecord(record.seq, id = record.id, description = '', name = ''))

    print("\nRemoved " + str(len(ids_removed)) + " sequences (< 200 nt)")
    print("Removed " + str(n_N) + " trailing Ns.")

    return polished_records, ids_removed


def clean_gff(annotation_gff, assembly_fasta, to_exclude, to_trim, split, ids_removed):
    """Take an annotation file in GFF format and removes the annotation overlapping
    with the removed / trimmed sequences in the assembly.
    """
    gff_contents = list()
    features_removed = list()

    with open(annotation_gff, 'r') as gff_handle:
        gff_contents = [line.rstrip() for line in gff_handle.readlines()]
        len_gff = len(gff_contents)

    if(len(to_exclude) > 0):
        gff_contents, features_removed = exclude_gff(gff_contents, to_exclude)

    if(len(to_trim) > 0):
        #assembly_records will be used to get the sequences length to check for split sequences
        assembly_records = list(SeqIO.parse(assembly_fasta, 'fasta'))
        gff_contents, features_trimmed = trim_gff(gff_contents, assembly_records, to_trim, split)
        if(len(features_trimmed) > 0):
            features_removed.append("\n".join(features_trimmed))

    # We remove the features corresponding to the sequences that have been polished (len < 200 nt)
    gff_contents, features_polished = exclude_gff(gff_contents, ids_removed)
    if(len(features_polished) > 0):
        features_removed.append("\n".join(features_polished))

    print("The cleaned gff contains " + str(len(gff_contents)) + " lines (We removed " + str(len_gff - len(gff_contents)) + " lines).")

    return gff_contents, features_removed


def exclude_gff(gff_contents, to_exclude):
    """Will exclude lines in the gff which have IDs belonging to the excluded sequences.
    """
    after_exclude_gff = list()
    features_removed = list()

    for gff_line in gff_contents:
        if(gff_line[0] == "#" or gff_line.split("\t")[0] not in to_exclude):
            after_exclude_gff.append(gff_line)
        else:
            features_removed.append(gff_line)
    
    return after_exclude_gff, features_removed


def trim_gff(gff_contents, assembly_records, to_trim, split):
    """Will remove from the GFF the lines which have features overlapping with the
    trimmed sequences. And rename the IDs of the features belonging to split sequences.
    """
    trimmed_gff = list()
    features_trimmed = list()
    assembly_dict = SeqIO.to_dict(assembly_records)

    for gff_line in gff_contents:
        gff_line_split = gff_line.split("\t")
        # If the line is a comment (#), meta-data (##) or not in the trimmed sequences
        # we just add it to the trimmed gff list.
        if(gff_line[0] == "#" or gff_line_split[0] not in to_trim.keys()):
            trimmed_gff.append(gff_line)
        else:
            # GFF is 1-based os no need for the -1 here.
            trim_start = int(to_trim.get(gff_line_split[0])[0])
            trim_stop = int(to_trim.get(gff_line_split[0])[1])
            seq_len = len(assembly_dict[gff_line_split[0]].seq)
            # We only keep the feature if it is outside the trimmed region:
            if(int(gff_line_split[3]) < trim_start and int(gff_line_split[4]) < trim_start):
                # If the trimmed region was internal, we rename the ID accordingly:
                if(trim_start != 1 and trim_stop != seq_len and split):
                    gff_line_split[0] = gff_line_split[0] + "_1"
                trimmed_gff.append("\t".join(gff_line_split))
            # We only keep the feature if it is outside the trimmed region:
            # /!\ if the sequence is after the trimming we need to update the coordinates accrodingly
            elif(int(gff_line_split[3]) > trim_stop and int(gff_line_split[4]) > trim_stop):
                # If the trimmed region was internal, we rename the ID accordingly:
                if(trim_start != 1 and trim_stop != seq_len and split):
                    gff_line_split[0] = gff_line_split[0] + "_2"
                    # We update the feature's start and stop, we shift by stop position,
                    # as the new 1st position of splitted contig_2 starts at "trim_stop"
                    gff_line_split[3] = str(int(gff_line_split[3]) - trim_stop)
                    gff_line_split[4] = str(int(gff_line_split[4]) - trim_stop)
                else:
                    # We update the feature's start and stop position, we shift by trimmed region's length:
                    gff_line_split[3] = str(int(gff_line_split[3]) - (trim_stop - trim_start + 1))
                    gff_line_split[4] = str(int(gff_line_split[4]) - (trim_stop - trim_start + 1))

                trimmed_gff.append("\t".join(gff_line_split))
            else:
                features_trimmed.append(gff_line)

    return trimmed_gff, features_trimmed


def main():
    assembly_fasta, annotation_gff, contamination_file, output_folder, split = get_arguments()
    check_arguments(assembly_fasta, annotation_gff, contamination_file, output_folder)

    os.makedirs(output_folder)
    output_fasta = str(output_folder+"/cleaned.fasta")
    output_gff = str(output_folder+"/cleaned.gff")
    removed_gff = str(output_folder+"/removed_features.gff")

    to_exclude, to_trim = parse_contamination(contamination_file)
    
    cleaned_records, ids_removed = clean_assembly(assembly_fasta, to_exclude, to_trim, split)
    print("\nWriting cleaned fasta file to:", output_fasta)
    SeqIO.write(cleaned_records, output_fasta, "fasta")

    if(annotation_gff != None):
        print("\nCleaning the GFF file:")
        gff_clean, features_removed = clean_gff(annotation_gff, assembly_fasta, to_exclude, to_trim, split, ids_removed)

        print("Writing the cleaned GFF to " + output_gff)
        with open(output_gff, 'w') as gff_out:
            gff_out.write('\n'.join(gff_clean) + "\n")

        print("Writing the removed features to " + removed_gff)
        with open(removed_gff, 'w') as removed_out:
            removed_out.write('\n'.join(features_removed) + "\n")

if __name__ == "__main__":
    main()
