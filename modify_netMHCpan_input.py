# This script modifies the fasta files so that they contain
# shorter headers for input into netMHCpan. It also produces a 
# dictionary file that maps the shorter header to the complete header

from Bio import SeqIO
import numpy as np

# Input and output file names
original_file = './tim_data/frameshiftPeptidesComplete.fasta'
corrected_file = './tim_data/MOD_frameshiftPeptidesComplete.fasta'

# Create a dictionary with id mapping
prefix = 'XID:'
ID_dict = {}

# Create a new FASTA file with changed named
with open(original_file) as original, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for it, record in enumerate(records):
        new_name = prefix+str(it)
        ID_dict[new_name] = str(record.id)
        record.id = new_name
        record.description = new_name
        SeqIO.write(record, corrected, 'fasta')
        if it%100000==0 and it !=0:
            print("Peptide number: {}".format(it))
            
# Save the look-up dictionary as well
np.save('./tim_data/ID_dict',ID_dict, allow_pickle=True)
