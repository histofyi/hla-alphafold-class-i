from typing import Tuple, List

from pdbfixer import PDBFixer
from openmm.app import PDBFile
from Bio.PDB import *
import csv


file_root = '../'

disordered_peptides = []
ignore_structures = ['4LCY']



class EndResidues(Select):
    def __init__(self, residue_ids):
        self.residue_ids = residue_ids

    def accept_residue(self, residue):
        if residue.id[1] in self.residue_ids:
            return 1
        else:
            return 0



def extract_peptide_ends(pdb_code:str, peptide_sequence:str) -> Tuple[bool, List]:
    """
        This function extracts the N-terminal 3 residues and the C-terminal 2 residues from the peptide structure in several ways

    """

    input_folder = 'fixed'
    output_folder = 'ends_plus'

    input_file_name = f'{file_root}/structures/peptides/{input_folder}/{pdb_code}.pdb'
    output_file_name = f'{file_root}/structures/peptides/{output_folder}/{pdb_code}.pdb'

    print ('------------------------')
    print (pdb_code)

    peptide_length = len(peptide_sequence)

    ends = [1,2,3, peptide_length -1, peptide_length]
    ends_plus = [1,2,3, 4, peptide_length -2, peptide_length -1, peptide_length]

    print (len(peptide_sequence))
    print (ends)
    print (ends_plus)

    parser = PDBParser(PERMISSIVE=0)
    structure = parser.get_structure(pdb_code, input_file_name)
    
    io=PDBIO()
    io.set_structure(structure)
    io.save(output_file_name, EndResidues(ends_plus))
    


errors = []

# open the CSV file
file = open(f'{file_root}human_class_i.csv')
# and read in the CSV 
csvreader = csv.reader(file)
# set the header to the first row
header = next(csvreader)

# iterate through the remaining rows
rows = []
for row in csvreader:
    # if the specified allele is present then add it to the curated rowset
    rows.append(row)

# next create a dictionary of the complexex relating to a specific allele
complex_set = []
for row in rows:
    complex_set.append({k:v for k,v in list(zip(header,row))})

# iterate through the complexes relating to that allele
for complex in complex_set:
    pdb_code = complex['pdb_code']
    if len(pdb_code) == 4:
        # and perform an action on them if it's not a structure with problems
        if pdb_code not in ignore_structures:
            extract_peptide_ends(pdb_code.lower(), complex['peptide'])

        