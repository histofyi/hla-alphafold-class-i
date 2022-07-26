from typing import Tuple, List
import argparse
import csv

from pdbfixer import PDBFixer
from openmm.app import PDBFile
from Bio.PDB import *


file_root = '../'

disordered_peptides = []
ignore_structures = ['4LCY']


class RemoveHydrogens(Select):
    def accept_atom(self, atom):
        return removable_atom(atom.get_name())     

def removable_atom(atom_name):
    is_removable = 1
    if 'H' in atom_name or 'OXT' in atom_name:
            if atom_name != 'OH':
                if atom_name != 'NH1':
                    if atom_name != 'NH2': 
                        if atom_name != 'CH2':
                            is_removable = 0
    return is_removable



def fix_peptide_structure(pdb_code:str, peptide_sequence:str) -> Tuple[bool, List]:
    """
        This function checks the peptide structure in several ways

        1. Using PDBFixer to clean up
        2. Using a readlines method on the cleaned .pdb file to make sure there is an empty alternate location indicator (index 16 on the line)
        3. Using BioPython to look for disordered residues or atoms

        if any of 2. or 3. find evidence of disorder, the 'not_disordered' is set to false

        Notes:

        PDBFixer sets the chain of a single chain .pdb file to 'A'. This is why it's set to 'C' in the BioPython block.

`       Args:

            pdb_code (str): the PDB code of the peptide to clean e.g. 3KYO
            peptide_sequence (str): the sequence of the peptide within the PDB file e.g. KLPAQFYIL
    
    #TODO - find which step removes all the B-factor information and sort it out.
    """

    input_folder = 'raw'
    output_folder = 'fixed'
    not_disordered = True
    reasons = []

    input_file_name = f'{file_root}structures/peptides/{input_folder}/{pdb_code}.pdb'
    output_file_name = f'{file_root}structures/peptides/{output_folder}/{pdb_code}.pdb'

    print ('------------------------')
    print (pdb_code)


    # 1. Use PDBFixer to do initial clean up

    fixer = PDBFixer(filename=input_file_name)
    fixer.findMissingResidues()
    fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # write out the cleaned file
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_file_name, 'w'))

    # 2. Check for a non-empty alternate location column in the .pdb file
    with open(output_file_name) as f:
        for line in f.readlines():
            if 'ATOM' in line:
                if line[16] != ' ':
                    print (f'{line[16]} - {line[17:20]}')
                    not_disordered = False
                    reasons.append('alternate_location')

    # 3. Use BioPython to check disordered state in residues or atoms 
    parser = PDBParser(PERMISSIVE=0)
    structure = parser.get_structure(pdb_code, output_file_name)
    i = 0
    for model in structure:
        for chain in model:
            last_residue_id = 0
            if chain.id != 'C':
                chain.id = 'C'
            for residue in chain:
                i += 1
                residue_id = residue.id[1]
                if i == 1:
                    if residue_id != 1:
                        print ('wrong_peptide_start_id')
                        not_disordered = False
                        reasons.append('wrong_peptide_start_id')
                if abs(last_residue_id - residue_id) > 1:
                    not_disordered = False
                    reasons.append(f'broken_chain_at_{last_residue_id}')
                last_residue_id = residue_id
                residue_name = residue.get_resname()
                if residue.is_disordered():
                    not_disordered = False
                    print (f'{residue_name}{residue_id} is disordered')
                    reasons.append(f'disordered_residue_{residue_id}')
                for atom in residue:
                    is_removable = removable_atom(atom.get_name())
                    if is_removable == 0:
                        print (f'{residue_name}{residue_id} - {atom.get_name()}')
                    if atom.is_disordered():
                        not_disordered = False
                        print (f'{residue_name}{residue_id} - {atom.id} is disordered')
                        reasons.append(f'disordered_atom_{atom.id}_at_residue_{residue_id}')
    if i != len(peptide_sequence):
        not_disordered = False
        reasons.append(f'unequal_chain_length_{i}_vs_{len(peptide_sequence)}_for_{peptide_sequence}')
    
    # 4. remove hydrogens
    # write out the cleaned file (with relabelled chain) the selector on the save removes hydrogens
    io=PDBIO()
    io.set_structure(structure)
    io.save(output_file_name, RemoveHydrogens())
    

    if not not_disordered:
        print ('disorder')
    print ('------------------------')
    return not_disordered, reasons


errors = []


def fix_peptides_in_csv(this_pdb_code):
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

    if this_pdb_code:
        complex_set = [complex for complex in complex_set if complex['pdb_code'] == this_pdb_code.upper()]
        if len(complex_set) == 0:
            errors.append('pdb_file_not_present')

    # iterate through the complexes relating to that allele
    for complex in complex_set:
        pdb_code = complex['pdb_code']
        if len(pdb_code) == 4:
            # and perform an action on them if it's not a structure with problems
            if pdb_code not in ignore_structures:                    
                success, reasons = fix_peptide_structure(pdb_code.lower(), complex['peptide'])
                if not success:
                    errors.append({'pdb_code':complex['pdb_code'], 'reasons':reasons})



parser = argparse.ArgumentParser()
parser.add_argument("--single_structure", help="specify a single file to run the function on")
args = parser.parse_args()
if args.single_structure:
    print(f'PDB code is {args.single_structure}')
    fix_peptides_in_csv(args.single_structure)
else:
    print ('working on complete set')
    fix_peptides_in_csv()


print (errors)

