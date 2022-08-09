import csv

from pymol import cmd

from common import act_on_set, locus_from_allele
from functions import split_peptide


file_root = '../'


def fetch_crystal_structure_peptide(pdb_code):
    success = False
    name = f'{complex["pdb_code"].lower()}_1_peptide'
    url = f'http://127.0.0.1:8080/structures/downloads/{name}.cif'
    print (pdb_code)
    try:
        cmd.load(url)
        chain = cmd.get_chains(name)[0]
        if chain != 'C':
            cmd.alter(f'chain {chain}', 'chain="C"')
            print (f'Change chain name from {chain} to C')
        file_key = f'{file_root}/structures/peptides/{pdb_code.lower()}.pdb'
        cmd.save(file_key)
        cmd.delete(name)
        success = True
    except:
        print (f'Cannot load {pdb_code}')
    return success




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
    if len(complex['pdb_code']) == 4:
        # and perform an action on them
        success = fetch_crystal_structure_peptide(complex['pdb_code'])
        if not success:
            errors.append(complex['pdb_code'])


print (errors)
        
