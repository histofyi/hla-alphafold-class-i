from pymol import cmd
import os
import csv


file_root = '../'



def build_filepath(directory_type, locus, allele, peptide):
    return f'{file_root}{locus.lower()}/{allele.lower()}/{directory_type}/{peptide}'


def build_runpath(run_folder, directory_type, locus, allele, peptide):
    return f'{build_filepath(directory_type, locus, allele, peptide)}/{run_folder}'


def align_prediction_set(locus, allele, peptide):
    runs = [folder for folder in os.listdir(build_filepath('unrelaxed', locus, allele, peptide)) if '.result' in folder]

    for run_folder in runs:
        original_filepath = build_runpath(run_folder, 'unrelaxed', locus, allele, peptide)
        print (original_filepath)
        try:
            aligned_filepath = build_runpath(run_folder, 'unrelaxed_aligned', locus, allele, peptide)
            os.makedirs(aligned_filepath)

            print('directory created')
        except:
            print('directory exists')
        pdb_files = [file for file in os.listdir(build_runpath(run_folder, 'unrelaxed', locus, allele, peptide)) if '.pdb' in file]
        for pdb_file in pdb_files:
            align_to_canonical(pdb_file, original_filepath, aligned_filepath)
            print (pdb_file)
        print ('-------')


def align_to_canonical(filename:str, original_filepath:str, aligned_filepath:str):
    print ('-----')
    print (filename)
    cmd.load('canonical/cannonical_class_i_1hhk.pdb', quiet=0)
    original_pdb_file = f'{original_filepath}/{filename}'
    print (original_pdb_file)
    cmd.load(original_pdb_file, quiet=0)
    name = filename.replace('.pdb','')
    print (name)
    align = cmd.cealign('cannonical_class_i_1hhk', name)
    align['rmsd'] = align['RMSD']
    del align['RMSD']
    cmd.delete('cannonical_class_i_1hhk')
    print (align['rmsd'])
    file_types = ['pdb']
    for file_type in file_types:
        file_key = f'{aligned_filepath}/{name}.{file_type}'
        print (file_key)
        cmd.save(file_key)
    cmd.delete(name)
    return align






file = open(f'{file_root}human_class_i.csv')
csvreader = csv.reader(file)
header = []
header = next(csvreader)

print (header)


allele = 'HLA-G*0104'

locus = allele.split('*')[0]
print (allele)
print (locus)

rows = []
for row in csvreader:
    if allele in row:
        rows.append(row)

print (rows)

alleleset = []
for row in rows:
    alleleset.append({k:v for k,v in list(zip(header,row))})

for complex_row in alleleset:
    print (complex_row['peptide'])
    align_prediction_set(locus, allele, complex_row['peptide'])

