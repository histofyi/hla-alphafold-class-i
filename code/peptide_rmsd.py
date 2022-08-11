from biopandas.pdb import PandasPdb
import csv
import os
import re
import json
from typing import List, Dict, Tuple

from common import locus_from_allele, build_filepath, build_runpath 


pdb_code = '6AT9'

file_root = '../'

real_structure = f'{file_root}structures/peptides/fixed/{pdb_code.lower()}.pdb'


real = PandasPdb().read_pdb(real_structure)
test = PandasPdb()

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

alpha_fold_row = None
# next create a dictionary of the complexex relating to a specific allele
complex_set = []
for row in rows:
    if pdb_code in row:
        alpha_fold_row = {k:v for k,v in list(zip(header,row))}


i = 0

rmsds = {
    'pdb_code':pdb_code,
}
recycles = [recycle for recycle in range(3,11)]


def parse_colabfold_filename(filename:str) -> Dict:
    run_details = {}
    labels = []
    components = filename.split('_')
    run_details['allele'] = '-'.join(components[0:2])
    run_details['peptide'] = components[2]
    run_details['uid'] = components[3]
    run_details['recycles'] = re.findall(r'\d+', components[4])[0]
    relaxed = lambda x : False if x == 'unrelaxed' else True
    run_details['relaxed'] = relaxed(components[5])
    run_details['rank'] = components[7]
    run_details['model'] = components[9].replace('.pdb','')
    return run_details


def process_colabfold_structure(run_file_path:str, real) -> Tuple[Dict,bool,List]:
    predicted_structure_data = {}
    test = PandasPdb().read_pdb(run_file_path)

    # Lamda function to select only side chain atoms in a dataframe (masks out backbone atoms)
    sel = lambda df: df['atom_name'].str.contains('[^(CA|HA|N|C|O|HN|H)]$')


    #TODO put in checking that the whole peptide dataframes are the same shape. They may not be in the case of missing atoms in residues
    backbone_rmsd = PandasPdb.rmsd(test.df['ATOM'], real.df['ATOM'], s='main chain') 
    all_atoms_rmsd = PandasPdb.rmsd(test.df['ATOM'], real.df['ATOM'], s='heavy')
    side_chain_rmsd = PandasPdb.rmsd(test.df['ATOM'][sel], real.df['ATOM'][sel], s='heavy')
    
    print ('----------')
    print (allele)
    print (peptide)
    print (run_file)
    print (f'pymol {real_structure} {run_file_path}')
    print ('----------')


    print (f'Full length - backbone rmsd: {backbone_rmsd}')
    print (f'Full length - all atom rmsd: {all_atoms_rmsd}')
    print ('----------')

    all_positions_plddt_sum = 0
    predicted_structure_data['all_positions'] = {
        'backbone_rmsd':backbone_rmsd,
        'all_atoms_rmsd':all_atoms_rmsd,
        'side_chain_rmsd':side_chain_rmsd
    }
    predicted_structure_data['individual_positions'] = {}
    for residue_position in residue_positions:
        
        test_residue = test.df['ATOM'][test.df['ATOM']['residue_number'] == residue_position]
        real_residue = real.df['ATOM'][real.df['ATOM']['residue_number'] == residue_position]
        predicted_structure_data['individual_positions'][str(residue_position)] = {
            'amino_acid':test_residue["residue_name"].max(), 
            'position':residue_position}
        print (f'P{residue_position}{test_residue["residue_name"].max()}')
        
        if test_residue.shape != real_residue.shape:
            print ('----------')
            print (f'P{residue_position} {real_residue.shape}')
            print (real_residue)
            print (f'P{residue_position} {test_residue.shape}')
            print (test_residue)
            print ('----------')
        else:
            residue_pldddt = round(test_residue['b_factor'].mean(), 2)
            all_positions_plddt_sum += residue_pldddt
            print (f'Residue - plddt: {residue_pldddt}')
            residue_backbone_rmsd = PandasPdb.rmsd(real_residue, test_residue, s='main chain')
            residue_all_atom_rmsd = PandasPdb.rmsd(real_residue, test_residue, s='heavy')
            if test_residue["residue_name"].max() != 'GLY':
                residue_side_chain_rmsd = PandasPdb.rmsd(real_residue[sel], test_residue[sel], s='heavy')
            else:
                residue_side_chain_rmsd = None
            print (f'Residue - backbone rmsd: {residue_backbone_rmsd}')
            print (f'Residue - all atom rmsd: {residue_all_atom_rmsd}')   
            print (f'Residue - side chain rmsd: {residue_side_chain_rmsd}')  
            predicted_structure_data['individual_positions'][str(residue_position)]['plddt'] = residue_pldddt
            predicted_structure_data['individual_positions'][str(residue_position)]['backbone_rmsd'] = residue_backbone_rmsd
            predicted_structure_data['individual_positions'][str(residue_position)]['all_atom_rmsd'] = residue_all_atom_rmsd
            predicted_structure_data['individual_positions'][str(residue_position)]['side_chain_rmsd'] = residue_side_chain_rmsd
        print ('----------')
    predicted_structure_data['all_positions']['plddt'] = round(all_positions_plddt_sum / len(peptide), 4)
    predicted_structure_data['pymol'] = f'pymol {real_structure.replace(file_root,"")} {run_file_path.replace(file_root,"")}'
    return predicted_structure_data, True, []






if alpha_fold_row:
    allele = alpha_fold_row['allele']
    locus = locus_from_allele(allele)
    peptide = alpha_fold_row['peptide']

    residue_positions = [position + 1 for position in range(0, len(peptide))]


    rmsds['locus'] = locus
    rmsds['allele'] = allele
    rmsds['peptide'] = peptide

    rmsds['recycles'] = {}

    for recycle in recycles:
        rmsds['recycles'][str(recycle)] = {'recycle_count': recycle, 'models':[]}


    folder_type = 'unrelaxed_split_peptide'
    folder_path = build_filepath(folder_type, locus, allele, peptide, file_root)
    run_folders = [folder for folder in os.listdir(folder_path) if '.result' in folder]

    for run_folder in run_folders:
        run_folder_path = f'{folder_path}/{run_folder}'
        run_files = [run_file for run_file in os.listdir(run_folder_path) if '.pdb' in run_file]
        for run_file in run_files:
            run_file_path = f'{run_folder_path}/{run_file}'
            run_details = parse_colabfold_filename(run_file)
            if i == 0:
                rmsds['relaxed'] = run_details['relaxed']
            predicted_structure_data, success, errors = process_colabfold_structure(run_file_path, real)
            rmsds['recycles'][run_details['recycles']]['models'].append({
                'model':run_details['model'], 
                'rank': run_details['rank'],
                'predicted_structure_data': predicted_structure_data
            })            
    i += 1

    try:
        statistics_directory = f'{file_root}{locus.lower()}/{allele.lower()}/statistics'
        os.makedirs(statistics_directory)
        print('directory created')
    except:
        print('directory already exists')
    file_name = f'{statistics_directory}/{peptide.lower()}.json'
    print (file_name)
    with open(file_name, 'w') as output_file:
        json.dump(rmsds, output_file, sort_keys = True, indent = 4, ensure_ascii = True)

