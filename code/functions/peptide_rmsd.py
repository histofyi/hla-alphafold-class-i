from typing import Dict, List, Tuple

from biopandas.pdb import PandasPdb
from Bio.Data.IUPACData import protein_letters_3to1

import re
import json


from common import get_run_file_dict, locus_from_allele, parse_colabfold_filename, build_filepath


file_root = '../'



# -- Functions for generating RMSDs between peptides -- #


def quality_bin_structure(rmsd):
    quality_bins = [5, 3, 2, 1.5, 1.25]
    quality_labels = ['poor', 'neutral', 'good', 'very_good', 'excellent']
    i = 0
    for quality_bin in quality_bins:
        if rmsd > quality_bin and i ==0:
            this_quality = 'very_poor'
        elif rmsd <= quality_bin:
            this_quality = quality_labels[i]
        i += 1
    return this_quality


def generate_run_rmsd_statistics(complex_path:str, run_file_dict:Dict, ground_truth_structure) -> Dict:
    best_rmsd = 100
    best_model = 0
    best_rank = 0
    best_plddt = 0
    best_recycle = 0
    best_filename = ''
    best_complex_path = ''
    rmsd_statistics = {
        'best':{},
        'recycles':{}
    }
    for recycle in run_file_dict:
        recycle_name_components = recycle.split('_')
        recycle_count = re.findall(r'\d+', recycle_name_components[len(recycle_name_components)-1])[0]
        recycle_rmsd_statistics = generate_recycle_rmsd_statistics(complex_path, recycle, run_file_dict[recycle], ground_truth_structure)
        rmsd_statistics['recycles'][recycle_count] = recycle_rmsd_statistics
        if recycle_rmsd_statistics['best']['allatom_rmsd'] < best_rmsd:
            best_rmsd = recycle_rmsd_statistics['best']['allatom_rmsd']
            best_sidechain_rmsd = recycle_rmsd_statistics['best']['sidechain_rmsd']
            best_backbone_rmsd = recycle_rmsd_statistics['best']['backbone_rmsd']
            best_model = recycle_rmsd_statistics['best']['model']
            best_rank = recycle_rmsd_statistics['best']['rank']
            best_plddt = recycle_rmsd_statistics['best']['plddt']
            best_filename = recycle_rmsd_statistics['best']['filename']
            best_recycle = recycle_count
            best_complex_path = f'{complex_path}/{recycle}'
    rmsd_statistics['best']['allatom_rmsd'] = best_rmsd
    rmsd_statistics['best']['sidechain_rmsd'] = best_sidechain_rmsd
    rmsd_statistics['best']['backbone_rmsd'] = best_backbone_rmsd
    rmsd_statistics['best']['model'] = best_model
    rmsd_statistics['best']['rank'] = best_rank
    rmsd_statistics['best']['recycle'] = best_recycle
    rmsd_statistics['best']['plddt'] = best_plddt
    rmsd_statistics['best']['quality'] = quality_bin_structure(best_rmsd)
    rmsd_statistics['best']['complex_path'] = best_complex_path
    rmsd_statistics['best']['filename'] = best_filename
    print (rmsd_statistics['best'])
    return rmsd_statistics


def generate_recycle_rmsd_statistics(complex_path:str, recycle:str, recycle_models:List, ground_truth_structure) -> Dict:
    best_rmsd = 100
    best_model = 0
    best_rank = 0
    best_plddt = 0
    best_filename = ''
    recycle_rmsd_statistics = {'best':{},'models':{}}
    for recycle_model in recycle_models:
        model_rmsd_statistics = generate_model_rmsd_statistics(complex_path, recycle, recycle_model, ground_truth_structure)
        recycle_rmsd_statistics['models'][model_rmsd_statistics['model']] = model_rmsd_statistics
        if model_rmsd_statistics['allatom_rmsd'] < best_rmsd:
            best_rmsd = model_rmsd_statistics['allatom_rmsd']
            best_sidechain_rmsd = model_rmsd_statistics['sidechain_rmsd']
            best_backbone_rmsd = model_rmsd_statistics['backbone_rmsd']
            best_model = model_rmsd_statistics['model']
            best_rank = model_rmsd_statistics['rank']
            best_plddt = model_rmsd_statistics['plddt']
            best_filename = recycle_model
    recycle_rmsd_statistics['best']['allatom_rmsd'] = best_rmsd
    recycle_rmsd_statistics['best']['sidechain_rmsd'] = best_sidechain_rmsd
    recycle_rmsd_statistics['best']['backbone_rmsd'] = best_backbone_rmsd
    recycle_rmsd_statistics['best']['model'] = best_model
    recycle_rmsd_statistics['best']['rank'] = best_rank
    recycle_rmsd_statistics['best']['plddt'] = best_plddt
    recycle_rmsd_statistics['best']['filename'] = best_filename
    #print (recycle_rmsd_statistics['best'])
    return recycle_rmsd_statistics


def generate_model_rmsd_statistics(complex_path:str, recycle:str, recycle_model:str, ground_truth_structure) -> Dict:
    # Lamda function to select only side chain atoms in a dataframe (masks out backbone atoms)
    sel = lambda df: df['atom_name'].str.contains('[^(CA|HA|N|C|O|HN|H)]$')
    
    model_details = parse_colabfold_filename(recycle_model) 
    model_filename = f'{complex_path}/{recycle}/{recycle_model}'
    model_filepath = f'{file_root}{model_filename}'

    model_rmsd_statistics = {
        'filename': model_filename,
        'model': model_details['model'],
        'rank': model_details['rank'],
        'residues':{}
    }

    model_structure = PandasPdb().read_pdb(model_filepath)

    if model_structure.df['ATOM'].shape == ground_truth_structure.df['ATOM'].shape:
        model_rmsd_statistics['backbone_rmsd'] = PandasPdb.rmsd(model_structure.df['ATOM'], ground_truth_structure.df['ATOM'], s='main chain') 
        model_rmsd_statistics['allatom_rmsd'] = PandasPdb.rmsd(model_structure.df['ATOM'], ground_truth_structure.df['ATOM'], s='heavy')
        model_rmsd_statistics['sidechain_rmsd'] = PandasPdb.rmsd(model_structure.df['ATOM'][sel], ground_truth_structure.df['ATOM'][sel], s='heavy')
    else:
        #TODO error handling
        print ('SHAPE MISMATCH')

    peptide_length = len(ground_truth_structure.amino3to1())
    
    residue_positions = [position + 1 for position in range(0, peptide_length)]

    peptide_plddt_sum = 0

    for residue_position in residue_positions:
        residue_rmsd_statistics = generate_residue_level_rmsds(model_structure, ground_truth_structure, residue_position)
        peptide_plddt_sum += residue_rmsd_statistics['plddt']
        model_rmsd_statistics['residues'][residue_position] = residue_rmsd_statistics
    model_rmsd_statistics['plddt'] = round(peptide_plddt_sum/peptide_length, 2)
    return model_rmsd_statistics


def generate_peptide_level_rmsds(model_structure, ground_truth_structure) -> Dict:
    pass


def generate_residue_level_rmsds(model_strcucture, ground_truth_structure, residue_position) -> Dict:
    sel = lambda df: df['atom_name'].str.contains('[^(CA|HA|N|C|O|HN|H)]$')

    model_residue = model_strcucture.df['ATOM'][model_strcucture.df['ATOM']['residue_number'] == residue_position]
    ground_truth_residue = ground_truth_structure.df['ATOM'][ground_truth_structure.df['ATOM']['residue_number'] == residue_position]
    residue_plddt = round(model_residue['b_factor'].mean(), 2)
    residue_backbone_rmsd = PandasPdb.rmsd(model_residue, ground_truth_residue, s='main chain')
    residue_allatom_rmsd = PandasPdb.rmsd(model_residue, ground_truth_residue, s='heavy')

    if model_residue["residue_name"].max() != 'GLY':
        residue_sidechain_rmsd = PandasPdb.rmsd(model_residue[sel], ground_truth_residue[sel], s='heavy')
    else:
        residue_sidechain_rmsd = None
    residue_rmsd_statistics = {
        'plddt': residue_plddt,
        'backbone_rmsd': residue_backbone_rmsd,
        'allatom_rmsd': residue_allatom_rmsd,
        'sidechain_rmsd': residue_sidechain_rmsd
    }    
    return residue_rmsd_statistics





def generate_rmsd_set(row):
    input_folder = 'unrelaxed_split_peptide'
    allele = row['allele']
    locus = locus_from_allele(allele)
    peptide = row['peptide']
    pdb_code = row['pdb_code'].upper()
    ground_truth_structure_filename = f'structures/peptides/fixed/{pdb_code.lower()}.pdb'

    rmsd_set = {
        'metadata':{
            'allele':allele,
            'locus':locus,
            'peptide':peptide,
            'pdb_code':pdb_code,
            'peptide_structure':ground_truth_structure_filename
        }
    }

    complex_path = build_filepath(input_folder, locus, allele, peptide, '')

    

    run_file_dict = get_run_file_dict(locus, allele, peptide, input_folder, file_root)
    
    print ('--------')
    print (pdb_code)

    ground_truth_structure_filepath = f'{file_root}{ground_truth_structure_filename}' 
    
    try:
        ground_truth_structure = PandasPdb().read_pdb(ground_truth_structure_filepath)
    except:
        ground_truth_structure = None
        print ('==========')
        print ('==========')
        print (f'CANNOT LOAD {pdb_code}')
        print ('==========')
        print ('==========')


    if ground_truth_structure:
        run_rmsd_statistics = generate_run_rmsd_statistics(complex_path, run_file_dict, ground_truth_structure)
    
        rmsd_set['best'] = run_rmsd_statistics['best']
        rmsd_set['recycles'] = run_rmsd_statistics['recycles']
    
        try:
            statistics_directory = f'{file_root}{locus.lower()}/{allele.lower()}/statistics'
            os.makedirs(statistics_directory)
            print('directory created')
        except:
            print('directory already exists')
        
        file_name = f'{statistics_directory}/{peptide.lower()}.json'

        print (file_name)
        with open(file_name, 'w') as output_file:
            json.dump(rmsd_set, output_file, sort_keys = True, indent = 4, ensure_ascii = True)

    print ('--------')

    return rmsd_set


