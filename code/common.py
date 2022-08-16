from typing import Dict, List, Tuple, Optional, Union, Callable

import os
import re
import json
import csv
import argparse


file_root = '../'


def locus_from_allele(allele:str) -> str:
    if '*' in allele:
        locus = allele.split('*')[0]
    else:
        locus = None
    return locus


def build_filepath(directory_type:str, locus:str, allele:str, peptide:str, file_root:str) -> str:
    """
    This function creates the filepath for the folder containing the AlphaFold prediction runs

    Args:
        directory_type (str): the type of directory, e.g. 'unrelaxed', 'unrelaxed_aligned'
        locus (str): the locus of the complex e.g. 'HLA-A'
        allele (str): the allele of the complex e.g. 'HLA-A*68:02'
        peptide (str): the peptide in the complex e.g. 'TLTSCNTSV'
        file_root (str): the root directory for all the files e.g. '../'

    Returns: 
        string: the returned filepath
    
    """
    return f'{file_root}{locus.lower()}/{allele.lower()}/{directory_type}/{peptide}'


def build_runpath(run_folder, directory_type, locus, allele, peptide, file_root:str) -> str:
    """
    This function creates the filepath for the folder containing a specific AlphaFold prediction run

    Args:
        run_folder (str): the folder name of a specific run, e.g. 'HLA_A6802_TLTSCNTSV_8d0ef_3un.result'
        directory_type (str): the type of directory, e.g. 'unrelaxed', 'unrelaxed_aligned'
        locus (str): the locus of the complex e.g. 'HLA-A'
        allele (str): the allele of the complex e.g. 'HLA-A*68:02'
        peptide (str): the peptide in the complex e.g. 'TLTSCNTSV'
        file_root (str): the root directory for all the files e.g. '../'

    Returns: 
        string: the returned filepath
    
    """
    return f'{build_filepath(directory_type, locus, allele, peptide, file_root)}/{run_folder}'



def parse_args() -> Tuple[str, str, List]:
    errors = []
    parser = argparse.ArgumentParser()
    parser.add_argument("--single_structure", help="specify a single file to run the function on")
    args = parser.parse_args()
    if args.single_structure:
        mode = 'single'
        pdb_code = args.single_structure.upper()
        if len(pdb_code) > 4:
            print (f'The pdb code specified is too long : {pdb_code}')
            errors.append('pdb_code_too_long')
    else:
        mode = 'multiple'
        pdb_code = None
    return pdb_code, mode, errors




def load_list_from_csv(listname:str) -> Tuple[List, bool, List]:
    # open the CSV  file
    file = open(f'{file_root}{listname}')
    # and read in the CSV 
    csvreader = csv.reader(file)

    # set the header to the first row
    header = next(csvreader)

    # iterate through the remaining rows
    all_rows = []
    success = True
    errors = []
    for row in csvreader:
        alpha_fold_row = {k:v for k,v in list(zip(header,row))}
        if len(alpha_fold_row['pdb_file_issues']) == 0:
            all_rows.append(alpha_fold_row)
        else:
            errors.append({
                'pdb_code':alpha_fold_row['pdb_code'],
                'error':alpha_fold_row['pdb_file_issues']})
    return all_rows, success, errors




def perform_action(action:Callable, listname:str) -> List:
    pdb_code, mode, errors = parse_args()
    if len(errors) == 0:
        actable_list, success, list_errors = load_list_from_csv(listname)
        if len(list_errors) > 0:
            for error in list_errors:
                errors.append(error)
        action_errors = act_on_list(action, actable_list, pdb_code)
        if len(action_errors) > 0:
            for error in action_errors:
                errors.append(error)
    return errors



def act_on_list(action:Callable, rowlist=List, pdb_code:Union[str,None] = None, ) -> List:
    errors = []
    if pdb_code:
        pdb_code_found = False 
        print (f'Act only on {pdb_code}')
        for row in rowlist:
            if row['pdb_code'] == pdb_code:
                action(row)
                pdb_code_found = True
        if not pdb_code_found:
            errors.append({
                'pdb_code': pdb_code,
                'error': 'not_present_in_list'
            })
    else:
        print ('Act on all in list')
        for row in rowlist:
            action(row)
    return errors


def parse_colabfold_filename(filename:str) -> Dict:
    #TODO check naming of all files before release
    #Ones to rename have either 9 or 11 components in .pdb filename
    run_details = {}
    components = filename.split('_')
    relaxed = lambda x : False if x == 'unrelaxed' else True
    if len(components) == 11:
        run_details['allele'] = '-'.join(components[0:2])+components[2]
        run_details['peptide'] = components[3]
        run_details['uid'] = components[4]
        run_details['recycles'] = re.findall(r'\d+', components[5])[0]
        run_details['relaxed'] = relaxed(components[6])
        run_details['rank'] = components[8]
        run_details['model'] = components[10].replace('.pdb','')
    else:
        run_details['allele'] = '-'.join(components[0:2])
        run_details['peptide'] = components[2]
        run_details['uid'] = components[3]
        run_details['recycles'] = re.findall(r'\d+', components[4])[0]
        run_details['relaxed'] = relaxed(components[5])
        run_details['rank'] = components[7]
        run_details['model'] = components[9].replace('.pdb','')
    return run_details



def get_run_file_dict(locus:str, allele:str, peptide:str, input_folder:str, file_root:str) -> Dict:
    """
    This function builds a set of run folders and run files for a specific allele and peptide and returns a dictionary of them
    
    Args:
        locus (str): the locus of the complex e.g. 'HLA-A'
        allele (str): the allele of the complex e.g. 'HLA-A*68:02'
        peptide (str): the peptide in the complex e.g. 'TLTSCNTSV'
        input_folder (str): the input folder for the set e.g. 'unrelaxed_split_peptide'
        file_root (str): the root directory for all the files e.g. '../'

    Returns:
        Dict (run_file_dict): the dictionary of directories and files for the run
    """
    run_file_dict = {}
    run_folders = [folder for folder in os.listdir(build_filepath(input_folder, locus, allele, peptide, file_root)) if '.result' in folder]
    for recycle_folder in run_folders:
        run_files = [file for file in os.listdir(build_runpath(recycle_folder, input_folder, locus, allele, peptide, file_root)) if '.pdb' in file]
        run_file_dict[recycle_folder] = run_files
    return run_file_dict



def act_on_set(locus:str, allele:str, peptide:str, action_dict:Dict, file_root:str):
    """
    This function builds a set of run folders for a specific allele and peptide and then creates new folders if required and runs a specific functionn
    
    Args:
        locus (str): the locus of the complex e.g. 'HLA-A'
        allele (str): the allele of the complex e.g. 'HLA-A*68:02'
        peptide (str): the peptide in the complex e.g. 'TLTSCNTSV'
        action_dict (Dict): the dictionary for a specific action. It will contain an input_folder, output_folder and an action (function) to perform
        file_root (str): the root directory for all the files e.g. '../'


    """
    # first build a set of run folders for a specific allele and peptide
    run_folders = [folder for folder in os.listdir(build_filepath(action_dict['input_folder'], locus, allele, peptide, file_root)) if '.result' in folder]
    print (run_folders)


    for run_folder in run_folders:
        original_filepath = build_runpath(run_folder, action_dict['input_folder'], locus, allele, peptide, file_root)
        print (original_filepath)
        try:
            transformed_filepath = build_runpath(run_folder, action_dict['output_folder'], locus, allele, peptide, file_root)
            os.makedirs(transformed_filepath)
            print('directory created')
        except:
            print('directory exists')
        pdb_files = [file for file in os.listdir(build_runpath(run_folder, action_dict['input_folder'], locus, allele, peptide, file_root)) if '.pdb' in file]
        for pdb_file in pdb_files:
            action_dict['action'](pdb_file, original_filepath, transformed_filepath)
            print (pdb_file)
        print ('-------')