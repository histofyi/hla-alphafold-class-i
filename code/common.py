from typing import Dict, List

import os

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