from biopandas.pdb import PandasPdb
import csv
import os


from common import locus_from_allele, build_filepath, build_runpath 


pdb_code = '3BP4'

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

if alpha_fold_row:
    allele = alpha_fold_row['allele']
    locus = locus_from_allele(allele)
    peptide = alpha_fold_row['peptide']

    residue_positions = [position + 1 for position in range(0, len(peptide))]


    rmsds['locus'] = locus
    rmsds['allele'] = allele
    rmsds['peptide'] = peptide
    rmsds['relaxed'] = 'unrelaxed'

    rmsds['recycles'] = {}

    for recycle in recycles:
        rmsds['recycles'][str(recycle)] = {'recycle_count': recycle}


    folder_type = 'unrelaxed_split_peptide'
    folder_path = build_filepath(folder_type, locus, allele, peptide, file_root)
    run_folders = [folder for folder in os.listdir(folder_path) if '.result' in folder]

    for run_folder in run_folders:
        run_folder_path = f'{folder_path}/{run_folder}'
        run_files = [run_file for run_file in os.listdir(run_folder_path) if '.pdb' in run_file]
        for run_file in run_files:
            run_file_path = f'{run_folder_path}/{run_file}'
            
            if i < 10:
                predicted_structure_data = {}
                
                test = None
                test = PandasPdb().read_pdb(run_file_path)
                #TODO put in checking that the whole peptide dataframes are the same shape. They may not be in the case of missing atoms in residues
                backbone_rmsd = PandasPdb.rmsd(test.df['ATOM'], real.df['ATOM'], s='main chain') 
                all_atom_rmsd = PandasPdb.rmsd(test.df['ATOM'], real.df['ATOM'], s='heavy')
                print ('----------')
                print (allele)
                print (peptide)
                print (run_file)
                print (f'pymol {real_structure} {run_file_path}')
                print ('----------')

                print (f'Full length - backbone rmsd: {backbone_rmsd}')
                print (f'Full length - all atom rmsd: {all_atom_rmsd}')
                print ('----------')
                for residue_position in residue_positions:

                    test_residue = test.df['ATOM'][test.df['ATOM']['residue_number'] == residue_position]
                    real_residue = real.df['ATOM'][real.df['ATOM']['residue_number'] == residue_position]
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
                        print (f'Residue - plddt: {residue_pldddt}')
                        residue_backbone_rmsd = PandasPdb.rmsd(real_residue, test_residue, s='main chain')
                        residue_all_atom_rmsd = PandasPdb.rmsd(real_residue, test_residue, s='heavy')
                        print (f'Residue - backbone rmsd: {residue_backbone_rmsd}')
                        print (f'Residue - all atom rmsd: {residue_all_atom_rmsd}')                        
                    print ('----------')
            i += 1


    print (rmsds)

