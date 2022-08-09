import csv

from common import act_on_set, locus_from_allele
from functions import split_peptide


file_root = '../'

actions = {
    'split_peptide':{
        'input_folder': 'unrelaxed_aligned',
        'output_folder': 'unrelaxed_split_peptide',
        'action': split_peptide
    },
}


allele = 'HLA-A*6802'
action = 'split_peptide'

locus = locus_from_allele(allele)




if locus and action in actions:
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
        if allele in row:
            rows.append(row)

    print (rows)

    # next create a dictionary of the complexex relating to a specific allele
    complex_set = []
    for row in rows:
        complex_set.append({k:v for k,v in list(zip(header,row))})

    # iterate through the complexes relating to that allele
    for complex in complex_set:
        print (complex['peptide'])
        # and perform an action on them
        act_on_set(locus, allele, complex['peptide'], actions[action], file_root)

