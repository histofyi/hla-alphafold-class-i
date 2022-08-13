import csv
import json

from common import locus_from_allele

file_root = '../'

pdb_code = '1HHG'


models = [model for model in range(1,6)]
recycles = [recycle_count for recycle_count in range(3,11)]


def get_row_labels():
    info_row_labels = ['allele_and_peptide','pdb_code', 'length']

    row_labels = info_row_labels

    data_row_labels = ['rank', 'lpddt', 'main_chain_rmsd', 'all_atom_rmsd', 'side_chain_rmsd']

    

    for model in models:
        for item in data_row_labels:
            row_labels.append(f'model_{model}_{item}')
    print (row_labels)
    return row_labels


recycle_counts = {}
for recycle_count in recycles:
    recycle_counts[str(recycle_count)] = {'pdb_codes':[], 'rmsds':[], 'alleles':[], 'peptides':[], 'ranks':[], 'plddts':[]}

all_runs = {}
all_statistics = {}

def process_row(dataset_row_labels, run_data, alpha_fold_row):
    min_rmsd = 100
    min_rmsd_recycle_count = 0
    min_model_number = 0
    min_rank = 0
    min_plddt = 0
    prev_best_model = None
    prev_min_rmsd = None
    rmsd_delta = None
    max_delta = 0
    first_rmsd = 0
    all_runs[alpha_fold_row['pdb_code']] = {'runs':{}, 'metadata':{}, 'max_delta':{'delta':0, 'recycle':3}}
    all_runs[alpha_fold_row['pdb_code']]['metadata'] = {
            'allele': alpha_fold_row['allele'],
            'peptide': alpha_fold_row['peptide'],
            'uid': alpha_fold_row['uid'],
            'pdb_code': alpha_fold_row['pdb_code'],
            'locus': locus_from_allele(alpha_fold_row['allele'])
        }
    
    for recycle_count in recycles:
        print (recycle_count)
        run_array = []
        row_min_rmsd = 100
        row_min_model_number = 0
        row_min_rank = 0
        row_min_plddt = 0
        prev_best_model = None
        best_model_changed = None
        this_run = run_data['recycles'][str(recycle_count)]
        for model in this_run['models']:
            model_number = model['model']
            if model['predicted_structure_data']['all_positions']['all_atoms_rmsd'] < min_rmsd:
                min_rmsd = model['predicted_structure_data']['all_positions']['all_atoms_rmsd']
                min_rmsd_recycle_count = recycle_count
                min_model_number = model_number
                min_rank = model['rank']
                min_plddt = model['predicted_structure_data']['all_positions']['plddt']
            if model['predicted_structure_data']['all_positions']['all_atoms_rmsd'] < row_min_rmsd:
                row_min_rmsd = model['predicted_structure_data']['all_positions']['all_atoms_rmsd']
                row_min_model_number = model_number
                row_min_rank = model['rank']
                row_min_plddt = model['predicted_structure_data']['all_positions']['plddt']
        if prev_min_rmsd:
            row_rmsd_delta = prev_min_rmsd - row_min_rmsd
            if prev_best_model is not None and prev_best_model != row_min_model_number:
                best_model_changed = True
            else:
                best_model_changed = False
        else:
            row_rmsd_delta = None
        prev_min_rmsd = row_min_rmsd
        all_runs[alpha_fold_row['pdb_code']]['runs'][str(recycle_count)] = {
            'rmsd': row_min_rmsd,
            'rank': row_min_rank,
            'model': row_min_model_number,
            'plddt' : row_min_plddt,
            'model_changed': best_model_changed,
            'rmsd_delta': row_rmsd_delta
        }
        print (row_rmsd_delta)
        if recycle_count == 3:
            first_rmsd = row_min_rmsd
            current_delta = 0
        else:
            current_delta = first_rmsd - row_min_rmsd
        if current_delta > max_delta:
            max_delta = current_delta
            all_runs[alpha_fold_row['pdb_code']]['max_delta']['delta'] = max_delta
            all_runs[alpha_fold_row['pdb_code']]['max_delta']['recycle'] = recycle_count

    print (min_rmsd)
    print (min_rmsd_recycle_count)
    print (min_model_number)
    print (min_rank)
    recycle_counts[str(min_rmsd_recycle_count)]['pdb_codes'].append(alpha_fold_row['pdb_code'])
    recycle_counts[str(min_rmsd_recycle_count)]['rmsds'].append(min_rmsd)
    recycle_counts[str(min_rmsd_recycle_count)]['plddts'].append(min_plddt)
    recycle_counts[str(min_rmsd_recycle_count)]['ranks'].append(min_rank)
    recycle_counts[str(min_rmsd_recycle_count)]['alleles'].append(alpha_fold_row['allele'])
    recycle_counts[str(min_rmsd_recycle_count)]['peptides'].append(alpha_fold_row['peptide'])




# open the CSV file
file = open(f'{file_root}human_class_i.csv')
# and read in the CSV 
csvreader = csv.reader(file)

# set the header to the first row
header = next(csvreader)

# iterate through the remaining rows
rows = []
dataset_row_labels = get_row_labels() 
for row in csvreader:
    alpha_fold_row = {k:v for k,v in list(zip(header,row))}
    if len(alpha_fold_row['pdb_file_issues']) == 0:
        pdb_code = alpha_fold_row['pdb_code']
        allele = alpha_fold_row['allele']
        peptide = alpha_fold_row['peptide']
        locus = locus_from_allele(allele)
        run_datafile_name = f'{file_root}{locus.lower()}/{allele.lower()}/statistics/{peptide.lower()}.json'
        print (run_datafile_name)
        try:
            with open(run_datafile_name, 'r')as run_datafile:
                run_data = json.load(run_datafile)
            
        except:
            print ('no_statistics')
        all_statistics[alpha_fold_row['pdb_code']] = run_data
        process_row(dataset_row_labels, run_data, alpha_fold_row)





quality_bins = [5, 3, 2, 1.5, 1.25]
quality_labels = ['poor', 'neutral', 'good', 'very_good', 'excellent']

quality_set = {
    'very_poor':[]
}

for quality_label in quality_labels:
    quality_set[quality_label] = []



#print (all_runs)

for recycle_count in recycle_counts:
    i = 0
    average = sum(recycle_counts[recycle_count]['rmsds']) / len(recycle_counts[recycle_count]['rmsds'])
    print (f'\nRecycle count {recycle_count} - [{len(recycle_counts[recycle_count]["rmsds"])} structures] - avg_best_rmsd = {round(average, 2)}')
    print ('Set members')
    for rmsd in recycle_counts[recycle_count]['rmsds']:
        this_pdb_code = recycle_counts[recycle_count]["pdb_codes"][i]
        print (f'{this_pdb_code} {recycle_counts[recycle_count]["alleles"][i]}/{recycle_counts[recycle_count]["peptides"][i]}: {rmsd} / {recycle_counts[recycle_count]["ranks"][i]} / {recycle_counts[recycle_count]["plddts"][i]}')
        i += 1
        j = 0
        for quality_bin in quality_bins:
            if rmsd > quality_bin and j ==0:
                this_quality = 'very_poor'
            elif rmsd <= quality_bin:
                this_quality = quality_labels[j]
            j += 1
        quality_set[this_quality].append(this_pdb_code)

        print (f'{this_pdb_code} - prediction quality : {this_quality}')            



    print ('Statistics')

good_to_excellent = []
for this_quality in ['good', 'very_good', 'excellent']:
    for item in quality_set[this_quality]:
        good_to_excellent.append(item)

print (good_to_excellent)


good_to_excellent_runs = []

k = 0
for row in all_runs:
    for this_quality in ['good', 'very_good', 'excellent']:
        if all_runs[row]['metadata']['pdb_code'] in quality_set[this_quality]:
            print ('')
            print (all_runs[row])
            all_runs[row]['metadata']['quality'] = this_quality
            good_to_excellent_runs.append(all_runs[row])
            print ('')
    k += 1


with open('../statistics/all_runs.json', 'w') as output_file:
    json.dump(all_runs, output_file, sort_keys = True, indent = 4, ensure_ascii = True)

with open('../statistics/all.json', 'w') as output_file:
    json.dump(all_statistics, output_file, sort_keys = True, indent = 4, ensure_ascii = True)

with open('../statistics/good_to_excellent.json', 'w') as output_file:
    json.dump(good_to_excellent_runs, output_file, sort_keys = True, indent = 4, ensure_ascii = True)

with open('../statistics/quality.json', 'w') as output_file:
    json.dump(quality_set, output_file, sort_keys = True, indent = 4, ensure_ascii = True)
