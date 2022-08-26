import json

from common import load_list_from_csv, locus_from_allele

file_root = '../'

all_stats = {}
diffs = {}

qualities = {
    'excellent':[],
    'very_good':[],
    'good':[],
    'neutral':[],
    'poor':[],
    'very_poor':[]
}

rmsds = {}


alpha_fold_list, success, errors =  load_list_from_csv('human_class_i.csv')

for alpha_fold_row in alpha_fold_list:
    pdb_code = alpha_fold_row['pdb_code']
    allele = alpha_fold_row['allele']
    locus = locus_from_allele(allele)
    peptide = alpha_fold_row['peptide']
    statistics_directory = f'{locus.lower()}/{allele.lower()}/statistics'
    file_name = f'{statistics_directory}/{peptide.lower()}.json'
    file_path = f'{file_root}{file_name}'
    try:
        with open(file_path, 'r') as input_file:
            run_data = json.load(input_file)
    except:
        run_data = None
    if run_data:
        all_stats[pdb_code] = {
            'best':run_data['best'],
            'metadata':run_data['metadata']
        }
        quality = run_data['best']['quality']
        rmsd = round(run_data['best']['allatom_rmsd'],1)
        qualities[quality].append(pdb_code)
        if rmsd not in rmsds:
            rmsds[rmsd] = []
        rmsds[rmsd].append(pdb_code)
        print ('------')
        print (pdb_code)
        diffs[pdb_code] = {'recycles':{}}
        for recycle in run_data['recycles']:
            diffs[pdb_code]['recycles'][recycle] = {}
            recycle_max_rmsd = 0
            recycle_min_rmsd = 100
            recycle_max_plddt = 0
            recycle_min_plddt = 100
            print (recycle)
            for model in run_data['recycles'][recycle]['models']:
                if run_data['recycles'][recycle]['models'][model]['allatom_rmsd'] < recycle_min_rmsd:
                    recycle_min_rmsd = run_data['recycles'][recycle]['models'][model]['allatom_rmsd']
                if run_data['recycles'][recycle]['models'][model]['allatom_rmsd'] > recycle_max_rmsd:
                    recycle_max_rmsd = run_data['recycles'][recycle]['models'][model]['allatom_rmsd']
                if run_data['recycles'][recycle]['models'][model]['plddt'] > recycle_max_plddt:
                    recycle_max_plddt = run_data['recycles'][recycle]['models'][model]['plddt']
                if run_data['recycles'][recycle]['models'][model]['plddt'] < recycle_min_plddt:
                    recycle_min_plddt = run_data['recycles'][recycle]['models'][model]['plddt']
            diffs[pdb_code]['recycles'][recycle]['plddt_min'] = recycle_min_plddt
            diffs[pdb_code]['recycles'][recycle]['plddt_max'] = recycle_max_plddt
            diffs[pdb_code]['recycles'][recycle]['plddt_diff'] = round(recycle_max_plddt - recycle_min_plddt, 2)
            diffs[pdb_code]['recycles'][recycle]['rmsd_min'] = recycle_min_rmsd
            diffs[pdb_code]['recycles'][recycle]['rmsd_max'] = recycle_max_rmsd
            diffs[pdb_code]['recycles'][recycle]['rmsd_diff'] = round(recycle_max_rmsd - recycle_min_rmsd, 2)
            
            print (quality)
            print (diffs[pdb_code]['recycles'][recycle])



with open('../statistics/all.json', 'w') as output_file:
    json.dump(all_stats, output_file, sort_keys = True, indent = 4, ensure_ascii = True)
with open('../statistics/rmsds.json', 'w') as output_file:
    json.dump(rmsds, output_file, sort_keys = True, indent = 4, ensure_ascii = True)
with open('../statistics/qualities.json', 'w') as output_file:
    json.dump(qualities, output_file, sort_keys = True, indent = 4, ensure_ascii = True)


    