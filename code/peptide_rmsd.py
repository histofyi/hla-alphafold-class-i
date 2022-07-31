from biopandas.pdb import PandasPdb
import requests

url = 'https://raw.githubusercontent.com/histofyi/hla-alphafold-class-i/main/hla-a/hla-a*0101/unrelaxed/AQDIYRASYY/HLA_A0101_AQDIYRASYY_ecc32_5un.result/HLA_A0101_AQDIYRASYY_ecc32_5un_unrelaxed_rank_1_model_1.pdb'

r = requests.get(url)
print (r.status_code)
print (r.headers['content-type'])
print (r.text)

ppdb = PandasPdb()
pdb_data = r.text.splitlines()
ppdb.read_pdb_from_list(pdb_data)

print('PDB Code: %s' % ppdb.code)
print('PDB Header Line: %s' % ppdb.header)

print (ppdb.df.keys())
print (ppdb.df['ATOM'].head())

print (ppdb.df['ATOM'].groupby(['chain_id']))

