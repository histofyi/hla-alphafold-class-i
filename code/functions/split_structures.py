from pymol import cmd


def split_peptide(filename:str, input_filepath:str, output_filepath:str):
    original_pdb_file = f'{input_filepath}/{filename}'
    print (original_pdb_file)
    name = filename.replace('.pdb','')
    print (name)
    cmd.load(original_pdb_file, quiet=0)
    cmd.alter('chain B', 'chain="A"')
    cmd.alter('chain C', 'chain="B"')
    cmd.alter('chain D', 'chain="C"')
    cmd.remove('chain A')
    cmd.remove('chain B')
    file_types = ['pdb']
    for file_type in file_types:
        file_key = f'{output_filepath}/{name}.{file_type}'
        print (file_key)
        cmd.save(file_key)
    cmd.delete(name)
