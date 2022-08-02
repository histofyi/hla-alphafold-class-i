# Processing pipeline for the HLA Class I ColabFold structures

1. Align to canonical Antigen Binding Domain (ABD) of 1hhk (same structure used for alignments on histo.fyi) This is performed by the 'align_structures.py' script. NOTE: run this outside of a pipenv shell

2. Re-letter chains, currently B (Class I alpha chain -> A), C (Beta-2 microglobulin ->B), D (peptide -> C).

3. Split into components for anlysis. ABD and peptide.

4. Energy minimisation of the complexes.

5. Split minimised complexes into components for anlysis. ABD and peptide.

6. Compare peptide main chain torsion angles with crystal structure and with the torsion angles from clustering

7. Compare peptide RMSD with other conformers in the prediction run and with crystal structure

8. Look for signals from AlphaFold which would give a sense of accuracy of prediction. See if they correlate to RMSD/torsion angles

