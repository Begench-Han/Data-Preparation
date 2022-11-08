import os, glob
import os.path
import glob
import pandas as pd
import sys
import inspect

def align_structures(pdb1, pdb2, file1, file2):
    """Take two structure and superimpose pdb1 on pdb2"""
    import Bio.PDB
    import subprocess

    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    # Get the structures
    ref_structure = pdb_parser.get_structure("ref", pdb1)
    sample_structure = pdb_parser.get_structure("sample", pdb2)

    aligner = Bio.PDB.CEAligner()
    aligner.set_reference(ref_structure)
    aligner.align(sample_structure)

    # file1name = 'outputs/' + file1
    # file2name = 'outputs/' + file2
    # io = Bio.PDB.PDBIO()
    # io.set_structure(ref_structure)
    # io.save(file1name)
    # io.set_structure(sample_structure)
    # io.save(file2name)

    return aligner.rms, "outputs/reference.pdb", f"outputs/af_prediction_aligned.pdb"
    
file_path = os.path.realpath(__file__)
k = file_path.rfind("/")
# Loop through all the files in /native with extension *.pdb
ext = ('.pdb')
pdb_path = "/align_data/native"
mpnn_path = "/align_data/mpnn_igfold"
ing_path = "/align_data/new_ing"
print(f'Path with PDBs to test: {pdb_path}')
# iterating over all files
pdbCounter = len(glob.glob1(pdb_path,"*.pdb"))
print(f'Total of PDBs to test: {pdbCounter}')
mpnn_results = []
ing_results = []
pdb_names = []
file_ids = []
csv_data = []
file_id = 1
sequence_id = 1
unique_id = 1
for file in os.listdir(pdb_path):
    if not file.endswith(ext):
            continue
    else:
        print(file_id)
        for i in range(50):
            # print(i)
            sequence_id = i+1
            split_id = file.rfind("_")
            mpnn_file = file[:split_id] + '_mpnn_igfold_' + str(sequence_id) + '.pdb'
            ing_file = file[:split_id] + '_igfold_' + str(sequence_id) + '.pdb'
            # print(ing_path + '/' + ing_file)
            # print(mpnn_path + '/' + mpnn_file)
            if os.path.isfile(mpnn_path + '/' + mpnn_file) and os.path.isfile(ing_path + '/' + ing_file):  
                rmsd_mpnn, input_pdb, aligned_pdb = align_structures(pdb_path + '/' + file, mpnn_path + '/' + mpnn_file, file, mpnn_file)
                rmsd_ing, input_pdb, aligned_pdb = align_structures(pdb_path + '/' + file, ing_path + '/' + ing_file, file, ing_file)
                pdb_names.append(file[:split_id])
                file_ids.append(file_id)
                mpnn_results.append(rmsd_mpnn)
                ing_results.append(rmsd_ing)
                csv_data.append([unique_id, file_id, file[:split_id], sequence_id, rmsd_mpnn, rmsd_ing])
                unique_id += 1
        file_id += 1

# save sequence generation data to dataframe for further analysis
df_seq = pd.DataFrame(csv_data, columns=['unique_id', 'file_id', 'pdb_name', 'sequence_id', 'mpnn_rmsd', 'ing_rmsd'])
csv_seqfile = 'results.csv'
df_seq.to_csv(csv_seqfile, index=False)