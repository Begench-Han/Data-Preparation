# from igfold import IgFoldRunner, init_pyrosetta
from igfold import IgFoldRunner
import glob
from anarci import anarci
from anarci import number
from numpy import integer
import torch

import pandas as pd
pdb_names = []
prmsd = []
method = []
csv_data = []
# init_pyrosetta()
# biounit_names = glob.glob('/Users/begenchhangeldiyev/Desktop/dataprep/all_structures_new/imgt/'+'*.pdb')
neww= glob.glob('/Users/begenchhangeldiyev/Desktop/dataprep/alignments/'+'*.fa')
sequences = {
    "H": "EVQLVQSGPEVKKPGTSVKVSCKASGFTFMSSAVQWVRQARGQRLEWIGWIVIGSGNTNYAQKFQERVTITRDMSTSTAYMELSSLRSEDTAVYYCAAPYCSSISCNDGFDIWGQGTMVTVS",
    "L": "DVVMTQTPFSLPVSLGDQASISCRSSQSLVHSNGNTYLHWYLQKPGQSPKLLIYKVSNRFSGVPDRFSGSGSGTDFTLKISRVEAEDLGVYFCSQSTHVPYTFGGGTKLEIK"
}
new_seq = {}
error = []
nna = ''
pred_pdb = "my_antibody.pdb"
igfold = IgFoldRunner()
for name in neww:
    # print(name)
    a = 0
    new_seq ={}
    b= 1
    print(name)
    for line in open(name,"rb"):
        line = line.decode("utf-8","ignore").rstrip()
        
        a+=1
        if a%2 == 0:
            final =(name[54:58] , line)
            ufo = []
            print(final)
            ufo.append(final)
            res = anarci(ufo, scheme="imgt", output=False)
            # print(res)
            numbering, alignment_details, hit_tables = res
            # print(numbering)
            if len(numbering[0]) == None:
                error.append([name[54:58],b ])
                continue
            for j in range(len(numbering[0])):
                domain_numbering, start_index, end_index = numbering[0][j]
                if alignment_details[0][j]["chain_type"] =="H":
                    new_seq["H"] = line[start_index:end_index]
                else:
                    new_seq["L"] = line[start_index:end_index]
            print(new_seq)
            print(b)
            if a ==2:
                bhh = name[54:58]+'_igfold_native' +'.pdb'
            else:
                bhh = name[54:58]+'_igfold_'+str(b)+".pdb"
            out = igfold.fold(
            bhh, # Output PDB file
            sequences=new_seq, # Antibody sequences
            do_refine=True, # Refine the antibody structure with PyRosetta
            do_renum=False, # Renumber predicted antibody structure (Chothia)
        )
            rmsd = out.prmsd
            rmsd = torch.flatten(rmsd)
            mean = torch.mean(rmsd)
            mean = torch.IntTensor.item(mean)
            print(mean)
            if a ==2:
                nna = name[54:58] + '_'+'native'
            else:
                nna = name[54:58] + '_'+str(b)
            pdb_names.append(nna)
            prmsd.append(rmsd)
            csv_data.append([nna ,mean, 'ingraham' ])
            if a !=2:
                b+=1
            print(b)
        if a==20:
            break
    break
print(nna)
print(prmsd)
error_seq = pd.DataFrame(error, columns=['name', 'number'])
error_seqfile = 'error_seq.csv'
error_seq.to_csv(error_seqfile, index=False)
df_seq = pd.DataFrame(csv_data, columns=['PDB_NAMES', 'RMSD', 'METHOD'])
csv_seqfile = 'df_seq.csv'
df_seq.to_csv(csv_seqfile, index=False)
# print(new_seq)