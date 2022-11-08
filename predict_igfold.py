# from igfold import IgFoldRunner, init_pyrosetta
from igfold import IgFoldRunner
import glob
# init_pyrosetta()
# biounit_names = glob.glob('/Users/begenchhangeldiyev/Desktop/dataprep/all_structures_new/imgt/'+'*.pdb')
neww= glob.glob('/Users/begenchhangeldiyev/Desktop/my_results/seqs/'+'*.fa')
sequences = {
    "H": "EVQLVQSGPEVKKPGTSVKVSCKASGFTFMSSAVQWVRQARGQRLEWIGWIVIGSGNTNYAQKFQERVTITRDMSTSTAYMELSSLRSEDTAVYYCAAPYCSSISCNDGFDIWGQGTMVTVS",
    "L": "DVVMTQTPFSLPVSLGDQASISCRSSQSLVHSNGNTYLHWYLQKPGQSPKLLIYKVSNRFSGVPDRFSGSGSGTDFTLKISRVEAEDLGVYFCSQSTHVPYTFGGGTKLEIK"
}
new_seq = {}
pred_pdb = "my_antibody.pdb"
igfold = IgFoldRunner()
# print(biounit_names)
print(neww)
# with open("/Users/begenchhangeldiyev/Desktop/my_results/seqs/7fgj1.fa", 'r') as hm:
#     print(hm)
for name in neww:
    a = 0
    print(name)
    for line in open(name,"rb"):
        line = line.decode("utf-8","ignore").rstrip()
        a+=1
        if a==4:
            if '/' in line:
                hm = line.find('/')
                print(hm)
                print(name)
                new_seq["H"] = line[:hm]
                new_seq["L"] = line[hm+1:]
                # print(new_seq)
            else:
                new_seq["H"] = line
                # print(new_seq)
            igfold.fold(
            name[50:54]+'_igfold'+".pdb", # Output PDB file
            sequences=new_seq, # Antibody sequences
            do_refine=True, # Refine the antibody structure with PyRosetta
            do_renum=False, # Renumber predicted antibody structure (Chothia)
        )
    
    
# print(new_seq)





