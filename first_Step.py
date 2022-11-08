import argparse
from dateutil import parser
import numpy as np
import os, time, gzip, json
import glob 
from iteration_utilities import unique_everseen
import difflib
import random
import json
# import difflib
import glob
biounit_names = glob.glob('/Users/begenchhangeldiyev/Desktop/dataprep/all_structures_new/imgt/'+'*.pdb')
neww= glob.glob('/Users/begenchhangeldiyev/Downloads/all_structures/imgt/'+'*.pdb')
post = []
for b in neww:
    post.append(b[56:60])
#print(post)
# print(neww)
pre = []
for i in biounit_names:
    pre.append(i[63:67])
last = list(set(post)-set(pre))
diff_ = []
for i in last:
    diff_.append('/Users/begenchhangeldiyev/Downloads/all_structures/imgt/' + i + '.pdb')
# print(hm)
#folder_with_pdbs_path = args.input_path
#save_path = args.output_path

alpha_1 = list("ARNDCQEGHILKMFPSTWYV")
states = len(alpha_1)
alpha_3 = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
            'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']

aa_1_N = {a:n for n,a in enumerate(alpha_1)}
aa_3_N = {a:n for n,a in enumerate(alpha_3)}
aa_N_1 = {n:a for n,a in enumerate(alpha_1)}
aa_1_3 = {a:b for a,b in zip(alpha_1,alpha_3)}
aa_3_1 = {b:a for a,b in zip(alpha_1,alpha_3)}

def AA_to_N(x):
    # ["ARND"] -> [[0,1,2,3]]
    x = np.array(x);
    if x.ndim == 0: x = x[None]
    return [[aa_1_N.get(a, states-1) for a in y] for y in x]

def N_to_AA(x):
    # [[0,1,2,3]] -> ["ARND"]
    x = np.array(x);
    if x.ndim == 1: x = x[None]
    return ["".join([aa_N_1.get(a,"-" ) for a in y]) for y in x]

def parse_PDB_biounits(x, atoms=['N','CA','C'], chain=None):
    '''
    input:  x = PDB filename
            atoms = atoms to extract (optional)
    output: (length, atoms, coords=(x,y,z)), sequence
    '''
    xyz,seq,min_resn,max_resn = {},{},1e6,-1e6
    for line in open(x,"rb"):
        line = line.decode("utf-8","ignore").rstrip()
        if line[:6] == "HETATM" and line[17:17+3] == "MSE":
            line = line.replace("HETATM","ATOM  ")
            line = line.replace("MSE","MET")

        if line[:4] == "ATOM":
            ch = line[21:22]
            if ch == chain or chain is None:
                atom = line[12:12+4].strip()# N, Ca, C , O
                resi = line[17:17+3]#name of amino acid
                resn = line[22:22+5].strip()#number of resiude
                x,y,z = [float(line[i:(i+8)]) for i in [30,38,46]]# coordinates

                if resn[-1].isalpha(): 
                    resa,resn = resn[-1],int(resn[:-1])-1
                else: 
                    resa,resn = "",int(resn)-1
        #         resn = int(resn)
                if resn < min_resn: 
                    min_resn = resn
                if resn > max_resn: 
                    max_resn = resn
                if resn not in xyz: 
                    xyz[resn] = {}
                if resa not in xyz[resn]: 
                    xyz[resn][resa] = {}
                if resn not in seq: 
                    seq[resn] = {}
                if resa not in seq[resn]: 
                    seq[resn][resa] = resi

                if atom not in xyz[resn][resa]:
                    xyz[resn][resa][atom] = np.array([x,y,z])
                

    # convert to numpy arrays, fill in missing values
    seq_,xyz_ = [],[]
    try:
            for resn in range(min_resn,max_resn+1):
                if resn in seq:
                    for k in sorted(seq[resn]): seq_.append(aa_3_N.get(seq[resn][k],20))
                # else: seq_.append(20)
                if resn in xyz:
                    for k in sorted(xyz[resn]):
                        for atom in atoms:
                            if atom in xyz[resn][k]: xyz_.append(xyz[resn][k][atom])
                            else: 
                                xyz_.append(np.full(3,np.nan))
                                # print(x)
                # else:
                #     for atom in atoms: xyz_.append(np.full(3,np.nan))
            return np.array(xyz_).reshape(-1,len(atoms),3), N_to_AA(np.array(seq_))
    except TypeError:
        return 'no_chain', 'no_chain'
init_alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G','H', 'I', 'J','K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T','U', 'V','W','X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g','h', 'i', 'j','k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't','u', 'v','w','x', 'y', 'z']
extra_alphabet = [str(item) for item in list(np.arange(300))]
chain_alphabet = init_alphabet + extra_alphabet
folder_with_pdbs_path = '/Users/begenchhangeldiyev/Downloads/all_structures/imgt/'
# biounit_names = glob.glob(folder_with_pdbs_path+'*.pdb')
# name = "/Users/begenchhangeldiyev/Desktop/dataprep/all_structures/imgt/3okk.pdb"
pdb_dict_list = []
aaaa=0
print(biounit_names)
for name in biounit_names:
    # print(name)
    print(aaaa)
    aaaa+=1
    d = 1
    sequence = []
    
    print(name)
    for line in open(name,"rb"):
        line = line.decode("utf-8","ignore").rstrip()
        concat_seq = ''
        dica = {}
        liiis = []
        hlist = []
        
        s = 0
        concat_N = []
        concat_CA = []
        concat_C = []
        concat_O = []
        concat_mask = []
        coords_dict = {}
        coords_dict_chain = {}
        my_dict = {}
        seq_list = ''
        if line[:6] == "REMARK":
            if line[11:20] == "PAIRED_HL":
                hlist=([line[line.find("HCHAIN")+ 7],line[line.find("LCHAIN")+ 7]])
            elif line[11:17] == "SINGLE":
                if line.find("HCHAIN") != -1:
                    hlist = line[line.find("HCHAIN")+ 7]
                else: 
                    hlist = line[line.find("LCHAIN")+ 7]
            for letter in hlist:
                xyz, seq = parse_PDB_biounits(name, atoms=['N','CA','C','O'], chain=letter)
                if type(xyz) != str:
                    s += 1
                    concat_seq +=seq[0]
                    concat_N += xyz[:,0,:].tolist()
                    concat_CA += xyz[:,1,:].tolist()
                    concat_C +=xyz[:,2,:].tolist()
                    concat_O += xyz[:,3,:].tolist()
                    if concat_seq in sequence:
                        hlist = []
                    if line[11:20] == "PAIRED_HL":
                        if '/' not in concat_seq:
                            concat_seq += "/"
            if concat_seq != '':
                sequence.append(concat_seq)
            if len(hlist)!=0:   
                if concat_seq != '':
                    coords_dict_chain['N'] = concat_N
                    coords_dict_chain['CA'] = concat_CA
                    coords_dict_chain['C'] = concat_C
                    coords_dict_chain['O'] = concat_O
                    my_dict['seq'] = concat_seq
                    my_dict['coords']=coords_dict_chain
                    my_dict['num_of_chains'] = s

                    fi = name.rfind("/")
                    # my_dict['name']=name[(fi+1):-4]+ str(d) ## For chains
                    my_dict['name']=name[(fi+1):-4]
                    
                    pdb_dict_list.append(my_dict)
                    d+=1
        if concat_seq != '':
            break
    # break

    # if aaaa == 50:
    #     break
# print((pdb_dict_list))
# print(pdb_dict_list)
name_list = []
for i in range(len(pdb_dict_list)):
    #  print(pdb_dict_list[i]['seq'])
    name_list.append(pdb_dict_list[i]['name'])
random.shuffle(name_list)
validation_name = name_list[:150]
test_name = name_list[150:300]
train_name = name_list[300:]
final = {}
final["test"] = test_name
final["train"] = train_name
final["validation"] = validation_name
# with open("/Users/begenchhangeldiyev/Desktop/dataprep/splits.json" , "w") as outf:
#     json.dump(final,outf)
# a= pdb_dict_list[72]["seq"]
# b= pdb_dict_list[73]["seq"]
# print(a, len(a))
# print(len(pdb_dict_list[72]['coords']['N']))
# print('')
# print(b, len(b))
# for i,s in enumerate(difflib.ndiff(a, b)):
#     if s[0]==' ': continue
#     elif s[0]=='-':
#         print(u'Delete "{}" from position {}'.format(s[-1],i))
#     elif s[0]=='+':
#         print(u'Add "{}" to position {}'.format(s[-1],i)) 
# m = [li for li in difflib.ndiff(a, b) if li[0] != ' ']

# for i in range(50):
#     for j in range(50):

# new = [i for n, i in enumerate(pdb_dict_list) if i not in pdb_dict_list[n + 1:]]
# pdb_dict_list = (unique_everseen(pdb_dict_list['seq']))
# print(pdb_dict_list)
# set_of_jsons = {json.dumps(d, sort_keys=True) for d in pdb_dict_list}
# X = [json.loads(t) for t in set_of_jsons]
# print(len(X))
# hm = set()
# new = []
# for d in pdb_dict_list: 
#     t = tuple(d.items())
#     if t not in hm:
#         hm.add(t)

#         new.append(d)
# print(new)
# for i in pdb_dict_list:
#     for a,b in i:
# pdb_dict_list = [dict(t) for t in {tuple(d.items()) for d in pdb_dict_list}]
#print((pdb_dict_list))
# print(my_dict['seq']) 
with open("/Users/begenchhangeldiyev/Desktop/dataprep/jsonl/new_hm.jsonl", 'w') as f:
        for entry in pdb_dict_list:
            f.write(json.dumps(entry) + '\n')
