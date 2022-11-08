from operator import index
from tkinter import E
import Bio.PDB as bpdb
import io
import json
from anarci import number
from anarci import anarci

class ResSelect(bpdb.Select):
    def accept_residue(self, res):
        if len(index_list) ==2:
            for i in index_list:
                if res.id[1] >= i[0] and res.id[1] <= i[1]:
                    if res.parent.id in key_list:
                        # print("dogry")
                        return True
                else:
                    # print("yalnysh")
                    return False
        elif len(index_list) ==1:
            if len(index_list[0])>2:
                if res.id[1] >= index_list[0][0] and res.id[1] <= index_list[0][1]: 
                    if res.parent.id in key_list:
                        return True
                    else:
                        return False
                if res.id[1] >= index_list[0][2] and res.id[1] <= index_list[0][3]:
                    if res.parent.id in key_list:
                        return True
                    else:
                        return False
            elif res.id[1] >= index_list[0][0] and res.id[1] <= index_list[0][1]:
                if res.parent.id in key_list:
                    # print("dogry-1")
                    return True
                else:
                    return False
                
                     
            else:
                # print("yalnysh-1")
                return False
# class ResSelect(bpdb.Select):
#     def accept_residue(self, res, start, end, id):
#         if res.id[1] >= start and res.id[1] <= end and res.parent.id == id:
#             return True
#         else:
#             return False
#     def accept_chain(self, chain):
#         if chain:
#             return 1
#         else:
#             return 0
io = bpdb.PDBIO()
s = bpdb.PDBParser().get_structure('temp', '/Users/begenchhangeldiyev/Desktop/dataprep/all_structures_new/imgt/1a5f.pdb')
diff = ["7vvo1", "7u9g1", "7t9n1", "7sbg1", "7sk51", "7wjy1", "7noz1", "7qtk1", "7vpy1", "7u0e1", "7wc01", "7um61", "7vgr1", "8dao1", "7vke1", "7skz1", "7upy1", "8dcs1", "7xox1", "7w0n1", "7fgl1", "8dt31", "7ura1", "7qtj1", "7um71", "7xke1", "7v241", "7fgr1", "7unb1", "7um51", "7sl51", "7u0a1", "8df51", "7sd21", "7vgs1", "7wo51", "8cw91", "7t9i1", "7y151", "7utz1", "7xms1", "7re71", "8djk1", "7sk31", "7sk91", "7z3a1", "7qvm1", "8d6z1", "7wp21", "7t9m1", "7wo41", "7xkd1", "7y121", "7urd1", "7woc1", "7y0g1", "7xov1", "7fgk1", "7t7b1", "7upw1", "7u091", "7sk81", "7w0l1", "7sk61", "7upx1", "7ure1", "7vvl1", "7vq01", "7sk71", "7sbd1", "7tcq1", "7vvn1", "7w0p1", "7xt81", "8dcr1", "7sue1", "7xmt1", "7xkf1", "7sk41", "7xou1", "7v231", "7str1", "7xow1", "7wob1", "7wog1", "7qe51", "7wp01", "7qti1", "7v271", "7vvj1", "7sts1", "7woa1", "7w0m1", "7xtp1", "7fjc1", "7urf1", "7xta1", "7w0o1", "8d361", "7u0p1", "7wcd1", "7vvk1", "7fh01", "7xtc1", "7urc1", "7fgj1", "8djm1", "7voa1", "7vvm1", "7re91", "7xtb1", "7wo71", "8dfl1"]
for i in range(len(diff)):
    diff[i] =diff[i][:4]
# print(diff)
diff_ = []
dicti = {}
for i in diff:
    diff_.append('/Users/begenchhangeldiyev/Desktop/dataprep/all_structures_new/imgt/' + i[:4] + '.pdb')
# print(diff_[10][67:71])
with open('/Users/begenchhangeldiyev/Desktop/dataprep/jsonl/new_hm.jsonl', 'r') as json_file:
    json_li = list(json_file)
    
    for l in range(len(json_li)):
        ghg = json.loads(json_li[l])
        # print(ghg['name'])
        if ghg['name'] in diff:
            # if '/' in ghg['seq']:
            #     list_seq = ghg['seq'].split('/')
            #     res_H = anarci(list_seq[0], scheme="imgt", output=False)
            #     numbering_H, alignment_detail_H, hit_tables_H = res_H
            #     res_L = anarci(list_seq[0], scheme="imgt", output=False)
            #     numbering_L, alignment_detail_L, hit_tables_L = res_L
            # else:
            # pla = ghg['seq'].find('/')
            # print(pla)
            # new_seq = ghg['seq'][:pla] + ghg['seq'][pla+1:]
            if '/' in ghg['seq']:
                new_seq = ghg['seq'].replace('/', '')
                new_seq = new_seq.replace('\n',"")
            else:
                new_seq = ghg['seq']
            # print(new_seq)
            # print(ghg['seq'])
            # print(ghg['name'])
            name=ghg['name']+ ':H'
            final = (ghg['name'], new_seq)
            ufo = []
            hlist = ''
            ufo.append(final)
            res = anarci(ufo, scheme="imgt", output=False)
            # print(res)
            numbering, alignment_details, hit_tables = res
            dicti[ghg['name']] = {}
            name_ = '/Users/begenchhangeldiyev/Desktop/dataprep/all_structures_new/imgt/' + ghg['name'] + '.pdb' 
            for line in open(name_,"rb"):
                line = line.decode("utf-8","ignore").rstrip()
                if line[:6] == "REMARK":
                    if line[11:20] == "PAIRED_HL":
                        hlist=line[line.find("HCHAIN")+ 7]
                        for j in range(len(numbering[0])):
                            if alignment_details[0][j]["chain_type"] == 'H':
                                domain_numbering, start_index, end_index = numbering[0][j]
                                dicti[ghg['name']][hlist] = [start_index, end_index ]
                        hlist=line[line.find("LCHAIN")+ 7]
                        for j in range(len(numbering[0])):
                            if alignment_details[0][j]["chain_type"] != 'H':
                                domain_numbering, start_index, end_index = numbering[0][j]
                                place = ghg['seq'].find('/')
                                dicti[ghg['name']][hlist] = [start_index-place, end_index-place ]
                    elif line[11:17] == "SINGLE":
                        if line.find("HCHAIN") != -1:
                            hlist = line[line.find("HCHAIN")+ 7]

                            for j in range(len(numbering[0])):
                                if alignment_details[0][j]["chain_type"] == 'H':
                                    domain_numbering, start_index, end_index = numbering[0][j]
                                    dicti[ghg['name']][hlist] = [start_index, end_index ]
                                else:
                                    domain_numbering, start_index, end_index = numbering[0][j]
                                    dicti[ghg['name']][hlist].append(start_index)
                                    dicti[ghg['name']][hlist].append(end_index)

                        else: 
                            hlist = line[line.find("LCHAIN")+ 7]
                            for j in range(len(numbering[0])):
                                if alignment_details[0][j]["chain_type"] != 'H':
                                    domain_numbering, start_index, end_index = numbering[0][j]
                                    dicti[ghg['name']][hlist] = [start_index, end_index ]
                                else:
                                    print("jenrivens")
                if hlist != '':
                    break
print(dicti)
start_res=20
end_res=30
chain_id = 'H'
hm = ResSelect()
for i in diff_:
    s = bpdb.PDBParser().get_structure('temp', i)
    io.set_structure(s)
    for k,v in dicti.items():
        if i[67:71] == k:
            key_list = list(v.keys())
            index_list = list(v.values())
            print(key_list)
            print(index_list)
            print(k)

            # print(list(v.values()))
            # for kl, vl in v.items():

            # io.save(i[67:71]+'_new.pdb', hm)