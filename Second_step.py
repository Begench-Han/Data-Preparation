from ctypes import alignment
from anarci import anarci
# Format the sequences that we want to number. 
sequences = [ ("7nxb1:H","VQLVESGGGVVQPGRSLRLSCAASAFTFSSYDMHWVRQAPGKGLEWVAVISYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDGGKLVWYYFDYWGQGTLVTVSSASTKGPSVFPLAPSGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKRVEPDIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTLALTFGGGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRG")]
res = anarci(sequences, scheme="imgt", output=False)
        # priint(res)
numbering, alignment_details, hit_tables = res
print(numbering)
print(alignment_details)
print(hit_tables)
import warnings
import re
from anarci import number
from Bio.PDB import PDBParser, PDBIO
from Bio.SeqUtils import seq1
import json
error = []
# ee = [176, 1099, 2034, 2149, 2532, 2644, 2657, 3383, 3465, 3776, 3845, 3917, 4467, 4701, 4891, 5686, 6135]
# ee.sort(reverse=True)
# print(ee)
diff = ["7vvo1", "7u9g1", "7t9n1", "7sbg1", "7sk51", "7wjy1", "7noz1", "7qtk1", "7vpy1", "7u0e1", "7wc01", "7um61", "7vgr1", "8dao1", "7vke1", "7skz1", "7upy1", "8dcs1", "7xox1", "7w0n1", "7fgl1", "8dt31", "7ura1", "7qtj1", "7um71", "7xke1", "7v241", "7fgr1", "7unb1", "7um51", "7sl51", "7u0a1", "8df51", "7sd21", "7vgs1", "7wo51", "8cw91", "7t9i1", "7y151", "7utz1", "7xms1", "7re71", "8djk1", "7sk31", "7sk91", "7z3a1", "7qvm1", "8d6z1", "7wp21", "7t9m1", "7wo41", "7xkd1", "7y121", "7urd1", "7woc1", "7y0g1", "7xov1", "7fgk1", "7t7b1", "7upw1", "7u091", "7sk81", "7w0l1", "7sk61", "7upx1", "7ure1", "7vvl1", "7vq01", "7sk71", "7sbd1", "7tcq1", "7vvn1", "7w0p1", "7xt81", "8dcr1", "7sue1", "7xmt1", "7xkf1", "7sk41", "7xou1", "7v231", "7str1", "7xow1", "7wob1", "7wog1", "7qe51", "7wp01", "7qti1", "7v271", "7vvj1", "7sts1", "7woa1", "7w0m1", "7xtp1", "7fjc1", "7urf1", "7xta1", "7w0o1", "8d361", "7u0p1", "7wcd1", "7vvk1", "7fh01", "7xtc1", "7urc1", "7fgj1", "8djm1", "7voa1", "7vvm1", "7re91", "7xtb1", "7wo71", "8dfl1"]
splits = []
names_to_delete = []
# results = anarci(sequences, scheme="imgt", output=False)
# f = open("/Users/begenchhangeldiyev/Desktop/dataprep/jsonl/hm.jsonl")
with open('/Users/begenchhangeldiyev/Desktop/dataprep/jsonl/hm.jsonl', 'r') as json_file:
    json_li = list(json_file)
    a=  0
    b = 0
    del json_li[4169]
    # length = len(json_li)
    # print(length)
    for l in range(len(json_li)):
        # if l==176:
        #     continue
        ghg = json.loads(json_li[l])
        name=ghg['name']+ ':H'
        # print(ghg['name'])
        # print(ghg['seq'])
        final = (ghg['name'], ghg['seq'])
        # print(ghg)
        ufo = []
        pr = ''
        coor_N = []
        coor_CA = []
        coor_C = []
        coor_O = []
        ufo.append(final)
        # print(ufo)
        # if l == 4169:
        #     print(json_li)
        #     print("fasdfadsasada")
        #     del json_li[l]
        #     continue
        res = anarci(ufo, scheme="imgt", output=False)
        # priint(res)
        numbering, alignment_details, hit_tables = res
        print(numbering)
        print(alignment_details)
        print(hit_tables)

        # print(res)
        if numbering[0] == None:
            b+=1
            error.append(l)
            if ghg['name'] in diff:
                diff.remove(ghg['name'])
            # print(json_li[l])
            # print(l)
            # del json_li[l]
            continue
        for j in range(len(numbering[0])):
            domain_numbering, start_index, end_index = numbering[0][j]
            pr += ghg['seq'][start_index:end_index+1]
            coor_N+=ghg['coords']["N"][start_index:end_index+1]
            coor_CA += ghg['coords']["CA"][start_index:end_index+1]
            coor_C += ghg['coords']["C"][start_index:end_index+1]
            coor_O += ghg['coords']["O"][start_index:end_index+1]
        ghg['seq'] = pr
        ghg['coords']["N"] = coor_N
        ghg['coords']["CA"] = coor_CA
        ghg['coords']["C"] = coor_C
        ghg['coords']["O"] = coor_O
        ghg['num_of_chains']= len(numbering[0])
        # print(len(ghg['seq']))
        # print(len(ghg['coords']["O"]))
        # print((coor))
        json_li[l] = ghg
        a+=1
        if ghg['name'] in diff:
            print("dsnvjnaskbvkjsdabkvbasj")
            continue
        splits.append(ghg['name'])
        print(a)

        # ghg['seq'] = pr
        # if l == 10:
# print(error)
        
# print(len(json_li))
error.sort(reverse=True)
for i in error:
    del json_li[i]
# print(len(json_li))
# print(len(splits))
# print(len(diff))
validation_name = splits[:150]
test_name = splits[150:300]
train_name = splits[300:]
final = {}
final["test"] = test_name
final["train"] = train_name
final["validation"] = validation_name
# print(json_li[:50])
# with open("/Users/begenchhangeldiyev/Desktop/dataprep/final_version.jsonl", 'w') as f:
#         for entry in json_li:
#             f.write(json.dumps(entry) + '\n')
# with open("/Users/begenchhangeldiyev/Desktop/dataprep/final_splits.json" , "w") as outf:
#     json.dump(final,outf)
# print(json_li)
# print(len(json_li))
# QLVESGGGLVQPGGSRKLSCSASGFAFSSFGMHWVRQAPEKGLEWVAYISSGSGTIYYADTVKGRFTISRDDPKNTLFLQMTSLRSEDTAMYYCVRSIYYYSGSPFDFWGQGTTLTVSSDIVMTQATSSVPVTPGESVSISCRSSKSLLHSNGNTYLYWFLQRPGQSPQLLIYRMSNLASGVPDRFSGSGSGTAFTLTISRLEAEDVGVYYCMQHLEYPLTFGAGTKL