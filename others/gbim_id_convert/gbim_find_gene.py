import sys
import csv

id_in = sys.argv[1]
# gff_in = "./data/Gbi_Genes.gff"
gff_in = sys.argv[2]

id_in_list = list()

with open(id_in, "r") as f:
    for line in f:
        id_in_list.append(line.rstrip())

id_note_dict = dict()

with open(gff_in, "r") as f:
    for line in f:
        if line[0] == "#":
            continue
        else:
            line = line.split("\t")
            if line[2] == "gene":
                feature_list = line[8].split(";")
                feature_list_tmp = [ i for i in feature_list if i[0:2] == "ID" or i[0:4] == "Note" ]
                feature_list_new = [ i.split("=")[1] for i in feature_list_tmp]
                id_note_dict[feature_list_new[0]] = feature_list_new[1]

# print(id_note_dict)
# print(id_in_list)

output = ''

for id in id_in_list:
    output += id + "," + id_note_dict[id] + "\n"

# print(output.rstrip())

with open("./output.csv", 'w') as f:
    f.write(output.rstrip())
f.close()
