import argparse

use = '\n'\
    '        python modify_cafe5_for_eggnog.py -t/--tsv [eggnog tsv file] -c/--cafe5 [cafe5-derived file] > output.tab'
des = ''

psr  = argparse.ArgumentParser(
    prog = 'modify_cafe5_for_eggnog',
    usage = use,
    description = des,
    formatter_class=argparse.RawTextHelpFormatter
    )

psr.add_argument('-t', '--tsv', required=True, help='Eggnog-mapper tsv file')
psr.add_argument('-c', '--cafe5', required=True, help='Tab-separated OG name, pval and increase/decrease number file')

args = psr.parse_args()

tsv_in = args.tsv
cafe5_in = args.cafe5

with open(tsv_in, 'r') as f:
    tsv_og_list = list()
    for line in f:
        if line[0] == "#":
            continue
        else:
            line = line.split("\t")
            tsv_og_list.append(line[0].split("_")[0])

# print(tsv_og_list)

with open(cafe5_in, 'r') as f:
    cafe5_og_and_pval = dict()
    cafe5_og_and_num = dict()
    cafe5_og = list()
    for line in f:
        line = line.split("\t")
        cafe5_og_and_pval[line[0]] = line[1]
        cafe5_og_and_num[line[0]] = line[2]
        cafe5_og.append(line[0])

# print(cafe5_og)

output = ''

for tsv_og in tsv_og_list:
    if tsv_og in cafe5_og:
        tmp = tsv_og + "\t" + cafe5_og_and_pval[tsv_og] + "\t" + cafe5_og_and_num[tsv_og].rstrip() + "\t" + "y"
        print(tmp)
    else:
        tmp = "-" + "\t" + "-" + "\t" + "-" + "\t" + "n"
        print(tmp)
