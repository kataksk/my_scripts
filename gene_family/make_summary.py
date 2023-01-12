import sys
import pandas as pd
import argparse

use = '\n''\n'\
    'python make_summary.py -p p-value -r Base_clade_results.txt -b Base_branch_probabilities.tab -c Base_change.tab'
des = ''

psr  = argparse.ArgumentParser(
    prog = 'nanocount_tpm_scatter_plot',
    usage = use,
    description = des,
    formatter_class=argparse.RawTextHelpFormatter
    )

psr.add_argument('-p', '--pval', required=True, help='p-value e.g. 0.05')
psr.add_argument('-r', '--clade_results', required=True, help='Base_clade_results.txt')
psr.add_argument('-b', '--branch_probabilities', required=True, help='Base_branch_probabilities.tab')
psr.add_argument('-c', '--change', required=True, help='Base_change.tab')

args = psr.parse_args()

pval = args.pval
clade_results = args.clade_results
branch_probabilities = args.branch_probabilities
change = args.change

def count_significant_og(taxon_name, df, pval):
    list_pval = df[taxon_name].to_list()
    cnt = sum(x < float(pval) for x in list_pval)
    return cnt

def output_significant_og_list(taxon_name, df, pval):
    df_sig = df.loc[df[taxon_name] < float(pval)]
    index_list = df_sig.index.to_list()
    return index_list

def count_changes(list_data):
    expansion = 0
    contraction = 0
    for num in list_data:
        if num > 0:
            expansion += 1
        elif num < 0:
            contraction += 1
    return expansion, contraction

clade_results = pd.read_table(
    clade_results,
    names=["Taxon", "Families Expansions", "Families Contractions"],
    skiprows=[0]
    )

clade_results = clade_results.set_index("Taxon")

# print(clade_results)

branch_probabilities = pd.read_table(
    branch_probabilities,
    index_col = 0
    )

# print(branch_probabilities)

clade_results["Total Significant Changes"] = ''

for taxon_name in clade_results.index:
    # print(taxon_name)
    clade_results["Total Significant Changes"][taxon_name] = count_significant_og(taxon_name, branch_probabilities, pval)

# print(clade_results)

change = pd.read_table(
    change,
    index_col = 0
    )

# print(change)

# print(len(clade_results.index))
# 22

clade_results["Signicant Expansions"] = ''
clade_results["Signicant Contractions"] = ''

for taxon_name in clade_results.index:
    sig_og_list = output_significant_og_list(taxon_name, branch_probabilities, pval)
    expansion_num = count_changes(change.loc[sig_og_list, : ][taxon_name].to_list())[0]
    contraction_num = count_changes(change.loc[sig_og_list, : ][taxon_name].to_list())[1]
    clade_results["Siginicant Expansions"][taxon_name] = expansion_num
    clade_results["Siginicant Contractions"][taxon_name] = contraction_num

# print(clade_results)

# output_significant_og_list(taxon_name, branch_probabilities, pval)

clade_results.to_csv('./out.csv')
