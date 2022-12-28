
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import math
from scipy import stats

def log10_tpm(df):
    return math.log10(df.iloc[1] + 1)

use = '\n'\
    '        python nanocount_tpm_scatter_plot.py -1 input1 -2 input2 -1n input1_label -2n input2_label'
des = ''

psr  = argparse.ArgumentParser(
    prog = 'nanocount_tpm_scatter_plot',
    usage = use,
    description = des,
    formatter_class=argparse.RawTextHelpFormatter
    )

psr.add_argument('-1', '--input_1', required=True, help='NanoCount input 1')
psr.add_argument('-2', '--input_2', required=True, help='NanoCount input 2')
psr.add_argument('-1n', '--input_1_label', required=True, help='NanoCount input 1 label')
psr.add_argument('-2n', '--input_2_label', required=True, help='NanoCount input 2 label')

args = psr.parse_args()

input_1 = args.input_1
input_2 = args.input_2
input_1_label = args.input_1_label
input_2_label = args.input_2_label

df1 = pd.read_table(input_1, usecols=[0, 3])
df2 = pd.read_table(input_2, usecols=[0, 3])

df1 = df1.rename(columns={'tpm': input_1_label})
df2 = df2.rename(columns={'tpm': input_2_label})

input_1_label_alt = input_1_label + "_alt"
input_2_label_alt = input_2_label + "_alt"

df1[input_1_label_alt] = df1.apply(log10_tpm, axis=1)
df2[input_2_label_alt] = df2.apply(log10_tpm, axis=1)

# print(df1)
# print(df2)

df_merge = pd.merge(df1, df2, on='transcript_name', how='outer')
df_merge = df_merge.fillna(0)

# print(df_merge)

correlation, pvalue = stats.pearsonr(df_merge[input_1_label_alt], df_merge[input_2_label_alt])
print(correlation, pvalue)

plt.figure()
df_merge.plot()
df_merge.plot.scatter(x = input_1_label_alt, y = input_2_label_alt, alpha=0.5)
plt.savefig('./scatter.png')
plt.close('all')