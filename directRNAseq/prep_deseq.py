import pandas as pd
import math, sys, csv

# Integrate Nanocount outputs into a matrix for DESeq2

# python prep_deseq.py tsv_list.txt

# tsv_list.txt
# control_1	transcript_counts.control_1.tsv
# control_2	transcript_counts.control_2.tsv
# thcv_1	transcript_counts.thcv_1.tsv
# thcv_2	transcript_counts.thcv_2.tsv

def split_transcript_name(transcript_name):
    transcript_name = transcript_name.split("|")
    return transcript_name[0]

def extract_est_count(tsv_in):
    output = pd.DataFrame()
    with open(tsv_in, "r") as f:
        reader = csv.reader(f, delimiter = "\t")
        header = next(reader)
        for cols in reader:
            col_tmp = pd.Series([split_transcript_name(cols[0]), cols[2]])
            output = output.append(col_tmp, ignore_index=True)
    output = output.rename(columns={0: 'transcript_name', 1: 'est_count'})
    output = output.set_index('transcript_name')
    output["est_count"] = output["est_count"].apply(float)
    output["est_count"] = output["est_count"].apply(math.floor)
    return output

def extract_tpm(tsv_in):
    output = pd.DataFrame()
    with open(tsv_in, "r") as f:
        reader = csv.reader(f, delimiter = "\t")
        header = next(reader)
        for cols in reader:
            col_tmp = pd.Series([split_transcript_name(cols[0]), cols[3]])
            output = output.append(col_tmp, ignore_index=True)
    output = output.rename(columns={0: 'transcript_name', 1: 'tpm'})
    output = output.set_index('transcript_name')
    return output

file_in = sys.argv[1]

df_est_count = pd.DataFrame()
df_tpm = pd.DataFrame()

with open(file_in, "r") as f:
    for line in f:
        line = line.split("\t")
        
        df_est_count_tmp = extract_est_count(line[1].rstrip())
        df_tpm_tmp = extract_tpm(line[1].rstrip())

        df_est_count_tmp = df_est_count_tmp.rename(columns={'est_count': line[0]})
        df_tpm_tmp = df_tpm_tmp.rename(columns={'tpm': line[0]})

        df_est_count = df_est_count.join([df_est_count_tmp], how = 'outer')
        df_tpm = df_tpm.join([df_tpm_tmp], how = 'outer')

df_est_count = df_est_count.fillna(0)
df_tpm = df_tpm.fillna(0)

df_est_count.to_csv("./transcript_count.csv")
df_tpm.to_csv("./transcript_tpm.csv")

