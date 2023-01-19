import pandas as pd
import csv

##########################

## BLASTP results ##

blastp_adp3 = "/work2/kataoka/amp/evigene230106/blastp/adp3.out.txt"
blastp_dbamp2 = "/work2/kataoka/amp/evigene230106/blastp/dbamp2.out.txt"
blastp_dramp2 = "/work2/kataoka/amp/evigene230106/blastp/dramp2.out.txt"
blastp_manual = "/work2/kataoka/amp/evigene230106/blastp/manual.out.txt"

## TBLASTN results ##

tblastn_adp3 = "/work2/kataoka/amp/evigene230106/tblastn/adp3.out.txt"
tblastn_dbamp2 = "/work2/kataoka/amp/evigene230106/tblastn/dbamp2.out.txt"
tblastn_dramp2 = "/work2/kataoka/amp/evigene230106/tblastn/dramp2.out.txt"
tblastn_manual = "/work2/kataoka/amp/evigene230106/tblastn/manual.out.txt"

## BLASTX results ##

blastx_adp3 = "/work2/kataoka/amp/evigene230106/blastx/adp3.out.txt"
blastx_dbamp2 = "/work2/kataoka/amp/evigene230106/blastx/dbamp2.out.txt"
blastx_dramp2 = "/work2/kataoka/amp/evigene230106/blastx/dramp2.out.txt"
blastx_manual = "/work2/kataoka/amp/evigene230106/blastx/manual.out.txt"

## ID - name reference ##

list_adp3 = "/work2/kataoka/amp/adp3_ID_name.txt"
list_dbamp2 = "/work2/kataoka/amp/dbAMP_pepinfo_20210922_dbAMP_Name.txt"
list_dramp2 = "/work2/kataoka/amp/general_amps_drampID_Name.txt"
list_manual = "/work2/kataoka/amp/target_id_name_final.txt"

##########################

outfmt_6 = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

def process_blastp_results(df, db_type):
    df = pd.read_table(df, names = [s + '_' + db_type for s in outfmt_6])
    df = df.set_index('qseqid' + '_' + db_type)
    df = df.iloc[:, [0, 9, 10]]
    return df

def process_tblastn_results(df, db_type):
    df = pd.read_table(df, names = [s + '_' + db_type for s in outfmt_6])
    df = df.set_index('sseqid' + '_' + db_type)
    df = df.iloc[:, [0, 9, 10]]
    return df

def process_blastx_results(df, db_type):
    df = pd.read_table(df, names = [s + '_' + db_type for s in outfmt_6])
    df = df.set_index('qseqid' + '_' + db_type)
    df = df.iloc[:, [0, 9, 10]]
    return df

def make_dict(list_file):
    out_dict = dict()
    with open(list_file, "r") as f:
        for line in f:
            line = line.split("\t")
            if len(line) == 2:
                out_dict[line[0]] = line[1].rstrip()
            else:
                out_dict[line[0]] = "unregistered"
    return out_dict

## Processing BLASTP ##

blastp_adp3 = process_blastp_results(blastp_adp3, "blastp_adp3")
blastp_dbamp2 = process_blastp_results(blastp_dbamp2, "blastp_dbamp2")
blastp_dramp2 = process_blastp_results(blastp_dramp2, "blastp_dramp2")
blastp_manual = process_blastp_results(blastp_manual, "blastp_manual")

merged_blastp = blastp_adp3.join([blastp_dbamp2, blastp_dramp2, blastp_manual], how = 'outer')

## Processing TBLASTN ##

tblastn_adp3 = process_tblastn_results(tblastn_adp3, "tblastn_adp3")
tblastn_dbamp2 = process_tblastn_results(tblastn_dbamp2, "tblastn_dbamp2")
tblastn_dramp2 = process_tblastn_results(tblastn_dramp2, "tblastn_dramp2")
tblastn_manual = process_tblastn_results(tblastn_manual, "tblastn_manual")

tblastn_adp3 = tblastn_adp3[tblastn_adp3.bitscore_tblastn_adp3 >= 90]
tblastn_dbamp2 = tblastn_dbamp2[tblastn_dbamp2.bitscore_tblastn_dbamp2 >= 90]
tblastn_dramp2 = tblastn_dramp2[tblastn_dramp2.bitscore_tblastn_dramp2 >= 90]
tblastn_manual = tblastn_manual[tblastn_manual.bitscore_tblastn_manual >= 90]

merged_tblastn = tblastn_adp3.join([tblastn_dbamp2, tblastn_dramp2, tblastn_manual], how = 'outer')

## Processing BLASTX ##

blastx_adp3 = process_blastx_results(blastx_adp3, "blastx_adp3")
blastx_dbamp2 = process_blastx_results(blastx_dbamp2, "blastx_dbamp2")
blastx_dramp2 = process_blastx_results(blastx_dramp2, "blastx_dramp2")
blastx_manual = process_blastx_results(blastx_manual, "blastx_manual")

merged_blastx = blastx_adp3.join([blastx_dbamp2, blastx_dramp2, blastx_manual], how = 'outer')

## Merge results ##

merged_blastp['blastp'] = "yes"
merged_tblastn['tblastn'] = "yes"
merged_blastx['blastx'] = "yes"

merged_df = merged_blastp.join([merged_tblastn, merged_blastx], how = 'outer')
merged_df = merged_df.sort_index(axis = 'index')
merged_df = merged_df.reindex(
    columns=[
        'sseqid_blastp_adp3',
        'sseqid_blastp_dbamp2', 
        'sseqid_blastp_dramp2',
        'sseqid_blastp_manual',
        'blastp',
        'qseqid_tblastn_adp3',
        'qseqid_tblastn_dbamp2',
        'qseqid_tblastn_dramp2', 
        'qseqid_tblastn_manual',
        'tblastn',
        'sseqid_blastx_adp3',
        'sseqid_blastx_dbamp2', 
        'sseqid_blastx_dramp2',
        'sseqid_blastx_manual',
        'blastx',
        'evalue_blastp_adp3',
        'evalue_blastp_dbamp2', 
        'evalue_blastp_dramp2',
        'evalue_blastp_manual',
        'bitscore_blastp_adp3',
        'bitscore_blastp_dbamp2',
        'bitscore_blastp_dramp2',
        'bitscore_blastp_manual', 
        'evalue_tblastn_adp3',
        'evalue_tblastn_dbamp2',
        'evalue_tblastn_dramp2',
        'evalue_tblastn_manual',
        'bitscore_tblastn_adp3',
        'bitscore_tblastn_dbamp2',
        'bitscore_tblastn_dramp2',
        'bitscore_tblastn_manual',
        'evalue_blastx_adp3',
        'evalue_blastx_dbamp2',
        'evalue_blastx_dramp2',
        'evalue_blastx_manual',
        'bitscore_blastx_adp3',
        'bitscore_blastx_dbamp2',
        'bitscore_blastx_dramp2',
        'bitscore_blastx_manual'
        ])

# print(merged_df)
# merged_df.to_csv("./merged_results.csv")

## Name editting ##

# merged_df.replace(
#     {
#         'sseqid_blastp_adp3': make_dict(list_adp3),
#         'sseqid_blastp_dbamp2': make_dict(list_dbamp2),
#         'sseqid_blastp_dramp2': make_dict(list_dramp2),
#         'sseqid_blastp_manual': make_dict(list_manual),
#         'qseqid_tblastn_adp3': make_dict(list_adp3),
#         'qseqid_tblastn_dbamp2': make_dict(list_dbamp2),
#         'qseqid_tblastn_dramp2': make_dict(list_dramp2),
#         'qseqid_tblastn_manual': make_dict(list_manual),
#         'sseqid_blastx_adp3': make_dict(list_adp3),
#         'sseqid_blastx_dbamp2': make_dict(list_dbamp2),
#         'sseqid_blastx_dramp2': make_dict(list_dramp2),
#         'sseqid_blastx_manual': make_dict(list_manual),
#         }
#         )

# print(make_dict(list_adp3)["02797|cOT2"])

merged_df = merged_df.replace({
    'sseqid_blastp_adp3': make_dict(list_adp3),
    'qseqid_tblastn_adp3': make_dict(list_adp3),
    'sseqid_blastx_adp3': make_dict(list_adp3),
    'sseqid_blastp_dbamp2': make_dict(list_dbamp2),
    'qseqid_tblastn_dbamp2': make_dict(list_dbamp2),
    'sseqid_blastx_dbamp2': make_dict(list_dbamp2),
    'sseqid_blastp_dramp2': make_dict(list_dramp2),
    'qseqid_tblastn_dramp2': make_dict(list_dramp2),
    'sseqid_blastx_dramp2': make_dict(list_dramp2),
    'sseqid_blastp_manual': make_dict(list_manual),
    'qseqid_tblastn_manual': make_dict(list_manual),
    'sseqid_blastx_manual': make_dict(list_manual)
    })

# dict_test = {'dbAMP_10728': 'test', 'dbAMP_08936' : 'test2'}

# merged_df = merged_df.replace({'sseqid_blastp_dbamp2': dict_test})

print(merged_df)

# print(make_dict(list_manual))

merged_df.to_csv("./merged_results.csv")
