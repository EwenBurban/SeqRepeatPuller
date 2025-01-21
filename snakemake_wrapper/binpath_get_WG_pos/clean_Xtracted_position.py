import numpy as np
import pandas as pd
import argparse # arguments management
import os # creating output directories
import pickle
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-b","--xbed",type=str,default=None,help="bed file containing TE of interest position")
    parser.add_argument("-o","--output",type=str,default=None,help="output name")
    parser.add_argument("-m","--metadata",type=str,default=None,help="metadata pickle file name")
    args = parser.parse_args()



    if not args.xbed:
        print("bed file is required!")
        exit(1)
    if not args.output:
        print("output name is required!")
        exit(1)


### load bed file
bed_file=pd.read_csv(args.xbed,sep='\t')


### remove duplicates
bed_file['unique_pos_name']=['{refName}-{refStart}-{refEnd}'.format(**bed_file.loc[x,:]) for x in range(bed_file.shape[0])]
duplicates=bed_file['unique_pos_name'].duplicated(keep=False)
if any(duplicates) :
    duplicated_rows=bed_file[duplicates]
    bed_file=bed_file[duplicates==False]
    

    duplicated_rows["AlignementScore"]= duplicated_rows["2Xtract_seq_tag"].str.extract(r"AS=(-?\d+)").astype(int)
    count=duplicated_rows['unique_pos_name'].value_counts()
    mean_duplication_rate=np.mean(count)

    idxMaxAS=duplicated_rows.groupby('unique_pos_name')["AlignementScore"].idxmax()
    rows_keept=duplicated_rows.loc[idxMaxAS,:].reset_index(drop=True)
    N_duplicated_position=rows_keept.shape[0]
    bed_file=pd.concat([bed_file,rows_keept],axis=0)

### clean gappy position
### gappy position is when the sequence of reference is shorter than the sequence i got (so there is gap in the alignement) 
# compare 2Xtract_ref_seq to diff between refStart and refEnd (+1) 
bed_file['true_seq_len']=bed_file['refEnd'] - bed_file['refStart'] + 1
bed_file['theo_seq_len']=bed_file['2Xtract_ref_seq'].str.len()
N_gappy=(bed_file['true_seq_len']!=bed_file['theo_seq_len']).sum()
bed_file=bed_file[bed_file['true_seq_len']==bed_file['theo_seq_len']]

bed_file= bed_file.drop(columns=['AlignementScore','unique_pos_name','true_seq_len','theo_seq_len'])

bed_file.to_csv(args.output,index=None,sep='\t')


### write in append on metadata pickle
with open(args.metadata,'wb') as f:
    pickle.dump({"Nb of gappy position removed":N_gappy,"Nb of unique duplicated position":N_duplicated_position,
                 "mean duplication rate":mean_duplication_rate},f)

