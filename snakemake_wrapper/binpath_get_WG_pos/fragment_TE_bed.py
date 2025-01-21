import numpy as np
import pandas as pd
import argparse # arguments management
import os # creating output directories

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-b","--bed",type=str,default=None,help="bed file containing TE of interest position")
    parser.add_argument("-o","--output",type=str,default=None,help="output name")
    parser.add_argument("--fragment_size",type=int,default=167,help="size  of the fragment produced")
    parser.add_argument("--fragment_step",type=int,default=20,help="size of the step between each fragment start. If fragment_size=fragment_step,then you produce non overlapping fragment")
    args = parser.parse_args()

    if not args.bed:
        print("bed file is required!")
        exit(1)
    if not args.output:
        print("output name is required!")
        exit(1)

bed_file=pd.read_csv(args.bed,sep='\t',header=None)




def generate_frag(TE_OI_data,fragment_size,fragment_step):
    start=TE_OI_data[1]
    end=TE_OI_data[2]

    windows = [(int(i), int(i + fragment_size - 1)) for i in range(start, end - fragment_size + 1, fragment_step)]
    df=pd.DataFrame(windows)
    df.insert(0,'chr',TE_OI_data[0])
    if len(TE_OI_data) > 3:

        new_colums=pd.DataFrame([TE_OI_data[3:] for i in range(df.shape[0])])
        new_colums.reset_index(inplace=True,drop=True)
        df=pd.concat([df,new_colums],axis=1)

    return df

df_list=[generate_frag(bed_file.loc[x,:],args.fragment_size,args.fragment_step) for x in range(bed_file.shape[0])]
fragmented_bed=pd.concat([df for df in df_list if not df.empty],axis=0)
# Replace NaN and convert to int64
fragmented_bed.iloc[:, 1] = fragmented_bed.iloc[:, 1].fillna(0).astype(int)
fragmented_bed.iloc[:, 2] = fragmented_bed.iloc[:, 2].fillna(0).astype(int)

# Check results
fragmented_bed.to_csv(args.output,header=None,index=None,sep='\t')


