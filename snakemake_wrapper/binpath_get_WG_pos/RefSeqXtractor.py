import pysam
import numpy as np
import pandas as pd
import argparse # arguments management
import os # creating output directories

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-b","--bam",type=str,default=None,help="bam file input")
    parser.add_argument("-p","--pos_file",type=str,default=None,help="position bed file containing the position to extract the sequences")
    parser.add_argument("-o","--output",type=str,default=None,help="output name")
    parser.add_argument("-s","--sep",type=str,default='\t',help='columns separator in output file')
    args = parser.parse_args()

    if not args.bam:
        print("bam file is required!")
        exit(1)
    if not args.pos_file:
        print("pos file required!")
        exit(1)
    if not args.output:
        print("output name is required!")
        exit(1)


# LOADING DATA
## load pos file
pos_bed=pd.read_csv(args.pos_file,sep='\t',header=None)
pos_bed.columns=['chr','start','end','name_tag']
## load bam file
bamfile=pysam.AlignmentFile(args.bam,'rb')



#PROCESS
## OUTPUT MANAGEMENT
output=open(args.output,'w')
header_content=['refName','refStart','refEnd','2Xtract_seq_tag','2Xtract_ref_seq']

output.write(args.sep.join(header_content)+'\n')
output.close()

output=open(args.output,'a')

## define extraction function
def extract_seq(row_pos):
    fetch_seq=bamfile.fetch(row_pos['chr'],row_pos['start'],row_pos['end'])
    res=list()
    for segment in fetch_seq:
        
        output_line=[]

        name=segment.query_name
        refName=name.split(':')[0]
        refStart=int(name.split(':')[1].split('-')[0])

        output_line.append(refName)## write the name of the position given in the fasta (chr:start-end)

        pairs=np.array(segment.get_aligned_pairs(matches_only=True),dtype=object)

        refpos=pairs[:,1]
        mask=(refpos != None) & (np.array(refpos,dtype=float)>=row_pos['start']) & (np.array(refpos,dtype=float)<=(row_pos['end']-1))
        IOI=pairs[mask,0].tolist()
        output_line.append(str(refStart+IOI[0]))## write the start
        output_line.append(str(refStart+IOI[len(IOI)-1]))## write the end

        output_line.append('{raw_tag}:NM={NM}:AS={AS}:TE_POS={TE_POS}'.format(raw_tag=row_pos['name_tag'],NM=segment.get_tag('NM'),AS=segment.get_tag('AS'),TE_POS=name))

        #IOI=[pos for pos, refpos in pairs if refpos is not None and row_pos['start'] <= int(refpos) <= row_pos['end']]
        seq=segment.get_forward_sequence()
        Xtracted_seq=''.join(seq[index] for index in IOI)
        output_line.append(Xtracted_seq)

        output.write(args.sep.join(output_line)+'\n')
    print('{}:{}-{} done'.format(*row_pos))





## apply the function to 
pos_bed.apply(extract_seq,axis=1)

    
output.close()
