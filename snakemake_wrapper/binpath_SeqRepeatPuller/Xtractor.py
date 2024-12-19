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
    parser.add_argument("-q","--mapping_quality_tag",type=str,default='NM',help='tag name used in bam input file to quantify the quality of the alignment. By default, it’ NM tag (Number of Mismatch)')
    parser.add_argument("-R","--right_ovelap_thr",type=int,default=10,help='minimun amount of bases that the read have after the end of seq of interest')
    parser.add_argument("-N","--name_tag",type=str,default=None,help='sample name tag')
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
sep=args.sep

# LOADING DATA
## load pos file
# CHANGE LOADING FOR PANDAS EFFICIENT LOADING AND SPECIFY THE NAME OF THE COLUNMS
pos_bed=pd.read_csv(args.pos_file,sep='\t',
                    usecols=['refName','refStart','refEnd','2Xtract_seq_tag','2Xtract_ref_seq'],
                    dtype={'refName': 'str','refStart':'int64','refEnd':'int64','2Xtract_seq_tag':'str','2Xtract_ref_seq':'str'})
## load bam file
bamfile=pysam.AlignmentFile(args.bam,'rb')



#PROCESS
## OUTPUT MANAGEMENT
output=open(args.output,'w')
header_content=['Read_ID','Read_Xtracted_seq','Xtracted_seq_tag','Xtracted_ref_seq','alignmentQual','Xtracted_seq_pos','Read_insert_start','Read_insert_end','Read_length','Name_tag']
output.write(sep.join(header_content)+'\n')

output.close()

output=open(args.output,'a')

## define extraction function
def extract_seq(row_pos):
    fetch_seq=bamfile.fetch(row_pos['refName'],row_pos['refStart'],row_pos['refEnd'])
    res=list()
    for segment in fetch_seq:
        # Segment right overlap checkpoint (purpose is that the segment must map args.right_overlap_thr nucleotide right away from the  
        if segment.reference_end-row_pos['refEnd'] < args.right_ovelap_thr:
            continue
        
        output_line=[]
        output_line.append(segment.query_name)

        pairs=np.array(segment.get_aligned_pairs(matches_only=True),dtype=object)

        refpos=pairs[:,1]
        mask=(refpos != None) & (np.array(refpos,dtype=float)>=row_pos['refStart']) & (np.array(refpos,dtype=float)<=(row_pos['refEnd']))
        IOI=pairs[mask,0].tolist()
        seq=segment.get_forward_sequence()
        Xtracted_seq=''.join(seq[index] for index in IOI)
        #Xtracted_seq=seq[(row_pos['refStart']-segment.reference_start):(row_pos['refEnd']-segment.reference_start+1)]
        output_line.append(Xtracted_seq)

        output_line.append(row_pos['2Xtract_seq_tag'])

        output_line.append(row_pos['2Xtract_ref_seq'])


        mapping_quality_score=segment.get_tag(args.mapping_quality_tag)
        output_line.append(str(mapping_quality_score))


        Xtracted_seq_pos='{chr}:{start}-{end}'.format(chr=row_pos['refName'],start=row_pos['refStart'],end=row_pos['refEnd'])
        output_line.append(Xtracted_seq_pos)

        output_line.append(str(segment.reference_start))
        output_line.append(str(segment.reference_end))
        output_line.append(str(segment.query_length))
        output_line.append(args.name_tag)

        output.write(sep.join(output_line)+'\n')
    print('{}:{}-{} done'.format(*row_pos))


## apply the function to 
pos_bed.apply(extract_seq,axis=1)

    
output.close()
