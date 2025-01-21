import argparse
import json
from datetime import datetime
import ast
import numpy as np
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input",type=str,required=True,default=None)
    parser.add_argument("-c","--config",type=str,required=True,default=None)
    parser.add_argument("-o","--output",type=str,required=True,default=None) 
    args = parser.parse_args()

#config=json.loads(args.config)
config=ast.literal_eval(args.config)
print(config)
metadata={"file_name" : args.input,
          "created_at": datetime.now().isoformat(),
          "version (git commit)of get_wg_pos": config['git_commit'],
          "aligner" : "bowtie2",
          "genome of reference" : config['genome_ref'],
          "transposable element position bed file" : config['TE_bed'],
          "transposable element reference sequenece" : config['TE_refseq'],
          "position of interest acrosse the te reference sequence" : config['seq_bed']}

if 'mode' in config.keys():
    metadata.update({"mode of alignement " : config['mode']})
    if config['mode'] == 'fragmented':
        supp_metadata={'fragment_size': config['fragment_size'], 'fragment_step': config['fragment_step']}
        metadata.update(supp_metadata)
else :
    metadata.update({"mode of alignement " : 'end-to-end'})

if 'clean_xtracted_results' in config.keys():
    metadata.update({"post process cleaning done" : config['clean_xtracted_results']})
    if config['clean_xtracted_results'] == 'True':

        import pickle # work
        with open("wg_xtracted_position_cleaned.xbed.pkl","rb") as file: #
            supp_metadata=pickle.load(file)
            metadata.update(supp_metadata)
else :

    metadata.update({"post process cleaning done" : 'False'})

# Ensure JSON compatibility
def convert_to_json_compatible(obj):
    if isinstance(obj, np.integer):  # Convert numpy int to Python int
        return int(obj)
    elif isinstance(obj, np.floating):  # Convert numpy float to Python float
        return float(obj)
    elif isinstance(obj, np.ndarray):  # Convert numpy arrays to lists
        return obj.tolist()
    else:
        return obj

# Convert the dictionary
metadata_cleaned = {k: convert_to_json_compatible(v) for k, v in metadata.items()}

          
with open(args.output,"w") as out_file:
    json.dump(metadata_cleaned,out_file,indent=4)
          

