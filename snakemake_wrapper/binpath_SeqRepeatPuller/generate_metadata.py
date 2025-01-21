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

if 'aligner' in config.keys():
    aligner=config['aligner']
else :
    aligner='bowtie2'

if 'genome_folder' in config.keys():
    genome_folder=config['genome_folder']
else :
    genome_folder= None

metadata={"file_name" : args.input,
          "created_at": datetime.now().isoformat(),
          "version (git commit)of SeqRepeatPuller": config['git_commit'],
          "aligner" : aligner,
          "genome of reference" : genome_folder,
          "pos_file" : config['pos_file']}



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
          

