import os
import yaml
from Bio import SeqIO
from ab1_to_fasta import convert_ab1_to_fasta
from sequence_matcher import find_matching_sequences
from sequence_aligner import align_sequences

def load_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def process_files(config_path):
    config = load_config(config_path)
    input_directory = config['input_directory']
    output_directory = config['output_directory']
    file_pattern = config['file_pattern']

    # Find matching sequences
    matching_pairs = find_matching_sequences(input_directory, file_pattern)
    
    # Convert, align sequences and generate consensus
    consensus_file = os.path.join(output_directory, "consensus.fasta")
    with open(consensus_file, "w") as outfile:
        for pair in matching_pairs:
            ab1_file_f = os.path.join(input_directory, pair['F'])
            ab1_file_r = os.path.join(input_directory, pair['R'])
            
            # Convert AB1 to FASTA in memory
            seq1 = SeqIO.read(ab1_file_f, "abi")
            seq2 = SeqIO.read(ab1_file_r, "abi")
            
            # Align sequences and generate consensus
            align_sequences(seq1, seq2, "temp_consensus.fasta")
            
            with open("temp_consensus.fasta", "r") as infile:
                outfile.write(infile.read())
    
    os.remove("temp_consensus.fasta")
    return consensus_file

if __name__ == "__main__":
    config_path = "config.yaml"
    final_consensus = process_files(config_path)
    print(f"Consensus file created: {final_consensus}")