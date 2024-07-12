import os
import sys
from Bio import SeqIO
from ab1_to_fasta import convert_ab1_to_fasta
from sequence_matcher import find_matching_sequences
from sequence_aligner import align_sequences

def process_files(input_files):
    # Find matching sequences
    matching_pairs = find_matching_sequences(input_files)
    
    results = {}
    for pair_number, pair in matching_pairs.items():
        ab1_file_f = pair['F']
        ab1_file_r = pair['R']
        
        # Convert AB1 to FASTA in memory
        seq1 = SeqIO.read(ab1_file_f, "abi")
        seq2 = SeqIO.read(ab1_file_r, "abi")
        
        # Align sequences and generate consensus
        temp_consensus = f"temp_consensus_{pair_number}.fasta"
        align_sequences(seq1, seq2, temp_consensus)
        
        # Read the consensus sequence
        with open(temp_consensus, "r") as infile:
            consensus = infile.read()
        
        # Store the result
        results[pair_number] = consensus
        
        # Clean up temporary file
        os.remove(temp_consensus)
    
    return results

def output_results(results):
    for pair_number, consensus in results.items():
        print(f"Pair {pair_number}:")
        print(consensus)
        print("---")  # Separator between pairs

if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_files = sys.argv[1:]
        results = process_files(input_files)
        output_results(results)
    else:
        print("Usage: python main.py input_file1.ab1 input_file2.ab1 ...")
