import os
import sys
from Bio import SeqIO
from ab1_to_fasta import convert_ab1_to_fasta
from sequence_matcher import find_matching_sequences
from sequence_aligner import align_sequences

def process_files(input_files):
    print(f"Processing files: {input_files}")
    matching_pairs = find_matching_sequences(input_files)
    print(f"Matching pairs found: {matching_pairs}")
    
    results = {}
    for pair_number, pair in matching_pairs.items():
        print(f"Processing pair {pair_number}: {pair}")
        ab1_file_f = pair['F']
        ab1_file_r = pair['R']
        
        # Convert AB1 to FASTA
        fasta_file_f = ab1_file_f.replace('.ab1', '.fasta')
        fasta_file_r = ab1_file_r.replace('.ab1', '.fasta')
        convert_ab1_to_fasta(ab1_file_f, fasta_file_f)
        convert_ab1_to_fasta(ab1_file_r, fasta_file_r)
        
        # Read sequences
        seq1 = SeqIO.read(fasta_file_f, "fasta")
        seq2 = SeqIO.read(fasta_file_r, "fasta")
        
        # Align sequences and generate consensus
        temp_consensus = f"temp_consensus_{pair_number}.fasta"
        consensus = align_sequences(seq1, seq2, temp_consensus)
        
        # Store the result
        results[pair_number] = str(consensus)
        
        # Clean up temporary files
        os.remove(fasta_file_f)
        os.remove(fasta_file_r)
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
