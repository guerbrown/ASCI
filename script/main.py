import os
from ab1_to_fasta import convert_ab1_to_fasta
from sequence_matcher import find_matching_sequences
from sequence_aligner import align_sequences
from Bio import SeqIO

def process_files(input_directory, output_directory):
    # Convert AB1 to FASTA
    for file in os.listdir(input_directory):
        if file.endswith(".ab1"):
            ab1_file = os.path.join(input_directory, file)
            fasta_file = os.path.join(output_directory, file.replace(".ab1", ".fasta"))
            convert_ab1_to_fasta(ab1_file, fasta_file)
    
    # Find matching sequences
    matching_pairs = find_matching_sequences(output_directory)
    
    # Align sequences and generate consensus
    consensus_file = os.path.join(output_directory, "consensus.fasta")
    with open(consensus_file, "w") as outfile:
        for pair in matching_pairs:
            seq1 = SeqIO.read(os.path.join(output_directory, pair['F']), "fasta")
            seq2 = SeqIO.read(os.path.join(output_directory, pair['R']), "fasta")
            align_sequences(seq1, seq2, "temp_consensus.fasta")
            
            with open("temp_consensus.fasta", "r") as infile:
                outfile.write(infile.read())
    
    os.remove("temp_consensus.fasta")
    return consensus_file

# Usage
# final_consensus = process_files("input_directory", "output_directory")
