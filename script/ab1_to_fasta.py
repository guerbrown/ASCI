from Bio import SeqIO

def convert_ab1_to_fasta(ab1_file, fasta_file):
    with open(fasta_file, "w") as output_handle:
        SeqIO.convert(ab1_file, "abi", output_handle, "fasta")

# Usage
# convert_ab1_to_fasta("input.ab1", "output.fasta")
