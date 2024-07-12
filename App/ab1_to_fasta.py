from Bio import SeqIO

def convert_ab1_to_fasta(ab1_file, fasta_file):
    with open(fasta_file, "w") as output_handle:
        SeqIO.convert(ab1_file, "abi", output_handle, "fasta")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2:
        convert_ab1_to_fasta(sys.argv[1], sys.argv[2])
    else:
        print("Usage: python ab1_to_fasta.py input.ab1 output.fasta")
