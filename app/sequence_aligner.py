from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from io import StringIO

def align_sequences(seq1, seq2, output_file):
    muscle_cline = MuscleCommandline(input="-")
    stdin = f">{seq1.id}\n{seq1.seq}\n>{seq2.id}\n{seq2.seq}"
    stdout, stderr = muscle_cline(stdin=stdin)
    align = AlignIO.read(StringIO(stdout), "fasta")
    consensus = align.consensus()
    with open(output_file, "w") as f:
        f.write(f">consensus\n{consensus}\n")

if __name__ == "__main__":
    from Bio import SeqIO
    import sys
    if len(sys.argv) > 3:
        seq1 = SeqIO.read(sys.argv[1], "fasta")
        seq2 = SeqIO.read(sys.argv[2], "fasta")
        align_sequences(seq1, seq2, sys.argv[3])
    else:
        print("Usage: python sequence_aligner.py seq1.fasta seq2.fasta output.fasta")
