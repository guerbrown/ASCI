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

# Usage
# align_sequences(seq1, seq2, "output_consensus.fasta")
