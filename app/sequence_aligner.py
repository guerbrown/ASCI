import os
import subprocess
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from io import StringIO
import sys
import shutil

def find_muscle():
    # First, check if MUSCLE is in the PATH
    muscle_path = shutil.which("muscle")
    if muscle_path:
        return muscle_path

    # If not in PATH, check common installation directories
    common_dirs = [
        "/usr/bin",
        "/usr/local/bin",
        "/opt/homebrew/bin",  # For macOS with Homebrew
        "C:\\Program Files\\MUSCLE",  # For Windows
        os.path.expanduser("~/bin"),  # User's home bin directory
    ]

    for directory in common_dirs:
        possible_path = os.path.join(directory, "muscle")
        if os.path.isfile(possible_path) and os.access(possible_path, os.X_OK):
            return possible_path

    # If MUSCLE is not found, raise an exception
    raise FileNotFoundError("MUSCLE executable not found. Please install MUSCLE or add it to your PATH.")

def calculate_consensus(alignment):
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus(threshold=0.5, ambiguous='N')
    return str(consensus)

def reverse_complement(seq):
    return SeqIO.SeqRecord(seq.seq.reverse_complement(), id=seq.id, description="")

def align_sequences(seq1, seq2, output_file):
    print(f"Aligning sequences: {seq1.id} and {seq2.id}")
    try:
        # Reverse complement the R sequence
        if seq2.id.endswith('R'):
            seq2 = reverse_complement(seq2)
            print(f"Reverse complemented: {seq2.id}")
        elif seq1.id.endswith('R'):
            seq1 = reverse_complement(seq1)
            print(f"Reverse complemented: {seq1.id}")

        muscle_path = find_muscle()
        muscle_cline = MuscleCommandline(
            muscle_path,
            input="-",
            gapopen=400,
            gapextend=0,
            cluster1="UPGMA"
        )
        stdin = f">{seq1.id}\n{seq1.seq}\n>{seq2.id}\n{seq2.seq}"
        
        # Use subprocess to run MUSCLE
        process = subprocess.Popen(str(muscle_cline), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
        stdout, stderr = process.communicate(stdin)
        
        if process.returncode != 0:
            print(f"MUSCLE error: {stderr}")
            raise Exception("MUSCLE alignment failed")

        align = AlignIO.read(StringIO(stdout), "fasta")
        consensus = calculate_consensus(align)
        
        with open(output_file, "w") as f:
            f.write(f">consensus\n{consensus}\n")
        
        print(f"Alignment saved to: {output_file}")
        return consensus
    except Exception as e:
        print(f"Error during alignment: {str(e)}")
        raise

if __name__ == "__main__":
    if len(sys.argv) > 3:
        seq1 = SeqIO.read(sys.argv[1], "fasta")
        seq2 = SeqIO.read(sys.argv[2], "fasta")
        output_file = sys.argv[3]
        align_sequences(seq1, seq2, output_file)
    else:
        print("Usage: python sequence_aligner.py seq1.fasta seq2.fasta output.fasta")
