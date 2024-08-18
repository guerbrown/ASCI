import os
import subprocess
from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from io import StringIO
import sys
import shutil
import json

def ab1_to_fasta(ab1_file):
    try:
        record = SeqIO.read(ab1_file, "abi")
        return record
    except FileNotFoundError:
        print(f"Error: File not found: {ab1_file}")
        return None

def find_muscle():
    muscle_path = shutil.which("muscle")
    if muscle_path:
        print(f"MUSCLE found in PATH: {muscle_path}")
        return muscle_path
    raise FileNotFoundError("MUSCLE executable not found. Please install MUSCLE or add it to your PATH.")

def get_muscle_version(muscle_path):
    try:
        result = subprocess.run([muscle_path, "-version"], capture_output=True, text=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        return "Unknown"

def align_sequences(seq_f, seq_r, output_file):
    print(f"Aligning sequences: {seq_f.id} and {seq_r.id}")
    try:
        # Always reverse complement the R sequence
        seq_r_rc = reverse_complement(seq_r)
        print(f"Reverse complemented: {seq_r.id}")

        muscle_path = find_muscle()
        muscle_version = get_muscle_version(muscle_path)
        print(f"MUSCLE version: {muscle_version}")

        with open("temp_input.fasta", "w") as temp_file:
            SeqIO.write([seq_f, seq_r_rc], temp_file, "fasta")
        
        if muscle_version.startswith("3"):
            command = [muscle_path, "-in", "temp_input.fasta", "-out", "temp_output.fasta"]
        else:  # MUSCLE v5 syntax
            command = [muscle_path, "-align", "temp_input.fasta", "-output", "temp_output.fasta"]

        print(f"Running command: {' '.join(command)}")
        process = subprocess.run(command, capture_output=True, text=True)
        
        print(f"MUSCLE stdout: {process.stdout}")
        print(f"MUSCLE stderr: {process.stderr}")
        
        if process.returncode != 0:
            raise Exception(f"MUSCLE alignment failed: {process.stderr}")

        align = AlignIO.read("temp_output.fasta", "fasta")
        consensus = calculate_consensus(align)
        
        os.remove("temp_input.fasta")
        os.remove("temp_output.fasta")
        
        return consensus
    except Exception as e:
        print(f"Error during alignment: {str(e)}")
        return None

def calculate_consensus(alignment):
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus(threshold=0.5, ambiguous='N')
    return str(consensus)

def reverse_complement(seq):
    return SeqIO.SeqRecord(seq.seq.reverse_complement(), id=seq.id, description="")

def process_pairs(pairs, output_file):
    all_consensuses = []
    for pair_id, pair in pairs.items():
        print(f"Processing pair {pair_id}: {pair}")
        seq_f = ab1_to_fasta(pair['F'])
        seq_r = ab1_to_fasta(pair['R'])
        if seq_f and seq_r:
            consensus = align_sequences(seq_f, seq_r, f"temp_consensus_{pair_id}.fasta")
            if consensus:
                # Use os.path.basename to get just the filename without the path
                header = os.path.basename(pair_id)
                all_consensuses.append((header, consensus))
        else:
            print(f"Error: Unable to read input files for pair {pair_id}")

    # Write all consensuses to a single file
    with open(output_file, "w") as f:
        for pair_id, consensus in all_consensuses:
            f.write(f">{pair_id}\n{consensus}\n")
    print(f"All consensus sequences saved to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) > 2:
        pairs = json.loads(sys.argv[1])
        output_file = sys.argv[2]
        process_pairs(pairs, output_file)
    else:
        print("Usage: python sequence_aligner.py '{\"CL023-1\": {\"F\": \"path/to/F.ab1\", \"R\": \"path/to/R.ab1\"}, ...}' output.fasta")
