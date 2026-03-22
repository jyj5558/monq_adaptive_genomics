#!/usr/bin/env python3
### usage ###
# python3 add_T_at_segment_ends.py concat.fa segments.tsv concat_with_T.fa
###

import sys

if len(sys.argv) != 4:
    sys.stderr.write("Usage: python3 add_T_at_segment_ends.py concat.fa segments.tsv out.fa\n")
    sys.exit(1)

fasta_in = sys.argv[1]
segments_in = sys.argv[2]
fasta_out = sys.argv[3]

def read_single_fasta(path):
    """
    Read a FASTA file that contains a single sequence.
    Returns (name, sequence_string).
    """
    name = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    # If multiple sequences exist, this will overwrite;
                    # here we assume only one sequence in the file.
                    pass
                name = line[1:].split()[0]
            else:
                chunks.append(line.strip())
    if name is None:
        raise ValueError("No FASTA header found in file: " + path)
    seq = "".join(chunks).upper()
    return name, seq

# Read the concat FASTA
seq_name, seq = read_single_fasta(fasta_in)
orig_len = len(seq)

# Read segment lengths from segments.tsv
segment_lengths = []
with open(segments_in) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            continue
        length = int(fields[1])  # length_in_concat column
        segment_lengths.append(length)

# Build new sequence by inserting 'T' after each segment
new_chunks = []
pos = 0  # 0-based index in the original concat sequence

for i, length in enumerate(segment_lengths):
    end = pos + length
    if end > orig_len:
        raise ValueError(
            f"Segment lengths exceed FASTA length: pos={pos}, length={length}, orig_len={orig_len}"
        )
    # Add the segment chunk
    new_chunks.append(seq[pos:end])
    # Add one 'T' as separator after the segment
    if i < len(segment_lengths) - 1:
        new_chunks.append("T")  # only between segments
    pos = end

# Optional: sanity check if we consumed the whole original sequence
if pos != orig_len:
    sys.stderr.write(
        f"[WARN] Sum of segment lengths ({pos}) != sequence length ({orig_len}); trailing bases will be ignored.\n"
    )

new_seq = "".join(new_chunks)

# Write new FASTA with 60 bp per line
with open(fasta_out, "w") as out:
    out.write(f">{seq_name}_with_T\n")
    for i in range(0, len(new_seq), 60):
        out.write(new_seq[i:i+60] + "\n")

sys.stderr.write(
    f"Original length: {orig_len}, new length (with Ts): {len(new_seq)}\n"
)
