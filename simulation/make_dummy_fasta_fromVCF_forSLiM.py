#!/usr/bin/env python3

"""
Make a single concatenated dummy FASTA for SLiM, plus segments.tsv.

- Chromosome names and lengths are taken from a .fai index file.
- At positions present in a VCF, the base is the REF allele from the VCF (A or G, etc.).
- At all other positions, the base is 'C'.
- All chromosomes are concatenated into a single sequence called 'concat',
  in natural chromosome order (chr1, chr2, chr3, ..., chr10, ...).
- A segments.tsv file is also written, describing where each original chromosome
  maps into the concatenated coordinate system (0-based).

Usage:
    python make_concat_dummy_from_vcf.py ref.fa.fai dummy2.vcf out_prefix
"""

import sys
import re
from collections import OrderedDict, defaultdict

if len(sys.argv) != 4:
    sys.stderr.write(
        "Usage: python make_concat_dummy_from_vcf.py ref.fa.fai dummy2.vcf out_prefix\n"
    )
    sys.exit(1)

fai_path = sys.argv[1]
vcf_path = sys.argv[2]
out_prefix = sys.argv[3]

# ---------------------------------------------------------------------
# helper: natural sort key for chromosome names
# ---------------------------------------------------------------------
def natural_chr_key(name: str):
    """
    Produce a sort key so that chr1, chr2, ..., chr10, ...
    come in numeric order rather than chr1, chr10, chr11, chr2, ...
    Works for names like 'chr1', 'chr2', 'chrX', etc.
    """
    # strip common prefix 'chr'
    m = re.match(r'^(chr|Chr|CHR)?(.+)$', name)
    if m:
        core = m.group(2)
    else:
        core = name

    # try to interpret the core as an integer (e.g. '1', '2', '10')
    if core.isdigit():
        return (0, int(core))   # numeric chromosomes come first, sorted by number
    else:
        return (1, core)        # non-numeric come after, lexicographically


# ---------------------------------------------------------------------
# 1) Read chrom lengths from .fai
# ---------------------------------------------------------------------
# store as dict so we can re-order later
chrom_len_dict = {}  # chrom -> length
with open(fai_path) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        # .fai format: name length offset line_bases line_width
        fields = line.split("\t")
        chrom = fields[0]
        length = int(fields[1])
        chrom_len_dict[chrom] = length

# build a sorted list of (chrom, length) using natural chromosome order
chrom_lengths = [
    (chrom, chrom_len_dict[chrom])
    for chrom in sorted(chrom_len_dict.keys(), key=natural_chr_key)
]

# ---------------------------------------------------------------------
# 2) Read VCF and store REF alleles for each chrom:pos
# ---------------------------------------------------------------------
# positions[chrom][pos] = ref_base  (pos is 1-based)
positions = defaultdict(dict)

with open(vcf_path) as f:
    for line in f:
        if line.startswith("#"):
            continue
        line = line.rstrip("\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue

        chrom = fields[0]
        try:
            pos = int(fields[1])  # 1-based genomic coordinate
        except ValueError:
            continue

        ref = fields[3].upper()  # REF allele

        # Optional sanity check: warn on weird REF
        if ref not in ("A", "C", "G", "T", "N"):
            sys.stderr.write(
                f"[WARN] Unexpected REF {ref} at {chrom}:{pos}; using it anyway.\n"
            )

        positions[chrom][pos] = ref

# ---------------------------------------------------------------------
# 3) Build concatenated dummy sequence (single chromosome "concat")
# ---------------------------------------------------------------------

concat_seq = []              # list of bases, will be joined into one big string
segments = OrderedDict()     # chrom -> length_in_concat (original chrom length)

for chrom, length in chrom_lengths:
    pos2ref = positions.get(chrom, {})  # dict: pos -> REF for this chrom

    # Record segment length: we use full chromosome length in concat
    segments[chrom] = length

    # Walk through positions 1..length for this chromosome
    for pos in range(1, length + 1):
        if pos in pos2ref:
            base = pos2ref[pos]
        else:
            base = "C"
        concat_seq.append(base)

# Join into a single string
concat_str = "".join(concat_seq)

# ---------------------------------------------------------------------
# 4) Write concatenated FASTA (single chromosome)
# ---------------------------------------------------------------------
fasta_out = out_prefix + ".fa"
with open(fasta_out, "w") as out:
    out.write(">concat\n")
    # Wrap lines at 60 bp per line
    for i in range(0, len(concat_str), 60):
        out.write(concat_str[i:i+60] + "\n")

# ---------------------------------------------------------------------
# 5) Write segments.tsv: original chrom blocks in concatenated coordinates
# ---------------------------------------------------------------------
seg_out = out_prefix + "_segments.tsv"
with open(seg_out, "w") as out:
    out.write("#chrom\tlength_in_concat\tsim_start\tsim_end\n")
    offset = 0  # 0-based coordinate in the concatenated sequence
    for chrom, length in segments.items():
        start = offset
        end = offset + length - 1
        out.write(f"{chrom}\t{length}\t{start}\t{end}\n")
        offset += length

sys.stderr.write(f"Wrote {fasta_out} and {seg_out}\n")