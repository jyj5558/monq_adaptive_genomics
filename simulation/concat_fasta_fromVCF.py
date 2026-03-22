#!/usr/bin/env python3

### usage ###
# python3 concat_fasta_fromVCF.py all.vcf outlier.vcf concat
### usage ###

import sys
from collections import OrderedDict

if len(sys.argv) != 4:
    sys.stderr.write("Usage: python concat_from_two_vcfs.py vcf1.vcf vcf2.vcf out_prefix\n")
    sys.exit(1)

vcf1_in = sys.argv[1]   # first VCF: defines all sites and their order (all start as 'A')
vcf2_in = sys.argv[2]   # second VCF: sites that should be changed from 'A' to 'G'
out_prefix = sys.argv[3]

# --- Read first VCF and build the concat sequence as all 'A' ---
concat_seq = []                  # list of characters, will be joined into final FASTA sequence
segments = OrderedDict()         # chr -> length in concat (number of sites from VCF1)
coord_to_index = {}              # "chrom:pos" -> index in concat_seq

index = 0  # 0-based index within concat_seq, increases for each site in VCF1

with open(vcf1_in) as f:
    for line in f:
        # Skip header lines
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2:
            continue

        chrom = fields[0]
        pos = fields[1]  # keep as string to avoid conversions back and forth

        key = f"{chrom}:{pos}"

        # For each site in the first VCF, we assign 'A' as the base
        concat_seq.append("A")

        # Record mapping from genomic coordinate to index in concat sequence
        coord_to_index[key] = index
        index += 1

        # Update segment length per chromosome
        segments.setdefault(chrom, 0)
        segments[chrom] += 1

# --- Read second VCF and update matching sites from 'A' to 'G' ---
updated_count = 0
skipped_count = 0

with open(vcf2_in) as f:
    for line in f:
        # Skip header lines
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2:
            continue

        chrom = fields[0]
        pos = fields[1]

        key = f"{chrom}:{pos}"

        # Only update if this site also exists in the first VCF (i.e. in concat_seq)
        if key in coord_to_index:
            idx = coord_to_index[key]
            concat_seq[idx] = "G"
            updated_count += 1
        else:
            # Site exists in second VCF but not in the first: skip and warn
            skipped_count += 1
            # Comment out the next line if you don't want warnings
            sys.stderr.write(f"[WARN] Site {key} in VCF2 not found in VCF1; skipped\n")

# Convert list of bases to a single string
concat_str = "".join(concat_seq)

# --- Write concat FASTA ---
fasta_out = out_prefix + ".fa"
with open(fasta_out, "w") as out:
    out.write(">concat\n")
    # Wrap lines at 60 bp per line (standard FASTA style)
    for i in range(0, len(concat_str), 60):
        out.write(concat_str[i:i+60] + "\n")

# --- Write segment information: original chromosome blocks in concat coordinates ---
seg_out = out_prefix + "_segments.tsv"
with open(seg_out, "w") as out:
    out.write("#chrom\tlength_in_concat\tsim_start\tsim_end\n")
    offset = 0
    for chrom, length in segments.items():
        # sim_start and sim_end are 0-based coordinates in the concatenated sequence
        start = offset
        end = offset + length - 1
        out.write(f"{chrom}\t{length}\t{start}\t{end}\n")
        offset += length

sys.stderr.write(f"Wrote {fasta_out} and {seg_out}\n")
sys.stderr.write(f"Updated {updated_count} sites to G; skipped {skipped_count} sites not in VCF1\n")
