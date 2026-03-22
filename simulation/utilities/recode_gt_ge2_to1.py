#!/usr/bin/env python3
### usage ###
# python3 recode_gt_ge2_to1.py input.vcf > output.vcf
###

import sys

if len(sys.argv) != 2:
    sys.stderr.write("Usage: python recode_gt_ge2_to1.py input.vcf > output.vcf\n")
    sys.exit(1)

vcf_in = sys.argv[1]

def recode_gt(gt_str: str) -> str:
    """
    Recode a GT string so that:
      - 0 stays 0
      - 1 stays 1
      - any allele index >= 2 becomes 1
      - '.' (missing) stays '.'
    Phasing (/ vs |) is preserved.
    """
    gt_str = gt_str.strip()
    if gt_str == "" or gt_str == ".":
        return gt_str

    # Detect separator: phased '|' or unphased '/'
    sep = "/"
    if "|" in gt_str:
        sep = "|"

    # Only split GT part (no subfields here)
    alleles = gt_str.replace("|", "/").split("/")

    new_alleles = []
    for a in alleles:
        a = a.strip()
        if a == ".":
            new_alleles.append(".")
            continue
        try:
            idx = int(a)
        except ValueError:
            # Unexpected token, keep as is
            new_alleles.append(a)
            continue

        if idx <= 1:
            new_alleles.append(str(idx))  # 0 or 1
        else:
            new_alleles.append("1")       # collapse >=2 to 1

    return sep.join(new_alleles)


with open(vcf_in) as f:
    for line in f:
        if line.startswith("#"):
            # Header lines: print unchanged
            sys.stdout.write(line)
            continue

        line = line.rstrip("\n")
        if not line:
            continue

        fields = line.split("\t")
        if len(fields) < 10:
            # No sample columns: print unchanged
            sys.stdout.write(line + "\n")
            continue

        # FORMAT column tells us where GT is
        format_field = fields[8]
        format_keys = format_field.split(":")
        try:
            gt_index = format_keys.index("GT")
        except ValueError:
            # No GT field: print unchanged
            sys.stdout.write(line + "\n")
            continue

        # Recode GT in each sample column
        for i in range(9, len(fields)):
            sample = fields[i]
            if sample == ".":
                continue
            parts = sample.split(":")
            if gt_index < len(parts):
                parts[gt_index] = recode_gt(parts[gt_index])
            fields[i] = ":".join(parts)

        sys.stdout.write("\t".join(fields) + "\n")
