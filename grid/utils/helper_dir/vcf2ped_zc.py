#!/usr/bin/env python3
import sys
import numpy as np
import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Convert a phased VCF into PED/MAP files for iLASH."
    )
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--snp-count",
        type=int,
        required=False,
        help="(ignored — recomputed internally for safety)",
    )
    args = parser.parse_args()

    vcf_addr = args.vcf
    output_addr = args.out

    # ── Pass 1: find header, count genuine variant lines ──────────────────────
    with open(vcf_addr, "r") as vcf_file:
        header_line = None
        while True:
            header_line = vcf_file.readline()
            if not header_line:
                sys.exit("ERROR: no #CHROM header line found in VCF")
            if header_line.startswith("#CHROM"):
                break

        id_list = header_line.strip().split("\t")[9:]
        n_samples = len(id_list)

        snp_count = 0
        skipped = 0
        for line in vcf_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                skipped += 1
                continue
            snp_count += 1

    print(
        f"Samples: {n_samples}  |  Variants: {snp_count}  "
        f"|  Skipped stray header/blank lines: {skipped}",
        file=sys.stderr,
    )

    # ── Pass 2: parse genotypes ────────────────────────────────────────────────
    genotype_array = np.zeros((snp_count * 2, n_samples), dtype="U1")

    with open(vcf_addr, "r") as vcf_file, open(output_addr + ".map", "w") as map_file:

        while True:
            line = vcf_file.readline()
            if line.startswith("#CHROM"):
                break

        row = 0
        for line in vcf_file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            data = line.split("\t")
            if len(data) != 9 + n_samples:
                sys.exit(
                    f"ERROR: variant row {row} has {len(data)} fields, "
                    f"expected {9 + n_samples}. VCF may be malformed."
                )

            map_file.write(f"{data[0]}\t{data[2]}\t0\t{data[1]}\n")

            alleles = [data[3], data[4]]
            calls = data[9:]

            for idx, call_pair in enumerate(calls):
                pair = call_pair.split("|")
                if len(pair) != 2:
                    sys.exit(
                        f"ERROR: sample {id_list[idx]} at variant row {row} "
                        f"is not phased (got '{call_pair}')"
                    )
                genotype_array[2 * row, idx] = alleles[int(pair[0])]
                genotype_array[2 * row + 1, idx] = alleles[int(pair[1])]

            row += 1

    assert (
        row == snp_count
    ), f"Row count mismatch after parsing: parsed {row}, expected {snp_count}"

    with open(output_addr + ".ped", "w") as ped_file:
        for idx, sample_id in enumerate(id_list):
            column = genotype_array[:, idx]
            if "" in column:
                sys.exit(
                    f"ERROR: sample {sample_id} has incomplete genotype data "
                    f"— refusing to write a malformed ped row."
                )
            ped_file.write(
                f"{sample_id} {sample_id} 0 0 0 -9 " + " ".join(column) + "\n"
            )

    print(
        f"Wrote {output_addr}.ped and {output_addr}.map "
        f"({n_samples} samples, {snp_count} variants)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
