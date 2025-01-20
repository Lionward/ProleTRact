import concurrent.futures
import pandas as pd 
import numpy as np
import importlib

import pysam

import matplotlib.pyplot as plt
# from upsetplot import UpSet

import parse_vcf_files
importlib.reload(parse_vcf_files)

import parse_fasta
importlib.reload(parse_fasta)

import exact_motif_counts
importlib.reload(exact_motif_counts)

import compare_with_asm
importlib.reload(compare_with_asm)

import plots
importlib.reload(plots)

import argparse

def main():
    parser = argparse.ArgumentParser(description="Process VCF, assembly, and output files.")
    
    # Add arguments
    parser.add_argument("--vcf", required=True, help="Path to the VCF file")
    parser.add_argument("--asm1", required=True, help="Path to the assembly file (haplotype 1)")
    parser.add_argument("--asm2", required=True, help="Path to the assembly file (haplotype 2)")
    parser.add_argument("--out", required=True, help="Path to the output file")
    parser.add_argument("--tool", required=True, help="Name of the TR genotyping tool used to produce the input VCF file")
    parser.add_argument("--log", required=True, help="Path to the log file")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Display parsed arguments
    print(f"VCF file: {args.vcf}")
    print(f"Assembly file: {args.asm1}")
    print(f"Assembly file: {args.asm2}")
    print(f"Output file: {args.out}")
    print(f"Log file: {args.log}")

    asm_h1 = parse_fasta.parse_fasta_reads_cut(fasta_file = args.asm1, seq_col_name = "asm_h1_seq")
    asm_h2 = parse_fasta.parse_fasta_reads_cut(fasta_file = args.asm2, seq_col_name = "asm_h2_seq")

    asm_df = pd.merge(asm_h1, asm_h2, on=['chrom', 'start', 'end'], how='inner')

    TR_df = pd.DataFrame()

    if args.tool == "TandemTwist":
        TR_df = parse_vcf_files.parse_vcf_tandemtwister(pysam.VariantFile(args.vcf))
    elif args.tool == "vamos":
        TR_df = parse_vcf_files.parse_vcf_vamos(pysam.VariantFile(args.vcf))
    elif args.tool == "TRGT":
        TR_df = TR_df = parse_vcf_files.parse_vcf_TRGT(pysam.VariantFile(args.vcf))
    
    TRs_df_compare = compare_with_asm.compare_sample_to_asm(TR_df.head(100000), asm_df)
    
    TRs_df_compare.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()