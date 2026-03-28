#!/usr/bin/env python3
import argparse, os, subprocess, sys, json
from pathlib import Path
import pysam
import pandas as pd
from collections import defaultdict
from Bio import SeqIO

def run(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        print("ERROR:", " ".join(cmd), file=sys.stderr)
        print(p.stderr, file=sys.stderr)
        sys.exit(p.returncode)
    return p

def bowtie2_index(fasta, prefix):
    exts = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
    if not all((Path(prefix + ext).exists() for ext in exts)):
        run(["bowtie2-build", fasta, prefix])

def map_pe(index_prefix, r1, r2, out_bam, threads=4):
    sam_tmp = out_bam.with_suffix(".sam")
    cmd = ["bowtie2", "-x", index_prefix, "-1", r1, "-2", r2, "-p", str(threads), "--very-sensitive"]
    with open(sam_tmp, "w") as sam:
        run(cmd + ["-S", str(sam_tmp)])
    run(["samtools", "sort", "-@","2", "-o", str(out_bam), str(sam_tmp)])
    sam_tmp.unlink(missing_ok=True)
    run(["samtools", "index", str(out_bam)])

def fasta_lengths(fasta):
    return {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta, "fasta")}

def per_target_coverage_identity(bam_path, target_lengths):
    bam = pysam.AlignmentFile(bam_path)
    covered = {t: set() for t in target_lengths}
    stats = defaultdict(lambda: {"mm":0, "alen":0})
    for aln in bam.fetch(until_eof=True):
        if aln.is_unmapped or aln.reference_name is None:
            continue
        t = aln.reference_name
        if t not in target_lengths:
            continue
        for qpos, rpos in aln.get_aligned_pairs(matches_only=True):
            if rpos is not None:
                covered[t].add(rpos)
        nm = aln.get_tag("NM") if aln.has_tag("NM") else 0
        alen = aln.query_alignment_length
        stats[t]["mm"] += nm
        stats[t]["alen"] += alen
    rows = []
    for t, L in target_lengths.items():
        cov_breadth = len(covered[t]) / L if L > 0 else 0.0
        mm = stats[t]["mm"]; alen = stats[t]["alen"]
        ident = max(0.0, 1.0 - (mm/alen)) if alen > 0 else 0.0
        rows.append({"target": t, "length": L, "breadth": cov_breadth, "identity": ident})
    return pd.DataFrame(rows)

def call_genes(df, min_breadth=0.90, min_ident=0.90):
    hits = df[(df["breadth"] >= min_breadth) & (df["identity"] >= min_ident)].copy()
    return hits.sort_values(["target"])

def parse_mlst_alleles(fasta):
    mapping = defaultdict(dict)
    for rec in SeqIO.parse(fasta, "fasta"):
        if "_" in rec.id:
            locus, allele = rec.id.split("_", 1)
            mapping[locus][rec.id] = len(rec.seq)
    return mapping

def best_allele_per_locus(df, mlst_loci):
    best = {}
    for locus in mlst_loci:
        sub = df[df["target"].str.startswith(locus + "_")]
        if sub.empty:
            best[locus] = None
            continue
        s = sub.sort_values(["identity","breadth"], ascending=[False, False]).iloc[0]
        best[locus] = s["target"]
    return best

def allele_name_to_number(allele_id):
    try:
        return int(allele_id.split("_",1)[1])
    except Exception:
        return None

def load_profiles(tsv):
    prof = pd.read_csv(tsv, sep="\t", dtype=str)
    prof.columns = [c.strip() for c in prof.columns]
    return prof

def match_ST(best_calls, profiles):
    row = {}
    for locus, allele in best_calls.items():
        row[locus] = str(allele_name_to_number(allele)) if allele else None
    loci = [c for c in profiles.columns if c != "ST"]
    candidates = profiles.copy()
    for locus in loci:
        if row.get(locus) is None:
            return {"ST":"ND", **row}
        candidates = candidates[candidates[locus] == row[locus]]
        if candidates.empty:
            return {"ST":"Novel/Partial", **row}
    st = candidates.iloc[0]["ST"]
    return {"ST": st, **row}

def main():
    import argparse
    ap = argparse.ArgumentParser(description="SRST2-like caller in Python 3")
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--genes_fasta", required=True)
    ap.add_argument("--mlst_fasta", required=True)
    ap.add_argument("--mlst_profiles", required=True)
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min_breadth", type=float, default=0.90)
    ap.add_argument("--min_ident", type=float, default=0.90)
    args = ap.parse_args()

    work = Path(args.out_prefix).parent
    work.mkdir(parents=True, exist_ok=True)

    gene_index = str(Path(args.genes_fasta).with_suffix(""))
    bowtie2_index(args.genes_fasta, gene_index)
    gene_bam = Path(f"{args.out_prefix}.genes.bam")
    map_pe(gene_index, args.r1, args.r2, gene_bam, threads=args.threads)
    genes_df = per_target_coverage_identity(gene_bam, fasta_lengths(args.genes_fasta))
    calls_df = call_genes(genes_df, args.min_breadth, args.min_ident)
    genes_df.to_csv(f"{args.out_prefix}.genes_metrics.tsv", sep="\t", index=False)
    calls_df.to_csv(f"{args.out_prefix}.genes_calls.tsv", sep="\t", index=False)

    mlst_index = str(Path(args.mlst_fasta).with_suffix(""))
    bowtie2_index(args.mlst_fasta, mlst_index)
    mlst_bam = Path(f"{args.out_prefix}.mlst.bam")
    map_pe(mlst_index, args.r1, args.r2, mlst_bam, threads=args.threads)
    mlst_df = per_target_coverage_identity(mlst_bam, fasta_lengths(args.mlst_fasta))
    mlst_df.to_csv(f"{args.out_prefix}.mlst_metrics.tsv", sep="\t", index=False)

    loci_map = parse_mlst_alleles(args.mlst_fasta)
    loci = sorted(loci_map.keys())
    best = best_allele_per_locus(mlst_df, loci)
    profiles = load_profiles(args.mlst_profiles)
    st_row = match_ST(best, profiles)

    with open(f"{args.out_prefix}.mlst_calls.json","w") as fh:
        json.dump(st_row, fh, indent=2)

    print(f"[OK] Wrote:\\n  {args.out_prefix}.genes_metrics.tsv\\n  {args.out_prefix}.genes_calls.tsv\\n  {args.out_prefix}.mlst_metrics.tsv\\n  {args.out_prefix}.mlst_calls.json")

if __name__ == "__main__":
    main()
