#!/usr/bin/env python3
"""
Build wild-type vs mutant protein sequences for any list of genes.

Inputs (all downloaded from the Galaxy history):
  --vcf     SnpEff-annotated SNP VCF            (Step 14 output, e.g. Galaxy241-...vcf)
  --wt      wild-type protein FASTA (gffread)  (Step 17 output, e.g. Galaxy246-...fasta)
  --mut     mutant protein FASTA (gffread)     (Step 16 output, e.g. Galaxy245-...fasta)
  --genes   text file, one gene symbol per line (e.g. "gene list.txt");
            optionally a transcript after the gene to force it, e.g. "Wwtr1 NM_133784.3"
  --out     output folder
Options:
  --all-transcripts   write every transcript of each gene, not just one per gene

How it works:
  1. For each gene, take the HIGH/MODERATE protein changes SnpEff wrote in the VCF ANN field
     (e.g. p.Arg245Leu on NM_011640.3).
  2. Choose one transcript per gene: the curated RefSeq NM_ transcript carrying the most changes
     (ties -> oldest accession, usually the canonical isoform). With --all-transcripts, every annotated transcript is used.
  3. Take the wild-type protein, put the new amino acid at each position, and check that
     the original residue matches and that the mutant FASTA differs only at those sites.
  4. Write one FASTA per gene (wild-type then mutant) + a combined FASTA + a summary TSV.
"""
import argparse, collections, csv, os, re, sys

AA = dict(Ala='A', Arg='R', Asn='N', Asp='D', Cys='C', Gln='Q', Glu='E', Gly='G', His='H', Ile='I',
          Leu='L', Lys='K', Met='M', Phe='F', Pro='P', Ser='S', Thr='T', Trp='W', Tyr='Y', Val='V',
          Ter='*', Sec='U')
HGVS = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)$")  # single amino-acid substitutions


def read_fasta(path):
    seqs, cur = collections.OrderedDict(), None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs.setdefault(cur, [])
            elif cur:
                seqs[cur].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def wrap(s, n=60):
    return "\n".join(s[i:i + n] for i in range(0, len(s), n))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--wt", required=True)
    ap.add_argument("--mut", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--all-transcripts", action="store_true")
    a = ap.parse_args()

    genes, forced = [], {}
    with open(a.genes, encoding="utf-8-sig") as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            genes.append(parts[0])
            if len(parts) > 1:          # optional 2nd column = transcript to use, e.g. "Wwtr1 NM_133784.3"
                forced[parts[0]] = parts[1]
    wanted = set(genes)
    wt, mut = read_fasta(a.wt), read_fasta(a.mut)
    os.makedirs(a.out, exist_ok=True)

    # variants[gene][transcript] -> list of variant dicts
    variants = collections.defaultdict(lambda: collections.defaultdict(list))
    other = collections.defaultdict(list)  # protein changes that are not simple substitutions
    with open(a.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
            fmt = dict(zip(f[8].split(":"), f[9].split(":"))) if len(f) > 9 else {}
            dp, ao = fmt.get("DP", ""), fmt.get("AO", "")
            vaf = round(int(ao) / int(dp), 2) if dp.isdigit() and ao.isdigit() and int(dp) else ""
            for ann in info.get("ANN", "").split(","):
                p = ann.split("|")
                if len(p) < 11 or p[3] not in wanted or p[2] not in ("HIGH", "MODERATE") or not p[10]:
                    continue
                v = dict(chrom=f[0], pos=f[1], ref=f[3], alt=f[4], effect=p[1], impact=p[2], tx=p[6],
                         c=p[9], p=p[10], gt=fmt.get("GT", ""), dp=dp, ao=ao, vaf=vaf, qual=f[5])
                m = HGVS.match(v["p"])
                if m and m.group(1) in AA and (m.group(3) in AA or m.group(3) == "*"):
                    v["aa_ref"], v["aa_pos"], v["aa_alt"] = AA[m.group(1)], int(m.group(2)), AA.get(m.group(3), "*")
                    variants[p[3]][p[6]].append(v)
                else:
                    other[p[3]].append(v)

    summary_rows, report = [], []
    combined = open(os.path.join(a.out, "selected_genes_WT_and_mutant_proteins.fasta"), "w")
    for g in genes:
        txs = {t: vs for t, vs in variants[g].items() if t in wt}
        if not txs:
            report.append(f"{g}: no protein-changing SNV found" + (" (other changes listed in summary)" if other[g] else ""))
            for v in other[g]:
                summary_rows.append([g, v["tx"], v["chrom"], v["pos"], v["ref"], v["alt"], v["effect"], v["c"], v["p"],
                                     "", v["gt"], v["dp"], v["ao"], v["vaf"], v["qual"], "", "", "not in FASTA (not a single-residue change)"])
            continue
        if a.all_transcripts:
            chosen = sorted(txs)
        else:
            def score(t):
                # curated NM_ first, then most changes, then the oldest accession (usually the canonical one)
                num = int(re.sub(r"\D", "", t.split(".")[0]) or 0)
                return (t.startswith("NM_"), len({x["aa_pos"] for x in txs[t]}), -num)
            if g in forced and forced[g] in txs:
                chosen = [forced[g]]
            else:
                if g in forced:
                    report.append(f"{g}: requested transcript {forced[g]} has no change, choosing automatically")
                chosen = [max(txs, key=score)]
        out_fh = open(os.path.join(a.out, f"{g}_WT_vs_mutant.fasta"), "w")
        for t in chosen:
            w = wt[t]
            m = list(w)
            uniq = {}
            for v in txs[t]:
                uniq[v["aa_pos"]] = v
            labels, notes = [], []
            for pos in sorted(uniq):
                v = uniq[pos]
                ok = pos <= len(w) and w[pos - 1] == v["aa_ref"]
                if not ok:
                    notes.append(f"{v['p']}: wild-type residue mismatch, skipped")
                    continue
                m[pos - 1] = v["aa_alt"]
                lab = f"{v['aa_ref']}{pos}{v['aa_alt']}"
                labels.append(lab)
                lo, hi = max(0, pos - 11), min(len(w), pos + 10)
                mut_check = mut.get(t, "")
                site = mut_check[pos - 1] if len(mut_check) >= pos else "?"
                check = "ok" if site in ("X", v["aa_alt"]) else f"mutant FASTA shows '{site}'"
                summary_rows.append([g, t, v["chrom"], v["pos"], v["ref"], v["alt"], v["effect"], v["c"], v["p"], lab,
                                     v["gt"], v["dp"], v["ao"], v["vaf"], v["qual"], w[lo:hi], "".join(m)[lo:hi], check])
            mseq = "".join(m)
            # stop-gained: cut the mutant protein at the new stop
            if "*" in mseq:
                mseq = mseq[:mseq.index("*")]
            lab = ",".join(labels) if labels else "none"
            for h, s in ((f">{g}|{t}|wild-type|{len(w)}aa", w), (f">{g}|{t}|mutant|{lab}|{len(mseq)}aa", mseq)):
                out_fh.write(h + "\n" + wrap(s) + "\n")
                combined.write(h + "\n" + wrap(s) + "\n")
            report.append(f"{g} {t}: {lab}" + (f"  [{'; '.join(notes)}]" if notes else ""))
        out_fh.close()
        for v in other[g]:
            summary_rows.append([g, v["tx"], v["chrom"], v["pos"], v["ref"], v["alt"], v["effect"], v["c"], v["p"],
                                 "", v["gt"], v["dp"], v["ao"], v["vaf"], v["qual"], "", "", "not in FASTA (not a single-residue change)"])
    combined.close()

    with open(os.path.join(a.out, "selected_genes_mutation_summary.tsv"), "w", newline="") as fh:
        cw = csv.writer(fh, delimiter="\t")
        cw.writerow(["Gene", "Transcript", "Chrom", "Position", "Ref_base", "Alt_base", "Effect", "DNA_change",
                     "Protein_change_HGVS", "Protein_change_short", "Genotype", "Read_depth", "Alt_reads",
                     "Allele_fraction", "QUAL", "WT_21mer", "Mutant_21mer", "Check"])
        cw.writerows(summary_rows)
    print("\n".join(report))


if __name__ == "__main__":
    main()
