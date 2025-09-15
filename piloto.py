#!/usr/bin/env python3
"""
microhap_windows.py - versão refinada

- Coluna alleles_freqs renomeia os alelos como A1, A2...
- Hardy-Weinberg: apenas uma coluna HWE (p-valor, 0 a 1)
- Remove coluna window_file
- Cria Score (0-1) para ranquear loci por utilidade em deconvolução de misturas
- Ordena resultado pelo Score
- Salva arquivos separados com alelos/frequências por microhaplótipo
- Mostra progresso no terminal
"""

import os
import argparse
from collections import Counter, defaultdict
import itertools
import math
import numpy as np
import pandas as pd
from scipy.stats import chi2


def parse_fasta_pairs(path):
    """Lê FASTA no formato SAMPLE_h1 / SAMPLE_h2."""
    records = {}
    cur_header = None
    cur_seq = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_header is not None:
                    records[cur_header] = "".join(cur_seq)
                cur_header = line[1:].strip()
                cur_seq = []
            else:
                if "_h" in line and set(line.lower()) <= set("abcdefghijklmnopqrstuvwxyz0123456789_"):
                    if cur_header is not None:
                        records[cur_header] = "".join(cur_seq)
                    cur_header = line.strip()
                    cur_seq = []
                else:
                    cur_seq.append(line.strip())
        if cur_header is not None:
            records[cur_header] = "".join(cur_seq)

    samples = {}
    for hdr, seq in records.items():
        if "_h" in hdr:
            base = hdr.split("_h")[0]
            hap = hdr.split("_h")[1]
            samples.setdefault(base, {})[hap] = seq
        else:
            if "_" in hdr:
                base, tail = hdr.rsplit("_", 1)
                samples.setdefault(base, {})[tail] = seq
            else:
                samples.setdefault(hdr, {})["1"] = seq

    final = {}
    for s, d in samples.items():
        h1 = d.get("1") or d.get("h1") or list(d.values())[0]
        h2 = d.get("2") or d.get("h2") or (list(d.values())[1] if len(d) > 1 else h1)
        final[s] = (h1.upper(), h2.upper())
    return final


def heterozygosity_expected(freqs):
    ps = np.array(list(freqs.values()))
    return float(1.0 - np.sum(ps * ps))


def heterozygosity_observed(genotypes):
    n = len(genotypes)
    if n == 0:
        return np.nan
    het = sum(1 for (a1, a2) in genotypes if a1 != a2)
    return het / n


def effective_number_of_alleles(freqs):
    ps = np.array(list(freqs.values()))
    s = np.sum(ps * ps)
    return float(1.0 / s) if s > 0 else np.nan


def pic_from_freqs(freqs):
    ps = list(freqs.values())
    s1 = sum(p * p for p in ps)
    s2 = 0.0
    for i in range(len(ps)):
        for j in range(i + 1, len(ps)):
            s2 += 2 * (ps[i] ** 2) * (ps[j] ** 2)
    return float(1.0 - s1 - s2)


def genotype_frequencies_under_hwe(freqs):
    alleles = sorted(freqs.keys())
    gf = {}
    for i, a in enumerate(alleles):
        gf[(a, a)] = freqs[a] * freqs[a]
        for j in range(i + 1, len(alleles)):
            b = alleles[j]
            gf[(a, b)] = 2 * freqs[a] * freqs[b]
    return gf


def match_probability_from_freqs(freqs):
    gf = genotype_frequencies_under_hwe(freqs)
    return float(sum(v * v for v in gf.values()))


def probability_all_different(freqs, m):
    from math import factorial

    alleles = list(freqs.keys())
    if m > len(alleles):
        return 0.0
    total = 0.0
    for combo in itertools.combinations(alleles, m):
        prod = 1.0
        for a in combo:
            prod *= freqs[a]
        total += prod
    return float(factorial(m) * total)


def hardy_weinberg_p(genotypes, freqs):
    """Retorna p-valor do teste chi2 de HWE."""
    from collections import Counter
    geno_obs = Counter()
    for a, b in genotypes:
        g = (a, a) if a == b else tuple(sorted((a, b)))
        geno_obs[g] += 1
    n = sum(geno_obs.values())
    if n == 0:
        return np.nan
    geno_exp_freq = genotype_frequencies_under_hwe(freqs)
    keys = sorted(set(geno_obs.keys()) | set(geno_exp_freq.keys()))
    obs = np.array([geno_obs.get(k, 0) for k in keys], dtype=float)
    exp = np.array([geno_exp_freq.get(k, 0) * n for k in keys], dtype=float)
    mask = exp > 0
    if not mask.any():
        return np.nan
    chi2_stat = np.sum(((obs[mask] - exp[mask]) ** 2) / exp[mask])
    df = np.sum(mask) - len(freqs)
    if df <= 0:
        return np.nan
    from scipy.stats import chi2
    pval = 1 - chi2.cdf(chi2_stat, df)
    return float(max(0.0, min(1.0, pval)))


def process_windows(samples, window_size, outdir):
    os.makedirs(outdir, exist_ok=True)
    sum_dir = os.path.join(outdir, "summaries")
    os.makedirs(sum_dir, exist_ok=True)
    allele_dir = os.path.join(outdir, "alleles")
    os.makedirs(allele_dir, exist_ok=True)

    lengths = [min(len(h1), len(h2)) for h1, h2 in samples.values()]
    if not lengths:
        raise ValueError("Nenhuma sequência válida encontrada.")
    L = min(lengths)

    results = []
    for start in range(0, L - window_size + 1):
        end = start + window_size
        genotypes = {}
        alleles_pool = []
        for s, (h1, h2) in samples.items():
            a1 = h1[start:end]
            a2 = h2[start:end]
            genotypes[s] = (a1, a2)
            alleles_pool.extend([a1, a2])

        allele_counter = Counter(alleles_pool)
        total = sum(allele_counter.values())
        freqs = {a: allele_counter[a] / total for a in allele_counter}
        num_alleles = len(freqs)

        # renomear para A1, A2...
        mapping = {allele: f"A{i+1}" for i, allele in enumerate(freqs.keys())}
        allele_freqs_str = ", ".join(
            [f"{mapping[a]}({freqs[a]:.4f})" for a in sorted(freqs, key=freqs.get, reverse=True)]
        )

        # salvar arquivo separado com frequências
        allele_out = os.path.join(allele_dir, f"window_{start+1}_{end}_alleles.tsv")
        with open(allele_out, "w") as fh:
            fh.write("Allele\tRenamed\tFrequency\n")
            for orig, renamed in mapping.items():
                fh.write(f"{orig}\t{renamed}\t{freqs[orig]:.6f}\n")

        # calcular métricas
        He = heterozygosity_expected(freqs)
        Ho = heterozygosity_observed(list(genotypes.values()))
        Neff = effective_number_of_alleles(freqs)
        PIC = pic_from_freqs(freqs)
        MP = match_probability_from_freqs(freqs)
        PD = 1.0 - MP
        prob_2_diff = probability_all_different(freqs, 2)
        prob_3_diff = probability_all_different(freqs, 3)
        prob_4_diff = probability_all_different(freqs, 4)
        HWE = hardy_weinberg_p(list(genotypes.values()), freqs)

        # criar score apenas com os selecionados
        metrics = [He, Ho, (Neff or 0) / 10, PIC, PD, prob_4_diff]
        metrics = [min(1, max(0, m)) for m in metrics if not math.isnan(m)]
        score = float(np.mean(metrics)) if metrics else 0

        results.append(
            {
                "start_1based": start + 1,
                "end_1based": end,
                "num_alleles": num_alleles,
                "alleles_freqs": allele_freqs_str,
                "He": He,
                "Ho": Ho,
                "Neff": Neff,
                "PIC": PIC,
                "HWE": HWE,
                "MP": MP,
                "PD": PD,
                "prob_2_diff": prob_2_diff,
                "prob_3_diff": prob_3_diff,
                "prob_4_diff": prob_4_diff,
                "Score": score,
            }
        )

        print(f"[OK] Janela {start+1}-{end} processada.")

    df = pd.DataFrame(results)
    df.sort_values("Score", ascending=False, inplace=True)
    df.to_csv(os.path.join(sum_dir, "all_windows_summary.tsv"), sep="\t", index=False)
    return df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--window", type=int, required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    samples = parse_fasta_pairs(args.fasta)
    df = process_windows(samples, args.window, args.outdir)
    print("Feito. Resumo salvo em:", os.path.join(args.outdir, "summaries", "all_windows_summary.tsv"))
    print(f"{len(df)} janelas processadas.")


if __name__ == "__main__":
    main()
