#!/usr/bin/env python3
import os
import argparse
from collections import Counter
import itertools
import math
import numpy as np
import pandas as pd

# ---------------------------
# Função para ler FASTA (pares de haplótipos)
# ---------------------------
def parse_fasta_pairs(path):
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

# ---------------------------
# Funções de estatísticas genéticas (sem HW/LD)
# ---------------------------
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
    s1 = sum(p*p for p in ps)
    s2 = sum(2*(ps[i]**2)*(ps[j]**2) for i in range(len(ps)) for j in range(i+1,len(ps)))
    return float(1.0 - s1 - s2)

def genotype_frequencies_under_hwe(freqs):
    alleles = sorted(freqs.keys())
    gf = {}
    for i, a in enumerate(alleles):
        gf[(a,a)] = freqs[a]*freqs[a]
        for j in range(i+1,len(alleles)):
            b = alleles[j]
            gf[(a,b)] = 2*freqs[a]*freqs[b]
    return gf

def match_probability_from_freqs(freqs):
    gf = genotype_frequencies_under_hwe(freqs)
    return float(sum(v*v for v in gf.values()))

def probability_all_different(freqs, m):
    alleles = list(freqs.keys())
    if m > len(alleles):
        return 0.0
    total = 0.0
    for combo in itertools.combinations(alleles, m):
        prod = 1.0
        for a in combo:
            prod *= freqs[a]
        total += prod
    from math import factorial
    return float(factorial(m)*total)

# ---------------------------
# Função otimizada: diversidade de alelos em combinações de indivíduos (amostragem)
# ---------------------------
def allele_diversity_combinations(genotypes_dict, n_individuals=2, n_samples=5000, seed=42):
    """
    Estima a % de combinações de indivíduos com 1..2*n alelos diferentes.
    Usa amostragem aleatória em vez de todas as combinações => rápido e leve.
    """
    rng = np.random.default_rng(seed)
    samples = list(genotypes_dict.keys())
    N = len(samples)
    if N < n_individuals:
        return {}

    counts = Counter()
    total = 0
    for _ in range(n_samples):
        chosen = rng.choice(samples, size=n_individuals, replace=False)
        alleles = set()
        for s in chosen:
            alleles.update(genotypes_dict[s])
        counts[len(alleles)] += 1
        total += 1

    result = {f"diff_{i}": counts[i]/total for i in range(1, 2*n_individuals+1)}
    return result

# ---------------------------
# Loop principal de janelas
# ---------------------------
def process_windows(samples, window_size, outdir):
    os.makedirs(outdir, exist_ok=True)
    sum_dir = os.path.join(outdir,"summaries")
    os.makedirs(sum_dir, exist_ok=True)
    allele_dir = os.path.join(outdir,"alleles")
    os.makedirs(allele_dir, exist_ok=True)

    lengths = [min(len(h1), len(h2)) for h1,h2 in samples.values()]
    if not lengths:
        raise ValueError("Nenhuma sequência válida encontrada.")
    L = min(lengths)

    results=[]
    for start in range(0,L-window_size+1):
        end = start+window_size
        genotypes={}
        alleles_pool=[]
        for s,(h1,h2) in samples.items():
            a1=h1[start:end]
            a2=h2[start:end]
            genotypes[s]=(a1,a2)
            alleles_pool.extend([a1,a2])
        allele_counter = Counter(alleles_pool)
        total = sum(allele_counter.values())
        freqs = {a: allele_counter[a]/total for a in allele_counter}
        num_alleles = len(freqs)
        mapping = {allele:f"A{i+1}" for i,allele in enumerate(freqs.keys())}
        allele_freqs_str = ", ".join([f"{mapping[a]}({freqs[a]:.6f})" for a in sorted(freqs,key=freqs.get,reverse=True)])

        # salva alelos
        allele_out = os.path.join(allele_dir,f"window_{start+1}_{end}_alleles.tsv")
        with open(allele_out,"w") as fh:
            fh.write("Allele\tRenamed\tFrequency\n")
            for orig,renamed in mapping.items():
                fh.write(f"{orig}\t{renamed}\t{freqs[orig]:.6f}\n")

        genotype_list = list(genotypes.values())
        He = heterozygosity_expected(freqs)
        Ho = heterozygosity_observed(genotype_list)
        Ae = effective_number_of_alleles(freqs)
        PIC = pic_from_freqs(freqs)
        MP = match_probability_from_freqs(freqs)
        PD = 1.0 - MP
        prob_2_diff = probability_all_different(freqs,2)
        prob_3_diff = probability_all_different(freqs,3)
        prob_4_diff = probability_all_different(freqs,4)

        # diversidade de alelos em combinações (amostragem)
        diversity_2 = allele_diversity_combinations(genotypes, n_individuals=2, n_samples=1000)
        diversity_3 = allele_diversity_combinations(genotypes, n_individuals=3, n_samples=1000)
        diversity_4 = allele_diversity_combinations(genotypes, n_individuals=4, n_samples=1000)

        # Score (sem LD/HW)
        metrics = {"He":He,"Ho":Ho,"Ae":Ae/10 if Ae else 0,"PIC":PIC,"PD":PD,
                   "prob_2_diff":prob_2_diff,"prob_3_diff":prob_3_diff,"prob_4_diff":prob_4_diff}
        metrics = {k: float(np.clip(v,0,1)) for k,v in metrics.items() if not math.isnan(v)}
        score = np.mean(list(metrics.values()))

        results.append({
            "start":start+1,"end":end,"num_alleles":num_alleles,"alleles_freqs":allele_freqs_str,
            "He":He,"Ho":Ho,"Ae":Ae,"PIC":PIC,
            "MP":MP,"PD":PD,"prob_2_diff":prob_2_diff,"prob_3_diff":prob_3_diff,"prob_4_diff":prob_4_diff,
            "Score":score,
            **{f"2ind_{k}":v for k,v in diversity_2.items()},
            **{f"3ind_{k}":v for k,v in diversity_3.items()},
            **{f"4ind_{k}":v for k,v in diversity_4.items()},
        })

        print(f"[OK] Janela {start+1}-{end} processada. He={He:.6f}, Ho={Ho:.6f}, Ae={Ae:.6f}, PIC={PIC:.6f}")

    df = pd.DataFrame(results)
    df.sort_values("Score",ascending=False,inplace=True)
    df.to_csv(os.path.join(sum_dir,"all_windows_summary.tsv"),sep="\t",index=False)
    return df

# ---------------------------
# Função main
# ---------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--window", type=int, required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    samples = parse_fasta_pairs(args.fasta)
    print(f"Carregadas {len(samples)} amostras.")
    
    df = process_windows(samples, args.window, args.outdir)
    
    print("\n"+"="*50)
    print("RESUMO DOS RESULTADOS:")
    print("="*50)
    print(f"Total de janelas processadas: {len(df)}")
    print(f"Arquivo principal salvo em: {os.path.join(args.outdir,'summaries','all_windows_summary.tsv')}")
    print(f"Arquivos de alelos salvos em: {os.path.join(args.outdir,'alleles/')}")

if __name__=="__main__":
    main()
