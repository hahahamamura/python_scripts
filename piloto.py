#!/usr/bin/env python3
import os
import argparse
from collections import Counter
import itertools
import math
import numpy as np
import pandas as pd
from scipy.stats import chi2
import tempfile

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
# Funções de estatísticas genéticas
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

def hardy_weinberg_p(genotypes, freqs):
    from collections import Counter
    geno_obs = Counter()
    for a,b in genotypes:
        g = (a,a) if a==b else tuple(sorted((a,b)))
        geno_obs[g] += 1
    n = sum(geno_obs.values())
    if n==0:
        return np.nan
    geno_exp_freq = genotype_frequencies_under_hwe(freqs)
    keys = sorted(set(geno_obs.keys())|set(geno_exp_freq.keys()))
    obs = np.array([geno_obs.get(k,0) for k in keys], dtype=float)
    exp = np.array([geno_exp_freq.get(k,0)*n for k in keys], dtype=float)
    mask = exp>0
    if not mask.any():
        return np.nan
    k = len(freqs)
    num_genotypes = np.sum(mask)
    df = num_genotypes - k
    if df <= 0:
        return np.nan
    chi2_stat = np.sum(((obs[mask]-exp[mask])**2)/exp[mask])
    pval = 1 - chi2.cdf(chi2_stat, df)
    return float(max(0.0, min(1.0,pval)))

def hardy_weinberg_exact(genotypes, freqs):
    if not genotypes:
        return np.nan
    n = len(genotypes)
    homo_obs = sum(1 for a,b in genotypes if a==b)
    hetero_obs = n - homo_obs
    He = heterozygosity_expected(freqs)
    hetero_exp = He*n
    homo_exp = n - hetero_exp
    if homo_exp>0 and hetero_exp>0:
        chi2_stat = ((homo_obs-homo_exp)**2/homo_exp + (hetero_obs-hetero_exp)**2/hetero_exp)
        pval = 1 - chi2.cdf(chi2_stat, 1)
        return float(max(0.0, min(1.0,pval)))
    else:
        return np.nan

def inbreeding_coefficient(genotypes, freqs):
    He = heterozygosity_expected(freqs)
    Ho = heterozygosity_observed(genotypes)
    if He==0:
        return 0
    return float((He-Ho)/He)

# ---------------------------
# Função LD
# ---------------------------
def calculate_linkage_disequilibrium(haplotypes, pos1, pos2):
    if pos1>=len(haplotypes[0]) or pos2>=len(haplotypes[0]):
        return {"r2": np.nan,"D_prime": np.nan,"D": np.nan}
    haplo_counts = Counter((hap[pos1], hap[pos2]) for hap in haplotypes)
    n = len(haplotypes)
    if n==0:
        return {"r2": np.nan,"D_prime": np.nan,"D": np.nan}
    allele1_counts = Counter(hap[pos1] for hap in haplotypes)
    allele2_counts = Counter(hap[pos2] for hap in haplotypes)
    if len(allele1_counts)==1 or len(allele2_counts)==1:
        return {"r2":0.0,"D_prime":0.0,"D":0.0}
    top1 = allele1_counts.most_common(2)
    top2 = allele2_counts.most_common(2)
    if len(top1)<2 or len(top2)<2:
        return {"r2":0.0,"D_prime":0.0,"D":0.0}
    a1,a2 = top1[0][0], top1[1][0]
    b1,b2 = top2[0][0], top2[1][0]
    pa1,pa2 = top1[0][1]/n, top1[1][1]/n
    pb1,pb2 = top2[0][1]/n, top2[1][1]/n
    pab11 = haplo_counts.get((a1,b1),0)/n
    pab12 = haplo_counts.get((a1,b2),0)/n
    pab21 = haplo_counts.get((a2,b1),0)/n
    pab22 = haplo_counts.get((a2,b2),0)/n
    D_vals = [pab11-pa1*pb1,pab12-pa1*pb2,pab21-pa2*pb1,pab22-pa2*pb2]
    D = max(D_vals,key=abs)
    if D>=0:
        Dmax = min(pa1*pb1,pa2*pb2) if D==D_vals[0] or D==D_vals[3] else min(pa1*pb2,pa2*pb1)
    else:
        Dmax = max(-pa1*pb1,-pa2*pb2) if D==D_vals[0] or D==D_vals[3] else max(-pa1*pb2,-pa2*pb1)
    D_prime = abs(D)/abs(Dmax) if Dmax!=0 else 0.0
    denom = pa1*pa2*pb1*pb2
    r2 = (D*D)/denom if denom>0 else 0.0
    return {"r2":float(max(0.0,min(1.0,r2))), "D_prime":float(max(0.0,min(1.0,D_prime))), "D":float(D)}

def calculate_mean_ld_in_window(haplotypes, window_size, sample_pairs=20):
    if not haplotypes or window_size<2:
        return {"r2_mean":np.nan,"D_prime_mean":np.nan}
    positions = list(range(window_size))
    if len(positions)>10:
        pairs=[]
        for _ in range(min(sample_pairs, len(positions)*(len(positions)-1)//2)):
            pos1,pos2 = np.random.choice(positions,2,replace=False)
            pairs.append((min(pos1,pos2),max(pos1,pos2)))
        pairs=list(set(pairs))
    else:
        pairs = [(i,j) for i in range(len(positions)) for j in range(i+1,len(positions))]
    if not pairs:
        return {"r2_mean":np.nan,"D_prime_mean":np.nan}
    r2_values=[]
    d_prime_values=[]
    for pos1,pos2 in pairs:
        ld_result = calculate_linkage_disequilibrium(haplotypes,pos1,pos2)
        if not np.isnan(ld_result["r2"]):
            r2_values.append(ld_result["r2"])
        if not np.isnan(ld_result["D_prime"]):
            d_prime_values.append(ld_result["D_prime"])
    r2_mean = np.mean(r2_values) if r2_values else np.nan
    d_prime_mean = np.mean(d_prime_values) if d_prime_values else np.nan
    return {"r2_mean":float(r2_mean) if not np.isnan(r2_mean) else np.nan,
            "D_prime_mean":float(d_prime_mean) if not np.isnan(d_prime_mean) else np.nan}

# ---------------------------
# Função nova: diversidade de alelos em combinações de indivíduos (streaming no disco)
# ---------------------------
def allele_diversity_combinations(genotypes_dict, n_individuals=2, tmpdir="./tmp_combos"):
    """
    Calcula a % de combinações de indivíduos com 1,2,...,2*n alelos diferentes.
    Usa arquivos temporários em disco para não explodir a RAM.
    """
    os.makedirs(tmpdir, exist_ok=True)
    samples = list(genotypes_dict.keys())
    if len(samples) < n_individuals:
        return {}

    # arquivo temporário com todas as combinações
    tmpfile = os.path.join(tmpdir, f"combos_{n_individuals}.txt")
    with open(tmpfile, "w") as fh:
        for combo in itertools.combinations(samples, n_individuals):
            fh.write(",".join(combo) + "\n")

    counts = Counter()
    total = 0

    # lê linha a linha e calcula diversidade
    with open(tmpfile) as fh:
        for line in fh:
            combo = line.strip().split(",")
            alleles = set()
            for s in combo:
                alleles.update(genotypes_dict[s])
            counts[len(alleles)] += 1
            total += 1

    # remove o arquivo temporário (se quiser manter, basta comentar a linha abaixo)
    os.remove(tmpfile)

    if total == 0:
        return {}

    result = {f"diff_{i}": counts[i] / total for i in counts}
    # garante todas as chaves até 2*n_individuals
    for i in range(1, 2*n_individuals + 1):
        if f"diff_{i}" not in result:
            result[f"diff_{i}"] = 0.0
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
        haplotypes=[]
        for s,(h1,h2) in samples.items():
            a1=h1[start:end]
            a2=h2[start:end]
            genotypes[s]=(a1,a2)
            alleles_pool.extend([a1,a2])
            haplotypes.extend([a1,a2])
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

        # LD
        np.random.seed(42)
        ld_results = calculate_mean_ld_in_window(haplotypes, window_size)
        r2_mean = ld_results["r2_mean"]
        d_prime_mean = ld_results["D_prime_mean"]

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
        HWE = hardy_weinberg_p(genotype_list,freqs)
        HWE_simple = hardy_weinberg_exact(genotype_list,freqs)
        Fis = inbreeding_coefficient(genotype_list,freqs)

        # diversidade de alelos em combinações de 2,3,4 indivíduos
        diversity_2 = allele_diversity_combinations(genotypes, n_individuals=2, tmpdir=os.path.join(outdir, "tmp"))
        diversity_3 = allele_diversity_combinations(genotypes, n_individuals=3, tmpdir=os.path.join(outdir, "tmp"))
        diversity_4 = allele_diversity_combinations(genotypes, n_individuals=4, tmpdir=os.path.join(outdir, "tmp"))


        # Score
        metrics = {"He":He,"Ho":Ho,"Ae":Ae/10 if Ae else 0,"PIC":PIC,"PD":PD,
                   "prob_2_diff":prob_2_diff,"prob_3_diff":prob_3_diff,"prob_4_diff":prob_4_diff,
                   "r2_mean":r2_mean if not np.isnan(r2_mean) else 0,
                   "D_prime_mean":d_prime_mean if not np.isnan(d_prime_mean) else 0}
        metrics = {k: float(np.clip(v,0,1)) for k,v in metrics.items() if not math.isnan(v)}
        ld_penalty = 1 - metrics.get("r2_mean",0)
        score = np.mean([metrics["He"],metrics["Ho"],metrics["Ae"],metrics["PIC"],metrics["PD"],
                         0.5*metrics.get("prob_2_diff",0),0.75*metrics.get("prob_3_diff",0),
                         metrics.get("prob_4_diff",0),0.3*ld_penalty])

        results.append({
            "start":start+1,"end":end,"num_alleles":num_alleles,"alleles_freqs":allele_freqs_str,
            "He":He,"Ho":Ho,"Ae":Ae,"PIC":PIC,"HWE":HWE,"HWE_simple":HWE_simple,"Fis":Fis,
            "MP":MP,"PD":PD,"prob_2_diff":prob_2_diff,"prob_3_diff":prob_3_diff,"prob_4_diff":prob_4_diff,
            "r2_mean":r2_mean,"D_prime_mean":d_prime_mean,"Score":score,
            **{f"2ind_{k}":v for k,v in diversity_2.items()},
            **{f"3ind_{k}":v for k,v in diversity_3.items()},
            **{f"4ind_{k}":v for k,v in diversity_4.items()},
        })

        r2_str = f"{r2_mean:.4f}" if not np.isnan(r2_mean) else "NA"
        d_prime_str = f"{d_prime_mean:.4f}" if not np.isnan(d_prime_mean) else "NA"
        print(f"[OK] Janela {start+1}-{end} processada. HWE={HWE:.6f}, Fis={Fis:.6f}, r²={r2_str}, D'={d_prime_str}")

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
