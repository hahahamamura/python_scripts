#!/usr/bin/env python3
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

# URL base fornecida
BASE_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"

# Cromossomos 1-22 + X
chromosomes = [str(i) for i in range(1, 23)] + ["X"]

# Pasta de saída
OUT_DIR = "/home/lab/Downloads/VCF_WHOLEGENOME"
os.makedirs(OUT_DIR, exist_ok=True)

def download_with_wget(url, out_path):
    """
    Executa o comando wget para baixar com opção -c (continue).
    Retorna True se sucesso, False caso contrário.
    """
    cmd = [
        "wget",
        "-c",   # continue downloads interrompidos
        "-O", out_path,
        url
    ]
    try:
        subprocess.run(cmd, check=True)
        print(f"[OK] {url} → {out_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[ERRO] wget falhou para {url}: {e}")
        return False

def download_chr(chrom):
    """
    Para um cromossomo, monta os nomes dos arquivos .vcf.gz e .tbi e
    baixa ambos usando wget.
    """
    if chrom == "X":
        vcf_name = f"1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    else:
        vcf_name = f"1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    tbi_name = vcf_name + ".tbi"

    urls = [
        f"{BASE_URL}/{vcf_name}",
        f"{BASE_URL}/{tbi_name}"
    ]
    out_paths = [
        os.path.join(OUT_DIR, vcf_name),
        os.path.join(OUT_DIR, tbi_name)
    ]

    results = []
    for url, out_path in zip(urls, out_paths):
        # Se já existe completamente (vcf + tbi), você pode pular, mas wget -c lida com isso.
        res = download_with_wget(url, out_path)
        results.append(res)
    return results

def main():
    # usar 2 downloads concorrentes ao mesmo tempo
    with ThreadPoolExecutor(max_workers=2) as executor:
        future_to_chr = {executor.submit(download_chr, chrom): chrom for chrom in chromosomes}
        for future in as_completed(future_to_chr):
            chrom = future_to_chr[future]
            try:
                res = future.result()
                print(f"[DONE] Cromossomo {chrom}: {res}")
            except Exception as exc:
                print(f"[ERRO] Cromossomo {chrom} gerou exceção: {exc}")

if __name__ == "__main__":
    main()
