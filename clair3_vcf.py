#!/usr/bin/env python3
"""
Script para executar Clair3 em todos os BAMs de um diretório.
Cada BAM será processado individualmente em regiões do BED,
gerando VCF e GVCF finais em uma pasta separada.
"""

import os
import sys
import subprocess
import glob
from pathlib import Path

# ================= CONFIGURAÇÕES =================
REF_FASTA = "/home/lab/Desktop/arq_joao/testes_freebayes_bams_marcel/REF/GRCh38_full_analysis_set_plus_decoy_hla.fa"
BAM_DIRECTORY = "/home/lab/Desktop/arq_joao/testes_freebayes_bams_marcel/BAMs_INDEX"
BED_FILE = "/home/lab/Desktop/arq_joao/testes_freebayes_bams_marcel/BED/livia.bed"
OUTPUT_DIR = "/home/lab/Desktop/arq_joao/testes_freebayes_bams_marcel/VCF_MERGED_CLAIR"
PLATFORM = "ont"   # opções: illumina, hifi, ont
THREADS = 4
MODEL = 'ont'#"r1041_e82_400bps_sup_v500"
# =================================================

def find_bam_files(bam_directory):
    """
    Encontra todos os arquivos BAM em um diretório
    """
    bam_pattern = os.path.join(bam_directory, "*.bam")
    bam_files = glob.glob(bam_pattern)
    
    if not bam_files:
        print(f"Nenhum arquivo BAM encontrado em: {bam_directory}")
        sys.exit(1)
    
    return sorted(bam_files)

def run_clair3(bam_file, ref_fasta, model, bed_file, output_dir, platform, threads):
    """
    Executa Clair3 para um BAM específico
    """
    bam_name = Path(bam_file).stem
    sample_outdir = Path(output_dir) / bam_name
    sample_outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "run_clair3.sh",
        f"--bam_fn={bam_file}",
        f"--ref_fn={ref_fasta}",
        f"--model_path=/home/lab/miniconda3/envs/clair3/bin/models/{model}",
        f"--threads={threads}",
        f"--platform={platform}",
        f"--output={sample_outdir}",
        f"--bed_fn={bed_file}",
        "--gvcf"
    ]
    
    print(f"\n[INFO] Executando Clair3 para {bam_name}")
    print(" ".join(cmd))
    
    try:
        subprocess.run(cmd, check=True)
        final_vcf = sample_outdir / "merge_output.vcf.gz"
        gvcf = sample_outdir / "merge_output.g.vcf.gz"

        if final_vcf.exists():
            print(f"[OK] VCF final gerado: {final_vcf}")
        else:
            print(f"[ERRO] VCF esperado não encontrado em {final_vcf}")

        if gvcf.exists():
            print(f"[OK] GVCF final gerado: {gvcf}")
        else:
            print(f"[ERRO] GVCF esperado não encontrado em {gvcf}")

        return str(final_vcf) if final_vcf.exists() else None

    except subprocess.CalledProcessError as e:
        print(f"[ERRO] Falha ao executar Clair3 para {bam_name}: {e}")
        return None
    except FileNotFoundError:
        print("[ERRO] run_clair3.sh não encontrado. Ative o ambiente conda do Clair3!")
        return None


def main():
    # Checar arquivos e diretórios
    if not os.path.exists(REF_FASTA):
        print(f"Erro: Arquivo de referência não encontrado: {REF_FASTA}")
        sys.exit(1)

    if not os.path.exists(BAM_DIRECTORY):
        print(f"Erro: Diretório de BAMs não encontrado: {BAM_DIRECTORY}")
        sys.exit(1)

    if not os.path.exists(BED_FILE):
        print(f"Erro: Arquivo BED não encontrado: {BED_FILE}")
        sys.exit(1)

    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    # Procurar BAMs
    print("Procurando arquivos BAM...")
    bam_files = find_bam_files(BAM_DIRECTORY)
    print(f"Encontrados {len(bam_files)} BAMs")

    # Executar Clair3 em cada BAM
    vcfs = []
    for bam_file in bam_files:
        vcf = run_clair3(bam_file, REF_FASTA, MODEL, BED_FILE, OUTPUT_DIR, PLATFORM, THREADS)
        if vcf:
            vcfs.append(vcf)

    # Resumo final
    print("\n===== RESUMO =====")
    if vcfs:
        print(f"Foram gerados {len(vcfs)} VCFs finais:")
        for vcf in vcfs:
            print(f" - {vcf}")
    else:
        print("Nenhum VCF foi gerado com sucesso.")

if __name__ == "__main__":
    main()
