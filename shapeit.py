#!/usr/bin/env python3
import os
import subprocess
import argparse

def run_cmd(cmd):
    print(f"\n>>> Rodando: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description="Pipeline SHAPEIT4 cromossomo por cromossomo com organização de pastas")
    parser.add_argument("--vcf", required=True, help="VCF de entrada (bgz + indexado)")
    parser.add_argument("--out_dir", required=True, help="Diretório de saída")
    parser.add_argument("--threads", default="4", help="Número de threads (default=4)")
    args = parser.parse_args()

    # Criar diretórios
    os.makedirs(args.out_dir, exist_ok=True)
    bial_dir = os.path.join(args.out_dir, "biallelic")
    multi_dir = os.path.join(args.out_dir, "multiallelic")
    log_dir = os.path.join(args.out_dir, "logs")
    os.makedirs(bial_dir, exist_ok=True)
    os.makedirs(multi_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    # === Passo 1: Normalizar para bialélico ===
    vcf_bial = os.path.join(bial_dir, "input_bial.vcf.gz")
    run_cmd(f"bcftools norm -m-any {args.vcf} -Oz -o {vcf_bial}")
    run_cmd(f"tabix -p vcf {vcf_bial}")

    # === Passo 2: Loop pelos cromossomos 1 a 22 ===
    for chr_num in range(1, 23):
        chrom_name = f"chr{chr_num}"  # Corrige o nome do cromossomo
        print(f"\n=== Faseando {chrom_name} ===")

        phased_vcf = os.path.join(bial_dir, f"{chrom_name}_phased.vcf.gz")  # salvar temporário na pasta biallelic
        phased_log = os.path.join(log_dir, f"{chrom_name}_phased.log")

        # Rodar SHAPEIT4 com parâmetros fixos
        run_cmd(
            f"shapeit4 "
            f"--input {vcf_bial} "
            f"--region {chrom_name} "
            f"--output {phased_vcf} "
            f"--thread {args.threads} "
            f"--log {phased_log}"
        )

        # === Passo 3: Converter de volta para multialélico ===
        phased_multi = os.path.join(multi_dir, f"{chrom_name}_phased_multial.vcf.gz")
        run_cmd(f"bcftools norm -m+any {phased_vcf} -Oz -o {phased_multi}")
        run_cmd(f"tabix -p vcf {phased_multi}")

    print("\n✅ Pipeline finalizado com sucesso!")
    print(f"Bialélicos -> {bial_dir}")
    print(f"Multialélicos -> {multi_dir}")
    print(f"Logs -> {log_dir}")

if __name__ == "__main__":
    main()
