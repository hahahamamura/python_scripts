#!/usr/bin/env python3
import os
import re
import subprocess

# Caminhos
vcf_dir = "/home/lab/Desktop/arq_joao/VCF_WHOLEGENOME"         # pasta onde estão os VCFs originais
bed_file = "/home/lab/Downloads/snps_positions.bed"            # arquivo BED com as regiões
out_dir = "/home/lab/Desktop/arq_joao/VCF_WHOLEGENOME_filtred" # pasta de saída
final_vcf = "/home/lab/Desktop/arq_joao/VCF_WHOLEGENOME_filtred/vcf_final_merged.vcf.gz"  # arquivo final concatenado

# Criar pasta de saída se não existir
os.makedirs(out_dir, exist_ok=True)

# Função para extrair número do cromossomo do nome do arquivo
def get_chr_number(filename):
    match = re.search(r"chr(\d+)", filename)
    if match:
        return int(match.group(1))
    else:
        return 9999  # se não encontrar, manda pro fim

# Listar todos os arquivos VCF da pasta (aceita .vcf e .vcf.gz)
vcf_files = [f for f in os.listdir(vcf_dir) if f.endswith(".vcf") or f.endswith(".vcf.gz")]

# Ordenar por número do cromossomo
vcf_files_sorted = sorted(vcf_files, key=get_chr_number)

filtrados = []  # lista de VCFs filtrados para concatenar

for vcf in vcf_files_sorted:
    input_path = os.path.join(vcf_dir, vcf)
    base_name = vcf.replace(".vcf.gz", "").replace(".vcf", "")
    output_name = base_name + ".filtered.vcf.gz"
    output_path = os.path.join(out_dir, output_name)
    
    # Comando bcftools view
    cmd = [
        "bcftools", "view",
        "-R", bed_file,
        "-Oz",
        "-o", output_path,
        input_path
    ]
    
    print(f"Filtrando {vcf} -> {output_name}")
    subprocess.run(cmd, check=True)
    
    # Indexar VCF filtrado
    subprocess.run(["bcftools", "index", output_path], check=True)

    filtrados.append(output_path)

# Concatenar todos os VCFs filtrados em um único
if filtrados:
    print("🔗 Concatenando todos os VCFs filtrados em", final_vcf)
    cmd_concat = ["bcftools", "concat", "-Oz", "-o", final_vcf] + filtrados
    subprocess.run(cmd_concat, check=True)
    subprocess.run(["bcftools", "index", final_vcf], check=True)

print("✅ Processo concluído! VCF final gerado em:", final_vcf)
