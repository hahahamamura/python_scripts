#!/usr/bin/env python3
"""
Script para filtrar VCFs usando arquivo BED e fazer merge por cromossomo
Requisitos: bcftools instalado no sistema
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path


def check_bcftools():
    """Verifica se bcftools está instalado"""
    try:
        subprocess.run(['bcftools', '--version'], 
                      capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("ERRO: bcftools não encontrado. Por favor, instale o bcftools.")
        print("Ubuntu/Debian: sudo apt-get install bcftools")
        print("conda: conda install -c bioconda bcftools")
        return False


def filter_vcf(vcf_file, bed_file, output_file):
    """
    Filtra um arquivo VCF usando um arquivo BED
    
    Args:
        vcf_file: Caminho para o arquivo VCF de entrada
        bed_file: Caminho para o arquivo BED com regiões
        output_file: Caminho para o arquivo VCF filtrado
    """
    cmd = [
        'bcftools', 'view',
        '-R', bed_file,  # Filtrar por regiões do BED
        '-O', 'z',       # Output comprimido (vcf.gz)
        '-o', output_file,
        vcf_file
    ]
    
    print(f"  Filtrando {os.path.basename(vcf_file)}...")
    subprocess.run(cmd, check=True)
    
    # Indexar o arquivo filtrado
    subprocess.run(['bcftools', 'index', '-t', output_file], check=True)


def merge_vcfs(vcf1, vcf2, output_file):
    """
    Faz merge de dois arquivos VCF
    
    Args:
        vcf1: Primeiro arquivo VCF
        vcf2: Segundo arquivo VCF
        output_file: Arquivo VCF de saída
    """
    cmd = [
        'bcftools', 'merge',
        '-O', 'z',       # Output comprimido
        '-o', output_file,
        vcf1, vcf2
    ]
    
    print(f"  Fazendo merge para {os.path.basename(output_file)}...")
    subprocess.run(cmd, check=True)
    
    # Indexar o arquivo final
    subprocess.run(['bcftools', 'index', '-t', output_file], check=True)


def process_chromosome(chr_num, dir1, dir2, bed_file, output_dir, temp_dir, 
                       pattern1="chr{}.vcf.gz", pattern2="chr{}.vcf.gz"):
    """
    Processa um cromossomo: filtra ambos VCFs e faz merge
    
    Args:
        chr_num: Número do cromossomo
        dir1: Diretório 1 com VCFs
        dir2: Diretório 2 com VCFs
        bed_file: Arquivo BED para filtro
        output_dir: Diretório de saída
        temp_dir: Diretório temporário
        pattern1: Padrão do nome dos arquivos VCF do diretório 1
        pattern2: Padrão do nome dos arquivos VCF do diretório 2
    """
    print(f"\nProcessando cromossomo {chr_num}...")
    
    # Caminhos dos arquivos
    vcf1_input = os.path.join(dir1, pattern1.format(chr_num))
    vcf2_input = os.path.join(dir2, pattern2.format(chr_num))
    
    # Verificar se os arquivos existem
    if not os.path.exists(vcf1_input):
        print(f"  AVISO: {vcf1_input} não encontrado. Pulando...")
        return False
    if not os.path.exists(vcf2_input):
        print(f"  AVISO: {vcf2_input} não encontrado. Pulando...")
        return False
    
    # Arquivos temporários filtrados
    vcf1_filtered = os.path.join(temp_dir, f"chr{chr_num}_dir1_filtered.vcf.gz")
    vcf2_filtered = os.path.join(temp_dir, f"chr{chr_num}_dir2_filtered.vcf.gz")
    
    # Arquivo de saída final
    output_file = os.path.join(output_dir, f"chr{chr_num}_merged.vcf.gz")
    
    try:
        # Filtrar ambos os VCFs
        filter_vcf(vcf1_input, bed_file, vcf1_filtered)
        filter_vcf(vcf2_input, bed_file, vcf2_filtered)
        
        # Fazer merge
        merge_vcfs(vcf1_filtered, vcf2_filtered, output_file)
        
        print(f"  ✓ Cromossomo {chr_num} concluído!")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"  ERRO ao processar cromossomo {chr_num}: {e}")
        return False
    
    finally:
        # Limpar arquivos temporários
        for temp_file in [vcf1_filtered, vcf2_filtered, 
                         f"{vcf1_filtered}.tbi", f"{vcf2_filtered}.tbi"]:
            if os.path.exists(temp_file):
                os.remove(temp_file)


def main():
    parser = argparse.ArgumentParser(
        description='Filtra e faz merge de VCFs por cromossomo usando arquivo BED'
    )
    parser.add_argument('--dir1', required=True, 
                       help='Diretório 1 com VCFs')
    parser.add_argument('--dir2', required=True, 
                       help='Diretório 2 com VCFs')
    parser.add_argument('--bed', required=True, 
                       help='Arquivo BED com regiões para filtro')
    parser.add_argument('--output', required=True, 
                       help='Diretório de saída (será criado se não existir)')
    parser.add_argument('--temp', default='./temp_vcf', 
                       help='Diretório temporário (padrão: ./temp_vcf)')
    parser.add_argument('--pattern1', 
                       default='1kGP_high_coverage_Illumina.chr{}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz',
                       help='Padrão do nome dos VCFs do diretório 1')
    parser.add_argument('--pattern2', 
                       default='chr{}_phased_multial.vcf.gz',
                       help='Padrão do nome dos VCFs do diretório 2')
    parser.add_argument('--chromosomes', default='1-22',
                       help='Cromossomos a processar (padrão: 1-22)')
    
    args = parser.parse_args()
    
    # Verificar bcftools
    if not check_bcftools():
        sys.exit(1)
    
    # Verificar se diretórios e arquivo BED existem
    if not os.path.exists(args.dir1):
        print(f"ERRO: Diretório 1 não encontrado: {args.dir1}")
        sys.exit(1)
    if not os.path.exists(args.dir2):
        print(f"ERRO: Diretório 2 não encontrado: {args.dir2}")
        sys.exit(1)
    if not os.path.exists(args.bed):
        print(f"ERRO: Arquivo BED não encontrado: {args.bed}")
        sys.exit(1)
    
    # Criar diretórios de saída e temporário
    os.makedirs(args.output, exist_ok=True)
    os.makedirs(args.temp, exist_ok=True)
    
    # Processar cromossomos
    chromosomes = range(1, 23)  # 1 a 22
    
    print("="*60)
    print("Iniciando processamento de VCFs")
    print("="*60)
    print(f"Diretório 1: {args.dir1}")
    print(f"Diretório 2: {args.dir2}")
    print(f"Arquivo BED: {args.bed}")
    print(f"Saída: {args.output}")
    print(f"Cromossomos: {list(chromosomes)}")
    
    success_count = 0
    fail_count = 0
    
    for chr_num in chromosomes:
        if process_chromosome(chr_num, args.dir1, args.dir2, args.bed, 
                             args.output, args.temp, args.pattern1, args.pattern2):
            success_count += 1
        else:
            fail_count += 1
    
    print("\n" + "="*60)
    print("RESUMO")
    print("="*60)
    print(f"Cromossomos processados com sucesso: {success_count}")
    print(f"Cromossomos com erro/pulados: {fail_count}")
    print(f"Arquivos de saída em: {args.output}")
    print("="*60)


if __name__ == "__main__":
    main()