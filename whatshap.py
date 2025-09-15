#!/usr/bin/env python3
"""
Script para processar VCF com WhatsHap separadamente por amostra
Autor: Script automatizado para WhatsHap
"""

import os
import subprocess
import argparse
import sys
from pathlib import Path
import tempfile
import shutil
import glob

def run_command(cmd, description=""):
    """Executa comando e verifica se foi bem-sucedido"""
    print(f"Executando: {description}")
    print(f"Comando: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"✓ Sucesso: {description}")
        return result
    except subprocess.CalledProcessError as e:
        print(f"✗ Erro em: {description}")
        print(f"Comando falhou: {' '.join(cmd)}")
        print(f"Stderr: {e.stderr}")
        sys.exit(1)

def get_samples_from_vcf(vcf_file):
    """Extrai nomes das amostras do arquivo VCF"""
    print(f"Extraindo amostras de: {vcf_file}")
    
    cmd = ["bcftools", "query", "-l", vcf_file]
    result = run_command(cmd, "Extraindo lista de amostras")
    
    samples = result.stdout.strip().split('\n')
    samples = [s for s in samples if s]  # Remove linhas vazias
    
    print(f"Encontradas {len(samples)} amostras: {', '.join(samples)}")
    return samples

def find_bam_for_sample(sample, bam_dir):
    """Encontra o arquivo BAM correspondente à amostra no diretório"""
    # Padrões possíveis para nomes de arquivos BAM
    patterns = [
        f"{sample}.bam",
        f"{sample}*.bam",
        f"*{sample}.bam",
        f"*{sample}*.bam"
    ]
    
    for pattern in patterns:
        bam_path = os.path.join(bam_dir, pattern)
        matches = glob.glob(bam_path)
        if matches:
            # Retorna o primeiro match encontrado
            bam_file = matches[0]
            print(f"BAM encontrado para amostra {sample}: {bam_file}")
            return bam_file
    
    # Se não encontrou, tenta buscar recursivamente
    for root, dirs, files in os.walk(bam_dir):
        for file in files:
            if file.endswith('.bam'):
                # Verifica se o nome da amostra está no nome do arquivo
                if sample in file:
                    bam_file = os.path.join(root, file)
                    print(f"BAM encontrado para amostra {sample}: {bam_file}")
                    return bam_file
    
    print(f"⚠️  Aviso: BAM não encontrado para amostra {sample} em {bam_dir}")
    return None

def split_vcf_by_sample(vcf_file, sample, output_dir):
    """Separa VCF por amostra usando bcftools"""
    output_file = os.path.join(output_dir, f"{sample}.vcf.gz")
    
    cmd = [
        "bcftools", "view",
        "-s", sample,
        "-O", "z",
        "-o", output_file,
        vcf_file
    ]
    
    run_command(cmd, f"Separando VCF para amostra {sample}")
    
    # Indexar o arquivo VCF
    cmd_index = ["bcftools", "index", "-t", output_file]
    run_command(cmd_index, f"Indexando VCF da amostra {sample}")
    
    return output_file

def run_whatshap_single_sample(vcf_file, bam_file, sample, output_dir, reference=None):
    """Executa WhatsHap para uma única amostra"""
    if not bam_file:
        print(f"✗ Pulando amostra {sample}: BAM não encontrado")
        return None
        
    output_file = os.path.join(output_dir, f"{sample}_phased.vcf.gz")
    
    cmd = [
        "whatshap", "phase",
        "--output", output_file,
        "--sample", sample
    ]
    
    if reference:
        cmd.extend(["--reference", reference])
    
    cmd.extend([vcf_file, bam_file])
    
    run_command(cmd, f"Executando WhatsHap para amostra {sample}")
    
    # Indexar o arquivo de saída
    cmd_index = ["bcftools", "index", "-t", output_file]
    run_command(cmd_index, f"Indexando VCF fasado da amostra {sample}")
    
    return output_file

def merge_phased_vcfs(phased_vcf_files, output_file):
    """Junta os arquivos VCF fasados usando bcftools merge"""
    if not phased_vcf_files:
        print("✗ Erro: Nenhum arquivo VCF fasado para juntar")
        sys.exit(1)
        
    print("Juntando arquivos VCF fasados...")
    
    cmd = [
        "bcftools", "merge",
        "-O", "z",
        "-o", output_file
    ] + phased_vcf_files
    
    run_command(cmd, "Juntando arquivos VCF fasados")
    
    # Indexar o arquivo final
    cmd_index = ["bcftools", "index", "-t", output_file]
    run_command(cmd_index, "Indexando VCF final")

def main():
    parser = argparse.ArgumentParser(
        description="Processa VCF com WhatsHap separadamente por amostra",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplo de uso:
    python whatshap_script.py -v input.vcf.gz -b /path/to/bam/directory -o output_phased.vcf.gz
    python whatshap_script.py -v input.vcf.gz -b ./bams/ -o output.vcf.gz -r reference.fasta
    python whatshap_script.py -v input.vcf.gz -b /data/bams --temp-dir /tmp/whatshap -o output.vcf.gz
        """
    )
    
    parser.add_argument("-v", "--vcf", required=True,
                      help="Arquivo VCF de entrada")
    parser.add_argument("-b", "--bam-dir", required=True,
                      help="Diretório contendo os arquivos BAM")
    parser.add_argument("-o", "--output", required=True,
                      help="Arquivo VCF de saída (fasado)")
    parser.add_argument("-r", "--reference",
                      help="Arquivo FASTA de referência (opcional)")
    parser.add_argument("--temp-dir",
                      help="Diretório temporário para arquivos intermediários")
    parser.add_argument("--keep-temp", action="store_true",
                      help="Manter arquivos temporários após processamento")
    
    args = parser.parse_args()
    
    # Verificar se arquivo VCF existe
    if not os.path.exists(args.vcf):
        print(f"Erro: Arquivo VCF não encontrado: {args.vcf}")
        sys.exit(1)
    
    # Verificar se diretório BAM existe
    if not os.path.isdir(args.bam_dir):
        print(f"Erro: Diretório BAM não encontrado: {args.bam_dir}")
        sys.exit(1)
    
    if args.reference and not os.path.exists(args.reference):
        print(f"Erro: Arquivo de referência não encontrado: {args.reference}")
        sys.exit(1)
    
    # Verificar se ferramentas estão disponíveis
    for tool in ["bcftools", "whatshap"]:
        try:
            subprocess.run([tool, "--version"], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print(f"Erro: {tool} não está instalado ou não está no PATH")
            sys.exit(1)
    
    # Criar diretório temporário
    if args.temp_dir:
        temp_dir = args.temp_dir
        os.makedirs(temp_dir, exist_ok=True)
    else:
        temp_dir = tempfile.mkdtemp(prefix="whatshap_")
    
    print(f"Usando diretório temporário: {temp_dir}")
    
    try:
        # 1. Extrair amostras do VCF
        samples = get_samples_from_vcf(args.vcf)
        
        if not samples:
            print("Erro: Nenhuma amostra encontrada no arquivo VCF")
            sys.exit(1)
        
        # 2. Encontrar arquivos BAM para cada amostra
        print("\n=== Mapeando BAMs para amostras ===")
        sample_bams = {}
        for sample in samples:
            bam_file = find_bam_for_sample(sample, args.bam_dir)
            sample_bams[sample] = bam_file
        
        # Verificar se pelo menos um BAM foi encontrado
        found_bams = [bam for bam in sample_bams.values() if bam is not None]
        if not found_bams:
            print("✗ Erro: Nenhum arquivo BAM encontrado para as amostras")
            sys.exit(1)
        
        print(f"Encontrados {len(found_bams)} de {len(samples)} arquivos BAM")
        
        # 3. Separar VCF por amostra
        print("\n=== Separando VCF por amostra ===")
        sample_vcfs = {}
        for sample in samples:
            sample_vcf = split_vcf_by_sample(args.vcf, sample, temp_dir)
            sample_vcfs[sample] = sample_vcf
        
        # 4. Executar WhatsHap para cada amostra
        print("\n=== Executando WhatsHap por amostra ===")
        phased_vcfs = []
        for sample in samples:
            print(f"\nProcessando amostra: {sample}")
            if sample_bams[sample]:
                phased_vcf = run_whatshap_single_sample(
                    sample_vcfs[sample], 
                    sample_bams[sample],
                    sample, 
                    temp_dir,
                    args.reference
                )
                if phased_vcf:
                    phased_vcfs.append(phased_vcf)
            else:
                print(f"✗ Pulando amostra {sample}: BAM não encontrado")
        
        # 5. Juntar arquivos VCF fasados
        print("\n=== Juntando resultados ===")
        merge_phased_vcfs(phased_vcfs, args.output)
        
        print(f"\n✓ Processamento concluído!")
        print(f"Arquivo de saída: {args.output}")
        print(f"Amostras processadas: {len(phased_vcfs)} de {len(samples)}")
        print(f"Arquivos temporários em: {temp_dir}")
        
    except KeyboardInterrupt:
        print("\n✗ Processamento interrompido pelo usuário")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ Erro durante processamento: {e}")
        sys.exit(1)
    finally:
        # Limpar arquivos temporários se solicitado
        if not args.keep_temp and not args.temp_dir:
            print(f"\nRemovendo arquivos temporários...")
            shutil.rmtree(temp_dir, ignore_errors=True)

if __name__ == "__main__":
    main()