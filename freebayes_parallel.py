#!/usr/bin/env python3
"""
Script simplificado para executar FreeBayes em todos os BAMs de um diretório,
usando regiões específicas de um arquivo BED, e concatenar os VCFs com bcftools.
"""

import os
import sys
import subprocess
import glob
from pathlib import Path

# ============= CONFIGURAÇÕES - EDITE AQUI =============
REF_FASTA = "/home/lab/Desktop/arq_joao/testes_freebayes_bams_marcel/REF/GRCh38_full_analysis_set_plus_decoy_hla.fa"
BAM_DIRECTORY = "/home/lab/Desktop/arq_joao/testes_freebayes_bams_marcel/BAMs_INDEX"
BED_FILE = "/home/lab/Desktop/arq_joao/testes_freebayes_bams_marcel/BED/livia.bed"
OUTPUT_DIR = "/home/lab/Desktop/arq_joao/testes_freebayes_bams_marcel/VCF_MERGED"
KEEP_INTERMEDIATE = False  # True para manter arquivos VCF intermediários
# =====================================================

def parse_bed_file(bed_file):
    """
    Parse do arquivo BED para extrair regiões no formato chr:start-end
    """
    regions = []
    try:
        with open(bed_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split('\t')
                if len(fields) < 3:
                    print(f"Aviso: Linha {line_num} do BED não tem pelo menos 3 campos: {line}")
                    continue
                
                chrom = fields[0]
                start = fields[1]
                end = fields[2]
                
                # Formato para freebayes: chr:start-end
                region = f"{chrom}:{start}-{end}"
                regions.append(region)
    
    except FileNotFoundError:
        print(f"Erro: Arquivo BED '{bed_file}' não encontrado.")
        sys.exit(1)
    except Exception as e:
        print(f"Erro ao ler arquivo BED: {e}")
        sys.exit(1)
    
    return regions

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

def run_freebayes(ref_fasta, bam_file, region, output_vcf):
    """
    Executa o FreeBayes para um BAM específico e região
    """
    cmd = [
        'freebayes',
        '-f', ref_fasta,
        '-r', region,
        bam_file
    ]
    
    print(f"Executando: {' '.join(cmd)} > {output_vcf}")
    
    try:
        with open(output_vcf, 'w') as outfile:
            result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, 
                                  text=True, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Erro ao executar FreeBayes: {e}")
        print(f"Stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print("Erro: FreeBayes não encontrado. Certifique-se de que está instalado e no PATH.")
        return False

def concatenate_with_bcftools(vcf_files, output_vcf):
    """
    Concatena múltiplos arquivos VCF usando bcftools concat
    """
    print(f"Concatenando {len(vcf_files)} arquivos VCF usando bcftools concat...")
    
    # Filtrar apenas arquivos que existem e não estão vazios
    valid_vcfs = []
    for vcf in vcf_files:
        if os.path.exists(vcf) and os.path.getsize(vcf) > 0:
            valid_vcfs.append(vcf)
    
    if not valid_vcfs:
        print("Nenhum arquivo VCF válido para concatenar")
        return False
    
    cmd = ['bcftools', 'concat'] + valid_vcfs + ['-O', 'v', '-o', output_vcf]
    
    try:
        print(f"Executando: {' '.join(cmd)}")
        result = subprocess.run(cmd, stderr=subprocess.PIPE, text=True, check=True)
        print(f"Concatenação concluída: {output_vcf}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Erro ao executar bcftools concat: {e}")
        print(f"Stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print("Erro: bcftools não encontrado. Certifique-se de que está instalado e no PATH.")
        return False

def main():
    # Verificar se arquivos e diretórios existem
    if not os.path.exists(REF_FASTA):
        print(f"Erro: Arquivo de referência não encontrado: {REF_FASTA}")
        sys.exit(1)
    
    if not os.path.exists(BAM_DIRECTORY):
        print(f"Erro: Diretório de BAMs não encontrado: {BAM_DIRECTORY}")
        sys.exit(1)
    
    if not os.path.exists(BED_FILE):
        print(f"Erro: Arquivo BED não encontrado: {BED_FILE}")
        sys.exit(1)
    
    # Criar diretório de saída
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(exist_ok=True)
    
    # Parse do arquivo BED
    print("Lendo arquivo BED...")
    regions = parse_bed_file(BED_FILE)
    print(f"Encontradas {len(regions)} regiões no arquivo BED")
    
    # Encontrar arquivos BAM
    print("Procurando arquivos BAM...")
    bam_files = find_bam_files(BAM_DIRECTORY)
    print(f"Encontrados {len(bam_files)} arquivos BAM")
    
    # Dicionário para armazenar VCFs por amostra
    sample_vcf_files = {}
    
    # Executar FreeBayes para cada combinação BAM x região
    total_jobs = len(bam_files) * len(regions)
    current_job = 0
    
    for bam_file in bam_files:
        bam_name = Path(bam_file).stem
        sample_vcf_files[bam_name] = []
        
        for i, region in enumerate(regions):
            current_job += 1
            print(f"\nProcessando job {current_job}/{total_jobs}")
            print(f"BAM: {bam_name}, Região: {region}")
            
            # Nome do arquivo VCF de saída
            region_safe = region.replace(':', '_').replace('-', '_')
            output_vcf = output_dir / f"{bam_name}_{region_safe}.vcf"
            
            # Executar FreeBayes
            success = run_freebayes(REF_FASTA, bam_file, region, str(output_vcf))
            
            if success and output_vcf.exists() and output_vcf.stat().st_size > 0:
                sample_vcf_files[bam_name].append(str(output_vcf))
            else:
                print(f"Aviso: Falha ou arquivo vazio para {output_vcf}")
    
    # Concatenar VCFs para cada amostra usando bcftools
    final_vcfs = []
    for bam_name, vcf_files in sample_vcf_files.items():
        if vcf_files:
            final_vcf = output_dir / f"{bam_name}_merged.vcf"
            print(f"\nConcatenando {len(vcf_files)} regiões para amostra {bam_name}...")
            success = concatenate_with_bcftools(vcf_files, str(final_vcf))
            
            if success:
                final_vcfs.append(str(final_vcf))
                print(f"VCF criado para {bam_name}: {final_vcf}")
            else:
                print(f"Erro na concatenação dos VCFs para {bam_name}")
        else:
            print(f"Nenhum VCF válido gerado para amostra {bam_name}")
    
    if final_vcfs:
        print(f"\nProcesso concluído!")
        print(f"VCFs finais criados: {len(final_vcfs)}")
        for vcf in final_vcfs:
            print(f"  - {vcf}")
        
        if not KEEP_INTERMEDIATE:
            print("\nRemoção de arquivos intermediários...")
            for sample_vcfs in sample_vcf_files.values():
                for vcf_file in sample_vcfs:
                    try:
                        os.remove(vcf_file)
                    except OSError:
                        pass
            print("Arquivos intermediários removidos")
        else:
            print(f"Arquivos VCF intermediários mantidos em: {output_dir}")
    else:
        print("Nenhum arquivo VCF final foi gerado com sucesso")

if __name__ == "__main__":
    main()