#!/usr/bin/env python3
"""
Script para executar FreeBayes em paralelo por regiões do arquivo BED
e concatenar os VCFs resultantes usando bcftools - VERSÃO CORRIGIDA PARA DUPLICATAS.
"""

import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import tempfile
import shutil

# =============================================================================
# CONFIGURAÇÕES - DEFINA OS CAMINHOS AQUI
# =============================================================================

# Caminhos dos arquivos necessários
BAM_DIRECTORY = "/home/lab/Desktop/arq_joao/NativoAmericanas/teste_bam"           # Arquivo com lista de BAMs
REFERENCE_GENOME = "/home/lab/Desktop/arq_joao/NativoAmericanas/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"      # Genoma de referência
BED_FILE = "/home/lab/Desktop/arq_joao/NativoAmericanas/livia.bed"                   # Arquivo BED com regiões
OUTPUT_DIR = "/home/lab/Desktop/arq_joao/NativoAmericanas/teste_vcf"                      # Diretório de saída
FINAL_VCF = "resultado_final.vcf.gz"                     # Nome do VCF final

# Configurações de paralelização
MAX_PARALLEL_JOBS = 4                                    # Número máximo de jobs simultâneos

# Parâmetros adicionais do FreeBayes (opcional)
FREEBAYES_EXTRA_PARAMS = ""                             # Ex: "--min-mapping-quality 30"

# =============================================================================

def validate_files():
    """Valida se todos os arquivos necessários existem e cria bamlist temporário."""
    files_to_check = [REFERENCE_GENOME, BED_FILE]
    
    for file_path in files_to_check:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Arquivo não encontrado: {file_path}")
    
    # Verifica se o diretório de BAMs existe
    if not os.path.exists(BAM_DIRECTORY):
        raise FileNotFoundError(f"Diretório de BAMs não encontrado: {BAM_DIRECTORY}")
    
    # Cria diretório de saída se não existir
    os.makedirs(OUTPUT_DIR, exist_ok=True)

def validate_bam_files(bam_files):
    """Valida se os arquivos BAM existem e têm índices."""
    valid_bams = []
    
    for bam_file in bam_files:
        # Verifica se o arquivo BAM existe
        if not os.path.exists(bam_file):
            print(f"Aviso: Arquivo BAM não existe: {bam_file}")
            continue
        
        # Verifica se há índice (.bai ou .bam.bai)
        index_files = [bam_file + ".bai", bam_file.replace(".bam", ".bai")]
        has_index = any(os.path.exists(idx) for idx in index_files)
        
        if not has_index:
            print(f"Aviso: Arquivo BAM sem índice: {bam_file}")
            # Tenta criar índice
            try:
                print(f"Tentando criar índice para {os.path.basename(bam_file)}...")
                subprocess.run(["samtools", "index", bam_file], check=True, 
                             capture_output=True, text=True)
                print(f"✓ Índice criado para {os.path.basename(bam_file)}")
            except subprocess.CalledProcessError as e:
                print(f"✗ Erro ao criar índice para {bam_file}: {e.stderr}")
                continue
        
        # Verifica se o BAM tem reads
        try:
            result = subprocess.run(
                ["samtools", "view", "-c", bam_file],
                capture_output=True, text=True, check=True
            )
            read_count = int(result.stdout.strip())
            if read_count == 0:
                print(f"Aviso: Arquivo BAM vazio: {bam_file}")
            else:
                print(f"✓ BAM válido: {os.path.basename(bam_file)} ({read_count:,} reads)")
        except Exception as e:
            print(f"Aviso: Não foi possível contar reads em {bam_file}: {e}")
        
        valid_bams.append(bam_file)
    
    return valid_bams

def create_bamlist():
    """Cria um arquivo temporário com a lista de BAMs do diretório."""
    bam_files = []
    
    # Procura por arquivos BAM no diretório
    bam_extensions = ['.bam']
    
    for file_name in os.listdir(BAM_DIRECTORY):
        if any(file_name.lower().endswith(ext) for ext in bam_extensions):
            full_path = os.path.join(BAM_DIRECTORY, file_name)
            bam_files.append(full_path)
    
    if not bam_files:
        raise Exception(f"Nenhum arquivo BAM encontrado em: {BAM_DIRECTORY}")
    
    # Ordena os arquivos BAM para ordem consistente
    bam_files.sort()
    
    # Valida arquivos BAM
    print("Validando arquivos BAM...")
    valid_bams = validate_bam_files(bam_files)
    
    if not valid_bams:
        raise Exception("Nenhum arquivo BAM válido encontrado!")
    
    # Cria arquivo temporário com lista de BAMs
    bamlist_file = os.path.join(OUTPUT_DIR, "bamlist_temp.txt")
    
    with open(bamlist_file, 'w') as f:
        for bam_file in valid_bams:
            f.write(f"{bam_file}\n")
    
    print(f"✓ {len(valid_bams)} arquivos BAM válidos encontrados")
    
    return bamlist_file

def check_overlapping_regions(regions):
    """Verifica e reporta regiões sobrepostas que podem causar duplicatas."""
    overlaps = []
    
    # Converte regiões para formato comparável
    parsed_regions = []
    for region_str, region_id in regions:
        # Format: chr:start-end
        chrom, coords = region_str.split(':')
        start, end = map(int, coords.split('-'))
        parsed_regions.append((chrom, start, end, region_str, region_id))
    
    # Verifica sobreposições
    for i, (chrom1, start1, end1, region1, id1) in enumerate(parsed_regions):
        for j, (chrom2, start2, end2, region2, id2) in enumerate(parsed_regions[i+1:], i+1):
            if chrom1 == chrom2:
                # Verifica se há sobreposição
                if not (end1 <= start2 or end2 <= start1):  # Há sobreposição
                    overlap_start = max(start1, start2)
                    overlap_end = min(end1, end2)
                    overlaps.append((region1, region2, f"{chrom1}:{overlap_start}-{overlap_end}"))
    
    if overlaps:
        print(f"\n⚠️  AVISO: {len(overlaps)} sobreposições detectadas entre regiões:")
        for reg1, reg2, overlap in overlaps:
            print(f"  {reg1} ↔ {reg2} (sobreposição: {overlap})")
        print("  Isso pode causar variantes duplicadas no resultado final!")
        
        # Pergunta se deve continuar
        response = input("\nContinuar mesmo assim? (y/N): ").lower()
        if response not in ['y', 'yes', 's', 'sim']:
            print("Operação cancelada pelo usuário.")
            sys.exit(1)
    
    return overlaps

def read_bed_regions(bed_file):
    """Lê o arquivo BED e retorna uma lista de regiões ordenadas."""
    regions = []
    
    with open(bed_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 3:
                print(f"Aviso: Linha {line_num} do BED inválida (menos de 3 colunas): {line}")
                continue
            
            chrom = parts[0]
            start = parts[1]
            end = parts[2]
            
            # Valida coordenadas
            try:
                start_pos = int(start)
                end_pos = int(end)
                if start_pos >= end_pos:
                    print(f"Aviso: Coordenadas inválidas na linha {line_num}: start >= end")
                    continue
            except ValueError:
                print(f"Aviso: Coordenadas não numéricas na linha {line_num}: {line}")
                continue
            
            # Formato para FreeBayes: chr:start-end
            region = f"{chrom}:{start}-{end}"
            regions.append((region, line_num, chrom, start_pos, end_pos))
    
    # Ordena as regiões por cromossomo e posição
    def sort_key(item):
        region, line_num, chrom, start_pos, end_pos = item
        # Tenta converter cromossomo para número, se falhar mantém como string
        try:
            if chrom.startswith('chr'):
                chrom_num = chrom[3:]
            else:
                chrom_num = chrom
            
            # Trata cromossomos especiais
            if chrom_num == 'X':
                return (23, start_pos)
            elif chrom_num == 'Y':
                return (24, start_pos)
            elif chrom_num == 'M' or chrom_num == 'MT':
                return (25, start_pos)
            else:
                return (int(chrom_num), start_pos)
        except ValueError:
            # Se não conseguir converter, ordena alfabeticamente
            return (1000, chrom, start_pos)
    
    regions.sort(key=sort_key)
    
    # Retorna apenas região e ID da linha (renumerado para ordem)
    return [(region, idx+1) for idx, (region, line_num, chrom, start_pos, end_pos) in enumerate(regions)]

def check_region_coverage(bamlist_file, region):
    """Verifica se há reads na região especificada."""
    try:
        with open(bamlist_file, 'r') as f:
            bam_files = [line.strip() for line in f if line.strip()]
        
        total_reads = 0
        for bam_file in bam_files:
            result = subprocess.run(
                ["samtools", "view", "-c", bam_file, region],
                capture_output=True, text=True, check=True
            )
            reads = int(result.stdout.strip())
            total_reads += reads
        
        return total_reads
    except Exception as e:
        print(f"Aviso: Erro ao verificar cobertura da região {region}: {e}")
        return -1

def run_freebayes_region(args):
    """Executa FreeBayes para uma região específica."""
    region, region_id, output_dir, bamlist_file, reference, extra_params = args
    
    # Nome do arquivo VCF para esta região
    vcf_file = os.path.join(output_dir, f"region_{region_id:04d}.vcf.gz")
    vcf_temp = vcf_file.replace('.gz', '')
    
    print(f"Iniciando FreeBayes para região {region} (ID: {region_id})")
    
    # Verifica cobertura da região
    read_count = check_region_coverage(bamlist_file, region)
    if read_count == 0:
        print(f"Aviso: Nenhum read encontrado na região {region}")
        # Ainda assim tenta executar o FreeBayes
    elif read_count > 0:
        print(f"✓ Região {region}: {read_count} reads encontrados")
    
    # Comando FreeBayes
    cmd = [
        "freebayes",
        "-L", bamlist_file,
        "-f", reference,
        "-r", region
    ]
    
    # Adiciona parâmetros extras se especificados
    if extra_params:
        cmd.extend(extra_params.split())
    
    print(f"Comando: {' '.join(cmd)}")
    
    try:
        # Executa FreeBayes
        with open(vcf_temp, 'w') as vcf_out:
            process = subprocess.run(
                cmd,
                stdout=vcf_out,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        
        # Verifica se o VCF tem conteúdo além do cabeçalho
        variant_count = 0
        try:
            with open(vcf_temp, 'r') as f:
                for line in f:
                    if not line.startswith('#') and line.strip():
                        variant_count += 1
        except Exception as e:
            print(f"Erro ao contar variantes em {vcf_temp}: {e}")
        
        print(f"VCF temporário gerado com {variant_count} variantes para região {region}")
        
        # Comprime o VCF usando bcftools mesmo se estiver vazio
        subprocess.run(
            ["bcftools", "view", "-O", "z", "-o", vcf_file, vcf_temp],
            check=True, capture_output=True
        )
        
        # Remove VCF descomprimido
        os.remove(vcf_temp)
        
        # Indexa o VCF comprimido usando bcftools
        subprocess.run(
            ["bcftools", "index", "-f", "-t", vcf_file],
            check=True, capture_output=True
        )
        
        print(f"✓ FreeBayes concluído para região {region} (ID: {region_id}) - {variant_count} variantes")
        return vcf_file, region, None, variant_count
        
    except subprocess.CalledProcessError as e:
        error_msg = f"Erro ao executar FreeBayes para região {region}: {e.stderr}"
        print(f"✗ {error_msg}")
        # Remove arquivo temporário se existir
        if os.path.exists(vcf_temp):
            os.remove(vcf_temp)
        return None, region, error_msg, 0

def concatenate_vcfs(vcf_files, output_file):
    """Concatena os VCFs usando bcftools concat e remove duplicatas."""
    print(f"\nConcatenando {len(vcf_files)} arquivos VCF...")
    
    # Ordena os arquivos VCF por nome (que corresponde à ordem das regiões)
    vcf_files.sort()
    
    # Cria arquivo temporário com lista de VCFs
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_file:
        for vcf_file in vcf_files:
            tmp_file.write(f"{vcf_file}\n")
        vcf_list_file = tmp_file.name
    
    # Arquivo temporário para VCF concatenado
    temp_concat = output_file.replace('.vcf.gz', '_concat.vcf.gz')
    # Arquivo temporário para VCF sem duplicatas
    temp_no_dup = output_file.replace('.vcf.gz', '_no_dup.vcf.gz')
    
    try:
        # Primeiro: Concatena os VCFs (permite sobreposições)
        print("Passo 1: Concatenando VCFs...")
        cmd_concat = [
            "bcftools", "concat",
            "-f", vcf_list_file,
            "-O", "z",
            "--allow-overlaps",  # Permite sobreposições
            "-o", temp_concat
        ]
        
        print(f"Executando: {' '.join(cmd_concat)}")
        subprocess.run(cmd_concat, check=True)
        
        # Segundo: Remove duplicatas usando bcftools norm
        print("Passo 2: Removendo variantes duplicadas...")
        cmd_norm = [
            "bcftools", "norm",
            "-d", "both",  # Remove duplicatas (both = posição e alelo)
            "-O", "z",
            "-o", temp_no_dup,
            temp_concat
        ]
        
        print(f"Executando: {' '.join(cmd_norm)}")
        result = subprocess.run(cmd_norm, check=True, capture_output=True, text=True)
        
        # Mostra estatísticas de duplicatas removidas
        if result.stderr:
            for line in result.stderr.split('\n'):
                if 'duplicate' in line.lower() or 'removed' in line.lower():
                    print(f"  {line.strip()}")
        
        # Remove arquivo temporário concatenado
        os.remove(temp_concat)
        
        # Terceiro: Ordena o VCF final
        print("Passo 3: Ordenando VCF...")
        cmd_sort = [
            "bcftools", "sort",
            "-O", "z",
            "-o", output_file,
            temp_no_dup
        ]
        
        print(f"Executando: {' '.join(cmd_sort)}")
        subprocess.run(cmd_sort, check=True)
        
        # Remove arquivo temporário sem duplicatas
        os.remove(temp_no_dup)
        
        # Quarto: Indexa o VCF final ordenado
        print("Passo 4: Indexando VCF final...")
        subprocess.run(["bcftools", "index", "-f", "-t", output_file], check=True)
        
        # Conta variantes no arquivo final
        try:
            result = subprocess.run(
                ["bcftools", "view", "-H", output_file],
                capture_output=True, text=True, check=True
            )
            final_variant_count = len([line for line in result.stdout.split('\n') if line.strip()])
            print(f"✓ VCF final criado: {output_file} ({final_variant_count} variantes únicas)")
        except:
            print(f"✓ VCF final criado: {output_file}")
        
    except subprocess.CalledProcessError as e:
        # Remove arquivos temporários em caso de erro
        for temp_file in [temp_concat, temp_no_dup]:
            if os.path.exists(temp_file):
                os.remove(temp_file)
        raise Exception(f"Erro ao processar VCFs: {e}")
    finally:
        # Remove arquivo temporário da lista
        os.unlink(vcf_list_file)

def cleanup_intermediate_files(vcf_files):
    """Remove arquivos VCF intermediários."""
    print("\nLimpando arquivos intermediários...")
    for vcf_file in vcf_files:
        try:
            if os.path.exists(vcf_file):
                os.remove(vcf_file)
            # Remove também o índice (.tbi ou .csi)
            for ext in [".tbi", ".csi"]:
                index_file = vcf_file + ext
                if os.path.exists(index_file):
                    os.remove(index_file)
        except Exception as e:
            print(f"Aviso: Não foi possível remover {vcf_file}: {e}")

def check_dependencies():
    """Verifica se as ferramentas necessárias estão disponíveis."""
    tools = ["freebayes", "bcftools", "samtools"]
    
    for tool in tools:
        if shutil.which(tool) is None:
            raise Exception(f"Ferramenta não encontrada: {tool}")
        
        # Verifica versão
        try:
            result = subprocess.run([tool, "--version"], capture_output=True, text=True)
            version_info = result.stdout.split('\n')[0] if result.stdout else result.stderr.split('\n')[0]
            print(f"✓ {tool}: {version_info}")
        except:
            print(f"✓ {tool}: disponível")

def main():
    """Função principal."""
    print("=== FreeBayes Paralelo - Versão Anti-Duplicatas ===")
    print(f"Diretório de BAMs: {BAM_DIRECTORY}")
    print(f"Genoma de referência: {REFERENCE_GENOME}")
    print(f"Arquivo BED: {BED_FILE}")
    print(f"Diretório de saída: {OUTPUT_DIR}")
    print(f"Jobs paralelos: {MAX_PARALLEL_JOBS}")
    print()
    
    bamlist_file = None
    
    try:
        # Verifica dependências
        print("Verificando dependências...")
        check_dependencies()
        
        # Valida arquivos
        print("Validando arquivos...")
        validate_files()
        print("✓ Todos os arquivos estão acessíveis")
        
        # Cria lista de BAMs
        print("Criando lista de arquivos BAM...")
        bamlist_file = create_bamlist()
        
        # Lê regiões do BED
        print("Lendo regiões do arquivo BED...")
        regions = read_bed_regions(BED_FILE)
        print(f"✓ {len(regions)} regiões encontradas")
        
        if not regions:
            print("Nenhuma região válida encontrada no arquivo BED!")
            return 1
        
        # Verifica sobreposições entre regiões
        print("Verificando sobreposições entre regiões...")
        overlaps = check_overlapping_regions(regions)
        if not overlaps:
            print("✓ Nenhuma sobreposição detectada entre regiões")
        
        # Teste com uma região primeiro (opcional)
        if len(regions) > 1:
            print(f"\nTestando primeira região: {regions[0][0]}")
            test_args = (regions[0][0], regions[0][1], OUTPUT_DIR, bamlist_file, 
                        REFERENCE_GENOME, FREEBAYES_EXTRA_PARAMS)
            test_result = run_freebayes_region(test_args)
            if test_result[0] is None:
                print("Teste falhou! Verifique configurações antes de continuar.")
                return 1
            print("✓ Teste bem-sucedido, continuando com todas as regiões...\n")
        
        # Prepara argumentos para o pool de processos
        args_list = [
            (region, region_id, OUTPUT_DIR, bamlist_file, REFERENCE_GENOME, FREEBAYES_EXTRA_PARAMS)
            for region, region_id in regions
        ]
        
        # Executa FreeBayes em paralelo
        print(f"Executando FreeBayes para {len(regions)} regiões com {MAX_PARALLEL_JOBS} jobs paralelos...")
        
        successful_vcfs = []
        failed_regions = []
        total_variants = 0
        
        with ProcessPoolExecutor(max_workers=MAX_PARALLEL_JOBS) as executor:
            # Submete todos os jobs
            future_to_region = {
                executor.submit(run_freebayes_region, args): args[0] 
                for args in args_list
            }
            
            # Coleta resultados conforme completam
            for future in as_completed(future_to_region):
                region = future_to_region[future]
                try:
                    vcf_file, region_name, error, variant_count = future.result()
                    if vcf_file:
                        successful_vcfs.append(vcf_file)
                        total_variants += variant_count
                    else:
                        failed_regions.append((region_name, error))
                except Exception as e:
                    failed_regions.append((region, str(e)))
        
        # Relatório de execução
        print(f"\n=== Relatório de Execução ===")
        print(f"Regiões processadas com sucesso: {len(successful_vcfs)}")
        print(f"Regiões com falha: {len(failed_regions)}")
        print(f"Total de variantes encontradas (com possíveis duplicatas): {total_variants}")
        
        if failed_regions:
            print("\nRegiões que falharam:")
            for region, error in failed_regions:
                print(f"  - {region}: {error}")
        
        if not successful_vcfs:
            print("Nenhum VCF foi gerado com sucesso!")
            return 1
        
        # Concatena VCFs e remove duplicatas
        final_vcf_path = os.path.join(OUTPUT_DIR, FINAL_VCF)
        concatenate_vcfs(successful_vcfs, final_vcf_path)
        
        # Limpeza (opcional - descomente se desejar remover arquivos intermediários)
        # cleanup_intermediate_files(successful_vcfs)
        
        print(f"\n✓ Processo concluído!")
        print(f"VCF final: {final_vcf_path}")
        print(f"Arquivos VCF intermediários mantidos em: {OUTPUT_DIR}")
        
        return 0
        
    except Exception as e:
        print(f"\n✗ Erro: {e}")
        return 1
    finally:
        # Remove arquivo bamlist temporário
        if bamlist_file and os.path.exists(bamlist_file):
            try:
                os.remove(bamlist_file)
                print(f"Arquivo temporário removido: {bamlist_file}")
            except:
                pass

if __name__ == "__main__":
    sys.exit(main())