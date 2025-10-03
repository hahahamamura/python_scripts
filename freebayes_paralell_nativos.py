#!/usr/bin/env python3

import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import tempfile
import resource
import math

# CONFIGURAÇÕES
##############################################################################################
BAM_DIRECTORY = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/BAMs" #Diretório dos bams
REFERENCE_GENOME = "/home/lab/Desktop/arq_joao/NativoAmericanas/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" #Genoma de referencia
BED_FILE = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/reference_panel.bed" #Arquivo .bed
OUTPUT_DIR = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/freebayes" #Diretório de saida dos vcfs
FINAL_VCF = "HGDP_51_NATIVO_AMERICANOS.vcf.gz" #Nome do vcf final
MAX_PARALLEL_JOBS = 4 #Número de threads
##############################################################################################

# Configurações para concatenação
BATCH_SIZE = 500  # Processar em lotes menores varias vezes pra concatenar
MAX_OPEN_FILES = 900  # Limite seguro de arquivos abertos

def increase_file_limits():
    """Aumenta os limites de arquivos abertos quando possível."""
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        # Tenta aumentar para o máximo permitido
        new_soft = min(hard, 8192)
        resource.setrlimit(resource.RLIMIT_NOFILE, (new_soft, hard))
        print(f"Limite de arquivos aumentado para: {new_soft}")
        return new_soft
    except:
        soft, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
        print(f"Mantendo limite atual de arquivos: {soft}")
        return soft

def create_bamlist():
    """Cria lista de BAMs do diretório."""
    bam_files = []
    for file_name in os.listdir(BAM_DIRECTORY):
        if file_name.lower().endswith('.bam'):
            bam_files.append(os.path.join(BAM_DIRECTORY, file_name))
    
    bam_files.sort()
    bamlist_file = os.path.join(OUTPUT_DIR, "bamlist_temp.txt")
    
    with open(bamlist_file, 'w') as f:
        for bam_file in bam_files:
            f.write(f"{bam_file}\n")
    
    print(f"Criada lista com {len(bam_files)} arquivos BAM")
    return bamlist_file

def read_bed_regions(bed_file):
    """Lê regiões do BED."""
    regions = []
    with open(bed_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            
            chrom, start, end = parts[0], parts[1], parts[2]
            region = f"{chrom}:{start}-{end}"
            regions.append((region, line_num))
    
    print(f"Lidas {len(regions)} regiões do arquivo BED")
    return regions

def check_existing_vcfs(regions, output_dir):
    """Verifica quais VCFs já existem e estão válidos."""
    existing_vcfs = []
    missing_regions = []
    invalid_vcfs = []
    
    print("Verificando VCFs existentes...")
    
    for region, region_id in regions:
        vcf_file = os.path.join(output_dir, f"region_{region_id:04d}.vcf.gz")
        idx_file = vcf_file + ".tbi"
        
        # Verifica se o arquivo existe e tem tamanho > 0
        if os.path.exists(vcf_file) and os.path.getsize(vcf_file) > 0:
            # Verifica se o índice existe
            if os.path.exists(idx_file):
                # Testa se o arquivo é válido
                try:
                    # Teste rápido: tenta ler o cabeçalho
                    result = subprocess.run([
                        "bcftools", "view", "-h", vcf_file
                    ], capture_output=True, check=True)
                    
                    existing_vcfs.append(vcf_file)
                    
                except subprocess.CalledProcessError:
                    print(f"  ⚠️  VCF corrompido: {vcf_file}")
                    invalid_vcfs.append((vcf_file, idx_file))
                    missing_regions.append((region, region_id))
            else:
                print(f"  ⚠️  Falta índice para: {vcf_file}")
                # Tenta criar o índice
                try:
                    subprocess.run(["bcftools", "index", "-f", "-t", vcf_file], 
                                  check=True, capture_output=True)
                    existing_vcfs.append(vcf_file)
                    print(f"  ✅ Índice criado para: {vcf_file}")
                except subprocess.CalledProcessError:
                    print(f"  ❌ Erro ao indexar: {vcf_file}")
                    invalid_vcfs.append((vcf_file, None))
                    missing_regions.append((region, region_id))
        else:
            missing_regions.append((region, region_id))
    
    # Remove arquivos corrompidos
    for vcf_file, idx_file in invalid_vcfs:
        try:
            if os.path.exists(vcf_file):
                os.remove(vcf_file)
            if idx_file and os.path.exists(idx_file):
                os.remove(idx_file)
            print(f"  🗑️  Removido arquivo corrompido: {vcf_file}")
        except OSError:
            print(f"  ⚠️  Não foi possível remover: {vcf_file}")
    
    print(f"\nResumo da verificação:")
    print(f"  VCFs existentes válidos: {len(existing_vcfs)}")
    print(f"  Regiões que precisam ser processadas: {len(missing_regions)}")
    print(f"  VCFs corrompidos removidos: {len(invalid_vcfs)}")
    
    return existing_vcfs, missing_regions

def run_freebayes_region(args):
    """Executa FreeBayes para uma região."""
    region, region_id, output_dir, bamlist_file, reference = args
    
    vcf_file = os.path.join(output_dir, f"region_{region_id:04d}.vcf.gz")
    vcf_temp = vcf_file.replace('.gz', '')
    
    cmd = [
        "freebayes",
        "-L", bamlist_file,
        "-f", reference,
        "-r", region
    ]
    
    try:
        print(f"  Processando região {region_id}: {region}")
        
        with open(vcf_temp, 'w') as vcf_out:
            result = subprocess.run(cmd, stdout=vcf_out, stderr=subprocess.PIPE, check=True)
        
        # Verifica se o arquivo foi criado e não está vazio
        if not os.path.exists(vcf_temp) or os.path.getsize(vcf_temp) == 0:
            if os.path.exists(vcf_temp):
                os.remove(vcf_temp)
            return None
        
        # Comprime
        subprocess.run(["bcftools", "view", "-O", "z", "-o", vcf_file, vcf_temp], 
                      check=True, capture_output=True)
        os.remove(vcf_temp)
        
        # Indexa
        subprocess.run(["bcftools", "index", "-f", "-t", vcf_file], 
                      check=True, capture_output=True)
        
        print(f"  ✅ Concluído região {region_id}: {region}")
        return vcf_file
        
    except subprocess.CalledProcessError as e:
        print(f"  ❌ Erro na região {region} (ID: {region_id}): {e.stderr.decode() if e.stderr else str(e)}")
        if os.path.exists(vcf_temp):
            os.remove(vcf_temp)
        return None

def concatenate_vcfs_in_batches(vcf_files, output_file, batch_size=BATCH_SIZE):
    """Concatena VCFs em lotes menores para evitar limite de arquivos abertos."""
    vcf_files.sort()
    num_files = len(vcf_files)
    num_batches = math.ceil(num_files / batch_size)
    
    print(f"Concatenando {num_files} arquivos VCF em {num_batches} lotes...")
    
    batch_files = []
    
    try:
        # Primeira fase: criar lotes intermediários
        for batch_num in range(num_batches):
            start_idx = batch_num * batch_size
            end_idx = min(start_idx + batch_size, num_files)
            batch_vcfs = vcf_files[start_idx:end_idx]
            
            batch_output = os.path.join(OUTPUT_DIR, f"batch_{batch_num:04d}.vcf.gz")
            
            print(f"Processando lote {batch_num + 1}/{num_batches} ({len(batch_vcfs)} arquivos)...")
            
            # Cria arquivo de lista para este lote
            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_file:
                for vcf_file in batch_vcfs:
                    tmp_file.write(f"{vcf_file}\n")
                batch_list_file = tmp_file.name
            
            try:
                # Concatena este lote
                subprocess.run([
                    "bcftools", "concat",
                    "-f", batch_list_file,
                    "-O", "z",
                    "--allow-overlaps",
                    "--remove-duplicates",
                    "-o", batch_output
                ], check=True, capture_output=True)
                
                # Indexa o lote
                subprocess.run(["bcftools", "index", "-f", "-t", batch_output], 
                              check=True, capture_output=True)
                
                batch_files.append(batch_output)
                
            finally:
                os.unlink(batch_list_file)
        
        # Segunda fase: concatenar os lotes
        if len(batch_files) == 1:
            # Se só há um lote, apenas renomeia
            os.rename(batch_files[0], output_file)
            os.rename(batch_files[0] + ".tbi", output_file + ".tbi")
        else:
            print(f"Concatenando {len(batch_files)} lotes finais...")
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_file:
                for batch_file in batch_files:
                    tmp_file.write(f"{batch_file}\n")
                final_list_file = tmp_file.name
            
            try:
                # Concatenação final
                subprocess.run([
                    "bcftools", "concat",
                    "-f", final_list_file,
                    "-O", "z",
                    "--allow-overlaps",
                    "--remove-duplicates",
                    "-o", output_file
                ], check=True, capture_output=True)
                
            finally:
                os.unlink(final_list_file)
        
        # Remove duplicatas e normaliza o arquivo final
        temp_final = output_file.replace('.vcf.gz', '_normalized.vcf.gz')
        
        print("Removendo duplicatas e normalizando...")
        subprocess.run([
            "bcftools", "norm",
            "-d", "all",  # Remove todas as duplicatas
            "-O", "z",
            "-o", temp_final,
            output_file
        ], check=True, capture_output=True)
        
        # Substitui pelo arquivo normalizado
        os.rename(temp_final, output_file)
        
        # Indexa o arquivo final
        subprocess.run(["bcftools", "index", "-f", "-t", output_file], 
                      check=True, capture_output=True)
        
        print(f"Arquivo final criado: {output_file}")
        
    finally:
        # Limpa arquivos temporários de lote
        for batch_file in batch_files:
            if os.path.exists(batch_file):
                os.remove(batch_file)
            if os.path.exists(batch_file + ".tbi"):
                os.remove(batch_file + ".tbi")

def cleanup_temp_files(successful_vcfs):
    """Remove arquivos VCF temporários após concatenação."""
    print("Limpando arquivos temporários...")
    removed_count = 0
    for vcf_file in successful_vcfs:
        try:
            if os.path.exists(vcf_file):
                os.remove(vcf_file)
                removed_count += 1
            
            # Remove índice também
            idx_file = vcf_file + ".tbi"
            if os.path.exists(idx_file):
                os.remove(idx_file)
        except OSError:
            pass  # Ignora erros de remoção
    
    print(f"Removidos {removed_count} arquivos temporários")

def get_vcf_stats(vcf_file):
    """Obtém estatísticas básicas do VCF final."""
    try:
        # Conta variantes
        result = subprocess.run([
            "bcftools", "view", "-H", vcf_file, "|", "wc", "-l"
        ], shell=True, capture_output=True, text=True, check=True)
        
        variant_count = result.stdout.strip()
        
        # Conta amostras
        result = subprocess.run([
            "bcftools", "query", "-l", vcf_file, "|", "wc", "-l"
        ], shell=True, capture_output=True, text=True, check=True)
        
        sample_count = result.stdout.strip()
        
        file_size = os.path.getsize(vcf_file) / (1024 * 1024)  # MB
        
        print(f"\nEstatísticas do arquivo final:")
        print(f"  Variantes: {variant_count}")
        print(f"  Amostras: {sample_count}")
        print(f"  Tamanho: {file_size:.2f} MB")
        
    except subprocess.CalledProcessError:
        print("Não foi possível obter estatísticas do arquivo final")

def main():
    print("=== FreeBayes Parallel - Processamento de Nativos Americanos ===")
    
    # Aumenta limites de arquivos
    file_limit = increase_file_limits()
    
    # Ajusta tamanho do lote baseado no limite de arquivos
    global BATCH_SIZE
    BATCH_SIZE = min(BATCH_SIZE, file_limit // 2)
    print(f"Tamanho do lote ajustado para: {BATCH_SIZE}")
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Lê regiões do BED
    regions = read_bed_regions(BED_FILE)
    
    if not regions:
        print("Erro: Nenhuma região válida encontrada no arquivo BED!")
        return
    
    # Verifica VCFs existentes
    existing_vcfs, missing_regions = check_existing_vcfs(regions, OUTPUT_DIR)
    
    # Se há regiões faltantes, processa elas
    newly_created_vcfs = []
    if missing_regions:
        print(f"\n🔧 Processando {len(missing_regions)} regiões faltantes com FreeBayes...")
        
        # Cria bamlist
        bamlist_file = create_bamlist()
        
        try:
            # Prepara argumentos apenas para regiões faltantes
            args_list = [
                (region, region_id, OUTPUT_DIR, bamlist_file, REFERENCE_GENOME)
                for region, region_id in missing_regions
            ]
            
            print(f"Iniciando processamento paralelo com {MAX_PARALLEL_JOBS} workers...")
            
            # Executa em paralelo
            failed_count = 0
            
            with ProcessPoolExecutor(max_workers=MAX_PARALLEL_JOBS) as executor:
                results = list(executor.map(run_freebayes_region, args_list))
                
                for result in results:
                    if result is not None:
                        newly_created_vcfs.append(result)
                    else:
                        failed_count += 1
            
            print(f"\nResultados do processamento das regiões faltantes:")
            print(f"  Sucessos: {len(newly_created_vcfs)}")
            print(f"  Falhas: {failed_count}")
        
        finally:
            # Limpa bamlist
            if os.path.exists(bamlist_file):
                os.remove(bamlist_file)
                print("Arquivo bamlist temporário removido")
    else:
        print("✅ Todos os VCFs já existem e são válidos!")
    
    # Combina VCFs existentes + novos
    all_vcfs = existing_vcfs + newly_created_vcfs
    
    if all_vcfs:
        print(f"\n📁 Total de VCFs para concatenação: {len(all_vcfs)}")
        print(f"  - VCFs já existentes: {len(existing_vcfs)}")
        print(f"  - VCFs recém-criados: {len(newly_created_vcfs)}")
        
        # Concatena todos os VCFs
        final_vcf_path = os.path.join(OUTPUT_DIR, FINAL_VCF)
        concatenate_vcfs_in_batches(all_vcfs, final_vcf_path)
        
        # Obtém estatísticas
        get_vcf_stats(final_vcf_path)
        
        # Remove apenas os VCFs temporários (não os que já existiam)
        if newly_created_vcfs:
            cleanup_temp_files(all_vcfs)  # Remove todos após concatenação
        
        print(f"\n✅ Processamento concluído com sucesso!")
        print(f"Arquivo final: {final_vcf_path}")
    else:
        print("\n❌ Erro: Nenhum arquivo VCF válido encontrado!")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n⚠️  Processamento interrompido pelo usuário")
    except Exception as e:
        print(f"\n❌ Erro inesperado: {str(e)}")
        raise