#!/usr/bin/env python3
import os
import subprocess
import concurrent.futures
from threading import Lock
import time
import shutil

# Lock para sincronizar prints
print_lock = Lock()

def safe_print(*args, **kwargs):
    """Print thread-safe"""
    with print_lock:
        print(*args, **kwargs)

def run_command(cmd, timeout=None, description="comando"):
    """
    Executa um comando e retorna o resultado
    
    Args:
        cmd: Lista com comando e argumentos
        timeout: Timeout em segundos
        description: Descrição para logs
    
    Returns:
        tuple: (success, stdout, stderr, returncode)
    """
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False
        )
        return True, result.stdout, result.stderr, result.returncode
    except subprocess.TimeoutExpired:
        return False, "", f"Timeout de {timeout}s excedido", -1
    except Exception as e:
        return False, "", str(e), -1

def download_cram_with_index(url, download_dir, max_retries=3):
    """
    Baixa arquivo CRAM e seu índice usando wget com retry
    
    Args:
        url: URL do arquivo CRAM
        download_dir: Diretório para download
        max_retries: Número máximo de tentativas
    
    Returns:
        tuple: (success, cram_path, crai_path, message)
    """
    filename = os.path.basename(url)
    sample_id = filename.split(".")[0]
    cram_path = os.path.join(download_dir, filename)
    crai_path = cram_path + ".crai"
    crai_url = url + ".crai"
    
    # Verificar se arquivos já existem
    if os.path.exists(cram_path) and os.path.exists(crai_path):
        # Verificar integridade básica (tamanho > 0)
        try:
            if os.path.getsize(cram_path) > 0 and os.path.getsize(crai_path) > 0:
                return True, cram_path, crai_path, "Arquivos já existem"
        except:
            pass
    
    # Tentar download com retry
    for attempt in range(max_retries):
        try:
            safe_print(f"[DOWNLOAD] {sample_id}: Tentativa {attempt + 1}/{max_retries}")
            
            # Download do arquivo CRAM
            safe_print(f"[DOWNLOAD] {sample_id}: Baixando CRAM...")
            cram_cmd = [
                "wget", 
                "--continue",           # Continuar downloads interrompidos
                "--timeout=300",        # Timeout de conexão: 5 min
                "--tries=3",           # Tentativas internas do wget
                "--waitretry=30",      # Espera entre tentativas
                "--progress=dot:giga", # Progress mais limpo
                "-O", cram_path,
                url
            ]
            
            success, stdout, stderr, returncode = run_command(
                cram_cmd, timeout=7200, description=f"download CRAM {sample_id}"  # 2 horas
            )
            
            if not success or returncode != 0:
                safe_print(f"[AVISO] {sample_id}: Falha no download CRAM: {stderr}")
                # Limpar arquivo parcial se existir
                if os.path.exists(cram_path):
                    try:
                        os.remove(cram_path)
                    except:
                        pass
                continue
            
            # Download do índice CRAI
            safe_print(f"[DOWNLOAD] {sample_id}: Baixando índice CRAI...")
            crai_cmd = [
                "wget", 
                "--continue",
                "--timeout=60",
                "--tries=3",
                "--waitretry=10",
                "--progress=dot:mega",
                "-O", crai_path,
                crai_url
            ]
            
            success, stdout, stderr, returncode = run_command(
                crai_cmd, timeout=600, description=f"download CRAI {sample_id}"  # 10 min
            )
            
            if not success or returncode != 0:
                safe_print(f"[AVISO] {sample_id}: Falha no download CRAI: {stderr}")
                # Limpar arquivos se CRAI falhou
                for f in [cram_path, crai_path]:
                    if os.path.exists(f):
                        try:
                            os.remove(f)
                        except:
                            pass
                continue
            
            # Verificar se ambos arquivos foram baixados com sucesso
            if os.path.exists(cram_path) and os.path.exists(crai_path):
                cram_size = os.path.getsize(cram_path)
                crai_size = os.path.getsize(crai_path)
                
                if cram_size > 1024 and crai_size > 0:  # Tamanhos mínimos razoáveis
                    return True, cram_path, crai_path, f"Download concluído (CRAM: {cram_size/1024/1024:.1f}MB, CRAI: {crai_size/1024:.1f}KB)"
                else:
                    safe_print(f"[AVISO] {sample_id}: Arquivos com tamanhos suspeitos (CRAM: {cram_size}, CRAI: {crai_size})")
            
        except Exception as e:
            safe_print(f"[AVISO] {sample_id}: Erro na tentativa {attempt + 1}: {str(e)}")
            # Limpar arquivos parciais
            for f in [cram_path, crai_path]:
                if os.path.exists(f):
                    try:
                        os.remove(f)
                    except:
                        pass
        
        # Esperar entre tentativas
        if attempt < max_retries - 1:
            delay = min(60 * (attempt + 1), 300)  # Até 5 min
            safe_print(f"[INFO] {sample_id}: Aguardando {delay}s antes da próxima tentativa...")
            time.sleep(delay)
    
    return False, None, None, f"Falha no download após {max_retries} tentativas"

def process_cram_to_bam(cram_path, crai_path, bed_file, output_bam):
    """
    Processa arquivo CRAM para BAM usando regiões do BED
    
    Args:
        cram_path: Caminho do arquivo CRAM
        crai_path: Caminho do arquivo de índice CRAI
        bed_file: Arquivo BED com regiões
        output_bam: Caminho de saída do BAM
    
    Returns:
        tuple: (success, message)
    """
    sample_id = os.path.basename(cram_path).split(".")[0]
    
    try:
        safe_print(f"[PROCESS] {sample_id}: Recortando regiões do CRAM...")
        
        # Comando samtools view
        cmd = [
            "samtools", "view",
            "-b",              # Saída em formato BAM
            "-ML", bed_file,   # Filtrar usando arquivo BED
            "-o", output_bam,  # Arquivo de saída
            cram_path          # Arquivo CRAM de entrada
        ]
        
        start_time = time.time()
        success, stdout, stderr, returncode = run_command(
            cmd, timeout=3600, description=f"samtools view {sample_id}"  # 1 hora
        )
        elapsed = time.time() - start_time
        
        if not success or returncode != 0:
            return False, f"Erro no samtools: {stderr}"
        
        # Verificar se o BAM foi gerado corretamente
        if not os.path.exists(output_bam):
            return False, "Arquivo BAM não foi criado"
        
        bam_size = os.path.getsize(output_bam)
        if bam_size < 1024:  # Menos de 1KB é suspeito
            return False, f"Arquivo BAM muito pequeno ({bam_size} bytes)"
        
        # Validar BAM usando samtools quickcheck
        cmd_check = ["samtools", "quickcheck", output_bam]
        success_check, _, stderr_check, returncode_check = run_command(
            cmd_check, timeout=60, description=f"quickcheck {sample_id}"
        )
        
        if not success_check or returncode_check != 0:
            return False, f"BAM corrompido (quickcheck falhou): {stderr_check}"
        
        return True, f"Processamento concluído em {elapsed:.1f}s (BAM: {bam_size/1024/1024:.2f}MB)"
        
    except Exception as e:
        return False, f"Erro durante processamento: {str(e)}"

def secure_delete_file(file_path):
    """
    Remove arquivo de forma segura e definitiva
    
    Args:
        file_path: Caminho do arquivo a ser removido
    
    Returns:
        bool: True se removido com sucesso
    """
    try:
        if os.path.exists(file_path):
            # Remover arquivo definitivamente
            os.remove(file_path)
            return True
        return True  # Se não existe, considerar como sucesso
    except Exception as e:
        safe_print(f"[ERRO] Falha ao remover {file_path}: {e}")
        return False

def process_single_sample(url, bed_file, output_dir, temp_dir):
    """
    Processa uma única amostra: download -> recorte -> limpeza
    
    Args:
        url: URL do arquivo CRAM
        bed_file: Arquivo BED com regiões
        output_dir: Diretório de saída dos BAMs
        temp_dir: Diretório temporário para CRAMs
    
    Returns:
        tuple: (success, sample_id, message)
    """
    if not url.strip():
        return False, "", "URL vazia"
    
    filename = os.path.basename(url.strip())
    sample_id = filename.split(".")[0]
    output_bam = os.path.join(output_dir, f"{sample_id}.bam")
    
    # Verificar se BAM já existe e está válido
    if os.path.exists(output_bam):
        cmd_check = ["samtools", "quickcheck", output_bam]
        success, _, _, returncode = run_command(cmd_check, timeout=30)
        if success and returncode == 0:
            return True, sample_id, "BAM válido já existe"
    
    try:
        # Etapa 1: Download do CRAM e índice
        success, cram_path, crai_path, download_msg = download_cram_with_index(
            url.strip(), temp_dir
        )
        
        if not success:
            return False, sample_id, f"Falha no download: {download_msg}"
        
        safe_print(f"[INFO] {sample_id}: {download_msg}")
        
        # Etapa 2: Processar CRAM para BAM
        success, process_msg = process_cram_to_bam(cram_path, crai_path, bed_file, output_bam)
        
        if not success:
            # Limpar arquivos temporários em caso de erro
            secure_delete_file(cram_path)
            secure_delete_file(crai_path)
            return False, sample_id, f"Falha no processamento: {process_msg}"
        
        safe_print(f"[INFO] {sample_id}: {process_msg}")
        
        # Etapa 3: Limpeza dos arquivos temporários
        safe_print(f"[CLEANUP] {sample_id}: Removendo arquivos temporários...")
        cram_deleted = secure_delete_file(cram_path)
        crai_deleted = secure_delete_file(crai_path)
        
        if cram_deleted and crai_deleted:
            cleanup_msg = "Arquivos temporários removidos"
        else:
            cleanup_msg = "Aviso: Alguns arquivos temporários não foram removidos"
        
        return True, sample_id, f"Concluído com sucesso. {cleanup_msg}"
        
    except Exception as e:
        # Tentar limpar em caso de exceção
        if 'cram_path' in locals():
            secure_delete_file(cram_path)
        if 'crai_path' in locals():
            secure_delete_file(crai_path)
        return False, sample_id, f"Erro inesperado: {str(e)}"

def main():
    # Configurações #nay
    bed_file = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/reference_panel.bed"
    output_dir = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/BAMs"
    temp_dir = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/temp_crams"
    links_file = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/1kg_hgdp_cram.txt"
    
    # Configurações de processamento
    MAX_WORKERS = 2  # Processar 2 amostras simultaneamente
    
    # Criar diretórios
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    
    # Ler URLs
    urls = []
    try:
        with open(links_file, "r") as f:
            urls = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"[ERRO] Arquivo não encontrado: {links_file}")
        return
    
    if not urls:
        print("[ERRO] Nenhuma URL encontrada no arquivo")
        return
    
    print(f"[INFO] Encontradas {len(urls)} amostras para processar")
    print(f"[INFO] Usando {MAX_WORKERS} processos paralelos")
    print(f"[INFO] Arquivos temporários em: {temp_dir}")
    print(f"[INFO] Arquivos finais em: {output_dir}")
    print("-" * 80)
    
    # Verificar dependências
    for tool in ["wget", "samtools"]:
        success, _, stderr, returncode = run_command(["which", tool], timeout=10)
        if not success or returncode != 0:
            print(f"[ERRO] {tool} não encontrado no sistema")
            return
    
    # Contadores
    successful = 0
    failed = 0
    already_exists = 0
    start_time = time.time()
    
    # Processar amostras
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # Submeter tarefas
        future_to_url = {
            executor.submit(process_single_sample, url, bed_file, output_dir, temp_dir): url
            for url in urls
        }
        
        # Processar resultados
        for future in concurrent.futures.as_completed(future_to_url):
            url = future_to_url[future]
            try:
                success, sample_id, message = future.result()
                if success:
                    if "já existe" in message:
                        safe_print(f"[SKIP] {sample_id}: {message}")
                        already_exists += 1
                    else:
                        safe_print(f"[OK] {sample_id}: {message}")
                        successful += 1
                else:
                    safe_print(f"[ERRO] {sample_id}: {message}")
                    failed += 1
            except Exception as exc:
                sample_id = os.path.basename(url).split(".")[0] if url else "UNKNOWN"
                safe_print(f"[ERRO] {sample_id}: Exceção: {exc}")
                failed += 1
    
    # Relatório final
    total_time = time.time() - start_time
    print("-" * 80)
    print(f"[RESUMO] Processamento concluído em {total_time:.1f}s ({total_time/3600:.1f}h)")
    print(f"[RESUMO] Processados com sucesso: {successful}")
    print(f"[RESUMO] Já existiam (válidos): {already_exists}")
    print(f"[RESUMO] Falharam: {failed}")
    print(f"[RESUMO] Total: {successful + failed + already_exists}")
    
    # Verificar diretório temporário
    try:
        temp_files = os.listdir(temp_dir)
        if temp_files:
            print(f"[AVISO] {len(temp_files)} arquivos restantes no diretório temporário")
            print(f"[INFO] Você pode limpar manualmente: rm -rf {temp_dir}/*")
        else:
            print(f"[INFO] Diretório temporário limpo")
    except:
        pass
    
    print(f"[INFO] Arquivos BAM salvos em: {output_dir}")

if __name__ == "__main__":
    main()