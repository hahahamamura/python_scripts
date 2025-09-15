import os
import subprocess
import concurrent.futures
from threading import Lock
import time
import random

# Lock para sincronizar prints
print_lock = Lock()

def safe_print(*args, **kwargs):
    """Print thread-safe"""
    with print_lock:
        print(*args, **kwargs)

def validate_bam_file(bam_path):
    """
    Valida se o arquivo BAM está íntegro usando samtools quickcheck
    
    Args:
        bam_path: Caminho para o arquivo BAM
    
    Returns:
        bool: True se o arquivo está íntegro, False caso contrário
    """
    if not os.path.exists(bam_path):
        return False
    
    try:
        # samtools quickcheck retorna 0 se o arquivo está OK
        result = subprocess.run(
            ["samtools", "quickcheck", bam_path],
            capture_output=True,
            check=False
        )
        return result.returncode == 0
    except Exception:
        return False

def remove_corrupted_file(file_path):
    """Remove arquivo corrompido de forma segura"""
    try:
        if os.path.exists(file_path):
            os.remove(file_path)
            safe_print(f"[CLEANUP] Arquivo corrompido removido: {file_path}")
    except Exception as e:
        safe_print(f"[ERRO] Falha ao remover arquivo corrompido {file_path}: {e}")

def process_sample_with_retry(url, bed_file, output_dir, max_retries=3):
    """
    Processa uma única amostra com sistema de retry
    
    Args:
        url: URL do arquivo CRAM
        bed_file: Caminho para o arquivo .bed
        output_dir: Diretório de saída
        max_retries: Número máximo de tentativas
    
    Returns:
        tuple: (success, sample_id, message)
    """
    if not url.strip():
        return False, "", "URL vazia"
    
    # Extrair o nome da amostra
    filename = os.path.basename(url.strip())
    sample_id = filename.split(".")[0]
    output_bam = os.path.join(output_dir, f"{sample_id}.bam")
    
    # Verificar se arquivo já existe e está íntegro
    if os.path.exists(output_bam):
        if validate_bam_file(output_bam):
            return True, sample_id, f"Arquivo válido já existe: {output_bam}"
        else:
            safe_print(f"[AVISO] {sample_id}: Arquivo existe mas está corrompido, removendo...")
            remove_corrupted_file(output_bam)
    
    # Tentar download com retry
    for attempt in range(max_retries):
        try:
            safe_print(f"[INFO] {sample_id}: Tentativa {attempt + 1}/{max_retries}")
            
            # Montar o comando samtools
            cmd = [
                "samtools", "view",
                "-b", "-ML", bed_file,
                "-o", output_bam,
                url.strip()
            ]
            
            # Executar comando
            start_time = time.time()
            result = subprocess.run(
                cmd, 
                check=True, 
                capture_output=True, 
                text=True,
                timeout=1800  # Timeout de 30 minutos por tentativa
            )
            elapsed = time.time() - start_time
            
            # Validar arquivo gerado
            if validate_bam_file(output_bam):
                return True, sample_id, f"Concluído em {elapsed:.2f}s (tentativa {attempt + 1}) - {output_bam}"
            else:
                safe_print(f"[AVISO] {sample_id}: Arquivo gerado está corrompido (tentativa {attempt + 1})")
                remove_corrupted_file(output_bam)
                raise Exception("Arquivo gerado está corrompido")
                
        except subprocess.TimeoutExpired:
            safe_print(f"[AVISO] {sample_id}: Timeout na tentativa {attempt + 1}")
            remove_corrupted_file(output_bam)
            if attempt < max_retries - 1:
                delay = min(60 * (2 ** attempt) + random.uniform(0, 10), 300)  # Max 5 min
                safe_print(f"[INFO] {sample_id}: Aguardando {delay:.1f}s antes da próxima tentativa...")
                time.sleep(delay)
            
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr if e.stderr else str(e)
            safe_print(f"[AVISO] {sample_id}: Erro samtools na tentativa {attempt + 1}: {error_msg}")
            remove_corrupted_file(output_bam)
            
            # Se for erro de conexão, aguardar antes de tentar novamente
            if any(keyword in error_msg.lower() for keyword in 
                   ['connection reset', 'connection refused', 'timeout', 'network', 'failed to open']):
                if attempt < max_retries - 1:
                    delay = min(30 * (2 ** attempt) + random.uniform(0, 5), 180)  # Max 3 min
                    safe_print(f"[INFO] {sample_id}: Erro de rede, aguardando {delay:.1f}s...")
                    time.sleep(delay)
            else:
                # Para outros tipos de erro, não vale a pena tentar novamente
                break
                
        except Exception as e:
            safe_print(f"[AVISO] {sample_id}: Erro geral na tentativa {attempt + 1}: {str(e)}")
            remove_corrupted_file(output_bam)
            if attempt < max_retries - 1:
                delay = min(20 * (2 ** attempt) + random.uniform(0, 5), 120)  # Max 2 min
                safe_print(f"[INFO] {sample_id}: Aguardando {delay:.1f}s antes da próxima tentativa...")
                time.sleep(delay)
    
    # Se chegou aqui, todas as tentativas falharam
    return False, sample_id, f"Falhou após {max_retries} tentativas"

def main():
    # Configurações
    bed_file = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/reference_panel.bed"
    output_dir = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/BAMs"
    links_file = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/1kg_hgdp_cram.txt"
    
    # Configurações de processamento
    MAX_WORKERS = 3  # Reduzido para ser mais gentil com o servidor
    MAX_RETRIES = 3  # Número de tentativas por amostra
    
    # Criar diretório de saída
    os.makedirs(output_dir, exist_ok=True)
    
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
    print(f"[INFO] Usando {MAX_WORKERS} threads paralelas")
    print(f"[INFO] Máximo de {MAX_RETRIES} tentativas por amostra")
    print("-" * 70)
    
    # Contadores
    successful = 0
    failed = 0
    already_exists = 0
    start_time = time.time()
    
    # Processar em paralelo
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # Submeter todas as tarefas
        future_to_url = {
            executor.submit(process_sample_with_retry, url, bed_file, output_dir, MAX_RETRIES): url 
            for url in urls
        }
        
        # Processar resultados conforme completam
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
                safe_print(f"[ERRO] {url}: Exceção não tratada: {exc}")
                failed += 1
    
    # Relatório final
    total_time = time.time() - start_time
    print("-" * 70)
    print(f"[RESUMO] Processamento concluído em {total_time:.2f}s ({total_time/3600:.1f}h)")
    print(f"[RESUMO] Downloads bem-sucedidos: {successful}")
    print(f"[RESUMO] Arquivos já existentes (válidos): {already_exists}")
    print(f"[RESUMO] Falhas definitivas: {failed}")
    print(f"[RESUMO] Total processado: {successful + failed + already_exists}")
    
    if failed > 0:
        print(f"[AVISO] {failed} amostras falharam após múltiplas tentativas.")
        print(f"[DICA] Você pode executar o script novamente para tentar as que falharam.")
    
    print(f"[INFO] Arquivos salvos em: {output_dir}")

if __name__ == "__main__":
    main()