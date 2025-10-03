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
    
    # Verificar se o arquivo tem tamanho mínimo razoável
    try:
        file_size = os.path.getsize(bam_path)
        if file_size < 1024:  # Menos de 1KB provavelmente está corrompido
            return False
    except:
        return False
    
    try:
        # samtools quickcheck retorna 0 se o arquivo está OK
        result = subprocess.run(
            ["samtools", "quickcheck", bam_path],
            capture_output=True,
            check=False,
            timeout=30  # Timeout para validação
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

def get_data_source(url):
    """Identifica a origem dos dados baseado na URL"""
    if "1000genomes.ebi.ac.uk" in url:
        return "1KG_EBI"
    elif "ftp.sra.ebi.ac.uk" in url:
        return "SRA_EBI"
    else:
        return "UNKNOWN"

def convert_ftp_to_https(url):
    """
    Converte URLs FTP para HTTPS baseado nos padrões dos servidores EBI
    
    Args:
        url: URL original (FTP)
    
    Returns:
        list: Lista de URLs alternativas para tentar (HTTPS primeiro, depois FTP)
    """
    alternatives = []
    
    if url.startswith("ftp://"):
        if "ftp.1000genomes.ebi.ac.uk" in url:
            # Converter 1000 Genomes FTP para HTTPS
            https_url = url.replace("ftp://ftp.1000genomes.ebi.ac.uk", 
                                  "https://ftp.1000genomes.ebi.ac.uk")
            alternatives.append(https_url)
            
        elif "ftp.sra.ebi.ac.uk" in url:
            # Converter SRA FTP para HTTPS
            https_url = url.replace("ftp://ftp.sra.ebi.ac.uk", 
                                  "https://ftp.sra.ebi.ac.uk")
            alternatives.append(https_url)
        
        # Sempre manter FTP como fallback
        alternatives.append(url)
    else:
        # Se já for HTTPS ou outro protocolo, usar como está
        alternatives.append(url)
    
    return alternatives

def process_sample_with_retry(url, bed_file, output_dir, max_retries=3):
    """
    Processa uma única amostra com sistema de retry, tentando HTTPS primeiro
    
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
    data_source = get_data_source(url)
    
    # Verificar se arquivo já existe e está íntegro
    if os.path.exists(output_bam):
        if validate_bam_file(output_bam):
            return True, sample_id, f"Arquivo válido já existe: {output_bam}"
        else:
            safe_print(f"[AVISO] {sample_id}: Arquivo existe mas está corrompido, removendo...")
            remove_corrupted_file(output_bam)
    
    # Configurar timeouts baseado na fonte
    if data_source == "SRA_EBI":
        timeout_per_attempt = 3600  # 1 hora para SRA (mais lento)
        base_delay = 60
    else:
        timeout_per_attempt = 1800  # 30 min para 1KG
        base_delay = 30
    
    # Obter URLs alternativas (HTTPS primeiro, depois FTP)
    url_alternatives = convert_ftp_to_https(url)
    
    # Tentar cada URL alternativa
    for url_idx, current_url in enumerate(url_alternatives):
        protocol = "HTTPS" if current_url.startswith("https://") else "FTP"
        
        # Se for a segunda tentativa (FTP), só tentar se HTTPS falhou por erro de protocolo
        if url_idx > 0:
            safe_print(f"[INFO] {sample_id}: Tentando fallback para {protocol}")
        
        # Tentar download com retry para esta URL
        for attempt in range(max_retries):
            temp_bam = output_bam + ".tmp"
            try:
                safe_print(f"[INFO] {sample_id} ({data_source}/{protocol}): Tentativa {attempt + 1}/{max_retries}")
                
                # Montar o comando samtools com configurações otimizadas
                cmd = [
                    "samtools", "view",
                    "-b", "-ML", bed_file,
                    "-o", temp_bam,  # Usar arquivo temporário primeiro
                    current_url.strip()
                ]
                
                # Executar comando
                start_time = time.time()
                result = subprocess.run(
                    cmd, 
                    check=True, 
                    capture_output=True, 
                    text=True,
                    timeout=timeout_per_attempt
                )
                elapsed = time.time() - start_time
                
                # Validar arquivo temporário
                if validate_bam_file(temp_bam):
                    # Mover arquivo temporário para destino final
                    os.rename(temp_bam, output_bam)
                    
                    return True, sample_id, f"Concluído via {protocol} em {elapsed:.2f}s (tentativa {attempt + 1}) - {output_bam}"
                else:
                    safe_print(f"[AVISO] {sample_id}: Arquivo gerado está corrompido (tentativa {attempt + 1})")
                    remove_corrupted_file(temp_bam)
                    raise Exception("Arquivo gerado está corrompido")
                    
            except subprocess.TimeoutExpired:
                safe_print(f"[AVISO] {sample_id}: Timeout de {timeout_per_attempt}s na tentativa {attempt + 1}")
                remove_corrupted_file(temp_bam)
                if attempt < max_retries - 1:
                    delay = min(base_delay * (2 ** attempt) + random.uniform(0, 10), 300)
                    safe_print(f"[INFO] {sample_id}: Aguardando {delay:.1f}s antes da próxima tentativa...")
                    time.sleep(delay)
                
            except subprocess.CalledProcessError as e:
                error_msg = e.stderr if e.stderr else str(e)
                safe_print(f"[AVISO] {sample_id}: Erro samtools ({protocol}) na tentativa {attempt + 1}: {error_msg.strip()}")
                remove_corrupted_file(temp_bam)
                
                # Analisar tipo de erro
                error_lower = error_msg.lower()
                is_network_error = any(keyword in error_lower for keyword in 
                       ['connection reset', 'connection refused', 'timeout', 'network', 
                        'failed to open', 'eof marker', 'truncated', 'seek at offset',
                        'container header crc32 failure', 'retrieval of region', 'error closing'])
                
                is_server_error = any(keyword in error_lower for keyword in 
                       ['no such file', '404', 'not found', 'access denied', 'permission'])
                
                is_protocol_error = any(keyword in error_lower for keyword in 
                       ['protocol not supported', 'unsupported protocol', 'ssl', 'certificate'])
                
                if is_server_error:
                    # Para erros de servidor, não vale a pena tentar novamente esta URL
                    safe_print(f"[INFO] {sample_id}: Erro de servidor com {protocol}, tentando próxima URL se disponível")
                    break
                
                elif is_protocol_error and protocol == "HTTPS":
                    # Se HTTPS falhou por problema de protocolo, pular para FTP
                    safe_print(f"[INFO] {sample_id}: Erro de protocolo HTTPS, pulando para FTP")
                    break
                
                elif is_network_error and attempt < max_retries - 1:
                    # Para erros de rede, aguardar mais tempo
                    delay = min(base_delay * (2 ** attempt) + random.uniform(0, 10), 600)  # Max 10 min
                    safe_print(f"[INFO] {sample_id}: Erro de rede/conectividade, aguardando {delay:.1f}s...")
                    time.sleep(delay)
                elif not is_network_error:
                    # Para outros tipos de erro, falhar rapidamente nesta URL
                    break
                    
            except Exception as e:
                safe_print(f"[AVISO] {sample_id}: Erro geral na tentativa {attempt + 1}: {str(e)}")
                remove_corrupted_file(temp_bam)
                if attempt < max_retries - 1:
                    delay = min(20 * (2 ** attempt) + random.uniform(0, 5), 120)
                    safe_print(f"[INFO] {sample_id}: Aguardando {delay:.1f}s antes da próxima tentativa...")
                    time.sleep(delay)
    
    # Se chegou aqui, todas as URLs e tentativas falharam
    return False, sample_id, f"Falhou com HTTPS e FTP após {max_retries} tentativas cada"

def main():
    # Configurações #nay
    bed_file = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/reference_panel.bed"
    output_dir = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/BAMs"
    links_file = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/1kg_hgdp_cram.txt"
    
    # Configurações de processamento (reduzidas para ser mais gentil)
    MAX_WORKERS = 2  # Reduzido ainda mais para evitar sobrecarga
    MAX_RETRIES = 4  # Aumentado número de tentativas
    
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
    
    # Separar por fonte de dados para estatísticas
    sources = {}
    for url in urls:
        source = get_data_source(url)
        sources[source] = sources.get(source, 0) + 1
    
    print(f"[INFO] Encontradas {len(urls)} amostras para processar")
    for source, count in sources.items():
        print(f"[INFO]   {source}: {count} amostras")
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
                sample_id = os.path.basename(url).split(".")[0]
                safe_print(f"[ERRO] {sample_id}: Exceção não tratada: {exc}")
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
        print(f"[DICA] Execute novamente para tentar as que falharam.")
        print(f"[DICA] Considere verificar conectividade de rede e disponibilidade dos servidores.")
    
    print(f"[INFO] Arquivos salvos em: {output_dir}")

if __name__ == "__main__":
    main()