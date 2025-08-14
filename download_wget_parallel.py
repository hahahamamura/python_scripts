import concurrent.futures
import subprocess
from pathlib import Path

# Caminho para o arquivo com os links
LINKS_FILE = "download.txt"
# Quantidade de downloads simultâneos
MAX_WORKERS = 3

def baixar_arquivo(link):
    """Executa o comando wget com --continue para baixar o arquivo."""
    try:
        subprocess.run(["wget", "--continue", link], check=True)
        print(f"[✔] Download concluído: {link}")
    except subprocess.CalledProcessError as e:
        print(f"[✘] Erro ao baixar: {link}\n{e}")

def main():
    path = Path(LINKS_FILE)
    if not path.is_file():
        print(f"Arquivo {LINKS_FILE} não encontrado.")
        return

    # Lê os links do arquivo, removendo linhas em branco
    with open(LINKS_FILE, "r") as f:
        links = [linha.strip() for linha in f if linha.strip()]

    print(f"Total de arquivos para baixar: {len(links)}")

    # Executa até 3 downloads paralelamente
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        executor.map(baixar_arquivo, links)

if __name__ == "__main__":
    main()

