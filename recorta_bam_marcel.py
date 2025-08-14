import os
import subprocess

# Caminhos base
base_dir = ""
bam_dir = os.path.join("/media/lab/Storage")
output_dir = os.path.join("/home/lab/Downloads", "BAMs_RECORTADOS")
bed_path = os.path.join("/home/lab/Desktop/arq_joao/python/data/BEDs", "optimized.bed")
amostras_path = os.path.join("/home/lab/Downloads", "amostras.txt")

# Garante que a pasta de saída exista
os.makedirs(output_dir, exist_ok=True)

# Lê a lista de amostras (linhas como "hg1")
with open(amostras_path) as f:
    amostras = [linha.strip().lower() for linha in f if linha.strip()]

for amostra in amostras:
    # Busca o arquivo BAM correspondente
    bam_nome = None
    for arquivo in os.listdir(bam_dir):
        if arquivo.lower().startswith(amostra) and arquivo.endswith(".bam"):
            bam_nome = arquivo
            break

    if not bam_nome:
        print(f"❌ BAM não encontrado para a amostra: {amostra}")
        continue

    bam_path = os.path.join(bam_dir, bam_nome)
    bam_saida = os.path.join(output_dir, f"{amostra}_recortado.bam")

    # Comando samtools view -L
    cmd = [
        "samtools", "view", "-b",
        "-L", bed_path,
        "-o", bam_saida,
        bam_path
    ]
    print(f"Comando: samtools view -b -L {bed_path} -o {bam_saida} {bam_path}")
    print(f"🔄 Processando {bam_nome} → {bam_saida}")
    subprocess.run(cmd)

print("✅ Recorte de todos os BAMs concluído.")
