import os
import subprocess

# Caminhos
bed_file = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/reference_panel.bed"   # <<== Arquivo .bed
output_dir = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/BAMs"         # <<== Pasta dos BAMs 
os.makedirs(output_dir, exist_ok=True)

# Arquivo com links
links_file = "/home/lab/Desktop/arq_joao/ANCESTRY_PANEL/1kg_hgdp_cram.txt"    # <<== Arquivo .txt com os links

with open(links_file, "r") as f:
    for line in f:
        url = line.strip()
        if not url:
            continue
        
        # Extrair o nome da amostra (HGDP0000099) do link
        filename = os.path.basename(url)
        sample_id = filename.split(".")[0]  # HGDP00704.alt_bwamem... -> HGDP00704

        output_bam = os.path.join(output_dir, f"{sample_id}.bam")
        
        # Montar o comando samtools
        cmd = [
            "samtools", "view",
            "-b", "-ML", bed_file,
            "-o", output_bam,
            url
        ]

        print(f"[INFO] Processando {sample_id} ...")
        try:
            subprocess.run(cmd, check=True)
            print(f"[OK] {sample_id} salvo em {output_bam}")
        except subprocess.CalledProcessError as e:
            print(f"[ERRO] Falhou para {sample_id}: {e}")
