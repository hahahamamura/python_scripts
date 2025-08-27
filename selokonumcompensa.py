import os

# caminho para o arquivo com as URLs
url_file = "/home/lab/Desktop/arq_joao/NativoAmericanas/amostras_nativos.txt"

# diretório onde estão os BAMs baixados
bam_dir = "/home/lab/Desktop/arq_joao/NativoAmericanas/bams"

# lê todas as linhas do arquivo
with open(url_file) as f:
    urls = [line.strip() for line in f if line.strip()]

# extrai o ID da amostra (ex: HGDP00995) de cada URL
ids = [url.split("/")[-1].split(".")[0] for url in urls]

# lista os BAMs já baixados
baixados = {f.split(".")[0] for f in os.listdir(bam_dir) if f.endswith(".bam")}

# verifica cada amostra
for sample_id in ids:
    if sample_id in baixados:
        print(f"{sample_id} ✅ já baixado")
    else:
        print(f"{sample_id} ❌ falta baixar")
