#!/usr/bin/env python3
import os
import glob
import subprocess
from pathlib import Path

# Caminhos
INPUT_DIR = "/home/lab/Desktop/arq_joao/testes_freebayes_bams_marcel/BAMs"
OUTPUT_DIR = "/home/lab/Desktop/arq_joao/testes_freebayes_bams_marcel/BAMs_RG"

# Criar diretório de saída se não existir
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

# Buscar todos os arquivos BAM que terminam com PIG706_sorted.bam
bam_files = glob.glob(os.path.join(INPUT_DIR, "*.bam"))

# Loop em todos os arquivos BAM
for bam in bam_files:
    filename = os.path.basename(bam)
    
    # Extrair amostra do nome do arquivo (ex: HG1 → hg1)
    amostra = filename.split('-')[0].lower()
    
    # Definir nome de saída
    output_bam = os.path.join(OUTPUT_DIR, filename.replace('.bam', '_rg.bam'))
    
    # Comando para adicionar read group
    addreplacerg_cmd = [
        'samtools', 'addreplacerg',
        '-r', f'ID:{amostra}',
        '-r', 'LB:lib1',
        '-r', 'PL:ONT',
        '-r', 'PU:unit1',
        '-r', f'SM:{amostra}',
        '-o', output_bam,
        bam
    ]
    
    print(f"Processando: {filename}")
    print(f"Amostra: {amostra}")
    print(f"Saída: {os.path.basename(output_bam)}")
    
    try:
        # Executar comando samtools addreplacerg
        subprocess.run(addreplacerg_cmd, check=True)
        print(f"✓ Read group adicionado com sucesso")
        
        # Indexar o arquivo BAM
        index_cmd = ['samtools', 'index', output_bam]
        subprocess.run(index_cmd, check=True)
        print(f"✓ Arquivo indexado com sucesso")
        
    except subprocess.CalledProcessError as e:
        print(f"✗ Erro ao processar {filename}: {e}")
    
    print("-" * 50)

print("Processamento concluído!")