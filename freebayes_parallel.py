import os
import subprocess
from glob import glob
from concurrent.futures import ThreadPoolExecutor
from threading import Lock

# Caminhos principais
base_dir = "/home/lab/Desktop/arq_joao/testes_freebayes_bams_marcel"
bam_dir = os.path.join(base_dir, "BAMs_RG")
ref_dir = os.path.join(base_dir, "REF")
bed_dir = os.path.join(base_dir, "BED")
output_dir = os.path.join(base_dir, "VCFs_regions")
temp_bed_dir = os.path.join(base_dir, "temp_beds")

# Criar diretórios
os.makedirs(output_dir, exist_ok=True)
os.makedirs(temp_bed_dir, exist_ok=True)

# Referência e BED
reference = glob(os.path.join(ref_dir, "*.fa"))[0]
bed_file = glob(os.path.join(bed_dir, "*output_fixed.bed"))[0]

# BAMs
bam_files = glob(os.path.join(bam_dir, "*.bam"))
bam_list = " ".join(bam_files)

# Lista de VCFs gerados (thread-safe)
vcf_files = []
vcf_lock = Lock()

def read_bed_lines(bed_file):
    """Lê todas as linhas do arquivo BED"""
    print(f"📋 Lendo arquivo BED: {bed_file}")
    with open(bed_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    print(f"📊 Total de regiões encontradas: {len(lines)}")
    return lines

def process_bed_region(region_data):
    """Processa uma região específica do BED"""
    region_index, bed_line = region_data
    
    # Criar arquivo BED temporário para esta região
    temp_bed = os.path.join(temp_bed_dir, f"region_{region_index:06d}.bed")
    with open(temp_bed, 'w') as f:
        f.write(bed_line + '\n')
    
    # Nome do VCF de saída para esta região
    output_vcf = os.path.join(output_dir, f"region_{region_index:06d}.vcf")
    
    # Extrair informações da região para o log
    parts = bed_line.split('\t')
    if len(parts) >= 3:
        chrom, start, end = parts[0], parts[1], parts[2]
        region_info = f"{chrom}:{start}-{end}"
    else:
        region_info = f"região_{region_index}"
    
    cmd = (
        f"~/freebayes -f {reference} "
        f"-t {temp_bed} "
        f"{bam_list} > {output_vcf}"
    )
    
    print(f"🔄 Processando {region_info} (região {region_index+1})...")
    
    try:
        # Executar FreeBayes
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        
        # Verificar se o VCF foi gerado e não está vazio
        if os.path.exists(output_vcf) and os.path.getsize(output_vcf) > 0:
            with vcf_lock:
                vcf_files.append(output_vcf)
            print(f"✅ {region_info} processado com sucesso")
        else:
            print(f"⚠️ VCF vazio para {region_info}, pulando.")
            
    except subprocess.CalledProcessError as e:
        print(f"❌ Erro ao processar {region_info}: {e}")
        if e.stderr:
            print(f"   Erro: {e.stderr}")
    
    finally:
        # Limpar arquivo BED temporário
        if os.path.exists(temp_bed):
            os.remove(temp_bed)

def clean_temp_files():
    """Remove arquivos temporários"""
    if os.path.exists(temp_bed_dir):
        for file in glob(os.path.join(temp_bed_dir, "*")):
            os.remove(file)
        os.rmdir(temp_bed_dir)

# Ler linhas do BED
bed_lines = read_bed_lines(bed_file)

if not bed_lines:
    print("❌ Nenhuma região válida encontrada no arquivo BED!")
    exit(1)

# Processar regiões em paralelo
print(f"🚀 Iniciando processamento paralelo de {len(bed_lines)} regiões com 4 threads...")

# Criar lista de tuplas (index, bed_line) para passar para as threads
region_data = list(enumerate(bed_lines))

try:
    with ThreadPoolExecutor(max_workers=4) as executor:
        executor.map(process_bed_region, region_data)
    
    print(f"📊 {len(vcf_files)} regiões processadas com sucesso de {len(bed_lines)} total")
    
    if not vcf_files:
        print("❌ Nenhum VCF foi gerado! Verifique os arquivos de entrada e comandos.")
        clean_temp_files()
        exit(1)
    
    # Ordenar VCFs por nome para manter ordem das regiões
    vcf_files.sort()
    
    # Comprimir e indexar todos os VCFs
    print("📦 Comprimindo e indexando VCFs...")
    compressed_vcfs = []
    
    for i, vcf_file in enumerate(vcf_files, 1):
        compressed_vcf = f"{vcf_file}.gz"
        
        print(f"🗜️ Comprimindo {i}/{len(vcf_files)}: {os.path.basename(vcf_file)}...")
        try:
            subprocess.run(f"bgzip -c {vcf_file} > {compressed_vcf}", shell=True, check=True)
            
            print(f"📇 Indexando {os.path.basename(compressed_vcf)}...")
            subprocess.run(f"tabix -p vcf {compressed_vcf}", shell=True, check=True)
            
            compressed_vcfs.append(compressed_vcf)
            
        except subprocess.CalledProcessError as e:
            print(f"❌ Erro ao comprimir/indexar {vcf_file}: {e}")
    
    if not compressed_vcfs:
        print("❌ Nenhum VCF foi comprimido com sucesso!")
        clean_temp_files()
        exit(1)
    
    # Concatenar todos os VCFs comprimidos em um final
    final_vcf = os.path.join(base_dir, "merged_regions.vcf.gz")
    compressed_vcf_input = " ".join(compressed_vcfs)
    
    print("🧬 Concatenando todos os VCFs comprimidos...")
    try:
        subprocess.run(
            f"bcftools concat --threads 14 --allow-overlaps --remove-duplicates -o {final_vcf} -Oz {compressed_vcf_input}",
            shell=True,
            check=True
        )
        
        # Indexar o VCF final
        print("📇 Indexando VCF final...")
        subprocess.run(f"tabix -p vcf {final_vcf}", shell=True, check=True)
        
        print(f"\n✅ Processo concluído! Arquivo final: {final_vcf}")
        print(f"📁 Regiões processadas: {len(vcf_files)}/{len(bed_lines)}")
        print(f"📦 Arquivos comprimidos: {len(compressed_vcfs)}")
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Erro na concatenação final: {e}")

except KeyboardInterrupt:
    print("\n⚠️ Processo interrompido pelo usuário")
    
finally:
    # Limpar arquivos temporários
    clean_temp_files()
    print("🧹 Arquivos temporários removidos")