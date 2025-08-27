import os
import subprocess
from glob import glob

# Caminhos
vcf_dir = "/home/lab/Desktop/arq_joao/testes_freebayes_bams_marcel/VCF_MERGED"         # pasta com seus VCFs
output_dir = "/home/lab/Desktop/arq_joao/testes_freebayes_bams_marcel/VCF_QC"  # pasta de saída
ref_fa = "/home/lab/Desktop/arq_joao/testes_freebayes_bams_marcel/REF/GRCh38_full_analysis_set_plus_decoy_hla.fa"         # referência hg38.fa

os.makedirs(output_dir, exist_ok=True)

# Procurar arquivos .vcf na pasta
vcf_files = glob(os.path.join(vcf_dir, "*.vcf"))

for vcf in vcf_files:
    sample_name = os.path.basename(vcf).replace(".vcf", "")
    output_file = os.path.join(output_dir, f"{sample_name}_filtrado_norm.vcf")

    print(f"[INFO] Processando {vcf} -> {output_file}")

    # Pipeline: bcftools view -> bcftools norm
    view_cmd = ["bcftools", "view", "--exclude", "QUAL<1", vcf]
    norm_cmd = ["bcftools", "norm", "-f", ref_fa]

    with open(output_file, "w") as out:
        p1 = subprocess.Popen(view_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(norm_cmd, stdin=p1.stdout, stdout=out)
        p1.stdout.close()  # permite p1 receber SIGPIPE se p2 morrer
        p2.communicate()

    print(f"[OK] Resultado salvo em {output_file}")
