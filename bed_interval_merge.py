#!/usr/bin/env python3
"""
Script para otimizar arquivo BED expandindo regiões e juntando intervalos próximos.
Adiciona 500pb antes do start e 500pb depois do end, e junta regiões próximas (≤500pb).
"""

import sys
import argparse
from typing import List, Tuple

def read_bed_file(filename: str) -> List[Tuple[str, int, int]]:
    """
    Lê arquivo BED e retorna lista de tuplas (chr, start, end).
    """
    intervals = []
    try:
        with open(filename, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 3:
                    print(f"Aviso: Linha {line_num} ignorada (formato inválido): {line}")
                    continue
                
                try:
                    chr_name = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    
                    if start >= end:
                        print(f"Aviso: Linha {line_num} ignorada (start >= end): {line}")
                        continue
                    
                    intervals.append((chr_name, start, end))
                    
                except ValueError:
                    print(f"Aviso: Linha {line_num} ignorada (coordenadas inválidas): {line}")
                    continue
                    
    except FileNotFoundError:
        print(f"Erro: Arquivo '{filename}' não encontrado.")
        sys.exit(1)
    except Exception as e:
        print(f"Erro ao ler arquivo: {e}")
        sys.exit(1)
    
    return intervals

def expand_intervals(intervals: List[Tuple[str, int, int]], expansion: int = 500) -> List[Tuple[str, int, int]]:
    """
    Expande cada intervalo adicionando 'expansion' pb antes do start e depois do end.
    """
    expanded = []
    for chr_name, start, end in intervals:
        new_start = max(0, start - expansion)  # Não permite coordenadas negativas
        new_end = end + expansion
        expanded.append((chr_name, new_start, new_end))
    
    return expanded

def merge_overlapping_intervals(intervals: List[Tuple[str, int, int]], max_gap: int = 500) -> List[Tuple[str, int, int]]:
    """
    Junta intervalos sobrepostos ou próximos (gap ≤ max_gap) no mesmo cromossomo.
    """
    if not intervals:
        return []
    
    # Agrupa por cromossomo
    chr_groups = {}
    for chr_name, start, end in intervals:
        if chr_name not in chr_groups:
            chr_groups[chr_name] = []
        chr_groups[chr_name].append((start, end))
    
    merged_intervals = []
    
    # Processa cada cromossomo separadamente
    for chr_name in sorted(chr_groups.keys()):
        chr_intervals = chr_groups[chr_name]
        
        # Ordena por posição de início
        chr_intervals.sort(key=lambda x: x[0])
        
        # Junta intervalos sobrepostos ou próximos
        merged_chr = []
        current_start, current_end = chr_intervals[0]
        
        for start, end in chr_intervals[1:]:
            # Se o intervalo atual está próximo ou sobreposto ao anterior
            if start <= current_end + max_gap:
                # Estende o intervalo atual
                current_end = max(current_end, end)
            else:
                # Adiciona o intervalo atual e inicia um novo
                merged_chr.append((chr_name, current_start, current_end))
                current_start, current_end = start, end
        
        # Adiciona o último intervalo
        merged_chr.append((chr_name, current_start, current_end))
        merged_intervals.extend(merged_chr)
    
    return merged_intervals

def write_bed_file(intervals: List[Tuple[str, int, int]], filename: str):
    """
    Escreve intervalos no formato BED.
    """
    try:
        with open(filename, 'w') as f:
            for chr_name, start, end in intervals:
                f.write(f"{chr_name}\t{start}\t{end}\n")
        print(f"Arquivo otimizado salvo: {filename}")
    except Exception as e:
        print(f"Erro ao escrever arquivo: {e}")
        sys.exit(1)

def print_statistics(original: List[Tuple[str, int, int]], 
                    optimized: List[Tuple[str, int, int]]):
    """
    Imprime estatísticas da otimização.
    """
    print("\n=== Estatísticas da Otimização ===")
    print(f"Intervalos originais: {len(original)}")
    print(f"Intervalos otimizados: {len(optimized)}")
    print(f"Redução: {len(original) - len(optimized)} intervalos ({((len(original) - len(optimized)) / len(original) * 100):.1f}%)")
    
    # Calcula cobertura total
    original_coverage = sum(end - start for _, start, end in original)
    optimized_coverage = sum(end - start for _, start, end in optimized)
    
    print(f"Cobertura original: {original_coverage:,} pb")
    print(f"Cobertura otimizada: {optimized_coverage:,} pb")
    print(f"Expansão da cobertura: {optimized_coverage - original_coverage:,} pb ({((optimized_coverage / original_coverage - 1) * 100):.1f}%)")

def main():
    parser = argparse.ArgumentParser(
        description="Otimiza arquivo BED expandindo regiões e juntando intervalos próximos",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos de uso:
  python bed_optimizer.py input.bed -o output.bed
  python bed_optimizer.py markers.bed -e 1000 -g 200 -o optimized.bed
  python bed_optimizer.py input.bed --expansion 500 --max-gap 500
        """
    )
    
    parser.add_argument('input_bed', 
                       help='Arquivo BED de entrada')
    parser.add_argument('-o', '--output', 
                       help='Arquivo BED de saída (default: input_optimized.bed)')
    parser.add_argument('-e', '--expansion', 
                       type=int, default=500,
                       help='Bases para expandir antes/depois (default: 500)')
    parser.add_argument('-g', '--max-gap', 
                       type=int, default=500,
                       help='Distância máxima para juntar intervalos (default: 500)')
    parser.add_argument('--stats', 
                       action='store_true',
                       help='Mostrar estatísticas da otimização')
    
    args = parser.parse_args()
    
    # Define arquivo de saída se não especificado
    if not args.output:
        if args.input_bed.endswith('.bed'):
            args.output = args.input_bed.replace('.bed', '_optimized.bed')
        else:
            args.output = args.input_bed + '_optimized.bed'
    
    print(f"Processando arquivo: {args.input_bed}")
    print(f"Expansão: {args.expansion} pb")
    print(f"Gap máximo para junção: {args.max_gap} pb")
    
    # Lê arquivo BED original
    original_intervals = read_bed_file(args.input_bed)
    
    if not original_intervals:
        print("Erro: Nenhum intervalo válido encontrado no arquivo.")
        sys.exit(1)
    
    print(f"Intervalos carregados: {len(original_intervals)}")
    
    # Expande intervalos
    expanded_intervals = expand_intervals(original_intervals, args.expansion)
    
    # Junta intervalos próximos/sobrepostos
    optimized_intervals = merge_overlapping_intervals(expanded_intervals, args.max_gap)
    
    # Escreve arquivo otimizado
    write_bed_file(optimized_intervals, args.output)
    
    # Mostra estatísticas se solicitado
    if args.stats:
        print_statistics(original_intervals, optimized_intervals)
    
    print("Otimização concluída!")

if __name__ == "__main__":
    main()