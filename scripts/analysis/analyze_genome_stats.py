#!/usr/bin/env python3
"""
Análisis de Estadísticas Genómicas de E. coli K-12 MG1655
=========================================================
Este script calcula estadísticas genómicas completas incluyendo:
- Contenido GC
- Densidad génica
- Regiones codificantes vs no codificantes
- Distribución de tamaños de genes

Autor: Proyecto Bioinformática
Fecha: 2024
"""

import os
import json
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import pandas as pd
import numpy as np

# ============================================================================
# CONFIGURACIÓN
# ============================================================================

INPUT_FILE = os.path.expanduser('~/data/raw/ecoli_k12_mg1655.gbk')
OUTPUT_DIR = os.path.expanduser('~/projects/bioinfo/results/tables')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# FUNCIONES DE ANÁLISIS
# ============================================================================

def calcular_contenido_gc(secuencia):
    """
    Calcula el contenido GC de una secuencia.
    
    Args:
        secuencia: Objeto Seq de BioPython
    
    Returns:
        float: Porcentaje de GC
    """
    gc_pct = gc_fraction(secuencia) * 100
    return gc_pct

def extraer_genes_completos(record):
    """
    Extrae información completa de todos los genes y CDS.
    
    Args:
        record: Objeto SeqRecord de BioPython
    
    Returns:
        dict: Información detallada de genes
    """
    genes_info = []
    cds_info = []
    
    for feature in record.features:
        if feature.type == 'gene':
            inicio = int(feature.location.start)
            fin = int(feature.location.end)
            longitud = fin - inicio
            
            genes_info.append({
                'tipo': 'gene',
                'inicio': inicio,
                'fin': fin,
                'longitud': longitud,
                'strand': '+' if feature.location.strand == 1 else '-',
                'locus_tag': feature.qualifiers.get('locus_tag', ['NA'])[0]
            })
        
        elif feature.type == 'CDS':
            inicio = int(feature.location.start)
            fin = int(feature.location.end)
            longitud = fin - inicio
            
            # Extraer secuencia del CDS
            cds_seq = feature.extract(record.seq)
            gc_cds = gc_fraction(cds_seq) * 100
            
            cds_info.append({
                'tipo': 'CDS',
                'inicio': inicio,
                'fin': fin,
                'longitud': longitud,
                'strand': '+' if feature.location.strand == 1 else '-',
                'producto': feature.qualifiers.get('product', ['Unknown'])[0],
                'locus_tag': feature.qualifiers.get('locus_tag', ['NA'])[0],
                'gc_content': gc_cds,
                'protein_id': feature.qualifiers.get('protein_id', ['NA'])[0]
            })
    
    return {
        'genes': genes_info,
        'cds': cds_info
    }

def calcular_regiones_codificantes(cds_list, tamano_genoma):
    """
    Calcula el total de bases en regiones codificantes.
    
    Args:
        cds_list: Lista de CDS
        tamano_genoma: Tamaño total del genoma
    
    Returns:
        dict: Estadísticas de regiones codificantes
    """
    total_bp_codificantes = sum(cds['longitud'] for cds in cds_list)
    porcentaje_codificante = (total_bp_codificantes / tamano_genoma) * 100
    porcentaje_no_codificante = 100 - porcentaje_codificante
    
    return {
        'bp_codificantes': total_bp_codificantes,
        'bp_no_codificantes': tamano_genoma - total_bp_codificantes,
        'porcentaje_codificante': porcentaje_codificante,
        'porcentaje_no_codificante': porcentaje_no_codificante
    }

def calcular_estadisticas_tamanos(longitudes):
    """
    Calcula estadísticas descriptivas de tamaños.
    
    Args:
        longitudes: Lista de longitudes
    
    Returns:
        dict: Estadísticas descriptivas
    """
    arr = np.array(longitudes)
    
    return {
        'total': len(arr),
        'media': float(np.mean(arr)),
        'mediana': float(np.median(arr)),
        'desviacion_std': float(np.std(arr)),
        'minimo': int(np.min(arr)),
        'maximo': int(np.max(arr)),
        'percentil_25': float(np.percentile(arr, 25)),
        'percentil_75': float(np.percentile(arr, 75)),
        'rango_intercuartil': float(np.percentile(arr, 75) - np.percentile(arr, 25))
    }

def calcular_densidad_genica(num_genes, tamano_genoma):
    """
    Calcula la densidad génica.
    
    Args:
        num_genes: Número de genes
        tamano_genoma: Tamaño del genoma en bp
    
    Returns:
        dict: Densidad génica en diferentes unidades
    """
    return {
        'genes_por_mb': (num_genes / tamano_genoma) * 1_000_000,
        'genes_por_kb': (num_genes / tamano_genoma) * 1000,
        'bp_por_gen': tamano_genoma / num_genes if num_genes > 0 else 0
    }

def analizar_distribucion_strands(genes_list):
    """
    Analiza la distribución de genes en hebras + y -.
    
    Args:
        genes_list: Lista de genes
    
    Returns:
        dict: Distribución por strand
    """
    strand_plus = sum(1 for g in genes_list if g['strand'] == '+')
    strand_minus = sum(1 for g in genes_list if g['strand'] == '-')
    total = len(genes_list)
    
    return {
        'strand_plus': strand_plus,
        'strand_minus': strand_minus,
        'porcentaje_plus': (strand_plus / total * 100) if total > 0 else 0,
        'porcentaje_minus': (strand_minus / total * 100) if total > 0 else 0
    }

# ============================================================================
# ANÁLISIS PRINCIPAL
# ============================================================================

def main():
    print("="*80)
    print("ANÁLISIS DE ESTADÍSTICAS GENÓMICAS - E. coli K-12 MG1655")
    print("="*80)
    
    # Cargar genoma
    print(f"\n[1/8] Cargando genoma desde: {INPUT_FILE}")
    record = SeqIO.read(INPUT_FILE, 'genbank')
    
    tamano_genoma = len(record.seq)
    
    print(f"✓ Genoma cargado: {record.id}")
    print(f"✓ Descripción: {record.description[:60]}...")
    print(f"✓ Tamaño: {tamano_genoma:,} bp")
    
    # ========================================================================
    # PASO 2.2: Calcular contenido GC
    # ========================================================================
    
    print(f"\n[2/8] Calculando contenido GC del genoma completo...")
    gc_total = calcular_contenido_gc(record.seq)
    
    # Contar nucleótidos individuales
    seq_str = str(record.seq).upper()
    count_a = seq_str.count('A')
    count_t = seq_str.count('T')
    count_g = seq_str.count('G')
    count_c = seq_str.count('C')
    count_n = seq_str.count('N')
    
    print(f"✓ Contenido GC: {gc_total:.2f}%")
    print(f"  → A: {count_a:,} ({count_a/tamano_genoma*100:.2f}%)")
    print(f"  → T: {count_t:,} ({count_t/tamano_genoma*100:.2f}%)")
    print(f"  → G: {count_g:,} ({count_g/tamano_genoma*100:.2f}%)")
    print(f"  → C: {count_c:,} ({count_c/tamano_genoma*100:.2f}%)")
    if count_n > 0:
        print(f"  → N: {count_n:,} ({count_n/tamano_genoma*100:.2f}%)")
    
    # ========================================================================
    # Extraer genes
    # ========================================================================
    
    print(f"\n[3/8] Extrayendo información de genes y CDS...")
    genes_data = extraer_genes_completos(record)
    
    num_genes = len(genes_data['genes'])
    num_cds = len(genes_data['cds'])
    
    print(f"✓ Total genes: {num_genes:,}")
    print(f"✓ Total CDS: {num_cds:,}")
    
    # ========================================================================
    # PASO 2.3: Calcular densidad génica
    # ========================================================================
    
    print(f"\n[4/8] Calculando densidad génica...")
    densidad = calcular_densidad_genica(num_cds, tamano_genoma)
    
    print(f"✓ Genes por Mb: {densidad['genes_por_mb']:.2f}")
    print(f"✓ Genes por kb: {densidad['genes_por_kb']:.4f}")
    print(f"✓ Promedio bp por gen: {densidad['bp_por_gen']:.2f}")
    
    # ========================================================================
    # Regiones codificantes vs no codificantes
    # ========================================================================
    
    print(f"\n[5/8] Analizando regiones codificantes...")
    regiones = calcular_regiones_codificantes(genes_data['cds'], tamano_genoma)
    
    print(f"✓ Bases codificantes: {regiones['bp_codificantes']:,} bp ({regiones['porcentaje_codificante']:.2f}%)")
    print(f"✓ Bases no codificantes: {regiones['bp_no_codificantes']:,} bp ({regiones['porcentaje_no_codificante']:.2f}%)")
    
    # ========================================================================
    # Estadísticas de tamaños de genes
    # ========================================================================
    
    print(f"\n[6/8] Analizando distribución de tamaños de genes...")
    
    longitudes_genes = [g['longitud'] for g in genes_data['genes']]
    longitudes_cds = [c['longitud'] for c in genes_data['cds']]
    
    stats_genes = calcular_estadisticas_tamanos(longitudes_genes)
    stats_cds = calcular_estadisticas_tamanos(longitudes_cds)
    
    print(f"\n  Estadísticas de GENES:")
    print(f"    • Media: {stats_genes['media']:.2f} bp")
    print(f"    • Mediana: {stats_genes['mediana']:.2f} bp")
    print(f"    • Desv. Std: {stats_genes['desviacion_std']:.2f} bp")
    print(f"    • Rango: {stats_genes['minimo']:,} - {stats_genes['maximo']:,} bp")
    
    print(f"\n  Estadísticas de CDS:")
    print(f"    • Media: {stats_cds['media']:.2f} bp")
    print(f"    • Mediana: {stats_cds['mediana']:.2f} bp")
    print(f"    • Desv. Std: {stats_cds['desviacion_std']:.2f} bp")
    print(f"    • Rango: {stats_cds['minimo']:,} - {stats_cds['maximo']:,} bp")
    
    # ========================================================================
    # Análisis de strands
    # ========================================================================
    
    print(f"\n[7/8] Analizando distribución en hebras (strands)...")
    
    strand_genes = analizar_distribucion_strands(genes_data['genes'])
    strand_cds = analizar_distribucion_strands(genes_data['cds'])
    
    print(f"\n  Distribución de GENES:")
    print(f"    • Hebra (+): {strand_genes['strand_plus']:,} ({strand_genes['porcentaje_plus']:.2f}%)")
    print(f"    • Hebra (-): {strand_genes['strand_minus']:,} ({strand_genes['porcentaje_minus']:.2f}%)")
    
    print(f"\n  Distribución de CDS:")
    print(f"    • Hebra (+): {strand_cds['strand_plus']:,} ({strand_cds['porcentaje_plus']:.2f}%)")
    print(f"    • Hebra (-): {strand_cds['strand_minus']:,} ({strand_cds['porcentaje_minus']:.2f}%)")
    
    # ========================================================================
    # GC en CDS
    # ========================================================================
    
    print(f"\n[8/8] Analizando contenido GC en CDS...")
    
    gc_values_cds = [c['gc_content'] for c in genes_data['cds']]
    stats_gc_cds = {
        'media': np.mean(gc_values_cds),
        'mediana': np.median(gc_values_cds),
        'desviacion_std': np.std(gc_values_cds),
        'minimo': np.min(gc_values_cds),
        'maximo': np.max(gc_values_cds)
    }
    
    print(f"✓ GC promedio en CDS: {stats_gc_cds['media']:.2f}%")
    print(f"✓ GC mediana en CDS: {stats_gc_cds['mediana']:.2f}%")
    print(f"✓ Rango GC en CDS: {stats_gc_cds['minimo']:.2f}% - {stats_gc_cds['maximo']:.2f}%")
    
    # ========================================================================
    # GUARDAR RESULTADOS
    # ========================================================================
    
    print("\n" + "="*80)
    print("GUARDANDO RESULTADOS")
    print("="*80)
    
    # JSON con estadísticas completas
    resultado_json = {
        'genoma': {
            'id': record.id,
            'descripcion': record.description,
            'tamano_bp': tamano_genoma,
            'contenido_gc_pct': round(gc_total, 2),
            'composicion_nucleotidos': {
                'A': {'total': count_a, 'porcentaje': round(count_a/tamano_genoma*100, 2)},
                'T': {'total': count_t, 'porcentaje': round(count_t/tamano_genoma*100, 2)},
                'G': {'total': count_g, 'porcentaje': round(count_g/tamano_genoma*100, 2)},
                'C': {'total': count_c, 'porcentaje': round(count_c/tamano_genoma*100, 2)}
            }
        },
        'genes': {
            'total': num_genes,
            'estadisticas_tamanos': stats_genes,
            'distribucion_strands': strand_genes
        },
        'cds': {
            'total': num_cds,
            'estadisticas_tamanos': stats_cds,
            'distribucion_strands': strand_cds,
            'contenido_gc': {
                'media': round(stats_gc_cds['media'], 2),
                'mediana': round(stats_gc_cds['mediana'], 2),
                'desviacion_std': round(stats_gc_cds['desviacion_std'], 2),
                'minimo': round(stats_gc_cds['minimo'], 2),
                'maximo': round(stats_gc_cds['maximo'], 2)
            }
        },
        'densidad_genica': {
            'genes_por_mb': round(densidad['genes_por_mb'], 2),
            'genes_por_kb': round(densidad['genes_por_kb'], 4),
            'bp_por_gen': round(densidad['bp_por_gen'], 2)
        },
        'regiones': {
            'codificantes': {
                'bp': regiones['bp_codificantes'],
                'porcentaje': round(regiones['porcentaje_codificante'], 2)
            },
            'no_codificantes': {
                'bp': regiones['bp_no_codificantes'],
                'porcentaje': round(regiones['porcentaje_no_codificante'], 2)
            }
        }
    }
    
    json_file = os.path.join(OUTPUT_DIR, 'genome_statistics.json')
    with open(json_file, 'w') as f:
        json.dump(resultado_json, f, indent=2)
    print(f"✓ JSON guardado: {json_file}")
    
    # CSV resumen general
    df_resumen = pd.DataFrame([{
        'Genoma_ID': record.id,
        'Tamaño_bp': tamano_genoma,
        'GC_contenido_%': round(gc_total, 2),
        'Total_Genes': num_genes,
        'Total_CDS': num_cds,
        'Genes_por_Mb': round(densidad['genes_por_mb'], 2),
        'Porcentaje_Codificante': round(regiones['porcentaje_codificante'], 2),
        'Tamaño_Medio_Gen_bp': round(stats_genes['media'], 2),
        'Tamaño_Medio_CDS_bp': round(stats_cds['media'], 2)
    }])
    
    csv_file = os.path.join(OUTPUT_DIR, 'genome_statistics_summary.csv')
    df_resumen.to_csv(csv_file, index=False)
    print(f"✓ CSV resumen guardado: {csv_file}")
    
    # CSV con detalles de todos los CDS
    df_cds = pd.DataFrame(genes_data['cds'])
    cds_file = os.path.join(OUTPUT_DIR, 'cds_detailed_info.csv')
    df_cds.to_csv(cds_file, index=False)
    print(f"✓ Detalles CDS guardados: {cds_file}")
    
    # CSV con distribución de tamaños
    df_tamanos = pd.DataFrame({
        'longitud_bp': longitudes_cds
    })
    tamanos_file = os.path.join(OUTPUT_DIR, 'cds_size_distribution.csv')
    df_tamanos.to_csv(tamanos_file, index=False)
    print(f"✓ Distribución tamaños guardada: {tamanos_file}")
    
    # ========================================================================
    # RESUMEN FINAL
    # ========================================================================
    
    print("\n" + "="*80)
    print("RESUMEN DE ESTADÍSTICAS GENÓMICAS")
    print("="*80)
    
    print(f"\n📊 ESTADÍSTICAS GENERALES:")
    print(f"  • Tamaño del genoma:         {tamano_genoma:,} bp")
    print(f"  • Contenido GC:              {gc_total:.2f}%")
    print(f"  • Total genes:               {num_genes:,}")
    print(f"  • Total CDS:                 {num_cds:,}")
    
    print(f"\n🧬 DENSIDAD GÉNICA:")
    print(f"  • Genes por Mb:              {densidad['genes_por_mb']:.2f}")
    print(f"  • Promedio bp por gen:       {densidad['bp_por_gen']:.2f}")
    
    print(f"\n📏 REGIONES GENÓMICAS:")
    print(f"  • Codificante:               {regiones['porcentaje_codificante']:.2f}% ({regiones['bp_codificantes']:,} bp)")
    print(f"  • No codificante:            {regiones['porcentaje_no_codificante']:.2f}% ({regiones['bp_no_codificantes']:,} bp)")
    
    print(f"\n📐 TAMAÑOS DE GENES:")
    print(f"  • Media:                     {stats_genes['media']:.2f} bp")
    print(f"  • Mediana:                   {stats_genes['mediana']:.2f} bp")
    print(f"  • Rango:                     {stats_genes['minimo']:,} - {stats_genes['maximo']:,} bp")
    
    print(f"\n➕➖ DISTRIBUCIÓN EN HEBRAS:")
    print(f"  • Hebra (+):                 {strand_cds['strand_plus']:,} ({strand_cds['porcentaje_plus']:.2f}%)")
    print(f"  • Hebra (-):                 {strand_cds['strand_minus']:,} ({strand_cds['porcentaje_minus']:.2f}%)")
    
    print("\n" + "="*80)
    print("✅ ANÁLISIS DE ESTADÍSTICAS GENÓMICAS COMPLETADO")
    print("="*80 + "\n")

if __name__ == '__main__':
    main()
