#!/usr/bin/env python3
"""
Análisis de Distribución de Tamaños de Genes - E. coli K-12 MG1655
===================================================================
Este script realiza un análisis estadístico detallado de la distribución
de tamaños de genes, incluyendo:
- Estadísticas descriptivas completas
- Categorización por tamaño
- Análisis de outliers
- Genes extremos (más pequeños y más grandes)

Autor: Proyecto Bioinformática
Fecha: 2024
"""

import os
import json
from Bio import SeqIO
import pandas as pd
import numpy as np
from collections import Counter

# ============================================================================
# CONFIGURACIÓN
# ============================================================================

INPUT_FILE = os.path.expanduser('~/data/raw/ecoli_k12_mg1655.gbk')
RESULTS_DIR = os.path.expanduser('~/projects/bioinfo/results/tables')
OUTPUT_DIR = os.path.expanduser('~/projects/bioinfo/results/tables')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# FUNCIONES DE ANÁLISIS
# ============================================================================

def extraer_info_genes_completa(record):
    """
    Extrae información completa de todos los CDS.
    
    Args:
        record: Objeto SeqRecord de BioPython
    
    Returns:
        list: Lista de diccionarios con información de cada CDS
    """
    cds_info = []
    
    for feature in record.features:
        if feature.type == 'CDS':
            inicio = int(feature.location.start)
            fin = int(feature.location.end)
            longitud = fin - inicio
            
            # Calcular longitud de proteína
            longitud_proteina = longitud // 3
            
            cds_info.append({
                'locus_tag': feature.qualifiers.get('locus_tag', ['NA'])[0],
                'gene': feature.qualifiers.get('gene', ['NA'])[0],
                'producto': feature.qualifiers.get('product', ['Unknown'])[0],
                'inicio': inicio,
                'fin': fin,
                'longitud_nt': longitud,
                'longitud_aa': longitud_proteina,
                'strand': '+' if feature.location.strand == 1 else '-',
                'protein_id': feature.qualifiers.get('protein_id', ['NA'])[0]
            })
    
    return cds_info

def calcular_estadisticas_avanzadas(longitudes):
    """
    Calcula estadísticas descriptivas avanzadas.
    
    Args:
        longitudes: Array de longitudes
    
    Returns:
        dict: Estadísticas completas
    """
    arr = np.array(longitudes)
    
    return {
        'total': len(arr),
        'media': float(np.mean(arr)),
        'mediana': float(np.median(arr)),
        'moda': float(Counter(arr).most_common(1)[0][0]) if len(arr) > 0 else 0,
        'desviacion_std': float(np.std(arr)),
        'varianza': float(np.var(arr)),
        'minimo': int(np.min(arr)),
        'maximo': int(np.max(arr)),
        'rango': int(np.max(arr) - np.min(arr)),
        'percentil_5': float(np.percentile(arr, 5)),
        'percentil_10': float(np.percentile(arr, 10)),
        'percentil_25': float(np.percentile(arr, 25)),
        'percentil_50': float(np.percentile(arr, 50)),
        'percentil_75': float(np.percentile(arr, 75)),
        'percentil_90': float(np.percentile(arr, 90)),
        'percentil_95': float(np.percentile(arr, 95)),
        'rango_intercuartil': float(np.percentile(arr, 75) - np.percentile(arr, 25)),
        'coef_variacion': float(np.std(arr) / np.mean(arr) * 100) if np.mean(arr) > 0 else 0
    }

def categorizar_por_tamano(cds_list):
    """
    Categoriza genes por tamaño.
    
    Args:
        cds_list: Lista de CDS
    
    Returns:
        dict: Genes categorizados
    """
    categorias = {
        'muy_pequenos': [],      # < 300 bp (< 100 aa)
        'pequenos': [],          # 300-600 bp (100-200 aa)
        'medianos': [],          # 600-1200 bp (200-400 aa)
        'grandes': [],           # 1200-2400 bp (400-800 aa)
        'muy_grandes': []        # > 2400 bp (> 800 aa)
    }
    
    for cds in cds_list:
        longitud = cds['longitud_nt']
        
        if longitud < 300:
            categorias['muy_pequenos'].append(cds)
        elif longitud < 600:
            categorias['pequenos'].append(cds)
        elif longitud < 1200:
            categorias['medianos'].append(cds)
        elif longitud < 2400:
            categorias['grandes'].append(cds)
        else:
            categorias['muy_grandes'].append(cds)
    
    return categorias

def identificar_outliers(longitudes):
    """
    Identifica outliers usando el método IQR.
    
    Args:
        longitudes: Array de longitudes
    
    Returns:
        dict: Información de outliers
    """
    arr = np.array(longitudes)
    q1 = np.percentile(arr, 25)
    q3 = np.percentile(arr, 75)
    iqr = q3 - q1
    
    limite_inferior = q1 - 1.5 * iqr
    limite_superior = q3 + 1.5 * iqr
    
    outliers_bajos = arr[arr < limite_inferior]
    outliers_altos = arr[arr > limite_superior]
    
    return {
        'limite_inferior': float(limite_inferior),
        'limite_superior': float(limite_superior),
        'iqr': float(iqr),
        'num_outliers_bajos': len(outliers_bajos),
        'num_outliers_altos': len(outliers_altos),
        'total_outliers': len(outliers_bajos) + len(outliers_altos),
        'porcentaje_outliers': ((len(outliers_bajos) + len(outliers_altos)) / len(arr) * 100) if len(arr) > 0 else 0
    }

def analizar_distribucion_por_multiplos_3(longitudes):
    """
    Analiza si las longitudes son múltiplos de 3 (codones completos).
    
    Args:
        longitudes: Lista de longitudes
    
    Returns:
        dict: Análisis de múltiplos de 3
    """
    multiplos_3 = sum(1 for l in longitudes if l % 3 == 0)
    no_multiplos_3 = len(longitudes) - multiplos_3
    
    return {
        'multiplos_de_3': multiplos_3,
        'no_multiplos_de_3': no_multiplos_3,
        'porcentaje_multiplos_3': (multiplos_3 / len(longitudes) * 100) if len(longitudes) > 0 else 0,
        'interpretacion': 'Todos los CDS deben ser múltiplos de 3 (codones completos)'
    }

# ============================================================================
# ANÁLISIS PRINCIPAL
# ============================================================================

def main():
    print("="*80)
    print("ANÁLISIS DE DISTRIBUCIÓN DE TAMAÑOS DE GENES - E. coli K-12 MG1655")
    print("="*80)
    
    # Cargar genoma
    print(f"\n[1/7] Cargando genoma desde: {INPUT_FILE}")
    record = SeqIO.read(INPUT_FILE, 'genbank')
    print(f"✓ Genoma cargado: {record.id}")
    
    # Extraer información de CDS
    print(f"\n[2/7] Extrayendo información completa de CDS...")
    cds_list = extraer_info_genes_completa(record)
    print(f"✓ Total CDS extraídos: {len(cds_list):,}")
    
    # Extraer longitudes
    longitudes_nt = [cds['longitud_nt'] for cds in cds_list]
    longitudes_aa = [cds['longitud_aa'] for cds in cds_list]
    
    # ========================================================================
    # PASO 3.2: Estadísticas Descriptivas
    # ========================================================================
    
    print(f"\n[3/7] Calculando estadísticas descriptivas...")
    stats_nt = calcular_estadisticas_avanzadas(longitudes_nt)
    stats_aa = calcular_estadisticas_avanzadas(longitudes_aa)
    
    print(f"\n  📊 Estadísticas en nucleótidos (nt):")
    print(f"    • Total:             {stats_nt['total']:,} genes")
    print(f"    • Media:             {stats_nt['media']:.2f} bp")
    print(f"    • Mediana:           {stats_nt['mediana']:.2f} bp")
    print(f"    • Moda:              {stats_nt['moda']:.0f} bp")
    print(f"    • Desv. Std:         {stats_nt['desviacion_std']:.2f} bp")
    print(f"    • Rango:             {stats_nt['minimo']:,} - {stats_nt['maximo']:,} bp")
    print(f"    • Coef. Variación:   {stats_nt['coef_variacion']:.2f}%")
    
    print(f"\n  🧬 Estadísticas en aminoácidos (aa):")
    print(f"    • Media:             {stats_aa['media']:.2f} aa")
    print(f"    • Mediana:           {stats_aa['mediana']:.2f} aa")
    print(f"    • Rango:             {stats_aa['minimo']:,} - {stats_aa['maximo']:,} aa")
    
    # ========================================================================
    # Percentiles detallados
    # ========================================================================
    
    print(f"\n[4/7] Análisis de percentiles...")
    print(f"\n  Percentiles (nucleótidos):")
    print(f"    •  5%:  {stats_nt['percentil_5']:.0f} bp")
    print(f"    • 10%:  {stats_nt['percentil_10']:.0f} bp")
    print(f"    • 25%:  {stats_nt['percentil_25']:.0f} bp")
    print(f"    • 50%:  {stats_nt['percentil_50']:.0f} bp (mediana)")
    print(f"    • 75%:  {stats_nt['percentil_75']:.0f} bp")
    print(f"    • 90%:  {stats_nt['percentil_90']:.0f} bp")
    print(f"    • 95%:  {stats_nt['percentil_95']:.0f} bp")
    print(f"    • IQR:  {stats_nt['rango_intercuartil']:.0f} bp")
    
    # ========================================================================
    # Categorización por tamaño
    # ========================================================================
    
    print(f"\n[5/7] Categorizando genes por tamaño...")
    categorias = categorizar_por_tamano(cds_list)
    
    print(f"\n  Categorías de tamaño:")
    print(f"    • Muy pequeños (<300 bp):     {len(categorias['muy_pequenos']):4,} ({len(categorias['muy_pequenos'])/len(cds_list)*100:5.2f}%)")
    print(f"    • Pequeños (300-600 bp):      {len(categorias['pequenos']):4,} ({len(categorias['pequenos'])/len(cds_list)*100:5.2f}%)")
    print(f"    • Medianos (600-1200 bp):     {len(categorias['medianos']):4,} ({len(categorias['medianos'])/len(cds_list)*100:5.2f}%)")
    print(f"    • Grandes (1200-2400 bp):     {len(categorias['grandes']):4,} ({len(categorias['grandes'])/len(cds_list)*100:5.2f}%)")
    print(f"    • Muy grandes (>2400 bp):     {len(categorias['muy_grandes']):4,} ({len(categorias['muy_grandes'])/len(cds_list)*100:5.2f}%)")
    
    # ========================================================================
    # Identificar outliers
    # ========================================================================
    
    print(f"\n[6/7] Identificando outliers...")
    outliers_info = identificar_outliers(longitudes_nt)
    
    print(f"\n  Análisis de outliers (método IQR):")
    print(f"    • Límite inferior:   {outliers_info['limite_inferior']:.0f} bp")
    print(f"    • Límite superior:   {outliers_info['limite_superior']:.0f} bp")
    print(f"    • Outliers bajos:    {outliers_info['num_outliers_bajos']:,}")
    print(f"    • Outliers altos:    {outliers_info['num_outliers_altos']:,}")
    print(f"    • Total outliers:    {outliers_info['total_outliers']:,} ({outliers_info['porcentaje_outliers']:.2f}%)")
    
    # ========================================================================
    # Genes extremos
    # ========================================================================
    
    print(f"\n[7/7] Identificando genes extremos...")
    
    # Ordenar por tamaño
    cds_sorted = sorted(cds_list, key=lambda x: x['longitud_nt'])
    
    genes_mas_pequenos = cds_sorted[:10]
    genes_mas_grandes = cds_sorted[-10:]
    
    print(f"\n  🔬 10 GENES MÁS PEQUEÑOS:")
    for i, gene in enumerate(genes_mas_pequenos, 1):
        print(f"    {i:2}. {gene['locus_tag']:12} | {gene['longitud_nt']:4} bp ({gene['longitud_aa']:3} aa) | {gene['producto'][:40]}")
    
    print(f"\n  🔬 10 GENES MÁS GRANDES:")
    for i, gene in enumerate(genes_mas_grandes, 1):
        print(f"    {i:2}. {gene['locus_tag']:12} | {gene['longitud_nt']:5} bp ({gene['longitud_aa']:4} aa) | {gene['producto'][:40]}")
    
    # ========================================================================
    # Análisis de múltiplos de 3
    # ========================================================================
    
    multiplos_3 = analizar_distribucion_por_multiplos_3(longitudes_nt)
    
    print(f"\n  ✓ Validación: {multiplos_3['multiplos_de_3']:,}/{len(longitudes_nt):,} genes son múltiplos de 3 ({multiplos_3['porcentaje_multiplos_3']:.2f}%)")
    
    # ========================================================================
    # GUARDAR RESULTADOS
    # ========================================================================
    
    print("\n" + "="*80)
    print("GUARDANDO RESULTADOS")
    print("="*80)
    
    # JSON con análisis completo
    resultado_json = {
        'genoma_id': record.id,
        'total_cds': len(cds_list),
        'estadisticas_nucleotidos': stats_nt,
        'estadisticas_aminoacidos': stats_aa,
        'categorias_tamano': {
            'muy_pequenos': {
                'total': len(categorias['muy_pequenos']),
                'porcentaje': round(len(categorias['muy_pequenos'])/len(cds_list)*100, 2),
                'rango': '< 300 bp'
            },
            'pequenos': {
                'total': len(categorias['pequenos']),
                'porcentaje': round(len(categorias['pequenos'])/len(cds_list)*100, 2),
                'rango': '300-600 bp'
            },
            'medianos': {
                'total': len(categorias['medianos']),
                'porcentaje': round(len(categorias['medianos'])/len(cds_list)*100, 2),
                'rango': '600-1200 bp'
            },
            'grandes': {
                'total': len(categorias['grandes']),
                'porcentaje': round(len(categorias['grandes'])/len(cds_list)*100, 2),
                'rango': '1200-2400 bp'
            },
            'muy_grandes': {
                'total': len(categorias['muy_grandes']),
                'porcentaje': round(len(categorias['muy_grandes'])/len(cds_list)*100, 2),
                'rango': '> 2400 bp'
            }
        },
        'outliers': outliers_info,
        'multiplos_de_3': multiplos_3,
        'genes_extremos': {
            'mas_pequenos': genes_mas_pequenos[:10],
            'mas_grandes': genes_mas_grandes[-10:]
        }
    }
    
    json_file = os.path.join(OUTPUT_DIR, 'gene_size_distribution_analysis.json')
    with open(json_file, 'w') as f:
        json.dump(resultado_json, f, indent=2)
    print(f"✓ JSON guardado: {json_file}")
    
    # CSV con todos los genes y sus tamaños
    df_genes = pd.DataFrame(cds_list)
    genes_file = os.path.join(OUTPUT_DIR, 'all_genes_with_sizes.csv')
    df_genes.to_csv(genes_file, index=False)
    print(f"✓ Todos los genes guardados: {genes_file}")
    
    # CSV con estadísticas por categoría
    df_categorias = pd.DataFrame([
        {
            'Categoria': 'Muy pequeños',
            'Rango': '< 300 bp',
            'Total': len(categorias['muy_pequenos']),
            'Porcentaje': round(len(categorias['muy_pequenos'])/len(cds_list)*100, 2)
        },
        {
            'Categoria': 'Pequeños',
            'Rango': '300-600 bp',
            'Total': len(categorias['pequenos']),
            'Porcentaje': round(len(categorias['pequenos'])/len(cds_list)*100, 2)
        },
        {
            'Categoria': 'Medianos',
            'Rango': '600-1200 bp',
            'Total': len(categorias['medianos']),
            'Porcentaje': round(len(categorias['medianos'])/len(cds_list)*100, 2)
        },
        {
            'Categoria': 'Grandes',
            'Rango': '1200-2400 bp',
            'Total': len(categorias['grandes']),
            'Porcentaje': round(len(categorias['grandes'])/len(cds_list)*100, 2)
        },
        {
            'Categoria': 'Muy grandes',
            'Rango': '> 2400 bp',
            'Total': len(categorias['muy_grandes']),
            'Porcentaje': round(len(categorias['muy_grandes'])/len(cds_list)*100, 2)
        }
    ])
    
    cat_file = os.path.join(OUTPUT_DIR, 'gene_size_categories.csv')
    df_categorias.to_csv(cat_file, index=False)
    print(f"✓ Categorías guardadas: {cat_file}")
    
    # CSV con genes extremos
    df_extremos = pd.concat([
        pd.DataFrame(genes_mas_pequenos).assign(tipo='Más pequeño'),
        pd.DataFrame(genes_mas_grandes).assign(tipo='Más grande')
    ])
    extremos_file = os.path.join(OUTPUT_DIR, 'extreme_genes.csv')
    df_extremos.to_csv(extremos_file, index=False)
    print(f"✓ Genes extremos guardados: {extremos_file}")
    
    # ========================================================================
    # RESUMEN FINAL
    # ========================================================================
    
    print("\n" + "="*80)
    print("RESUMEN DE DISTRIBUCIÓN DE TAMAÑOS")
    print("="*80)
    
    print(f"\n📊 ESTADÍSTICAS PRINCIPALES:")
    print(f"  • Total genes analizados:    {len(cds_list):,}")
    print(f"  • Tamaño medio:              {stats_nt['media']:.2f} bp ({stats_aa['media']:.0f} aa)")
    print(f"  • Tamaño mediano:            {stats_nt['mediana']:.0f} bp ({stats_aa['mediana']:.0f} aa)")
    print(f"  • Rango:                     {stats_nt['minimo']:,} - {stats_nt['maximo']:,} bp")
    
    print(f"\n📏 DISTRIBUCIÓN POR CATEGORÍAS:")
    print(f"  • La mayoría de genes son de tamaño mediano (600-1200 bp)")
    print(f"  • {len(categorias['medianos']):,} genes ({len(categorias['medianos'])/len(cds_list)*100:.1f}%) en categoría mediana")
    
    print(f"\n🔍 OUTLIERS:")
    print(f"  • {outliers_info['total_outliers']:,} genes outliers ({outliers_info['porcentaje_outliers']:.2f}%)")
    
    print("\n" + "="*80)
    print("✅ ANÁLISIS DE DISTRIBUCIÓN DE GENES COMPLETADO")
    print("="*80 + "\n")

if __name__ == '__main__':
    main()
