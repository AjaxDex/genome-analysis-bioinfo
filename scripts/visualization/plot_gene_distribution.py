#!/usr/bin/env python3
"""
Visualización de Distribución de Tamaños de Genes
=================================================
Genera histogramas y gráficos de distribución de tamaños de genes.

Autor: Proyecto Bioinformática
Fecha: 2024
"""

import os
import json
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# ============================================================================
# CONFIGURACIÓN
# ============================================================================

RESULTS_DIR = os.path.expanduser('~/projects/bioinfo/results/tables')
OUTPUT_DIR = os.path.expanduser('~/projects/bioinfo/results/figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Configuración de estilo
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 11

# ============================================================================
# FUNCIONES DE VISUALIZACIÓN
# ============================================================================

def plot_gene_size_histogram():
    """Histograma de distribución de tamaños de genes."""
    
    # Cargar datos
    df = pd.read_csv(os.path.join(RESULTS_DIR, 'cds_size_distribution.csv'))
    longitudes = df['longitud_bp'].values
    
    with open(os.path.join(RESULTS_DIR, 'gene_size_distribution_analysis.json'), 'r') as f:
        stats = json.load(f)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # ========================================================================
    # Subplot 1: Histograma principal
    # ========================================================================
    
    axes[0, 0].hist(longitudes, bins=50, color='#3498db', alpha=0.7, edgecolor='black')
    axes[0, 0].axvline(stats['estadisticas_nucleotidos']['media'], 
                       color='red', linestyle='--', linewidth=2, label=f'Media: {stats["estadisticas_nucleotidos"]["media"]:.0f} bp')
    axes[0, 0].axvline(stats['estadisticas_nucleotidos']['mediana'], 
                       color='green', linestyle='--', linewidth=2, label=f'Mediana: {stats["estadisticas_nucleotidos"]["mediana"]:.0f} bp')
    
    axes[0, 0].set_xlabel('Longitud del gen (bp)', fontsize=12, fontweight='bold')
    axes[0, 0].set_ylabel('Frecuencia', fontsize=12, fontweight='bold')
    axes[0, 0].set_title('Distribución de Tamaños de Genes (Todos los genes)',
                         fontsize=14, fontweight='bold', pad=15)
    axes[0, 0].legend(fontsize=11)
    axes[0, 0].grid(axis='y', alpha=0.3)
    
    # ========================================================================
    # Subplot 2: Histograma sin outliers (mejor visualización)
    # ========================================================================
    
    # Filtrar outliers extremos para mejor visualización
    q1 = stats['estadisticas_nucleotidos']['percentil_25']
    q3 = stats['estadisticas_nucleotidos']['percentil_75']
    iqr = q3 - q1
    limite_superior = q3 + 3 * iqr
    
    longitudes_filtradas = longitudes[longitudes <= limite_superior]
    
    axes[0, 1].hist(longitudes_filtradas, bins=50, color='#2ecc71', alpha=0.7, edgecolor='black')
    axes[0, 1].axvline(stats['estadisticas_nucleotidos']['media'], 
                       color='red', linestyle='--', linewidth=2, label=f'Media: {stats["estadisticas_nucleotidos"]["media"]:.0f} bp')
    axes[0, 1].axvline(stats['estadisticas_nucleotidos']['mediana'], 
                       color='darkgreen', linestyle='--', linewidth=2, label=f'Mediana: {stats["estadisticas_nucleotidos"]["mediana"]:.0f} bp')
    
    axes[0, 1].set_xlabel('Longitud del gen (bp)', fontsize=12, fontweight='bold')
    axes[0, 1].set_ylabel('Frecuencia', fontsize=12, fontweight='bold')
    axes[0, 1].set_title('Distribución de Tamaños (Sin outliers extremos)',
                         fontsize=14, fontweight='bold', pad=15)
    axes[0, 1].legend(fontsize=11)
    axes[0, 1].grid(axis='y', alpha=0.3)
    
    # ========================================================================
    # Subplot 3: Box plot
    # ========================================================================
    
    bp = axes[1, 0].boxplot([longitudes], vert=False, widths=0.5, patch_artist=True,
                             boxprops=dict(facecolor='#3498db', alpha=0.7),
                             medianprops=dict(color='red', linewidth=2),
                             whiskerprops=dict(linewidth=1.5),
                             capprops=dict(linewidth=1.5))
    
    axes[1, 0].set_xlabel('Longitud del gen (bp)', fontsize=12, fontweight='bold')
    axes[1, 0].set_title('Box Plot de Tamaños de Genes',
                         fontsize=14, fontweight='bold', pad=15)
    axes[1, 0].grid(axis='x', alpha=0.3)
    axes[1, 0].set_yticks([])
    
    # Añadir estadísticas
    stats_text = f"""Estadísticas:
Min: {stats['estadisticas_nucleotidos']['minimo']:,} bp
Q1: {stats['estadisticas_nucleotidos']['percentil_25']:.0f} bp
Mediana: {stats['estadisticas_nucleotidos']['mediana']:.0f} bp
Q3: {stats['estadisticas_nucleotidos']['percentil_75']:.0f} bp
Max: {stats['estadisticas_nucleotidos']['maximo']:,} bp
IQR: {stats['estadisticas_nucleotidos']['rango_intercuartil']:.0f} bp"""
    
    axes[1, 0].text(0.98, 0.97, stats_text,
                    transform=axes[1, 0].transAxes,
                    fontsize=10,
                    verticalalignment='top',
                    horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # ========================================================================
    # Subplot 4: Distribución por categorías
    # ========================================================================
    
    categorias_data = stats['categorias_tamano']
    categorias = ['Muy\npequeños\n(<300)', 'Pequeños\n(300-600)', 
                  'Medianos\n(600-1200)', 'Grandes\n(1200-2400)', 
                  'Muy\ngrandes\n(>2400)']
    valores = [
        categorias_data['muy_pequenos']['total'],
        categorias_data['pequenos']['total'],
        categorias_data['medianos']['total'],
        categorias_data['grandes']['total'],
        categorias_data['muy_grandes']['total']
    ]
    porcentajes = [
        categorias_data['muy_pequenos']['porcentaje'],
        categorias_data['pequenos']['porcentaje'],
        categorias_data['medianos']['porcentaje'],
        categorias_data['grandes']['porcentaje'],
        categorias_data['muy_grandes']['porcentaje']
    ]
    
    colors = ['#e74c3c', '#f39c12', '#2ecc71', '#3498db', '#9b59b6']
    bars = axes[1, 1].bar(categorias, valores, color=colors, alpha=0.8, edgecolor='black')
    
    axes[1, 1].set_ylabel('Número de genes', fontsize=12, fontweight='bold')
    axes[1, 1].set_title('Distribución por Categorías de Tamaño',
                         fontsize=14, fontweight='bold', pad=15)
    axes[1, 1].grid(axis='y', alpha=0.3)
    
    # Añadir valores y porcentajes
    for bar, val, pct in zip(bars, valores, porcentajes):
        height = bar.get_height()
        axes[1, 1].text(bar.get_x() + bar.get_width()/2., height,
                       f'{val:,}\n({pct:.1f}%)',
                       ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # Título general
    fig.suptitle('Análisis de Distribución de Tamaños de Genes - E. coli K-12 MG1655',
                 fontsize=16, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    
    # Guardar
    output_file = os.path.join(OUTPUT_DIR, 'gene_size_distribution.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Gráfico guardado: {output_file}")
    plt.close()

def plot_gene_size_violin():
    """Violin plot de tamaños de genes."""
    
    # Cargar datos
    df = pd.read_csv(os.path.join(RESULTS_DIR, 'all_genes_with_sizes.csv'))
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # ========================================================================
    # Subplot 1: Violin plot - nucleótidos
    # ========================================================================
    
    parts = axes[0].violinplot([df['longitud_nt']], positions=[1], widths=0.7,
                                showmeans=True, showmedians=True)
    
    # Colorear
    for pc in parts['bodies']:
        pc.set_facecolor('#3498db')
        pc.set_alpha(0.7)
    
    axes[0].set_ylabel('Longitud (bp)', fontsize=12, fontweight='bold')
    axes[0].set_title('Distribución de Tamaños de Genes (Nucleótidos)',
                      fontsize=14, fontweight='bold', pad=15)
    axes[0].set_xticks([1])
    axes[0].set_xticklabels(['Todos los genes'])
    axes[0].grid(axis='y', alpha=0.3)
    
    # ========================================================================
    # Subplot 2: Violin plot - aminoácidos
    # ========================================================================
    
    parts2 = axes[1].violinplot([df['longitud_aa']], positions=[1], widths=0.7,
                                 showmeans=True, showmedians=True)
    
    # Colorear
    for pc in parts2['bodies']:
        pc.set_facecolor('#2ecc71')
        pc.set_alpha(0.7)
    
    axes[1].set_ylabel('Longitud (aminoácidos)', fontsize=12, fontweight='bold')
    axes[1].set_title('Distribución de Tamaños de Proteínas (Aminoácidos)',
                      fontsize=14, fontweight='bold', pad=15)
    axes[1].set_xticks([1])
    axes[1].set_xticklabels(['Todas las proteínas'])
    axes[1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    # Guardar
    output_file = os.path.join(OUTPUT_DIR, 'gene_size_violin.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Gráfico guardado: {output_file}")
    plt.close()

def plot_extreme_genes():
    """Gráfico de genes extremos."""
    
    # Cargar datos
    df_extremos = pd.read_csv(os.path.join(RESULTS_DIR, 'extreme_genes.csv'))
    
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    # ========================================================================
    # Subplot 1: Genes más pequeños
    # ========================================================================
    
    pequenos = df_extremos[df_extremos['tipo'] == 'Más pequeño'].head(10)
    
    y_pos = np.arange(len(pequenos))
    axes[0].barh(y_pos, pequenos['longitud_nt'], color='#e74c3c', alpha=0.8, edgecolor='black')
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels([f"{row['locus_tag']}\n({row['longitud_aa']} aa)" 
                             for _, row in pequenos.iterrows()], fontsize=9)
    axes[0].set_xlabel('Longitud (bp)', fontsize=12, fontweight='bold')
    axes[0].set_title('10 Genes Más Pequeños',
                      fontsize=14, fontweight='bold', pad=15)
    axes[0].grid(axis='x', alpha=0.3)
    axes[0].invert_yaxis()
    
    # Añadir valores
    for i, (_, row) in enumerate(pequenos.iterrows()):
        axes[0].text(row['longitud_nt'], i,
                    f" {row['longitud_nt']} bp",
                    va='center', fontsize=9, fontweight='bold')
    
    # ========================================================================
    # Subplot 2: Genes más grandes
    # ========================================================================
    
    grandes = df_extremos[df_extremos['tipo'] == 'Más grande'].tail(10)
    
    y_pos = np.arange(len(grandes))
    axes[1].barh(y_pos, grandes['longitud_nt'], color='#3498db', alpha=0.8, edgecolor='black')
    axes[1].set_yticks(y_pos)
    axes[1].set_yticklabels([f"{row['locus_tag']}\n({row['longitud_aa']} aa)" 
                             for _, row in grandes.iterrows()], fontsize=9)
    axes[1].set_xlabel('Longitud (bp)', fontsize=12, fontweight='bold')
    axes[1].set_title('10 Genes Más Grandes',
                      fontsize=14, fontweight='bold', pad=15)
    axes[1].grid(axis='x', alpha=0.3)
    axes[1].invert_yaxis()
    
    # Añadir valores
    for i, (_, row) in enumerate(grandes.iterrows()):
        axes[1].text(row['longitud_nt'], i,
                    f" {row['longitud_nt']:,} bp",
                    va='center', fontsize=9, fontweight='bold')
    
    # Título general
    fig.suptitle('Genes Extremos en E. coli K-12 MG1655',
                 fontsize=16, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    
    # Guardar
    output_file = os.path.join(OUTPUT_DIR, 'extreme_genes.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Gráfico guardado: {output_file}")
    plt.close()

# ============================================================================
# MAIN
# ============================================================================

def main():
    print("="*70)
    print("GENERANDO VISUALIZACIONES - DISTRIBUCIÓN DE GENES")
    print("="*70)
    
    print("\n[1/3] Generando histogramas de distribución...")
    plot_gene_size_histogram()
    
    print("[2/3] Generando violin plots...")
    plot_gene_size_violin()
    
    print("[3/3] Generando gráficos de genes extremos...")
    plot_extreme_genes()
    
    print("\n" + "="*70)
    print("✅ VISUALIZACIONES DE DISTRIBUCIÓN COMPLETADAS")
    print("="*70)
    print(f"\nGráficos guardados en: {OUTPUT_DIR}\n")

if __name__ == '__main__':
    main()
