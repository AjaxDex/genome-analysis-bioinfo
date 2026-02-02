#!/usr/bin/env python3
"""
Visualización de Distribución de Codones ATG
============================================
Genera gráficos profesionales del análisis de codones ATG.

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
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10

# ============================================================================
# FUNCIONES DE VISUALIZACIÓN
# ============================================================================

def plot_atg_comparison():
    """Gráfico de comparación ATG vs CDS."""
    
    # Cargar datos
    with open(os.path.join(RESULTS_DIR, 'atg_analysis.json'), 'r') as f:
        data = json.load(f)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # ========================================================================
    # Subplot 1: Barras de comparación
    # ========================================================================
    
    categorias = ['Total ATG\nen genoma', 'CDS\nanotados', 'ATG no\ncodificantes']
    valores = [
        data['analisis_atg']['total_atg'],
        data['genes_anotados']['total_cds'],
        data['comparacion']['atg_no_codificantes_estimados']
    ]
    colores = ['#3498db', '#2ecc71', '#e74c3c']
    
    bars = axes[0].bar(categorias, valores, color=colores, alpha=0.8, edgecolor='black')
    axes[0].set_ylabel('Cantidad', fontsize=12, fontweight='bold')
    axes[0].set_title('Comparación: Codones ATG vs Genes Funcionales', 
                      fontsize=14, fontweight='bold', pad=20)
    axes[0].grid(axis='y', alpha=0.3)
    
    # Añadir valores en las barras
    for bar in bars:
        height = bar.get_height()
        axes[0].text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height):,}',
                    ha='center', va='bottom', fontweight='bold')
    
    # ========================================================================
    # Subplot 2: Gráfico de pastel
    # ========================================================================
    
    labels = ['CDS funcionales\n(inicio real)', 'ATG no codificantes\n(falsos positivos)']
    sizes = [
        data['genes_anotados']['total_cds'],
        data['comparacion']['atg_no_codificantes_estimados']
    ]
    colors = ['#2ecc71', '#e74c3c']
    explode = (0.05, 0)
    
    wedges, texts, autotexts = axes[1].pie(sizes, explode=explode, labels=labels,
                                             colors=colors, autopct='%1.1f%%',
                                             startangle=90, textprops={'fontsize': 11})
    
    # Mejorar el texto
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(12)
    
    axes[1].set_title('Proporción de ATG Funcionales vs No Funcionales',
                      fontsize=14, fontweight='bold', pad=20)
    
    # Título general
    fig.suptitle('Análisis de Codones ATG en E. coli K-12 MG1655',
                 fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    
    # Guardar
    output_file = os.path.join(OUTPUT_DIR, 'atg_distribution.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Gráfico guardado: {output_file}")
    plt.close()

def plot_atg_density():
    """Gráfico de densidad de ATG."""
    
    # Cargar datos
    with open(os.path.join(RESULTS_DIR, 'atg_analysis.json'), 'r') as f:
        data = json.load(f)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Datos
    densidad_atg = data['analisis_atg']['densidad_por_kb']
    tamano_genoma_mb = data['genoma']['tamano_bp'] / 1_000_000
    total_atg = data['analisis_atg']['total_atg']
    
    # Crear visualización
    info_text = f"""
    Tamaño del genoma: {tamano_genoma_mb:.2f} Mb
    Total codones ATG: {total_atg:,}
    Densidad: {densidad_atg:.2f} ATG/kb
    
    Interpretación:
    • Se encuentra ~1 ATG cada {1000/densidad_atg:.0f} nucleótidos
    • Solo ~{data['genes_anotados']['total_cds']/total_atg*100:.1f}% son inicios reales de genes
    """
    
    # Gráfico de barras horizontal
    categories = ['Densidad ATG\n(por kilobase)']
    values = [densidad_atg]
    
    bars = ax.barh(categories, values, color='#3498db', alpha=0.8, edgecolor='black', height=0.5)
    ax.set_xlabel('Codones ATG por kb', fontsize=12, fontweight='bold')
    ax.set_title('Densidad de Codones ATG en el Genoma', fontsize=14, fontweight='bold', pad=20)
    ax.grid(axis='x', alpha=0.3)
    
    # Añadir valor
    for bar in bars:
        width = bar.get_width()
        ax.text(width, bar.get_y() + bar.get_height()/2.,
                f'{width:.2f}',
                ha='left', va='center', fontweight='bold', fontsize=14,
                bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='black'))
    
    # Añadir texto informativo
    ax.text(0.98, 0.02, info_text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment='bottom',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    # Guardar
    output_file = os.path.join(OUTPUT_DIR, 'atg_density.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Gráfico guardado: {output_file}")
    plt.close()

# ============================================================================
# MAIN
# ============================================================================

def main():
    print("="*70)
    print("GENERANDO VISUALIZACIONES - CODONES ATG")
    print("="*70)
    
    print("\n[1/2] Generando gráfico de comparación ATG...")
    plot_atg_comparison()
    
    print("[2/2] Generando gráfico de densidad ATG...")
    plot_atg_density()
    
    print("\n" + "="*70)
    print("✅ VISUALIZACIONES DE ATG COMPLETADAS")
    print("="*70)
    print(f"\nGráficos guardados en: {OUTPUT_DIR}\n")

if __name__ == '__main__':
    main()
