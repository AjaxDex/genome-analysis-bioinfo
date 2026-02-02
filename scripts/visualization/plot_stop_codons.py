#!/usr/bin/env python3
"""
Visualización de Distribución de Stop Codons
============================================
Genera gráficos profesionales del análisis de codones de terminación.

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

# ============================================================================
# FUNCIONES DE VISUALIZACIÓN
# ============================================================================

def plot_stop_codons_comparison():
    """Comparación de proporciones de stop codons: Genoma vs CDS."""
    
    # Cargar datos
    with open(os.path.join(RESULTS_DIR, 'stop_codons_analysis.json'), 'r') as f:
        data = json.load(f)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Preparar datos
    stop_codons = ['TAA', 'TAG', 'TGA']
    genoma_pct = [data['analisis_stop_codons']['por_tipo'][sc]['proporcion_pct'] for sc in stop_codons]
    cds_pct = [data['analisis_cds']['stop_codons_por_tipo'][sc]['proporcion_pct'] for sc in stop_codons]
    
    # ========================================================================
    # Subplot 1: Barras agrupadas
    # ========================================================================
    
    x = np.arange(len(stop_codons))
    width = 0.35
    
    bars1 = axes[0].bar(x - width/2, genoma_pct, width, label='Genoma completo',
                        color='#3498db', alpha=0.8, edgecolor='black')
    bars2 = axes[0].bar(x + width/2, cds_pct, width, label='CDS anotados',
                        color='#2ecc71', alpha=0.8, edgecolor='black')
    
    axes[0].set_xlabel('Codón de Terminación', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('Proporción (%)', fontsize=12, fontweight='bold')
    axes[0].set_title('Proporciones de Stop Codons: Genoma vs CDS',
                      fontsize=14, fontweight='bold', pad=20)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(stop_codons, fontsize=12, fontweight='bold')
    axes[0].legend(fontsize=11)
    axes[0].grid(axis='y', alpha=0.3)
    
    # Añadir valores
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            axes[0].text(bar.get_x() + bar.get_width()/2., height,
                        f'{height:.1f}%',
                        ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # ========================================================================
    # Subplot 2: Preferencias evolutivas (enriquecimiento)
    # ========================================================================
    
    # Calcular enriquecimiento
    with open(os.path.join(RESULTS_DIR, 'validation_results.json'), 'r') as f:
        val_data = json.load(f)
    
    preferencias = val_data['preferencias_evolutivas']
    enriquecimientos = [p['enriquecimiento'] for p in preferencias]
    clasificaciones = [p['preferencia'] for p in preferencias]
    
    # Colores según clasificación
    colors = []
    for c in clasificaciones:
        if 'PREFERIDO' in c:
            colors.append('#2ecc71')  # Verde
        elif 'EVITADO' in c:
            colors.append('#e74c3c')  # Rojo
        else:
            colors.append('#95a5a6')  # Gris
    
    bars = axes[1].bar(stop_codons, enriquecimientos, color=colors, alpha=0.8, edgecolor='black')
    axes[1].axhline(y=1.0, color='black', linestyle='--', linewidth=2, label='Sin preferencia (1.0x)')
    axes[1].set_xlabel('Codón de Terminación', fontsize=12, fontweight='bold')
    axes[1].set_ylabel('Factor de Enriquecimiento', fontsize=12, fontweight='bold')
    axes[1].set_title('Preferencias Evolutivas de Stop Codons en CDS',
                      fontsize=14, fontweight='bold', pad=20)
    axes[1].set_xticklabels(stop_codons, fontsize=12, fontweight='bold')
    axes[1].legend(fontsize=11)
    axes[1].grid(axis='y', alpha=0.3)
    
    # Añadir valores y clasificación
    for i, (bar, pref) in enumerate(zip(bars, preferencias)):
        height = bar.get_height()
        axes[1].text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.2f}x',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        # Añadir clasificación
        clasificacion = pref['preferencia'].split()[1] if len(pref['preferencia'].split()) > 1 else pref['preferencia']
        axes[1].text(bar.get_x() + bar.get_width()/2., -0.15,
                    clasificacion,
                    ha='center', va='top', fontsize=9, style='italic')
    
    # Ajustar límites
    axes[1].set_ylim(0, max(enriquecimientos) * 1.2)
    
    # Título general
    fig.suptitle('Análisis de Codones de Terminación en E. coli K-12 MG1655',
                 fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    
    # Guardar
    output_file = os.path.join(OUTPUT_DIR, 'stop_codons_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Gráfico guardado: {output_file}")
    plt.close()

def plot_stop_codons_pie():
    """Gráficos de pastel para stop codons."""
    
    # Cargar datos
    with open(os.path.join(RESULTS_DIR, 'stop_codons_analysis.json'), 'r') as f:
        data = json.load(f)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    stop_codons = ['TAA', 'TAG', 'TGA']
    colors = ['#e74c3c', '#3498db', '#f39c12']
    
    # ========================================================================
    # Subplot 1: Genoma completo
    # ========================================================================
    
    genoma_pct = [data['analisis_stop_codons']['por_tipo'][sc]['proporcion_pct'] for sc in stop_codons]
    
    wedges1, texts1, autotexts1 = axes[0].pie(genoma_pct, labels=stop_codons,
                                                colors=colors, autopct='%1.1f%%',
                                                startangle=90, textprops={'fontsize': 12})
    
    for autotext in autotexts1:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(14)
    
    for text in texts1:
        text.set_fontweight('bold')
        text.set_fontsize(13)
    
    axes[0].set_title('Stop Codons en\nGenoma Completo',
                      fontsize=14, fontweight='bold', pad=20)
    
    # ========================================================================
    # Subplot 2: CDS anotados
    # ========================================================================
    
    cds_pct = [data['analisis_cds']['stop_codons_por_tipo'][sc]['proporcion_pct'] for sc in stop_codons]
    explode = (0.1, 0, 0)  # Destacar TAA
    
    wedges2, texts2, autotexts2 = axes[1].pie(cds_pct, labels=stop_codons,
                                                colors=colors, autopct='%1.1f%%',
                                                startangle=90, explode=explode,
                                                textprops={'fontsize': 12})
    
    for autotext in autotexts2:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(14)
    
    for text in texts2:
        text.set_fontweight('bold')
        text.set_fontsize(13)
    
    axes[1].set_title('Stop Codons en\nCDS Anotados\n(⭐ TAA preferido)',
                      fontsize=14, fontweight='bold', pad=20)
    
    # Título general
    fig.suptitle('Distribución de Codones de Terminación',
                 fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    
    # Guardar
    output_file = os.path.join(OUTPUT_DIR, 'stop_codons_pie.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Gráfico guardado: {output_file}")
    plt.close()

# ============================================================================
# MAIN
# ============================================================================

def main():
    print("="*70)
    print("GENERANDO VISUALIZACIONES - STOP CODONS")
    print("="*70)
    
    print("\n[1/2] Generando gráfico de comparación...")
    plot_stop_codons_comparison()
    
    print("[2/2] Generando gráficos de pastel...")
    plot_stop_codons_pie()
    
    print("\n" + "="*70)
    print("✅ VISUALIZACIONES DE STOP CODONS COMPLETADAS")
    print("="*70)
    print(f"\nGráficos guardados en: {OUTPUT_DIR}\n")

if __name__ == '__main__':
    main()
