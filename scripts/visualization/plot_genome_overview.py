#!/usr/bin/env python3
"""
Visualización Overview del Genoma
==================================
Genera un dashboard completo con las estadísticas principales del genoma.

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

# ============================================================================
# FUNCIONES DE VISUALIZACIÓN
# ============================================================================

def plot_genome_overview():
    """Dashboard completo del genoma."""
    
    # Cargar todos los datos
    with open(os.path.join(RESULTS_DIR, 'genome_statistics.json'), 'r') as f:
        genome_data = json.load(f)
    
    with open(os.path.join(RESULTS_DIR, 'atg_analysis.json'), 'r') as f:
        atg_data = json.load(f)
    
    with open(os.path.join(RESULTS_DIR, 'stop_codons_analysis.json'), 'r') as f:
        stop_data = json.load(f)
    
    # Crear figura con subplots
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # ========================================================================
    # Panel 1: Información básica del genoma
    # ========================================================================
    
    ax1 = fig.add_subplot(gs[0, :])
    ax1.axis('off')
    
    info_text = f"""
    GENOMA: {genome_data['genoma']['id']}
    {genome_data['genoma']['descripcion']}
    
    Tamaño: {genome_data['genoma']['tamano_bp']:,} bp ({genome_data['genoma']['tamano_bp']/1e6:.2f} Mb)
    Contenido GC: {genome_data['genoma']['contenido_gc_pct']:.2f}%
    Total de genes: {genome_data['genes']['total']:,}
    Total de CDS: {genome_data['cds']['total']:,}
    Densidad génica: {genome_data['densidad_genica']['genes_por_mb']:.2f} genes/Mb
    Región codificante: {genome_data['regiones']['codificantes']['porcentaje']:.2f}%
    """
    
    ax1.text(0.5, 0.5, info_text, transform=ax1.transAxes,
             fontsize=14, verticalalignment='center', horizontalalignment='center',
             bbox=dict(boxstyle='round,pad=1', facecolor='lightblue', alpha=0.8),
             family='monospace', fontweight='bold')
    
    ax1.set_title('RESUMEN GENERAL DEL GENOMA',
                  fontsize=18, fontweight='bold', pad=20)
    
    # ========================================================================
    # Panel 2: Composición de nucleótidos
    # ========================================================================
    
    ax2 = fig.add_subplot(gs[1, 0])
    
    nucleotidos = ['A', 'T', 'G', 'C']
    valores = [
        genome_data['genoma']['composicion_nucleotidos']['A']['porcentaje'],
        genome_data['genoma']['composicion_nucleotidos']['T']['porcentaje'],
        genome_data['genoma']['composicion_nucleotidos']['G']['porcentaje'],
        genome_data['genoma']['composicion_nucleotidos']['C']['porcentaje']
    ]
    colors = ['#e74c3c', '#f39c12', '#2ecc71', '#3498db']
    
    bars = ax2.bar(nucleotidos, valores, color=colors, alpha=0.8, edgecolor='black')
    ax2.set_ylabel('Porcentaje (%)', fontsize=11, fontweight='bold')
    ax2.set_title('Composición de Nucleótidos', fontsize=12, fontweight='bold')
    ax2.grid(axis='y', alpha=0.3)
    
    for bar in bars:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.2f}%',
                ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # ========================================================================
    # Panel 3: Regiones codificantes vs no codificantes
    # ========================================================================
    
    ax3 = fig.add_subplot(gs[1, 1])
    
    labels = ['Codificante', 'No codificante']
    sizes = [
        genome_data['regiones']['codificantes']['porcentaje'],
        genome_data['regiones']['no_codificantes']['porcentaje']
    ]
    colors_regions = ['#2ecc71', '#e74c3c']
    explode = (0.05, 0)
    
    wedges, texts, autotexts = ax3.pie(sizes, explode=explode, labels=labels,
                                         colors=colors_regions, autopct='%1.1f%%',
                                         startangle=90, textprops={'fontsize': 10})
    
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    ax3.set_title('Regiones Genómicas', fontsize=12, fontweight='bold')
    
    # ========================================================================
    # Panel 4: Distribución de strands
    # ========================================================================
    
    ax4 = fig.add_subplot(gs[1, 2])
    
    strand_labels = ['Hebra (+)', 'Hebra (-)']
    strand_values = [
        genome_data['cds']['distribucion_strands']['porcentaje_plus'],
        genome_data['cds']['distribucion_strands']['porcentaje_minus']
    ]
    colors_strand = ['#3498db', '#9b59b6']
    
    wedges2, texts2, autotexts2 = ax4.pie(strand_values, labels=strand_labels,
                                            colors=colors_strand, autopct='%1.1f%%',
                                            startangle=90, textprops={'fontsize': 10})
    
    for autotext in autotexts2:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    ax4.set_title('Distribución en Hebras', fontsize=12, fontweight='bold')
    
    # ========================================================================
    # Panel 5: Codones ATG
    # ========================================================================
    
    ax5 = fig.add_subplot(gs[2, 0])
    
    categorias_atg = ['Total ATG', 'CDS', 'No funcionales']
    valores_atg = [
        atg_data['analisis_atg']['total_atg'],
        atg_data['genes_anotados']['total_cds'],
        atg_data['comparacion']['atg_no_codificantes_estimados']
    ]
    colors_atg = ['#3498db', '#2ecc71', '#e74c3c']
    
    bars_atg = ax5.bar(categorias_atg, valores_atg, color=colors_atg, alpha=0.8, edgecolor='black')
    ax5.set_ylabel('Cantidad', fontsize=11, fontweight='bold')
    ax5.set_title('Codones ATG', fontsize=12, fontweight='bold')
    ax5.grid(axis='y', alpha=0.3)
    ax5.tick_params(axis='x', labelsize=9)
    
    for bar in bars_atg:
        height = bar.get_height()
        ax5.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height):,}',
                ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # ========================================================================
    # Panel 6: Stop Codons en CDS
    # ========================================================================
    
    ax6 = fig.add_subplot(gs[2, 1])
    
    stop_codons = ['TAA', 'TAG', 'TGA']
    stop_values = [
        stop_data['analisis_cds']['stop_codons_por_tipo']['TAA']['proporcion_pct'],
        stop_data['analisis_cds']['stop_codons_por_tipo']['TAG']['proporcion_pct'],
        stop_data['analisis_cds']['stop_codons_por_tipo']['TGA']['proporcion_pct']
    ]
    
    colors_stop = ['#2ecc71', '#e74c3c', '#f39c12']
    
    bars_stop = ax6.bar(stop_codons, stop_values, color=colors_stop, alpha=0.8, edgecolor='black')
    ax6.set_ylabel('Porcentaje (%)', fontsize=11, fontweight='bold')
    ax6.set_title('Stop Codons en CDS', fontsize=12, fontweight='bold')
    ax6.grid(axis='y', alpha=0.3)
    
    for bar in bars_stop:
        height = bar.get_height()
        ax6.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}%',
                ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # ========================================================================
    # Panel 7: Estadísticas de tamaño de genes
    # ========================================================================
    
    ax7 = fig.add_subplot(gs[2, 2])
    
    stats_labels = ['Mínimo', 'Q1', 'Mediana', 'Q3', 'Máximo']
    stats_values = [
        genome_data['cds']['estadisticas_tamanos']['minimo'],
        genome_data['cds']['estadisticas_tamanos']['percentil_25'],
        genome_data['cds']['estadisticas_tamanos']['mediana'],
        genome_data['cds']['estadisticas_tamanos']['percentil_75'],
        genome_data['cds']['estadisticas_tamanos']['maximo']
    ]
    
    ax7.plot(stats_labels, stats_values, marker='o', linewidth=2, 
             markersize=8, color='#3498db')
    ax7.fill_between(range(len(stats_labels)), stats_values, alpha=0.3, color='#3498db')
    ax7.set_ylabel('Tamaño (bp)', fontsize=11, fontweight='bold')
    ax7.set_title('Distribución Tamaño Genes', fontsize=12, fontweight='bold')
    ax7.grid(True, alpha=0.3)
    ax7.tick_params(axis='x', labelsize=9, rotation=45)
    
    # Añadir valores
    for i, (label, value) in enumerate(zip(stats_labels, stats_values)):
        ax7.text(i, value, f'{int(value):,}',
                ha='center', va='bottom', fontsize=8, fontweight='bold')
    
    # Título general
    fig.suptitle('OVERVIEW COMPLETO DEL GENOMA - E. coli K-12 MG1655',
                 fontsize=20, fontweight='bold', y=0.98)
    
    # Guardar
    output_file = os.path.join(OUTPUT_DIR, 'genome_overview.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Gráfico guardado: {output_file}")
    plt.close()

# ============================================================================
# MAIN
# ============================================================================

def main():
    print("="*70)
    print("GENERANDO VISUALIZACIÓN - OVERVIEW DEL GENOMA")
    print("="*70)
    
    print("\n[1/1] Generando dashboard completo...")
    plot_genome_overview()
    
    print("\n" + "="*70)
    print("✅ OVERVIEW DEL GENOMA COMPLETADO")
    print("="*70)
    print(f"\nGráfico guardado en: {OUTPUT_DIR}\n")

if __name__ == '__main__':
    main()
