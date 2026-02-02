#!/usr/bin/env python3
"""
Validación de Resultados con Literatura Científica
==================================================
Este script compara los resultados obtenidos con valores
reportados en publicaciones científicas sobre E. coli K-12 MG1655.

Autor: Proyecto Bioinformática
Fecha: 2024
"""

import os
import json
import pandas as pd
from datetime import datetime

# ============================================================================
# CONFIGURACIÓN
# ============================================================================

RESULTS_DIR = os.path.expanduser('~/projects/bioinfo/results/tables')
OUTPUT_DIR = os.path.expanduser('~/projects/bioinfo/results/tables')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# VALORES ESPERADOS DE LA LITERATURA
# ============================================================================

VALORES_ESPERADOS = {
    'genoma': {
        'tamano_bp': 4_641_652,  # Valor exacto de NCBI
        'tolerancia_bp': 1000,
        'total_genes': 4300,  # Aproximado según literatura
        'tolerancia_genes': 200,
        'contenido_gc_pct': 50.8,
        'tolerancia_gc': 0.5,
        'fuente': 'Blattner et al. (1997) Science; Riley et al. (2006) Nucleic Acids Res'
    },
    'stop_codons_cds': {
        'TAA': {
            'proporcion_pct': 61.0,  # Basado en estudios de E. coli
            'tolerancia_pct': 5.0,
            'rango_esperado': [55, 70],
            'interpretacion': 'Codón de terminación preferido en E. coli'
        },
        'TAG': {
            'proporcion_pct': 9.0,
            'tolerancia_pct': 3.0,
            'rango_esperado': [5, 15],
            'interpretacion': 'Codón de terminación menos usado'
        },
        'TGA': {
            'proporcion_pct': 30.0,
            'tolerancia_pct': 5.0,
            'rango_esperado': [25, 35],
            'interpretacion': 'Codón de terminación intermedio'
        },
        'fuente': 'Nakamura et al. (2000) Nucleic Acids Res; Sharp et al. (2010)'
    },
    'ratio_atg_cds': {
        'valor_esperado': 17.5,
        'tolerancia': 2.0,
        'rango_esperado': [15, 20],
        'interpretacion': 'Mayoría de ATG no son inicios de genes funcionales',
        'fuente': 'Análisis genómicos de E. coli K-12'
    }
}

# ============================================================================
# FUNCIONES DE VALIDACIÓN
# ============================================================================

def cargar_resultados():
    """Carga los resultados de los análisis anteriores."""
    
    # Cargar ATG
    with open(os.path.join(RESULTS_DIR, 'atg_analysis.json'), 'r') as f:
        atg_data = json.load(f)
    
    # Cargar Stop Codons
    with open(os.path.join(RESULTS_DIR, 'stop_codons_analysis.json'), 'r') as f:
        stop_data = json.load(f)
    
    return atg_data, stop_data

def validar_valor(obtenido, esperado, tolerancia):
    """
    Valida si un valor obtenido está dentro del rango esperado.
    
    Returns:
        dict: Resultado de validación
    """
    diferencia = abs(obtenido - esperado)
    en_rango = diferencia <= tolerancia
    desviacion_pct = (diferencia / esperado) * 100 if esperado != 0 else 0
    
    return {
        'obtenido': obtenido,
        'esperado': esperado,
        'diferencia': round(diferencia, 2),
        'desviacion_pct': round(desviacion_pct, 2),
        'en_rango': en_rango,
        'tolerancia': tolerancia,
        'status': '✅ VÁLIDO' if en_rango else '⚠️ FUERA DE RANGO'
    }

def validar_rango(obtenido, rango_min, rango_max):
    """Valida si un valor está dentro de un rango."""
    en_rango = rango_min <= obtenido <= rango_max
    
    return {
        'obtenido': obtenido,
        'rango_min': rango_min,
        'rango_max': rango_max,
        'en_rango': en_rango,
        'status': '✅ DENTRO DEL RANGO' if en_rango else '⚠️ FUERA DEL RANGO'
    }

# ============================================================================
# ANÁLISIS DE VALIDACIÓN
# ============================================================================

def main():
    print("="*80)
    print("VALIDACIÓN DE RESULTADOS CON LITERATURA CIENTÍFICA")
    print("="*80)
    
    # Cargar resultados
    print("\n[1/5] Cargando resultados de análisis previos...")
    atg_data, stop_data = cargar_resultados()
    print("✓ Datos cargados correctamente")
    
    # Preparar estructura de validación
    validaciones = {
        'metadata': {
            'fecha_validacion': datetime.now().isoformat(),
            'genoma_id': atg_data['genoma']['id']
        },
        'validaciones': {}
    }
    
    # ========================================================================
    # VALIDACIÓN 1: Tamaño del Genoma
    # ========================================================================
    
    print("\n[2/5] VALIDACIÓN 1: Tamaño del Genoma")
    print("-" * 80)
    
    tamano_obtenido = atg_data['genoma']['tamano_bp']
    val_tamano = validar_valor(
        tamano_obtenido,
        VALORES_ESPERADOS['genoma']['tamano_bp'],
        VALORES_ESPERADOS['genoma']['tolerancia_bp']
    )
    
    print(f"Tamaño obtenido:     {val_tamano['obtenido']:,} bp")
    print(f"Tamaño esperado:     {val_tamano['esperado']:,} bp")
    print(f"Diferencia:          {val_tamano['diferencia']:,} bp")
    print(f"Estado:              {val_tamano['status']}")
    
    validaciones['validaciones']['tamano_genoma'] = val_tamano
    
    # ========================================================================
    # VALIDACIÓN 2: Número de Genes
    # ========================================================================
    
    print(f"\n[3/5] VALIDACIÓN 2: Número Total de Genes")
    print("-" * 80)
    
    genes_obtenidos = atg_data['genes_anotados']['total_cds']
    val_genes = validar_valor(
        genes_obtenidos,
        VALORES_ESPERADOS['genoma']['total_genes'],
        VALORES_ESPERADOS['genoma']['tolerancia_genes']
    )
    
    print(f"Genes obtenidos:     {val_genes['obtenido']:,}")
    print(f"Genes esperados:     ~{val_genes['esperado']:,}")
    print(f"Diferencia:          {val_genes['diferencia']:,}")
    print(f"Desviación:          {val_genes['desviacion_pct']:.2f}%")
    print(f"Estado:              {val_genes['status']}")
    
    validaciones['validaciones']['total_genes'] = val_genes
    
    # ========================================================================
    # VALIDACIÓN 3: Ratio ATG/CDS
    # ========================================================================
    
    print(f"\n[4/5] VALIDACIÓN 3: Ratio ATG/CDS")
    print("-" * 80)
    
    ratio_obtenido = atg_data['comparacion']['ratio_atg_vs_cds']
    val_ratio = validar_rango(
        ratio_obtenido,
        VALORES_ESPERADOS['ratio_atg_cds']['rango_esperado'][0],
        VALORES_ESPERADOS['ratio_atg_cds']['rango_esperado'][1]
    )
    
    print(f"Ratio obtenido:      {val_ratio['obtenido']:.2f}x")
    print(f"Rango esperado:      {val_ratio['rango_min']}-{val_ratio['rango_max']}x")
    print(f"Estado:              {val_ratio['status']}")
    print(f"Interpretación:      {VALORES_ESPERADOS['ratio_atg_cds']['interpretacion']}")
    
    validaciones['validaciones']['ratio_atg_cds'] = val_ratio
    validaciones['validaciones']['ratio_atg_cds']['interpretacion'] = VALORES_ESPERADOS['ratio_atg_cds']['interpretacion']
    
    # ========================================================================
    # VALIDACIÓN 4: Proporciones de Stop Codons en CDS
    # ========================================================================
    
    print(f"\n[5/5] VALIDACIÓN 4: Proporciones de Stop Codons en CDS")
    print("-" * 80)
    
    validaciones['validaciones']['stop_codons'] = {}
    
    for stop_codon in ['TAA', 'TAG', 'TGA']:
        print(f"\n→ Validando {stop_codon}:")
        
        prop_obtenida = stop_data['analisis_cds']['stop_codons_por_tipo'][stop_codon]['proporcion_pct']
        esperado = VALORES_ESPERADOS['stop_codons_cds'][stop_codon]
        
        val_stop = validar_rango(
            prop_obtenida,
            esperado['rango_esperado'][0],
            esperado['rango_esperado'][1]
        )
        
        print(f"  Proporción obtenida:  {val_stop['obtenido']:.2f}%")
        print(f"  Rango esperado:       {val_stop['rango_min']}-{val_stop['rango_max']}%")
        print(f"  Estado:               {val_stop['status']}")
        print(f"  Interpretación:       {esperado['interpretacion']}")
        
        validaciones['validaciones']['stop_codons'][stop_codon] = {
            **val_stop,
            'interpretacion': esperado['interpretacion']
        }
    
    # ========================================================================
    # ANÁLISIS DE PREFERENCIAS EVOLUTIVAS
    # ========================================================================
    
    print("\n" + "="*80)
    print("ANÁLISIS DE PREFERENCIAS EVOLUTIVAS DE STOP CODONS")
    print("="*80)
    
    preferencias = []
    for stop_codon in ['TAA', 'TAG', 'TGA']:
        cds_pct = stop_data['analisis_cds']['stop_codons_por_tipo'][stop_codon]['proporcion_pct']
        genoma_pct = stop_data['analisis_stop_codons']['por_tipo'][stop_codon]['proporcion_pct']
        enriquecimiento = cds_pct / genoma_pct if genoma_pct > 0 else 0
        
        if enriquecimiento > 1.2:
            preferencia = "⭐ FUERTEMENTE PREFERIDO"
        elif enriquecimiento > 0.8:
            preferencia = "○ NEUTRAL"
        else:
            preferencia = "❌ EVITADO"
        
        preferencias.append({
            'stop_codon': stop_codon,
            'cds_pct': cds_pct,
            'genoma_pct': genoma_pct,
            'enriquecimiento': enriquecimiento,
            'preferencia': preferencia
        })
        
        print(f"\n{stop_codon}:")
        print(f"  En CDS:              {cds_pct:.2f}%")
        print(f"  En genoma completo:  {genoma_pct:.2f}%")
        print(f"  Enriquecimiento:     {enriquecimiento:.2f}x")
        print(f"  Clasificación:       {preferencia}")
    
    validaciones['preferencias_evolutivas'] = preferencias
    
    # ========================================================================
    # GUARDAR RESULTADOS
    # ========================================================================
    
    print("\n" + "="*80)
    print("GUARDANDO RESULTADOS DE VALIDACIÓN")
    print("="*80)
    
    # JSON con validaciones completas
    json_file = os.path.join(OUTPUT_DIR, 'validation_results.json')
    with open(json_file, 'w') as f:
        json.dump(validaciones, f, indent=2)
    print(f"✓ JSON guardado: {json_file}")
    
    # CSV con resumen de validaciones
    df_validacion = pd.DataFrame([
        {
            'Parámetro': 'Tamaño Genoma (bp)',
            'Obtenido': f"{val_tamano['obtenido']:,}",
            'Esperado': f"{val_tamano['esperado']:,}",
            'Status': val_tamano['status']
        },
        {
            'Parámetro': 'Total Genes',
            'Obtenido': val_genes['obtenido'],
            'Esperado': f"~{val_genes['esperado']}",
            'Status': val_genes['status']
        },
        {
            'Parámetro': 'Ratio ATG/CDS',
            'Obtenido': f"{val_ratio['obtenido']:.2f}x",
            'Esperado': f"{val_ratio['rango_min']}-{val_ratio['rango_max']}x",
            'Status': val_ratio['status']
        }
    ])
    
    csv_file = os.path.join(OUTPUT_DIR, 'validation_summary.csv')
    df_validacion.to_csv(csv_file, index=False)
    print(f"✓ CSV guardado: {csv_file}")
    
    # CSV con análisis de stop codons
    df_stops = pd.DataFrame(preferencias)
    stops_file = os.path.join(OUTPUT_DIR, 'stop_codons_preferences.csv')
    df_stops.to_csv(stops_file, index=False)
    print(f"✓ Preferencias guardadas: {stops_file}")
    
    # ========================================================================
    # RESUMEN FINAL
    # ========================================================================
    
    print("\n" + "="*80)
    print("RESUMEN DE VALIDACIÓN")
    print("="*80)
    
    print("\n📊 ESTADÍSTICAS GENERALES:")
    print(f"  • Tamaño genoma:  {val_tamano['status']}")
    print(f"  • Total genes:    {val_genes['status']}")
    print(f"  • Ratio ATG/CDS:  {val_ratio['status']}")
    
    print("\n🧬 STOP CODONS EN CDS:")
    for stop_codon in ['TAA', 'TAG', 'TGA']:
        status = validaciones['validaciones']['stop_codons'][stop_codon]['status']
        print(f"  • {stop_codon}: {status}")
    
    print("\n⭐ PREFERENCIAS EVOLUTIVAS:")
    for pref in preferencias:
        print(f"  • {pref['stop_codon']}: {pref['preferencia']} ({pref['enriquecimiento']:.2f}x)")
    
    print("\n📚 FUENTES:")
    print(f"  • Genoma: {VALORES_ESPERADOS['genoma']['fuente']}")
    print(f"  • Stop codons: {VALORES_ESPERADOS['stop_codons_cds']['fuente']}")
    
    print("\n" + "="*80)
    print("✅ VALIDACIÓN COMPLETADA EXITOSAMENTE")
    print("="*80 + "\n")

if __name__ == '__main__':
    main()
