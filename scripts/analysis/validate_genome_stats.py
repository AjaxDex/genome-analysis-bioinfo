#!/usr/bin/env python3
"""
Validación de Estadísticas Genómicas con Literatura
===================================================
Valida los resultados de estadísticas genómicas contra
valores reportados en publicaciones científicas.

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

# ============================================================================
# VALORES ESPERADOS DE LA LITERATURA
# ============================================================================

VALORES_ESPERADOS = {
    'contenido_gc': {
        'valor': 50.8,
        'tolerancia': 0.2,
        'fuente': 'Blattner et al. (1997) Science',
        'interpretacion': 'Equilibrio AT/GC característico de E. coli'
    },
    'tamano_genoma': {
        'valor': 4_641_652,
        'tolerancia': 100,
        'fuente': 'NCBI RefSeq NC_000913.3',
        'interpretacion': 'Secuencia de referencia completa'
    },
    'porcentaje_codificante': {
        'valor': 87.8,
        'tolerancia': 2.0,
        'rango': [85, 90],
        'fuente': 'Riley et al. (2006) Nucleic Acids Res',
        'interpretacion': 'Genoma bacteriano altamente compacto'
    },
    'densidad_genica': {
        'genes_por_mb': 930,
        'tolerancia': 50,
        'fuente': 'Análisis genómico de E. coli K-12',
        'interpretacion': 'Alta densidad génica típica de procariotas'
    },
    'tamano_medio_gen': {
        'valor': 950,
        'tolerancia': 100,
        'rango': [800, 1100],
        'fuente': 'Estudios genómicos bacterianos',
        'interpretacion': 'Tamaño típico de genes bacterianos'
    },
    'distribucion_strands': {
        'equilibrio_esperado': 50,
        'tolerancia': 5,
        'fuente': 'Análisis de genomas bacterianos',
        'interpretacion': 'Balance esperado entre hebras'
    }
}

# ============================================================================
# FUNCIONES DE VALIDACIÓN
# ============================================================================

def validar_valor(obtenido, esperado, tolerancia, nombre):
    """Valida un valor contra el esperado."""
    diferencia = abs(obtenido - esperado)
    en_rango = diferencia <= tolerancia
    desviacion_pct = (diferencia / esperado) * 100 if esperado != 0 else 0
    
    return {
        'parametro': nombre,
        'obtenido': round(obtenido, 2),
        'esperado': esperado,
        'diferencia': round(diferencia, 2),
        'desviacion_pct': round(desviacion_pct, 2),
        'tolerancia': tolerancia,
        'en_rango': en_rango,
        'status': '✅ VÁLIDO' if en_rango else '⚠️ FUERA DE RANGO'
    }

def validar_rango(obtenido, rango_min, rango_max, nombre):
    """Valida si un valor está dentro de un rango."""
    en_rango = rango_min <= obtenido <= rango_max
    
    return {
        'parametro': nombre,
        'obtenido': round(obtenido, 2),
        'rango_min': rango_min,
        'rango_max': rango_max,
        'en_rango': en_rango,
        'status': '✅ DENTRO DEL RANGO' if en_rango else '⚠️ FUERA DEL RANGO'
    }

# ============================================================================
# ANÁLISIS PRINCIPAL
# ============================================================================

def main():
    print("="*80)
    print("VALIDACIÓN DE ESTADÍSTICAS GENÓMICAS")
    print("="*80)
    
    # Cargar resultados
    print("\n[1/6] Cargando resultados...")
    with open(os.path.join(RESULTS_DIR, 'genome_statistics.json'), 'r') as f:
        genome_data = json.load(f)
    print("✓ Datos cargados")
    
    validaciones = {
        'metadata': {
            'fecha_validacion': datetime.now().isoformat(),
            'genoma_id': genome_data['genoma']['id']
        },
        'validaciones': {}
    }
    
    # ========================================================================
    # VALIDACIÓN 1: Contenido GC
    # ========================================================================
    
    print("\n[2/6] VALIDACIÓN 1: Contenido GC")
    print("-" * 80)
    
    gc_obtenido = genome_data['genoma']['contenido_gc_pct']
    val_gc = validar_valor(
        gc_obtenido,
        VALORES_ESPERADOS['contenido_gc']['valor'],
        VALORES_ESPERADOS['contenido_gc']['tolerancia'],
        'Contenido GC'
    )
    
    print(f"Contenido GC obtenido:   {val_gc['obtenido']:.2f}%")
    print(f"Contenido GC esperado:   {val_gc['esperado']:.2f}%")
    print(f"Diferencia:              {val_gc['diferencia']:.2f}%")
    print(f"Estado:                  {val_gc['status']}")
    print(f"Interpretación:          {VALORES_ESPERADOS['contenido_gc']['interpretacion']}")
    
    validaciones['validaciones']['contenido_gc'] = {
        **val_gc,
        'interpretacion': VALORES_ESPERADOS['contenido_gc']['interpretacion']
    }
    
    # ========================================================================
    # VALIDACIÓN 2: Tamaño del Genoma
    # ========================================================================
    
    print("\n[3/6] VALIDACIÓN 2: Tamaño del Genoma")
    print("-" * 80)
    
    tamano_obtenido = genome_data['genoma']['tamano_bp']
    val_tamano = validar_valor(
        tamano_obtenido,
        VALORES_ESPERADOS['tamano_genoma']['valor'],
        VALORES_ESPERADOS['tamano_genoma']['tolerancia'],
        'Tamaño Genoma'
    )
    
    print(f"Tamaño obtenido:         {val_tamano['obtenido']:,.0f} bp")
    print(f"Tamaño esperado:         {val_tamano['esperado']:,} bp")
    print(f"Diferencia:              {val_tamano['diferencia']:.0f} bp")
    print(f"Estado:                  {val_tamano['status']}")
    
    validaciones['validaciones']['tamano_genoma'] = val_tamano
    
    # ========================================================================
    # VALIDACIÓN 3: Porcentaje Codificante
    # ========================================================================
    
    print("\n[4/6] VALIDACIÓN 3: Porcentaje de Región Codificante")
    print("-" * 80)
    
    codificante_obtenido = genome_data['regiones']['codificantes']['porcentaje']
    val_codif = validar_rango(
        codificante_obtenido,
        VALORES_ESPERADOS['porcentaje_codificante']['rango'][0],
        VALORES_ESPERADOS['porcentaje_codificante']['rango'][1],
        'Porcentaje Codificante'
    )
    
    print(f"Porcentaje obtenido:     {val_codif['obtenido']:.2f}%")
    print(f"Rango esperado:          {val_codif['rango_min']}-{val_codif['rango_max']}%")
    print(f"Estado:                  {val_codif['status']}")
    print(f"Interpretación:          {VALORES_ESPERADOS['porcentaje_codificante']['interpretacion']}")
    
    validaciones['validaciones']['porcentaje_codificante'] = {
        **val_codif,
        'interpretacion': VALORES_ESPERADOS['porcentaje_codificante']['interpretacion']
    }
    
    # ========================================================================
    # VALIDACIÓN 4: Densidad Génica
    # ========================================================================
    
    print("\n[5/6] VALIDACIÓN 4: Densidad Génica")
    print("-" * 80)
    
    densidad_obtenida = genome_data['densidad_genica']['genes_por_mb']
    val_densidad = validar_valor(
        densidad_obtenida,
        VALORES_ESPERADOS['densidad_genica']['genes_por_mb'],
        VALORES_ESPERADOS['densidad_genica']['tolerancia'],
        'Densidad Génica'
    )
    
    print(f"Genes/Mb obtenido:       {val_densidad['obtenido']:.2f}")
    print(f"Genes/Mb esperado:       ~{val_densidad['esperado']:.0f}")
    print(f"Diferencia:              {val_densidad['diferencia']:.2f}")
    print(f"Estado:                  {val_densidad['status']}")
    print(f"Interpretación:          {VALORES_ESPERADOS['densidad_genica']['interpretacion']}")
    
    validaciones['validaciones']['densidad_genica'] = {
        **val_densidad,
        'interpretacion': VALORES_ESPERADOS['densidad_genica']['interpretacion']
    }
    
    # ========================================================================
    # VALIDACIÓN 5: Tamaño Medio de Genes
    # ========================================================================
    
    print("\n[6/6] VALIDACIÓN 5: Tamaño Medio de Genes")
    print("-" * 80)
    
    tamano_medio_obtenido = genome_data['cds']['estadisticas_tamanos']['media']
    val_tamano_medio = validar_rango(
        tamano_medio_obtenido,
        VALORES_ESPERADOS['tamano_medio_gen']['rango'][0],
        VALORES_ESPERADOS['tamano_medio_gen']['rango'][1],
        'Tamaño Medio Gen'
    )
    
    print(f"Tamaño medio obtenido:   {val_tamano_medio['obtenido']:.2f} bp")
    print(f"Rango esperado:          {val_tamano_medio['rango_min']}-{val_tamano_medio['rango_max']} bp")
    print(f"Estado:                  {val_tamano_medio['status']}")
    print(f"Interpretación:          {VALORES_ESPERADOS['tamano_medio_gen']['interpretacion']}")
    
    validaciones['validaciones']['tamano_medio_gen'] = {
        **val_tamano_medio,
        'interpretacion': VALORES_ESPERADOS['tamano_medio_gen']['interpretacion']
    }
    
    # ========================================================================
    # ANÁLISIS ADICIONAL: Balance de Hebras
    # ========================================================================
    
    print("\n" + "="*80)
    print("ANÁLISIS ADICIONAL: Balance de Hebras")
    print("="*80)
    
    strand_plus_pct = genome_data['cds']['distribucion_strands']['porcentaje_plus']
    strand_minus_pct = genome_data['cds']['distribucion_strands']['porcentaje_minus']
    desbalance = abs(strand_plus_pct - 50)
    
    print(f"\nDistribución de CDS:")
    print(f"  • Hebra (+):             {strand_plus_pct:.2f}%")
    print(f"  • Hebra (-):             {strand_minus_pct:.2f}%")
    print(f"  • Desbalance:            {desbalance:.2f}%")
    
    if desbalance <= VALORES_ESPERADOS['distribucion_strands']['tolerancia']:
        print(f"  • Estado:                ✅ EQUILIBRADO")
        balance_status = "EQUILIBRADO"
    else:
        print(f"  • Estado:                ⚠️ DESBALANCEADO")
        balance_status = "DESBALANCEADO"
    
    validaciones['analisis_adicional'] = {
        'balance_hebras': {
            'strand_plus_pct': round(strand_plus_pct, 2),
            'strand_minus_pct': round(strand_minus_pct, 2),
            'desbalance': round(desbalance, 2),
            'status': balance_status,
            'interpretacion': VALORES_ESPERADOS['distribucion_strands']['interpretacion']
        }
    }
    
    # ========================================================================
    # COMPARACIÓN CON OTRAS BACTERIAS
    # ========================================================================
    
    print("\n" + "="*80)
    print("COMPARACIÓN CON OTRAS BACTERIAS")
    print("="*80)
    
    comparaciones = [
        {'bacteria': 'E. coli K-12', 'gc_pct': gc_obtenido, 'genes': genome_data['cds']['total'], 'tamano_mb': tamano_obtenido/1e6},
        {'bacteria': 'Salmonella (típico)', 'gc_pct': 52.2, 'genes': 4500, 'tamano_mb': 4.9},
        {'bacteria': 'B. subtilis (típico)', 'gc_pct': 43.5, 'genes': 4100, 'tamano_mb': 4.2},
        {'bacteria': 'P. aeruginosa (típico)', 'gc_pct': 66.6, 'genes': 5570, 'tamano_mb': 6.3}
    ]
    
    print("\n{:<25} {:>10} {:>10} {:>12}".format("Bacteria", "GC%", "Genes", "Tamaño (Mb)"))
    print("-" * 60)
    for comp in comparaciones:
        print("{:<25} {:>10.2f} {:>10,} {:>12.2f}".format(
            comp['bacteria'], comp['gc_pct'], comp['genes'], comp['tamano_mb']
        ))
    
    validaciones['comparacion_bacterias'] = comparaciones
    
    # ========================================================================
    # GUARDAR RESULTADOS
    # ========================================================================
    
    print("\n" + "="*80)
    print("GUARDANDO RESULTADOS")
    print("="*80)
    
    # JSON completo
    json_file = os.path.join(OUTPUT_DIR, 'genome_validation_results.json')
    with open(json_file, 'w') as f:
        json.dump(validaciones, f, indent=2)
    print(f"✓ JSON guardado: {json_file}")
    
    # CSV resumen
    df_validacion = pd.DataFrame([
        {
            'Parámetro': val['parametro'],
            'Obtenido': val.get('obtenido', 'N/A'),
            'Esperado': val.get('esperado', f"{val.get('rango_min', 'N/A')}-{val.get('rango_max', 'N/A')}"),
            'Status': val['status']
        }
        for val in validaciones['validaciones'].values()
    ])
    
    csv_file = os.path.join(OUTPUT_DIR, 'genome_validation_summary.csv')
    df_validacion.to_csv(csv_file, index=False)
    print(f"✓ CSV guardado: {csv_file}")
    
    # ========================================================================
    # RESUMEN FINAL
    # ========================================================================
    
    print("\n" + "="*80)
    print("RESUMEN DE VALIDACIÓN")
    print("="*80)
    
    print("\n✅ VALIDACIONES EXITOSAS:")
    validaciones_exitosas = sum(1 for v in validaciones['validaciones'].values() if '✅' in v['status'])
    total_validaciones = len(validaciones['validaciones'])
    print(f"  • {validaciones_exitosas}/{total_validaciones} validaciones pasaron")
    
    print("\n📊 PARÁMETROS VALIDADOS:")
    for nombre, val in validaciones['validaciones'].items():
        print(f"  • {val['parametro']}: {val['status']}")
    
    print("\n📚 FUENTES BIBLIOGRÁFICAS:")
    fuentes_unicas = set()
    for param, config in VALORES_ESPERADOS.items():
        if 'fuente' in config:
            fuentes_unicas.add(config['fuente'])
    for i, fuente in enumerate(sorted(fuentes_unicas), 1):
        print(f"  {i}. {fuente}")
    
    print("\n" + "="*80)
    print("✅ VALIDACIÓN DE ESTADÍSTICAS GENÓMICAS COMPLETADA")
    print("="*80 + "\n")

if __name__ == '__main__':
    main()
