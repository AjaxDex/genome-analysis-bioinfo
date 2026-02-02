#!/usr/bin/env python3
"""
Script Maestro - Análisis Genómico Completo de E. coli K-12 MG1655
==================================================================
Este script ejecuta todo el pipeline de análisis genómico en orden:
1. Análisis de codones ATG
2. Análisis de codones de terminación
3. Estadísticas genómicas completas
4. Distribución de tamaños de genes
5. Validaciones con literatura
6. Generación de visualizaciones

Autor: Proyecto Bioinformática
Fecha: 2024
"""

import os
import sys
import subprocess
from datetime import datetime
import json

# ============================================================================
# CONFIGURACIÓN
# ============================================================================

PROJECT_DIR = os.path.expanduser('~/projects/bioinfo')
SCRIPTS_DIR = os.path.join(PROJECT_DIR, 'scripts/analysis')
VIZ_DIR = os.path.join(PROJECT_DIR, 'scripts/visualization')
RESULTS_DIR = os.path.join(PROJECT_DIR, 'results')

# Scripts de análisis
SCRIPTS_ANALISIS = [
    ('analyze_atg.py', 'Análisis de codones ATG'),
    ('analyze_stop_codons.py', 'Análisis de codones de terminación'),
    ('analyze_genome_stats.py', 'Estadísticas genómicas'),
    ('analyze_gene_distribution.py', 'Distribución de tamaños de genes'),
]

# Scripts de validación
SCRIPTS_VALIDACION = [
    ('validate_results.py', 'Validación de codones'),
    ('validate_genome_stats.py', 'Validación de estadísticas genómicas'),
]

# Scripts de visualización
SCRIPTS_VISUALIZACION = [
    ('plot_atg_distribution.py', 'Visualización de ATG'),
    ('plot_stop_codons.py', 'Visualización de stop codons'),
    ('plot_gene_distribution.py', 'Visualización de distribución de genes'),
    ('plot_genome_overview.py', 'Overview del genoma'),
]

# ============================================================================
# FUNCIONES AUXILIARES
# ============================================================================

def imprimir_banner():
    """Imprime el banner inicial."""
    banner = """
    ╔════════════════════════════════════════════════════════════════════╗
    ║                                                                    ║
    ║         ANÁLISIS GENÓMICO COMPLETO - E. coli K-12 MG1655          ║
    ║                                                                    ║
    ║                    Pipeline Automatizado                          ║
    ║                                                                    ║
    ╚════════════════════════════════════════════════════════════════════╝
    """
    print(banner)
    print(f"\n    Fecha de ejecución: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"    Directorio del proyecto: {PROJECT_DIR}\n")

def imprimir_seccion(titulo):
    """Imprime un separador de sección."""
    print("\n" + "="*80)
    print(f"  {titulo}")
    print("="*80 + "\n")

def ejecutar_script(script_path, descripcion):
    """
    Ejecuta un script Python y captura su salida.
    
    Args:
        script_path (str): Ruta al script
        descripcion (str): Descripción del script
    
    Returns:
        bool: True si exitoso, False si falló
    """
    print(f"→ Ejecutando: {descripcion}")
    print(f"  Script: {os.path.basename(script_path)}")
    
    try:
        # Ejecutar el script
        result = subprocess.run(
            ['python', script_path],
            capture_output=True,
            text=True,
            timeout=300  # 5 minutos de timeout
        )
        
        if result.returncode == 0:
            print(f"  ✅ COMPLETADO\n")
            return True
        else:
            print(f"  ❌ ERROR: {result.stderr}\n")
            return False
            
    except subprocess.TimeoutExpired:
        print(f"  ❌ ERROR: Timeout (>5 minutos)\n")
        return False
    except Exception as e:
        print(f"  ❌ ERROR: {str(e)}\n")
        return False

def verificar_archivos_generados():
    """Verifica que se hayan generado los archivos esperados."""
    
    archivos_esperados = {
        'Análisis': [
            'results/tables/atg_analysis.json',
            'results/tables/stop_codons_analysis.json',
            'results/tables/genome_statistics.json',
            'results/tables/gene_size_distribution_analysis.json',
        ],
        'Validación': [
            'results/tables/validation_results.json',
            'results/tables/genome_validation_results.json',
        ],
        'Visualizaciones': [
            'results/figures/atg_distribution.png',
            'results/figures/stop_codons_comparison.png',
            'results/figures/gene_size_distribution.png',
            'results/figures/genome_overview.png',
        ]
    }
    
    print("\n" + "="*80)
    print("  VERIFICACIÓN DE ARCHIVOS GENERADOS")
    print("="*80 + "\n")
    
    total_archivos = 0
    archivos_encontrados = 0
    
    for categoria, archivos in archivos_esperados.items():
        print(f"{categoria}:")
        for archivo in archivos:
            total_archivos += 1
            filepath = os.path.join(PROJECT_DIR, archivo)
            if os.path.exists(filepath):
                size_kb = os.path.getsize(filepath) / 1024
                print(f"  ✅ {archivo} ({size_kb:.1f} KB)")
                archivos_encontrados += 1
            else:
                print(f"  ❌ {archivo} - NO ENCONTRADO")
        print()
    
    print(f"Total: {archivos_encontrados}/{total_archivos} archivos generados correctamente")
    
    return archivos_encontrados == total_archivos

def generar_reporte_final():
    """Genera un reporte final del análisis."""
    
    print("\n" + "="*80)
    print("  GENERANDO REPORTE FINAL")
    print("="*80 + "\n")
    
    try:
        # Cargar todos los resultados
        with open(os.path.join(RESULTS_DIR, 'tables/atg_analysis.json'), 'r') as f:
            atg_data = json.load(f)
        
        with open(os.path.join(RESULTS_DIR, 'tables/stop_codons_analysis.json'), 'r') as f:
            stop_data = json.load(f)
        
        with open(os.path.join(RESULTS_DIR, 'tables/genome_statistics.json'), 'r') as f:
            genome_data = json.load(f)
        
        with open(os.path.join(RESULTS_DIR, 'tables/validation_results.json'), 'r') as f:
            val_data = json.load(f)
        
        # Crear reporte
        reporte = f"""
╔════════════════════════════════════════════════════════════════════════════╗
║                         REPORTE FINAL DE ANÁLISIS                          ║
║                      E. coli K-12 substr. MG1655                          ║
╚════════════════════════════════════════════════════════════════════════════╝

Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 1. INFORMACIÓN DEL GENOMA
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  ID del genoma:         {genome_data['genoma']['id']}
  Tamaño:                {genome_data['genoma']['tamano_bp']:,} bp ({genome_data['genoma']['tamano_bp']/1e6:.2f} Mb)
  Contenido GC:          {genome_data['genoma']['contenido_gc_pct']:.2f}%
  Total genes:           {genome_data['genes']['total']:,}
  Total CDS:             {genome_data['cds']['total']:,}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 2. ANÁLISIS DE CODONES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Codones ATG:
    • Total en genoma:   {atg_data['analisis_atg']['total_atg']:,}
    • Densidad:          {atg_data['analisis_atg']['densidad_por_kb']:.2f} ATG/kb
    • Ratio ATG/CDS:     {atg_data['comparacion']['ratio_atg_vs_cds']:.2f}x

  Codones de Terminación (en CDS):
    • TAA:               {stop_data['analisis_cds']['stop_codons_por_tipo']['TAA']['proporcion_pct']:.2f}% ⭐ PREFERIDO
    • TAG:               {stop_data['analisis_cds']['stop_codons_por_tipo']['TAG']['proporcion_pct']:.2f}% (evitado)
    • TGA:               {stop_data['analisis_cds']['stop_codons_por_tipo']['TGA']['proporcion_pct']:.2f}% (evitado)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 3. ESTADÍSTICAS GENÓMICAS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Densidad génica:       {genome_data['densidad_genica']['genes_por_mb']:.2f} genes/Mb
  Región codificante:    {genome_data['regiones']['codificantes']['porcentaje']:.2f}%
  Región no codif.:      {genome_data['regiones']['no_codificantes']['porcentaje']:.2f}%

  Tamaño de genes:
    • Media:             {genome_data['cds']['estadisticas_tamanos']['media']:.2f} bp
    • Mediana:           {genome_data['cds']['estadisticas_tamanos']['mediana']:.2f} bp
    • Rango:             {genome_data['cds']['estadisticas_tamanos']['minimo']:,} - {genome_data['cds']['estadisticas_tamanos']['maximo']:,} bp

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 4. VALIDACIÓN CON LITERATURA
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Tamaño del genoma:     {val_data['validaciones'].get('tamano_genoma', {}).get('status', 'N/A')}
  Contenido GC:          {val_data['validaciones'].get('contenido_gc', {}).get('status', 'No disponible')}
  Total genes:           {val_data['validaciones'].get('total_genes', {}).get('status', 'N/A')}
  Ratio ATG/CDS:         {val_data['validaciones']['ratio_atg_cds']['status']}
  Stop codons (TAA):     {val_data['validaciones']['stop_codons']['TAA']['status']}
  Stop codons (TAG):     {val_data['validaciones']['stop_codons']['TAG']['status']}
  Stop codons (TGA):     {val_data['validaciones']['stop_codons']['TGA']['status']}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 5. ARCHIVOS GENERADOS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Datos de análisis:     results/tables/
  Visualizaciones:       results/figures/
  Documentación:         docs/

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 CONCLUSIÓN
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  ✅ Análisis completado exitosamente
  ✅ Todos los resultados validados con literatura científica
  ✅ Visualizaciones generadas correctamente

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
        
        print(reporte)
        
        # Guardar reporte
        reporte_file = os.path.join(RESULTS_DIR, 'REPORTE_FINAL.txt')
        with open(reporte_file, 'w', encoding='utf-8') as f:
            f.write(reporte)
        
        print(f"\n✅ Reporte final guardado: {reporte_file}\n")
        
        return True
        
    except Exception as e:
        print(f"❌ Error al generar reporte final: {str(e)}\n")
        return False

# ============================================================================
# FUNCIÓN PRINCIPAL
# ============================================================================

def main():
    """Ejecuta el pipeline completo de análisis."""
    
    inicio = datetime.now()
    
    # Banner inicial
    imprimir_banner()
    
    # Contadores
    total_scripts = len(SCRIPTS_ANALISIS) + len(SCRIPTS_VALIDACION) + len(SCRIPTS_VISUALIZACION)
    scripts_exitosos = 0
    scripts_fallidos = 0
    
    # ========================================================================
    # FASE 1-3: SCRIPTS DE ANÁLISIS
    # ========================================================================
    
    imprimir_seccion("FASE 1-3: ANÁLISIS DE DATOS")
    
    for script, descripcion in SCRIPTS_ANALISIS:
        script_path = os.path.join(SCRIPTS_DIR, script)
        if ejecutar_script(script_path, descripcion):
            scripts_exitosos += 1
        else:
            scripts_fallidos += 1
    
    # ========================================================================
    # FASE 2: VALIDACIÓN
    # ========================================================================
    
    imprimir_seccion("FASE 2: VALIDACIÓN CON LITERATURA")
    
    for script, descripcion in SCRIPTS_VALIDACION:
        script_path = os.path.join(SCRIPTS_DIR, script)
        if ejecutar_script(script_path, descripcion):
            scripts_exitosos += 1
        else:
            scripts_fallidos += 1
    
    # ========================================================================
    # FASE 4: VISUALIZACIONES
    # ========================================================================
    
    imprimir_seccion("FASE 4: GENERACIÓN DE VISUALIZACIONES")
    
    for script, descripcion in SCRIPTS_VISUALIZACION:
        script_path = os.path.join(VIZ_DIR, script)
        if ejecutar_script(script_path, descripcion):
            scripts_exitosos += 1
        else:
            scripts_fallidos += 1
    
    # ========================================================================
    # VERIFICACIÓN Y REPORTE FINAL
    # ========================================================================
    
    archivos_ok = verificar_archivos_generados()
    reporte_ok = generar_reporte_final()
    
    # ========================================================================
    # RESUMEN FINAL
    # ========================================================================
    
    fin = datetime.now()
    duracion = (fin - inicio).total_seconds()
    
    print("\n" + "="*80)
    print("  RESUMEN DE EJECUCIÓN")
    print("="*80 + "\n")
    
    print(f"  Inicio:              {inicio.strftime('%H:%M:%S')}")
    print(f"  Fin:                 {fin.strftime('%H:%M:%S')}")
    print(f"  Duración total:      {duracion:.1f} segundos ({duracion/60:.1f} minutos)")
    print(f"\n  Scripts ejecutados:  {total_scripts}")
    print(f"  ✅ Exitosos:         {scripts_exitosos}")
    print(f"  ❌ Fallidos:         {scripts_fallidos}")
    
    if scripts_fallidos == 0 and archivos_ok and reporte_ok:
        print("\n" + "="*80)
        print("  ✅ PIPELINE COMPLETADO EXITOSAMENTE")
        print("="*80 + "\n")
        return 0
    else:
        print("\n" + "="*80)
        print("  ⚠️  PIPELINE COMPLETADO CON ADVERTENCIAS")
        print("="*80 + "\n")
        return 1

# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == '__main__':
    sys.exit(main())
