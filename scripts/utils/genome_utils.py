#!/usr/bin/env python3
"""
Utilidades para Análisis Genómico
==================================
Funciones reutilizables para análisis de genomas bacterianos.

Autor: Proyecto Bioinformática
Fecha: 2024
"""

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import numpy as np
from collections import Counter

# ============================================================================
# FUNCIONES DE ANÁLISIS DE SECUENCIAS
# ============================================================================

def contar_codon(secuencia, codon):
    """
    Cuenta todas las apariciones de un codón en una secuencia.
    
    Args:
        secuencia (str): Secuencia de ADN
        codon (str): Codón a buscar (ej: 'ATG')
    
    Returns:
        list: Lista de posiciones donde aparece el codón (0-indexed)
    """
    secuencia = str(secuencia).upper()
    codon = codon.upper()
    posiciones = []
    
    for i in range(len(secuencia) - len(codon) + 1):
        if secuencia[i:i+len(codon)] == codon:
            posiciones.append(i)
    
    return posiciones

def calcular_contenido_gc(secuencia):
    """
    Calcula el contenido GC de una secuencia.
    
    Args:
        secuencia: Objeto Seq de BioPython o string
    
    Returns:
        float: Porcentaje de GC
    """
    try:
        return gc_fraction(secuencia) * 100
    except:
        seq_str = str(secuencia).upper()
        g_count = seq_str.count('G')
        c_count = seq_str.count('C')
        total = len(seq_str)
        return ((g_count + c_count) / total * 100) if total > 0 else 0

def extraer_features_por_tipo(record, feature_type):
    """
    Extrae todos los features de un tipo específico.
    
    Args:
        record: Objeto SeqRecord de BioPython
        feature_type (str): Tipo de feature ('gene', 'CDS', etc.)
    
    Returns:
        list: Lista de features del tipo especificado
    """
    return [f for f in record.features if f.type == feature_type]

def obtener_secuencia_feature(record, feature):
    """
    Obtiene la secuencia de un feature.
    
    Args:
        record: Objeto SeqRecord de BioPython
        feature: Feature de BioPython
    
    Returns:
        Seq: Secuencia extraída
    """
    return feature.extract(record.seq)

# ============================================================================
# FUNCIONES DE ESTADÍSTICAS
# ============================================================================

def calcular_estadisticas_descriptivas(valores):
    """
    Calcula estadísticas descriptivas completas.
    
    Args:
        valores (list): Lista de valores numéricos
    
    Returns:
        dict: Diccionario con estadísticas
    """
    arr = np.array(valores)
    
    if len(arr) == 0:
        return {
            'total': 0,
            'media': 0,
            'mediana': 0,
            'desviacion_std': 0,
            'minimo': 0,
            'maximo': 0
        }
    
    return {
        'total': len(arr),
        'media': float(np.mean(arr)),
        'mediana': float(np.median(arr)),
        'moda': float(Counter(arr).most_common(1)[0][0]) if len(arr) > 0 else 0,
        'desviacion_std': float(np.std(arr)),
        'varianza': float(np.var(arr)),
        'minimo': float(np.min(arr)),
        'maximo': float(np.max(arr)),
        'rango': float(np.max(arr) - np.min(arr)),
        'percentil_25': float(np.percentile(arr, 25)),
        'percentil_50': float(np.percentile(arr, 50)),
        'percentil_75': float(np.percentile(arr, 75)),
        'percentil_95': float(np.percentile(arr, 95)),
        'rango_intercuartil': float(np.percentile(arr, 75) - np.percentile(arr, 25))
    }

def calcular_densidad(count, tamano_genoma, unidad='kb'):
    """
    Calcula densidad por unidad genómica.
    
    Args:
        count (int): Número de elementos
        tamano_genoma (int): Tamaño del genoma en bp
        unidad (str): 'kb' o 'mb'
    
    Returns:
        float: Densidad calculada
    """
    if unidad == 'kb':
        return (count / tamano_genoma) * 1000
    elif unidad == 'mb':
        return (count / tamano_genoma) * 1_000_000
    else:
        raise ValueError("Unidad debe ser 'kb' o 'mb'")

def identificar_outliers_iqr(valores, multiplicador=1.5):
    """
    Identifica outliers usando el método IQR.
    
    Args:
        valores (list): Lista de valores
        multiplicador (float): Multiplicador del IQR (default 1.5)
    
    Returns:
        dict: Información de outliers
    """
    arr = np.array(valores)
    q1 = np.percentile(arr, 25)
    q3 = np.percentile(arr, 75)
    iqr = q3 - q1
    
    limite_inferior = q1 - multiplicador * iqr
    limite_superior = q3 + multiplicador * iqr
    
    outliers_bajos = arr[arr < limite_inferior]
    outliers_altos = arr[arr > limite_superior]
    
    return {
        'q1': float(q1),
        'q3': float(q3),
        'iqr': float(iqr),
        'limite_inferior': float(limite_inferior),
        'limite_superior': float(limite_superior),
        'num_outliers_bajos': len(outliers_bajos),
        'num_outliers_altos': len(outliers_altos),
        'total_outliers': len(outliers_bajos) + len(outliers_altos),
        'porcentaje_outliers': ((len(outliers_bajos) + len(outliers_altos)) / len(arr) * 100) if len(arr) > 0 else 0
    }

# ============================================================================
# FUNCIONES DE VALIDACIÓN
# ============================================================================

def validar_valor_esperado(obtenido, esperado, tolerancia):
    """
    Valida si un valor está dentro de la tolerancia esperada.
    
    Args:
        obtenido (float): Valor obtenido
        esperado (float): Valor esperado
        tolerancia (float): Tolerancia aceptable
    
    Returns:
        dict: Resultado de validación
    """
    diferencia = abs(obtenido - esperado)
    en_rango = diferencia <= tolerancia
    desviacion_pct = (diferencia / esperado) * 100 if esperado != 0 else 0
    
    return {
        'obtenido': round(obtenido, 2),
        'esperado': esperado,
        'diferencia': round(diferencia, 2),
        'desviacion_pct': round(desviacion_pct, 2),
        'tolerancia': tolerancia,
        'valido': en_rango,
        'status': '✅ VÁLIDO' if en_rango else '⚠️ FUERA DE RANGO'
    }

def validar_rango(obtenido, rango_min, rango_max):
    """
    Valida si un valor está dentro de un rango.
    
    Args:
        obtenido (float): Valor obtenido
        rango_min (float): Límite inferior
        rango_max (float): Límite superior
    
    Returns:
        dict: Resultado de validación
    """
    en_rango = rango_min <= obtenido <= rango_max
    
    return {
        'obtenido': round(obtenido, 2),
        'rango_min': rango_min,
        'rango_max': rango_max,
        'valido': en_rango,
        'status': '✅ DENTRO DEL RANGO' if en_rango else '⚠️ FUERA DEL RANGO'
    }

# ============================================================================
# FUNCIONES DE ARCHIVOS
# ============================================================================

def cargar_genoma_genbank(filepath):
    """
    Carga un archivo GenBank.
    
    Args:
        filepath (str): Ruta al archivo GenBank
    
    Returns:
        SeqRecord: Objeto record de BioPython
    """
    try:
        return SeqIO.read(filepath, 'genbank')
    except Exception as e:
        raise Exception(f"Error al cargar genoma: {str(e)}")

def guardar_json(data, filepath, indent=2):
    """
    Guarda datos en formato JSON.
    
    Args:
        data (dict): Datos a guardar
        filepath (str): Ruta del archivo de salida
        indent (int): Indentación del JSON
    """
    import json
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=indent)

def cargar_json(filepath):
    """
    Carga datos desde un archivo JSON.
    
    Args:
        filepath (str): Ruta del archivo JSON
    
    Returns:
        dict: Datos cargados
    """
    import json
    with open(filepath, 'r') as f:
        return json.load(f)

# ============================================================================
# FUNCIONES DE FORMATEO
# ============================================================================

def formatear_numero(numero, decimales=2):
    """
    Formatea un número con separadores de miles.
    
    Args:
        numero (float): Número a formatear
        decimales (int): Número de decimales
    
    Returns:
        str: Número formateado
    """
    if isinstance(numero, int):
        return f"{numero:,}"
    else:
        return f"{numero:,.{decimales}f}"

def crear_tabla_ascii(headers, rows):
    """
    Crea una tabla ASCII formateada.
    
    Args:
        headers (list): Lista de encabezados
        rows (list): Lista de filas (cada fila es una lista)
    
    Returns:
        str: Tabla formateada
    """
    # Calcular anchos de columnas
    col_widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            col_widths[i] = max(col_widths[i], len(str(cell)))
    
    # Crear separador
    separator = "+" + "+".join("-" * (w + 2) for w in col_widths) + "+"
    
    # Crear tabla
    table = [separator]
    
    # Encabezados
    header_row = "|" + "|".join(f" {h:<{col_widths[i]}} " for i, h in enumerate(headers)) + "|"
    table.append(header_row)
    table.append(separator)
    
    # Filas
    for row in rows:
        row_str = "|" + "|".join(f" {str(cell):<{col_widths[i]}} " for i, cell in enumerate(row)) + "|"
        table.append(row_str)
    
    table.append(separator)
    
    return "\n".join(table)

# ============================================================================
# FUNCIONES DE LOGGING
# ============================================================================

def imprimir_seccion(titulo, ancho=70):
    """
    Imprime un título de sección formateado.
    
    Args:
        titulo (str): Título de la sección
        ancho (int): Ancho de la línea
    """
    print("\n" + "="*ancho)
    print(titulo.center(ancho))
    print("="*ancho)

def imprimir_subseccion(titulo, ancho=70):
    """
    Imprime un título de subsección formateado.
    
    Args:
        titulo (str): Título de la subsección
        ancho (int): Ancho de la línea
    """
    print("\n" + "-"*ancho)
    print(titulo)
    print("-"*ancho)

def imprimir_progreso(paso_actual, total_pasos, descripcion):
    """
    Imprime el progreso de un proceso.
    
    Args:
        paso_actual (int): Paso actual
        total_pasos (int): Total de pasos
        descripcion (str): Descripción del paso
    """
    print(f"\n[{paso_actual}/{total_pasos}] {descripcion}")

# ============================================================================
# INFORMACIÓN DEL MÓDULO
# ============================================================================

__version__ = "1.0.0"
__author__ = "Proyecto Bioinformática"
__description__ = "Utilidades para análisis genómico de bacterias"

if __name__ == '__main__':
    print(f"Módulo de utilidades genómicas v{__version__}")
    print(f"Autor: {__author__}")
    print(f"\nEste módulo contiene funciones reutilizables para análisis genómico.")
    print("Importa este módulo en tus scripts para usar las funciones.")
