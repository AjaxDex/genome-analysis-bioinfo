#!/usr/bin/env python3
"""
Descarga E. coli K-12 MG1655 con SECUENCIA COMPLETA
"""

from Bio import Entrez, SeqIO
import os
import time

# Configuración
Entrez.email = '215788@unsaac.edu.pe'
OUTPUT_DIR = os.path.expanduser('~/data/raw')
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("="*70)
print("DESCARGA E. coli K-12 MG1655 CON BIOPYTHON (SECUENCIA COMPLETA)")
print("="*70)

# Accession número oficial de K-12 MG1655
accession = 'U00096.3'  # O NC_000913.3 (son equivalentes)

print(f"\n[1/3] Descargando {accession} desde NCBI...")
print("Esto puede tardar 1-2 minutos para el genoma completo...\n")

try:
    # Descargar usando Entrez.efetch
    handle = Entrez.efetch(
        db='nucleotide',
        id=accession,
        rettype='gbwithparts',  # GenBank con secuencia completa
        retmode='text'
    )
    
    # Leer el registro
    record = SeqIO.read(handle, 'genbank')
    handle.close()
    
    print(f"✓ Descarga completada")
    print(f"  ID: {record.id}")
    print(f"  Descripción: {record.description[:60]}...")
    print(f"  Tamaño: {len(record.seq):,} bp")
    print(f"  Features: {len(record.features)}")
    
    # Verificar secuencia
    print("\n[2/3] Verificando secuencia...")
    seq_str = str(record.seq)
    print(f"✓ Secuencia accesible")
    print(f"  Primeros 60 bp: {seq_str[:60]}")
    print(f"  Últimos 60 bp: {seq_str[-60:]}")
    
    # Guardar archivo
    print(f"\n[3/3] Guardando archivo...")
    output_file = os.path.join(OUTPUT_DIR, 'ecoli_k12_mg1655.gbk')
    
    # Hacer backup si ya existe
    if os.path.exists(output_file):
        backup = output_file + '.backup'
        os.rename(output_file, backup)
        print(f"  Backup creado: {backup}")
    
    SeqIO.write(record, output_file, 'genbank')
    print(f"✓ Archivo guardado: {output_file}")
    
    # Estadísticas finales
    print("\n" + "="*70)
    print("DESCARGA COMPLETADA EXITOSAMENTE")
    print("="*70)
    print(f"Archivo: {output_file}")
    print(f"Tamaño genoma: {len(record.seq):,} bp")
    print(f"Features anotados: {len(record.features):,}")
    
    # Contar tipos de features
    feature_types = {}
    for feat in record.features:
        feature_types[feat.type] = feature_types.get(feat.type, 0) + 1
    
    print("\nTipos de anotaciones:")
    for ftype, count in sorted(feature_types.items()):
        print(f"  {ftype:15} : {count:5,}")
    
    print("="*70)
    
except Exception as e:
    print(f"\n❌ ERROR: {e}")
    print("\nIntentando con NC_000913.3...")
    
    try:
        handle = Entrez.efetch(
            db='nucleotide',
            id='NC_000913.3',
            rettype='gbwithparts',
            retmode='text'
        )
        record = SeqIO.read(handle, 'genbank')
        handle.close()
        
        output_file = os.path.join(OUTPUT_DIR, 'ecoli_k12_mg1655.gbk')
        SeqIO.write(record, output_file, 'genbank')
        
        print(f"✓ Descarga exitosa con NC_000913.3")
        print(f"✓ Tamaño: {len(record.seq):,} bp")
        print(f"✓ Guardado en: {output_file}")
        
    except Exception as e2:
        print(f"❌ ERROR CRÍTICO: {e2}")

time.sleep(1)
