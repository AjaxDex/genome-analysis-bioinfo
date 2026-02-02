#!/usr/bin/env python3
"""
Verificación de librerías instaladas para el proyecto E. coli
"""

print("="*60)
print("VERIFICACIÓN DE LIBRERÍAS INSTALADAS")
print("="*60)

libs = {
    'Bio': 'BioPython',
    'pandas': 'Pandas',
    'matplotlib': 'Matplotlib',
    'seaborn': 'Seaborn',
    'numpy': 'NumPy'
}

for module, nombre in libs.items():
    try:
        mod = __import__(module)
        version = getattr(mod, '__version__', 'desconocida')
        print(f"✅ {nombre:15} - versión {version}")
    except ImportError:
        print(f"❌ {nombre:15} - NO INSTALADO")

print("="*60)

# Test rápido de BioPython con el archivo descargado
print("\nTest de BioPython con el genoma descargado...")
try:
    from Bio import SeqIO
    import os
    
    gbk_file = os.path.expanduser('~/data/raw/ecoli_k12_mg1655.gbk')
    
    if os.path.exists(gbk_file):
        record = SeqIO.read(gbk_file, 'genbank')
        print(f"✅ Archivo GenBank leído correctamente")
        print(f"   ID: {record.id}")
        print(f"   Descripción: {record.description[:60]}...")
        print(f"   Tamaño: {len(record.seq):,} bp")
        print(f"   Features anotados: {len(record.features)}")
    else:
        print(f"❌ Archivo no encontrado: {gbk_file}")
        
except Exception as e:
    print(f"❌ Error al leer GenBank: {e}")

print("\n✅ VERIFICACIÓN COMPLETADA")
