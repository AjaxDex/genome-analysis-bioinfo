import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
import os
import time

# Configuración
email = '215788@unsaac.edu.pe'
OUTPUT_DIR = os.path.expanduser('~/data/raw')

# Crear directorio si no existe
os.makedirs(OUTPUT_DIR, exist_ok=True)

def buscar_ids(db, term, retmax=5):
    """Buscar IDs en NCBI Entrez"""
    params = {
        'db': db,
        'term': term,
        'retmax': str(retmax),
        'retmode': 'xml',
        'email': email
    }
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?' + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url) as response:
        xml_data = response.read()
    root = ET.fromstring(xml_data)
    ids = [id_elem.text for id_elem in root.findall('./IdList/Id')]
    return ids

def descargar_genbank(db, id, output_file):
    """Descargar archivo GenBank completo"""
    params = {
        'db': db,
        'id': id,
        'rettype': 'gb',
        'retmode': 'text',
        'email': email
    }
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urllib.parse.urlencode(params)
    
    print(f"Descargando desde: {url}")
    with urllib.request.urlopen(url) as response:
        data = response.read().decode('utf-8')
    
    # Guardar archivo
    with open(output_file, 'w') as f:
        f.write(data)
    
    return data

def validar_genbank(filepath):
    """Validar que el archivo GenBank se descargó correctamente"""
    with open(filepath, 'r') as f:
        contenido = f.read()
    
    # Extraer información básica
    lineas = contenido.split('\n')
    info = {}
    
    for linea in lineas[:50]:  # Revisar primeras 50 líneas
        if linea.startswith('DEFINITION'):
            info['definicion'] = linea.replace('DEFINITION', '').strip()
        if linea.startswith('ACCESSION'):
            info['accession'] = linea.replace('ACCESSION', '').strip()
        if 'source' in linea.lower() and 'bp' in linea.lower():
            # Intentar extraer tamaño
            partes = linea.split()
            for parte in partes:
                if 'bp' in parte:
                    info['tamano'] = parte
    
    return info, len(contenido)

# ====== EJECUCIÓN PRINCIPAL ======

print("="*60)
print("DESCARGA DE E. coli K-12 MG1655 - FORMATO GENBANK")
print("="*60)

# Buscar específicamente K-12 MG1655
print("\n[1/4] Buscando E. coli K-12 MG1655 en NCBI...")
term = 'Escherichia coli str. K-12 substr. MG1655[Organism] AND complete genome'
ids = buscar_ids('nucleotide', term, retmax=1)

if not ids:
    print("❌ No se encontraron resultados. Intentando con accession directo...")
    # Buscar por accession conocido
    ids = ['U00096.3']  # Accession number oficial de K-12 MG1655

print(f"✓ ID encontrado: {ids[0]}")

# Descargar GenBank
print("\n[2/4] Descargando archivo GenBank...")
output_file = os.path.join(OUTPUT_DIR, 'ecoli_k12_mg1655.gbk')
gb_data = descargar_genbank('nucleotide', ids[0], output_file)

print(f"✓ Archivo guardado en: {output_file}")
print(f"✓ Tamaño del archivo: {len(gb_data)} caracteres")

# Validar descarga
print("\n[3/4] Validando archivo descargado...")
info, tamano = validar_genbank(output_file)

print("\n" + "="*60)
print("INFORMACIÓN DEL GENOMA DESCARGADO")
print("="*60)
for clave, valor in info.items():
    print(f"{clave.upper()}: {valor}")
print(f"TAMAÑO ARCHIVO: {tamano:,} caracteres")

# Mostrar primeras líneas
print("\n[4/4] Primeras líneas del archivo GenBank:")
print("-"*60)
with open(output_file, 'r') as f:
    for i, linea in enumerate(f):
        if i < 20:
            print(linea.rstrip())
        else:
            break
print("-"*60)

print("\n✅ DESCARGA COMPLETADA EXITOSAMENTE")
print(f"📁 Archivo ubicado en: {output_file}")
print("\nPróximos pasos:")
print("  1. Verificar el contenido con: micro ~/data/raw/ecoli_k12_mg1655.gbk")
print("  2. Instalar librerías Python: pip install biopython pandas matplotlib")
print("  3. Desarrollar scripts de análisis")

time.sleep(1)
