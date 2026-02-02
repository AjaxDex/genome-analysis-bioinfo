from Bio import SeqIO
import os

gbk_file = os.path.expanduser('~/data/raw/ecoli_k12_mg1655.gbk')
record = SeqIO.read(gbk_file, 'genbank')

print(f"ID: {record.id}")
print(f"Descripción: {record.description}")
print(f"Número de features: {len(record.features)}")
print(f"Tipo de secuencia: {type(record.seq)}")
print(f"¿Secuencia definida?: {hasattr(record.seq, '_data')}")

try:
    seq_str = str(record.seq)
    print(f"Tamaño secuencia: {len(seq_str)} bp")
    print(f"Primeros 100 bp: {seq_str[:100]}")
except Exception as e:
    print(f"ERROR al acceder a secuencia: {e}")
