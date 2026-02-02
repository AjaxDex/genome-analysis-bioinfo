# Análisis Genómico de *E. coli* K-12 MG1655

## 📋 Descripción del Proyecto

Proyecto de bioinformática que realiza un análisis computacional completo del genoma de *Escherichia coli* K-12 cepa MG1655, incluyendo análisis de codones, estadísticas genómicas, y visualizaciones profesionales.

**Institución:** Proyecto Académico de Bioinformática  
**Fecha:** 2024  
**Genoma Analizado:** NC_000913.3 (*E. coli* K-12 substr. MG1655)

---

## 🎯 Objetivos

1. ✅ Análisis completo de codones de inicio (ATG) y terminación (TAA, TAG, TGA)
2. ✅ Cálculo de estadísticas genómicas (GC%, densidad génica, distribución de genes)
3. ✅ Validación de resultados con literatura científica
4. ✅ Generación de visualizaciones profesionales
5. ✅ Documentación reproducible del análisis

---

## 🧬 Resultados Principales

### Estadísticas del Genoma

| Parámetro | Valor | Validación |
|-----------|-------|------------|
| **Tamaño** | 4,641,652 bp (4.64 Mb) | ✅ Exacto |
| **Contenido GC** | 50.79% | ✅ Esperado: ~50.8% |
| **Total de genes** | 4,651 | ✅ Válido |
| **Total CDS** | 4,318 | ✅ Esperado: ~4,300 |
| **Densidad génica** | 930.27 genes/Mb | ✅ Válido |
| **Región codificante** | 87.23% | ✅ Típico de bacterias |

### Análisis de Codones

**Codones ATG (Inicio):**
- Total en genoma: **76,282**
- Densidad: **16.43 ATG/kb**
- Ratio ATG/CDS: **17.67x** (mayoría no son inicios funcionales)

**Codones de Terminación (en CDS):**
- **TAA:** 63.97% ⭐ **Fuertemente preferido** (1.67x enriquecimiento)
- **TAG:** 7.09% ❌ **Evitado** (0.47x enriquecimiento)
- **TGA:** 28.94% ⚠️ **Menos preferido** (0.62x enriquecimiento)

### Hallazgos Científicos

1. 🧬 **Preferencia evolutiva clara por TAA** como codón de terminación
2. 📊 **Genoma altamente compacto** (87.23% codificante)
3. ⚖️ **Balance perfecto entre hebras** (48.70% vs 51.30%)
4. 📏 **Tamaño medio de genes:** 937 bp (~312 aminoácidos)

---

## 📁 Estructura del Proyecto

```
bioinfo/
├── data/
│   └── raw/
│       └── ecoli_k12_mg1655.gbk          # Genoma descargado
│
├── scripts/
│   ├── analysis/
│   │   ├── analyze_atg.py                # Análisis de ATG
│   │   ├── analyze_stop_codons.py        # Análisis de stop codons
│   │   ├── analyze_genome_stats.py       # Estadísticas genómicas
│   │   ├── analyze_gene_distribution.py  # Distribución de genes
│   │   ├── validate_results.py           # Validación de codones
│   │   └── validate_genome_stats.py      # Validación de estadísticas
│   │
│   ├── visualization/
│   │   ├── plot_atg_distribution.py      # Gráficos de ATG
│   │   ├── plot_stop_codons.py           # Gráficos de stop codons
│   │   ├── plot_gene_distribution.py     # Gráficos de distribución
│   │   └── plot_genome_overview.py       # Dashboard general
│   │
│   ├── utils/
│   │   ├── __init__.py
│   │   └── genome_utils.py               # Funciones reutilizables
│   │
│   ├── download_ecoli.py                 # Descarga del genoma
│   └── run_full_analysis.py              # Pipeline completo
│
├── results/
│   ├── tables/                           # Resultados en CSV/JSON
│   ├── figures/                          # Gráficos generados
│   └── REPORTE_FINAL.txt                 # Reporte consolidado
│
├── docs/
│   ├── FASE1_Analisis_Codones.md         # Documentación Fase 1
│   └── metodologia.md                    # Metodología detallada
│
├── README.md                              # Este archivo
└── requirements.txt                       # Dependencias Python
```

---

## 🚀 Instalación y Uso

### Requisitos Previos

- Python 3.8+
- pip (gestor de paquetes)
- Conexión a internet (para descargar el genoma)

### 1. Clonar el Repositorio

```bash
git clone https://github.com/AjaxDex/genome-analysis-bioinfo.git
cd ecoli-genome-analysis
```

### 2. Crear Entorno Virtual

```bash
python3 -m venv bioinfo
source bioinfo/bin/activate  # En Linux/Mac
# bioinfo\Scripts\activate   # En Windows
```

### 3. Instalar Dependencias

```bash
pip install -r requirements.txt
```

### 4. Descargar el Genoma

```bash
cd scripts
python download_ecoli.py
```

### 5. Ejecutar Análisis Completo

```bash
python run_full_analysis.py
```

O ejecutar análisis individuales:

```bash
# Análisis de codones ATG
python scripts/analysis/analyze_atg.py

# Análisis de stop codons
python scripts/analysis/analyze_stop_codons.py

# Estadísticas genómicas
python scripts/analysis/analyze_genome_stats.py

# Generar visualizaciones
python scripts/visualization/plot_genome_overview.py
```

---

## 📊 Visualizaciones Generadas

El proyecto genera 8 gráficos profesionales:

1. **atg_distribution.png** - Comparación de codones ATG funcionales vs no funcionales
2. **atg_density.png** - Densidad de ATG en el genoma
3. **stop_codons_comparison.png** - Comparación de proporciones de stop codons
4. **stop_codons_pie.png** - Distribución en gráficos de pastel
5. **gene_size_distribution.png** - Histogramas de tamaños de genes
6. **gene_size_violin.png** - Violin plots de distribución
7. **extreme_genes.png** - Genes más pequeños y más grandes
8. **genome_overview.png** - Dashboard completo del genoma

---

## 🔬 Metodología

### Herramientas Utilizadas

- **BioPython (1.86):** Parsing del genoma GenBank
- **Pandas (3.0.0):** Análisis de datos tabulares
- **NumPy (2.4.2):** Cálculos estadísticos
- **Matplotlib (3.10.8):** Visualizaciones base
- **Seaborn (0.13.2):** Visualizaciones avanzadas

### Pipeline de Análisis

1. **Descarga del genoma** desde NCBI usando API Entrez
2. **Extracción de features** (genes, CDS) del archivo GenBank
3. **Análisis de codones** mediante búsqueda de patrones
4. **Cálculos estadísticos** (media, mediana, percentiles, IQR)
5. **Validación** comparando con valores de literatura
6. **Visualización** con gráficos profesionales

### Control de Calidad

- ✅ Validación con valores publicados en literatura científica
- ✅ Verificación de integridad de datos
- ✅ Código documentado y reproducible
- ✅ Manejo de errores y excepciones

---

## 📚 Referencias Científicas

1. **Blattner et al. (1997)** - "The Complete Genome Sequence of *Escherichia coli* K-12" - *Science* 277(5331):1453-1462

2. **Riley et al. (2006)** - "*Escherichia coli* K-12: a cooperatively developed annotation snapshot" - *Nucleic Acids Research* 34(1):1-9

3. **Nakamura et al. (2000)** - "Codon usage tabulated from international DNA sequence databases" - *Nucleic Acids Research* 28(1):292

4. **Sharp et al. (2010)** - "Variation in the strength of selected codon usage bias among bacteria" - *Nucleic Acids Research* 33(4):1141-1153

5. **NCBI RefSeq** - NC_000913.3 (*E. coli* str. K-12 substr. MG1655)

---

## 👥 Autores

**Proyecto Bioinformática 2024**

---

## 📄 Licencia

Este proyecto es de código abierto para fines educativos y de investigación.

---

## 🙏 Agradecimientos

- NCBI por proporcionar acceso a las bases de datos genómicas
- Comunidad de BioPython por las herramientas de análisis
- Literatura científica que permitió la validación de resultados

---

## 📧 Contacto

Para preguntas o sugerencias sobre este proyecto, por favor abre un issue en GitHub.

---

**Última actualización:** 2026
