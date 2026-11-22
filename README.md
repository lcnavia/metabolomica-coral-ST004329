# Respuesta metabolómica de *Acropora hyacinthus* al estrés térmico agudo (27 °C vs 34 °C)

Proyecto final del curso **Metabolómica** (Maestría en Biotecnología, 2026-II).
Reanálisis del estudio público **ST004329** de Metabolomics Workbench.

**Autores:** [Laura Carolina Navia Alvarez] ([laura.naviaalvarez@u.icesi.edu.co])

---

## 1. Pregunta de investigación

¿El estrés térmico agudo severo (34 °C, 7 días) altera el perfil de metabolitos
identificados del holobionte de *Acropora hyacinthus* respecto al control (27 °C),
y qué rutas metabólicas se asocian a ese cambio?

## 2. Datos

| Campo | Valor |
|---|---|
| Repositorio | Metabolomics Workbench — Study **ST004329**, Project PR002743 |
| DOI | [10.21228/M83G2D](https://doi.org/10.21228/M83G2D) |
| Análisis | AN007214 (ESI+, HILIC) y AN007215 (ESI−, RPLC) |
| Muestra | Holobionte coralino completo (hospedero + Symbiodiniaceae), *n* = 6 por grupo |
| Diseño original | 3 especies × 3 temperaturas (27, 31, 34 °C) × 6 réplicas + 4 QC pool |
| Subconjunto analizado | *A. hyacinthus* 27 °C (AH27_1–6) vs 34 °C (AH34_1–6) |
| Variables | 20 metabolitos identificados (10 ESI+ y 10 ESI−), áreas de pico normalizadas por TIC |
| Licencia de los datos | Datos públicos de Metabolomics Workbench (NIH Common Fund); citar el DOI del estudio |

## 3. Estructura del repositorio

```
.
├── README.md
├── requirements.txt
├── LICENSE                        # MIT (aplica al código)
├── data/
│   ├── raw/                       # tablas "named_metabolite_data_working" (ESI+ / ESI−) del Workbench
│   ├── processed/
│   │   └── ST004329_for_metaboanalyst.csv   # matriz 12 muestras × 20 metabolitos (entrada a MetaboAnalyst)
│   └── metadata/                  # metadatos de muestras y diseño del estudio
├── scripts/
│   ├── utils.py                   # rutas, normalización log10 + Pareto, FDR
│   ├── 01_preprocesamiento.py     # unión ESI+/ESI−, QC (RSD), subconjunto, exportación
│   ├── 02_analisis_univariado.py  # normalización, t de Welch + FDR, fold change, volcano
│   ├── 03_analisis_multivariado.py# PCA, PLS-DA (LOOCV, permutaciones, VIP), heatmap
│   ├── 04_rutas_metabolicas.py    # resumen y figura del análisis de rutas KEGG
│   └── run_all.py                 # ejecuta todo el pipeline
├── metaboanalyst/
│   ├── figuras/                   # figuras exportadas de MetaboAnalyst 6.0 (usadas en el informe)
│   └── tablas/                    # name_map.csv y pathway_results.csv exportados de MetaboAnalyst
├── docs/
│   └── protocolo_metaboanalyst.md # parámetros exactos usados en la interfaz web
└── results/                       # se genera al ejecutar los scripts
    ├── tables/
    └── figures/
```

## 4. Cómo reproducir el análisis

Requisitos: Python ≥ 3.10.

```bash
git clone https://github.com/<usuario>/<repositorio>.git
cd <repositorio>
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python scripts/run_all.py          # ~30 s
```

Cada script también se puede correr por separado (`python scripts/02_analisis_univariado.py`).
Los resultados quedan en `results/tables/` (CSV) y `results/figures/` (PNG, 300 dpi).
La semilla aleatoria está fijada (`SEED = 2026` en `utils.py`), así que la validación
cruzada y las permutaciones dan los mismos números en cada ejecución.

La parte web (MetaboAnalyst 6.0) se reproduce siguiendo `docs/protocolo_metaboanalyst.md`
con el archivo `data/processed/ST004329_for_metaboanalyst.csv`.

## 5. Pipeline analítico

| Paso | Qué se hace | Por qué |
|---|---|---|
| 1. Preprocesamiento | Unión de ESI+ y ESI−; armonización de nombres; selección AH27 vs AH34 | Los datos depositados ya pasaron por XCMS (picos), CAMERA (aductos), QC-RLSC + ComBat (deriva) y normalización por TIC |
| 2. Control de calidad | RSD de cada metabolito en 4 QC pool | RSD ≤ 30 % indica medición reproducible (mediana 7.2 %, 20/20 aceptados) |
| 3. Faltantes | 0 faltantes / 0 ceros | No se requiere imputación ni filtrado |
| 4. Normalización | log10 + escalado de Pareto (sin normalización por muestra) | Reduce la asimetría y la dominancia de metabolitos muy abundantes sin inflar el ruido como el autoescalado |
| 5. Univariado | t de Welch, corrección FDR (Benjamini-Hochberg), fold change ≥ 1.5, volcano | Varianzas desiguales y *n* pequeño; control de falsos positivos |
| 6. Multivariado | PCA; PLS-DA con LOOCV, 1000 permutaciones (Q2) y VIP; heatmap (Ward) | Patrones globales, discriminación supervisada validada y metabolitos clave |
| 7. Interpretación | Análisis de rutas KEGG en MetaboAnalyst (global test + topología) | Contexto biológico de los cambios |

## 6. Resultados principales

* **QC:** RSD mediano 7.2 % (máx. 26.7 %).
* **PCA:** PC1 = 69.9 %, PC2 = 11.2 %. Tres réplicas de 34 °C (AH34_1–3) se separan
  claramente; AH34_4–6 se superponen con los controles.
* **Univariado:** ningún metabolito con FDR < 0.05. Con *p* crudo < 0.05: ácido
  β-muricólico (×10.5 a 34 °C, *p* = 0.009) y gluconolactona (×2.1, *p* = 0.035).
  14/20 metabolitos con FC ≥ 1.5 (10 aumentan y 4 disminuyen a 34 °C).
* **PLS-DA:** Q2 (LOOCV) = 0.46 con 2 componentes (permutación *p* = 0.010) y 0.75
  con 5 componentes (*p* = 0.004). VIP ≥ 1: ácido β-muricólico, agmatina, glutamato,
  treonina, histidina y ácido taurocólico (idénticos a MetaboAnalyst).
* **Rutas:** 29 rutas evaluadas; vía de las pentosas fosfato con *p* = 0.030
  (FDR = 0.27). Mayor impacto topológico: biosíntesis de Phe/Tyr/Trp (1.0),
  metabolismo de fenilalanina (0.44) y de alanina, aspartato y glutamato (0.43).

## 7. Correspondencia con MetaboAnalyst

Los scripts reproducen los números de la interfaz web: varianza de PCA (69.9 / 11.2 / 5.4 %),
p-valores de la prueba t de Welch y VIP del componente 1. Diferencias esperadas:

* **Validación cruzada:** MetaboAnalyst usa particiones aleatorias (k-fold); aquí se usa LOOCV,
  que es determinista. Las tendencias coinciden (Q2 sube con el número de componentes).
* **Permutaciones:** MetaboAnalyst usó como estadístico la separación en el entrenamiento, que
  con 12 muestras se satura (casi todas las permutaciones alcanzan el máximo y *p* = 0.26–0.93).
  Aquí el estadístico es Q2 de LOOCV, que no se satura; por eso es la prueba más informativa.

## 8. Limitaciones

* *n* = 6 por grupo y solo 20 metabolitos identificados (no la matriz completa de features).
* Señal del holobionte: no separa hospedero y simbionte.
* Las réplicas provienen de 3 colonias por especie; el efecto colonia no está en los metadatos.
* El análisis de rutas usa bibliotecas KEGG de organismos modelo; no hay biblioteca de corales.

## 9. Referencias de software

Python (pandas, NumPy, SciPy, scikit-learn, matplotlib, seaborn); MetaboAnalyst 6.0
(Pang et al., 2024, *Nucleic Acids Research*).

## 10. Cita del dataset

Tang, S. (2025). *Acute heat stress redirects coral carbon budgets: integrative evidence from
physiology and metabolomics*. Metabolomics Workbench, Study ST004329. https://doi.org/10.21228/M83G2D
