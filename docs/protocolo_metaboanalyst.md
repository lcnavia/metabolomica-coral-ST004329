# Protocolo en MetaboAnalyst 6.0 (https://www.metaboanalyst.ca)

Archivo de entrada: `data/processed/ST004329_for_metaboanalyst.csv`

## A. Statistical Analysis [one factor]

1. **Upload:** Data type = *Concentrations*; Format = *Samples in rows (unpaired)*.
2. **Data check:** 12 muestras, 20 variables, 2 grupos (Control_27C, Severe_34C); 0 valores faltantes.
3. **Missing values / filtering:** no aplica (sin faltantes; 20 variables, no se filtra).
4. **Normalization:**
   * Sample normalization: *None* (los datos ya están normalizados por TIC en el Workbench).
   * Data transformation: *Log transformation (base 10)*.
   * Data scaling: *Pareto scaling*.
   * Figura: `metaboanalyst/figuras/01_normalizacion_log10_pareto.png`.
5. **PCA:** scree plot y 2D scores (PC1 69.9 %, PC2 11.2 %). Figuras 02 y 03.
6. **t-test:** *Welch's t-test*, unpaired, umbral FDR 0.05 → 0 significativos. Figura 04.
7. **Volcano plot:** FC threshold 1.5 (dirección Control/Severe). Figura 05.
8. **PLS-DA:**
   * Cross validation: 2 y 5 componentes máximos (performance measure Q2). Figuras 06 y 07.
   * Permutation: 1000 permutaciones. Figuras 08 (p = 0.932) y 09 (p = 0.262).
   * VIP scores, componente 1. Figura 10.
9. **Heatmap:** datos normalizados, distancia euclidiana, clustering de Ward (parámetros por defecto). Figura 11.

## B. Pathway Analysis

1. Entrada: lista de 20 metabolitos con concentraciones (mismo CSV).
2. ID type: *Compound name*; mapeo 20/20 (ver `metaboanalyst/tablas/name_map.csv`).
3. Normalización: log10 + Pareto.
4. Parámetros: Enrichment = *Global test*; Topology = *Relative-betweenness centrality*;
   Pathway library = KEGG.
5. Resultados exportados: `metaboanalyst/tablas/pathway_results.csv`; figuras 12 y 13.
