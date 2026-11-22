"""
04_rutas_metabolicas.py - Interpretación biológica: análisis de rutas (KEGG).

El análisis de rutas se ejecutó en MetaboAnalyst 6.0 (módulo 'Pathway Analysis'):
  * Entrada: los 20 metabolitos con su tabla de concentraciones (log10 + Pareto).
  * Mapeo de nombres -> HMDB / KEGG / PubChem (tabla metaboanalyst/tablas/name_map.csv;
    20/20 metabolitos mapeados).
  * Biblioteca de rutas KEGG; prueba de enriquecimiento 'Global test';
    topología: centralidad de intermediación relativa (relative-betweenness).
  * Resultados exportados a metaboanalyst/tablas/pathway_results.csv.

Este script NO recalcula el enriquecimiento (requiere las bibliotecas KEGG de
MetaboAnalyst); lee el resultado exportado, lo resume, asigna los metabolitos de
la lista a cada ruta y genera la figura de impacto vs significancia.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utils import METABOANALYST, TABLES, savefig

# Metabolitos de la lista presentes en las rutas de mayor impacto o significancia
# (KEGG IDs del name_map; asignación manual según los mapas KEGG de referencia)
PATHWAY_HITS = {
    "Pentose phosphate pathway": "Gluconolactona (C00198)",
    "Biosynthesis of unsaturated fatty acids": "1 ácido graso insaturado (palmitoleico o DGLA)",
    "Histidine metabolism": "L-Histidina, L-Aspartato",
    "Arginine and proline metabolism": "L-Arginina, Agmatina",
    "Arginine biosynthesis": "L-Arginina, L-Aspartato, L-Glutamato",
    "Phenylalanine, tyrosine and tryptophan biosynthesis": "L-Fenilalanina, L-Tirosina",
    "Phenylalanine metabolism": "L-Fenilalanina, L-Tirosina",
    "Alanine, aspartate and glutamate metabolism": "L-Aspartato, L-Asparagina, L-Glutamato",
    "Tyrosine metabolism": "L-Tirosina",
    "Glycine, serine and threonine metabolism": "L-Treonina, Colina",
    "Citrate cycle (TCA cycle)": "L-Malato",
}


def main():
    print("[04] Análisis de rutas (resultados MetaboAnalyst)")
    pw = pd.read_csv(METABOANALYST / "tablas" / "pathway_results.csv", index_col=0)
    pw.index.name = "Ruta"
    pw["metabolitos_lista"] = pw.index.map(PATHWAY_HITS).fillna("")
    pw.to_csv(TABLES / "04_rutas_resumen.csv")

    print(f"  rutas evaluadas: {len(pw)}")
    print(f"  rutas con p crudo < 0.05: {(pw['Raw p'] < 0.05).sum()} "
          f"({', '.join(pw.index[pw['Raw p'] < 0.05])})")
    print(f"  rutas con FDR < 0.05: {(pw['FDR'] < 0.05).sum()}")
    top_imp = pw.sort_values("Impact", ascending=False).head(5)
    print("  mayor impacto topológico:\n" +
          top_imp[["Hits", "Raw p", "FDR", "Impact"]].to_string())

    # ---- Figura: impacto vs -log10(p) -------------------------------------- #
    fig, ax = plt.subplots(figsize=(8, 6))
    size = 60 + 900 * pw["Impact"]
    sc = ax.scatter(pw["Impact"], pw["-log10(p)"], s=size, c=pw["-log10(p)"],
                    cmap="YlOrRd", vmin=0, edgecolor="#52514e", lw=0.6, alpha=0.9)
    lab = pw[(pw["Raw p"] < 0.11) | (pw["Impact"] > 0.2)]
    for name, r in lab.iterrows():
        right = r.Impact > 0.8                       # evita salir del panel
        dy = -12 if name.startswith("Arginine and proline") else 4
        ax.annotate(name, (r.Impact, r["-log10(p)"]), fontsize=7,
                    xytext=(-8 if right else 6, dy), textcoords="offset points",
                    ha="right" if right else "left")
    ax.axhline(-np.log10(0.05), ls="--", c="#52514e", lw=1)
    ax.set_xlabel("Impacto de la ruta (topología)")
    ax.set_ylabel("-log10(p)")
    ax.set_title("Análisis de rutas KEGG (MetaboAnalyst 6.0)")
    fig.colorbar(sc, ax=ax, label="-log10(p)")
    ax.spines[["top", "right"]].set_visible(False)
    savefig(fig, "04_rutas_impacto.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
