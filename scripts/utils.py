"""
utils.py - Funciones y rutas compartidas por todos los scripts del proyecto.

Proyecto: Respuesta metabolómica de Acropora hyacinthus al estrés térmico agudo
Dataset : Metabolomics Workbench ST004329 (PR002743), AN007214 (ESI+) y AN007215 (ESI-)

Todas las rutas son relativas a la raíz del repositorio, de modo que los scripts
funcionan sin importar desde qué carpeta se ejecuten.
"""
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Rutas del proyecto
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
METADATA = ROOT / "data" / "metadata"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
METABOANALYST = ROOT / "metaboanalyst"

for _d in (PROCESSED, TABLES, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)

RAW_POS = RAW / "ST004329_named_metabolite_data_working_positive.txt"
RAW_NEG = RAW / "ST004329_named_metabolite_data_working_negative.txt"
MATRIX_CSV = PROCESSED / "ST004329_for_metaboanalyst.csv"  # formato MetaboAnalyst

# Semilla global: hace reproducibles la validación cruzada y las permutaciones
SEED = 2026

# Grupos comparados (etiquetas idénticas a las usadas en MetaboAnalyst)
CONTROL = "Control_27C"
HEAT = "Severe_34C"
GROUP_COLORS = {CONTROL: "#2a78d6", HEAT: "#eb6834"}  # azul / naranja (paleta validada CVD)

# Armonización de nombres: Workbench -> nombres usados en MetaboAnalyst
NAME_FIX = {
    "Beta-Muricholic acid": "beta-Muricholic acid",
    "(S)-MALATE": "L-Malic acid",
}


# --------------------------------------------------------------------------- #
# Lectura de la matriz procesada
# --------------------------------------------------------------------------- #
def load_matrix():
    """Devuelve (X, y): X = intensidades crudas (muestras x metabolitos),
    y = etiqueta de grupo de cada muestra."""
    df = pd.read_csv(MATRIX_CSV, index_col=0)
    y = df.pop("Label")
    return df.astype(float), y


# --------------------------------------------------------------------------- #
# Normalización (idéntica a la elegida en MetaboAnalyst)
#   1) Sin normalización por muestra (los datos ya vienen normalizados por TIC)
#   2) Transformación log10
#   3) Escalado de Pareto: (x - media) / sqrt(desv. estándar)
# --------------------------------------------------------------------------- #
def log10_transform(X):
    """Transformación log10. Todos los valores del dataset son > 0;
    si hubiera ceros se reemplazan por 1/5 del mínimo positivo de la variable."""
    X = X.copy()
    for c in X.columns:
        pos_min = X.loc[X[c] > 0, c].min()
        X.loc[X[c] <= 0, c] = pos_min / 5
    return np.log10(X)


def pareto_scale(X):
    """Centrado en la media y división por la raíz de la desviación estándar."""
    return (X - X.mean()) / np.sqrt(X.std(ddof=1))


def normalize(X):
    """Pipeline completo de normalización: log10 + Pareto."""
    return pareto_scale(log10_transform(X))


def bh_fdr(p):
    """Corrección de Benjamini-Hochberg (FDR) para un vector de p-valores."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def savefig(fig, name):
    """Guarda una figura en results/figures en PNG (300 dpi)."""
    path = FIGURES / name
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"  figura -> {path.relative_to(ROOT)}")
