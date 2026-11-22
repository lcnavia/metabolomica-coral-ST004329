"""
01_preprocesamiento.py - De las tablas del Workbench a la matriz de análisis.

Pasos
-----
1. Lee las matrices "working" de metabolitos identificados (ESI+ y ESI-) del
   estudio ST004329 descargadas de Metabolomics Workbench.
2. Une ambos modos de ionización (10 + 10 metabolitos = 20) y armoniza nombres.
3. Control de calidad: calcula el RSD (%) de cada metabolito en las 4 muestras QC
   (pool). Criterio habitual en metabolómica no dirigida: RSD_QC <= 30 %.
4. Revisa valores faltantes o ceros.
5. Selecciona la comparación de interés: Acropora hyacinthus 27 °C (control)
   vs 34 °C (estrés severo), n = 6 por grupo.
6. Exporta la matriz en formato MetaboAnalyst (muestras en filas,
   2.ª columna = 'Label') y verifica que coincide con la cargada en MetaboAnalyst.
"""
import numpy as np
import pandas as pd

from utils import (CONTROL, HEAT, MATRIX_CSV, NAME_FIX, RAW_NEG, RAW_POS,
                   TABLES)


def read_workbench(path, mode):
    """Lee una tabla 'named_metabolite_data_working' del Workbench.
    Fila 1 = nombres de muestra, fila 2 = factores, siguientes = metabolitos."""
    df = pd.read_csv(path, sep="\t", index_col=0, skiprows=[1])
    df = df.dropna(how="all")              # elimina la fila vacía final
    df = df.loc[:, ~df.columns.str.startswith("Unnamed")]
    df.index = df.index.str.strip()
    df = df.astype(float)
    df.insert(0, "ion_mode", mode)
    return df


def main():
    print("[01] Preprocesamiento")
    pos = read_workbench(RAW_POS, "ESI+")
    neg = read_workbench(RAW_NEG, "ESI-")
    data = pd.concat([pos, neg])
    data.index = [NAME_FIX.get(i, i) for i in data.index]
    modes = data.pop("ion_mode")
    print(f"  metabolitos: {data.shape[0]} (ESI+ {len(pos)}, ESI- {len(neg)}); "
          f"muestras: {data.shape[1]}")

    # ---- Control de calidad con el pool QC -------------------------------- #
    qc = data[[c for c in data.columns if c.startswith("QC")]]
    rsd = (qc.std(axis=1, ddof=1) / qc.mean(axis=1) * 100).round(1)
    qc_tab = pd.DataFrame({"ion_mode": modes, "RSD_QC_%": rsd,
                           "pasa_RSD<=30%": rsd <= 30}).sort_values("RSD_QC_%")
    qc_tab.to_csv(TABLES / "01_QC_RSD_por_metabolito.csv")
    print(f"  RSD en QC: mediana {rsd.median():.1f} %, máx {rsd.max():.1f} % "
          f"({(rsd <= 30).sum()}/{len(rsd)} metabolitos con RSD <= 30 %)")

    # ---- Faltantes / ceros ------------------------------------------------ #
    n_missing = int(data.isna().sum().sum())
    n_zero = int((data == 0).sum().sum())
    print(f"  valores faltantes: {n_missing}; ceros: {n_zero} -> no se requiere imputación")

    # ---- Subconjunto A. hyacinthus 27 vs 34 °C ---------------------------- #
    ctrl = [f"AH27_{i}" for i in range(1, 7)]
    heat = [f"AH34_{i}" for i in range(1, 7)]
    sub = data[ctrl + heat].T
    sub.insert(0, "Label", [CONTROL] * 6 + [HEAT] * 6)
    sub.index.name = "Sample"

    # Orden de columnas igual al archivo cargado en MetaboAnalyst
    if MATRIX_CSV.exists():
        ref = pd.read_csv(MATRIX_CSV, index_col=0)
        sub = sub[ref.columns]
        same = np.allclose(sub.drop(columns="Label").values,
                           ref.drop(columns="Label").values, rtol=1e-9)
        print(f"  verificación contra la matriz usada en MetaboAnalyst: "
              f"{'IDÉNTICA' if same else 'DIFERENTE - revisar'}")
    else:
        sub.to_csv(MATRIX_CSV)
        print(f"  matriz exportada -> {MATRIX_CSV}")

    # Tabla descriptiva por grupo
    desc = sub.groupby("Label").agg(["mean", "std"]).T
    desc.to_csv(TABLES / "01_descriptivos_por_grupo.csv")


if __name__ == "__main__":
    main()
