"""
02_analisis_univariado.py - Normalización, prueba t de Welch, FDR, fold change y volcano.

Equivalente a los módulos 'Normalization', 't-test' y 'Volcano plot' de
MetaboAnalyst 6.0 (Statistical Analysis [one factor]).

Decisiones
----------
* Normalización: log10 + escalado de Pareto (sin normalización por muestra,
  porque los datos depositados ya están normalizados por TIC).
* Prueba t de Welch (no asume varianzas iguales; n = 6 vs 6) sobre datos log10.
  El escalado de Pareto no cambia el estadístico t (es lineal por variable).
* Corrección por comparaciones múltiples: Benjamini-Hochberg (FDR < 0.05).
* Fold change (FC) = media(Control) / media(34 °C) sobre datos crudos, como en
  MetaboAnalyst; log2FC < 0 => metabolito MÁS abundante a 34 °C.
* Umbral de FC para el volcano: FC >= 1.5 (|log2FC| >= 0.585).
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from utils import (CONTROL, GROUP_COLORS, HEAT, TABLES, bh_fdr, load_matrix,
                   log10_transform, normalize, savefig)

FC_THR = 1.5
P_THR = 0.05


def plot_normalization(X, Xn):
    """Densidad y boxplots antes/después de normalizar (como en MetaboAnalyst)."""
    fig, ax = plt.subplots(1, 2, figsize=(11, 6), gridspec_kw={"wspace": 0.45})
    order = X.columns[::-1]
    ax[0].boxplot([X[c] for c in order], vert=False, tick_labels=order,
                  patch_artist=True, boxprops=dict(facecolor="#cde2fb"))
    ax[0].set_title("Antes de normalizar (área de pico)")
    ax[0].set_xlabel("Intensidad")
    ax[1].boxplot([Xn[c] for c in order], vert=False, tick_labels=order,
                  patch_artist=True, boxprops=dict(facecolor="#cde2fb"))
    ax[1].set_title("Después: log10 + Pareto")
    ax[1].set_xlabel("Intensidad normalizada")
    for a in ax:
        a.tick_params(labelsize=7)
        a.spines[["top", "right"]].set_visible(False)
    savefig(fig, "02_normalizacion_boxplots.png")
    plt.close(fig)


def main():
    print("[02] Análisis univariado")
    X, y = load_matrix()
    Xlog = log10_transform(X)
    Xn = normalize(X)
    Xn.assign(Label=y).to_csv(TABLES / "02_matriz_normalizada_log10_pareto.csv")
    plot_normalization(X, Xn)

    c, h = (y == CONTROL), (y == HEAT)
    t, p = stats.ttest_ind(Xlog[c], Xlog[h], equal_var=False)
    fc = X[c].mean() / X[h].mean()
    res = pd.DataFrame({
        "media_control_27C": X[c].mean(),
        "media_34C": X[h].mean(),
        "FC_control/34C": fc,
        "log2FC": np.log2(fc),
        "t_Welch": t,
        "p_raw": p,
        "FDR_BH": bh_fdr(p),
    })
    res["direccion"] = np.where(res["log2FC"] < 0, "mayor a 34 °C", "menor a 34 °C")
    res["pasa_FC>=1.5"] = (fc >= FC_THR) | (fc <= 1 / FC_THR)
    res["pasa_p<0.05"] = res["p_raw"] < P_THR
    res = res.sort_values("p_raw")
    res.round(5).to_csv(TABLES / "02_ttest_welch_fold_change.csv")

    print(f"  FDR < 0.05: {(res.FDR_BH < 0.05).sum()} metabolitos")
    print(f"  p crudo < 0.05: {(res.p_raw < 0.05).sum()} -> "
          f"{', '.join(res.index[res.p_raw < 0.05])}")
    print(f"  |FC| >= 1.5: {res['pasa_FC>=1.5'].sum()} "
          f"({((res['pasa_FC>=1.5']) & (res.log2FC < 0)).sum()} mayores a 34 °C, "
          f"{((res['pasa_FC>=1.5']) & (res.log2FC > 0)).sum()} menores a 34 °C)")
    print(res[["log2FC", "p_raw", "FDR_BH"]].round(4).head(6).to_string())

    # ---- Volcano plot ------------------------------------------------------- #
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    up34 = res["pasa_FC>=1.5"] & (res.log2FC < 0)
    dn34 = res["pasa_FC>=1.5"] & (res.log2FC > 0)
    ns = ~res["pasa_FC>=1.5"]
    ax.scatter(res.log2FC[ns], -np.log10(res.p_raw[ns]), s=60, c="#b4b2a9",
               edgecolor="white", linewidth=1.5, label=f"FC < 1.5 ({ns.sum()})")
    ax.scatter(res.log2FC[up34], -np.log10(res.p_raw[up34]), s=60,
               c=GROUP_COLORS[HEAT], edgecolor="white", linewidth=1.5,
               label=f"Mayor a 34 °C ({up34.sum()})")
    ax.scatter(res.log2FC[dn34], -np.log10(res.p_raw[dn34]), s=60,
               c=GROUP_COLORS[CONTROL], edgecolor="white", linewidth=1.5,
               label=f"Menor a 34 °C ({dn34.sum()})")
    ax.axhline(-np.log10(P_THR), ls="--", lw=1, c="#52514e")
    for v in (np.log2(FC_THR), -np.log2(FC_THR)):
        ax.axvline(v, ls=":", lw=1, c="#52514e")
    for name, r in res.iterrows():
        if r.p_raw < 0.1 or abs(r.log2FC) > 2:
            ax.annotate(name, (r.log2FC, -np.log10(r.p_raw)), fontsize=7,
                        xytext=(4, 3), textcoords="offset points", color="#0b0b0b")
    ax.set_xlabel("log2(FC)  [Control 27 °C / 34 °C]")
    ax.set_ylabel("-log10(p)  (t de Welch)")
    ax.set_title("Volcano plot - A. hyacinthus 27 °C vs 34 °C")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    savefig(fig, "02_volcano_plot.png")
    plt.close(fig)

    # ---- Boxplots de los metabolitos más discriminantes -------------------- #
    top = res.index[:6]
    fig, axes = plt.subplots(2, 3, figsize=(10, 6))
    for a, m in zip(axes.flat, top):
        data = [Xlog.loc[c, m], Xlog.loc[h, m]]
        bp = a.boxplot(data, tick_labels=["27 °C", "34 °C"], patch_artist=True,
                       widths=0.5)
        for patch, g in zip(bp["boxes"], (CONTROL, HEAT)):
            patch.set_facecolor(GROUP_COLORS[g])
            patch.set_alpha(0.55)
        for i, d in enumerate(data, start=1):
            a.scatter(np.full(len(d), i), d, s=18, c="#0b0b0b", zorder=3)
        a.set_title(f"{m}\np = {res.loc[m, 'p_raw']:.3f}", fontsize=9)
        a.set_ylabel("log10 intensidad", fontsize=8)
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    savefig(fig, "02_boxplots_top6.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
