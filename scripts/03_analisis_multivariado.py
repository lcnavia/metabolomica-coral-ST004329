"""
03_analisis_multivariado.py - PCA, PLS-DA (VIP, validación cruzada, permutaciones)
y heatmap con agrupamiento jerárquico.

Todos los análisis se hacen sobre la matriz normalizada log10 + Pareto.

* PCA (no supervisado): explora la estructura global y la presencia de outliers.
  Con log10 + Pareto se reproduce el resultado de MetaboAnalyst
  (PC1 = 69.9 %, PC2 = 11.2 %).
* PLS-DA (supervisado): Y codificada 0/1. Se evalúa con validación cruzada
  leave-one-out (LOOCV) para 1-5 componentes (exactitud, R2, Q2) y con una prueba
  de 1000 permutaciones de las etiquetas usando Q2 como estadístico.
  Con n = 12 y 20 variables el riesgo de sobreajuste es alto: la permutación es
  la prueba decisiva de si la separación es real.
* VIP (Variable Importance in Projection) del primer componente; VIP >= 1 se
  considera discriminante.
* Heatmap: z-score por metabolito, distancia euclidiana, enlace de Ward.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Ellipse
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA

from utils import (CONTROL, GROUP_COLORS, HEAT, SEED, TABLES, load_matrix,
                   normalize, savefig)

N_PERM = 1000
MAX_COMP = 5


# --------------------------------------------------------------------------- #
# Funciones auxiliares PLS-DA
# --------------------------------------------------------------------------- #
def vip_scores(pls, n_comp=None):
    """VIP según Wold: sqrt(p * sum_a(SS_a * w_ja^2 / ||w_a||^2) / sum_a SS_a)."""
    T, W, Q = pls.x_scores_, pls.x_weights_, pls.y_loadings_
    a = n_comp or T.shape[1]
    T, W, Q = T[:, :a], W[:, :a], Q[:, :a]
    p = W.shape[0]
    ss = (Q[0] ** 2) * (T ** 2).sum(axis=0)
    wn = (W / np.linalg.norm(W, axis=0)) ** 2
    return np.sqrt(p * (wn @ ss) / ss.sum())


def loocv(Xv, yv, k):
    """LOOCV de PLS-DA con k componentes: devuelve (exactitud, Q2)."""
    n = len(yv)
    pred = np.empty(n)
    for i in range(n):
        tr = np.arange(n) != i
        # se re-centra con la media del set de entrenamiento (sin fuga de info)
        mu = Xv[tr].mean(axis=0)
        m = PLSRegression(n_components=k, scale=False).fit(Xv[tr] - mu, yv[tr])
        pred[i] = m.predict((Xv[i] - mu)[None, :]).ravel()[0]
    acc = np.mean((pred > 0.5) == (yv == 1))
    q2 = 1 - np.sum((yv - pred) ** 2) / np.sum((yv - yv.mean()) ** 2)
    return acc, q2


def r2_fit(Xv, yv, k):
    m = PLSRegression(n_components=k, scale=False).fit(Xv, yv)
    pred = m.predict(Xv).ravel()
    return 1 - np.sum((yv - pred) ** 2) / np.sum((yv - yv.mean()) ** 2)


def conf_ellipse(ax, pts, color):
    """Elipse de confianza del 95 % para un grupo en el plano de scores."""
    cov = np.cov(pts.T)
    vals, vecs = np.linalg.eigh(cov)
    ang = np.degrees(np.arctan2(*vecs[:, -1][::-1]))
    w, h = 2 * np.sqrt(vals[::-1] * 5.991)  # chi2(0.95, 2 gl)
    ax.add_patch(Ellipse(pts.mean(axis=0), w, h, angle=ang, color=color,
                         alpha=0.15, lw=0))


def main():
    print("[03] Análisis multivariado")
    rng = np.random.default_rng(SEED)
    X, y = load_matrix()
    Xn = normalize(X)
    Xv = Xn.values
    yv = (y == HEAT).astype(float).values

    # ------------------------------ PCA -------------------------------------- #
    pca = PCA().fit(Xv)
    ev = pca.explained_variance_ratio_ * 100
    scores = pca.transform(Xv)
    print(f"  PCA: PC1 = {ev[0]:.1f} %, PC2 = {ev[1]:.1f} %, PC3 = {ev[2]:.1f} % "
          f"(acumulado PC1-2 = {ev[:2].sum():.1f} %)")
    pd.DataFrame(scores[:, :3], index=Xn.index, columns=["PC1", "PC2", "PC3"]) \
        .assign(Label=y).round(4).to_csv(TABLES / "03_PCA_scores.csv")
    pd.DataFrame(pca.components_[:2].T, index=Xn.columns, columns=["PC1", "PC2"]) \
        .round(4).to_csv(TABLES / "03_PCA_loadings.csv")

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    k = np.arange(1, 8)
    ax[0].plot(k, ev[:7], "o-", c="#2a78d6", lw=2, label="Individual")
    ax[0].plot(k, np.cumsum(ev[:7]), "s-", c="#52514e", lw=2, label="Acumulada")
    ax[0].set_xlabel("Componente principal")
    ax[0].set_ylabel("Varianza explicada (%)")
    ax[0].set_title("Scree plot")
    ax[0].legend(frameon=False)
    for g in (CONTROL, HEAT):
        sel = (y == g).values
        conf_ellipse(ax[1], scores[sel, :2], GROUP_COLORS[g])
        ax[1].scatter(scores[sel, 0], scores[sel, 1], s=70, c=GROUP_COLORS[g],
                      edgecolor="white", lw=1.5, label=g)
    for i, s in enumerate(Xn.index):
        ax[1].annotate(s, scores[i, :2], fontsize=7, xytext=(4, 3),
                       textcoords="offset points")
    ax[1].set_xlabel(f"PC1 ({ev[0]:.1f} %)")
    ax[1].set_ylabel(f"PC2 ({ev[1]:.1f} %)")
    ax[1].set_title("PCA scores (elipses 95 %)")
    ax[1].legend(frameon=False)
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    savefig(fig, "03_PCA_scree_scores.png")
    plt.close(fig)

    # --------------------------- PLS-DA: CV ---------------------------------- #
    cv = []
    for kc in range(1, MAX_COMP + 1):
        acc, q2 = loocv(Xv, yv, kc)
        cv.append({"componentes": kc, "exactitud_LOOCV": acc,
                   "R2": r2_fit(Xv, yv, kc), "Q2": q2})
    cv = pd.DataFrame(cv).round(3)
    cv.to_csv(TABLES / "03_PLSDA_validacion_cruzada.csv", index=False)
    print("  PLS-DA LOOCV:\n" + cv.to_string(index=False))

    best_k = int(cv.loc[cv.Q2.idxmax(), "componentes"])
    q2_obs = cv.loc[cv.componentes == best_k, "Q2"].item()

    # ----------------------- PLS-DA: permutaciones --------------------------- #
    # Se evalúa el modelo parsimonioso (2 componentes) y el de mayor Q2.
    # Estadístico: Q2 de LOOCV (no satura, a diferencia de la exactitud de
    # entrenamiento o la distancia B/W, que con n pequeño casi siempre es máxima).
    perm_res = []
    perm_q2 = None
    for kc in sorted({2, best_k}):
        q2_k = cv.loc[cv.componentes == kc, "Q2"].item()
        null = np.array([loocv(Xv, rng.permutation(yv), kc)[1]
                         for _ in range(N_PERM)])
        p_k = (np.sum(null >= q2_k) + 1) / (N_PERM + 1)
        perm_res.append({"componentes": kc, "Q2_observado": q2_k,
                         "Q2_perm_mediana": np.median(null),
                         "Q2_perm_p95": np.percentile(null, 95), "p_perm": p_k})
        print(f"  Permutación ({N_PERM}x, {kc} comp.): Q2 obs = {q2_k:.3f}, "
              f"p = {p_k:.3f}")
        if kc == best_k:
            perm_q2, p_perm = null, p_k
    pd.DataFrame(perm_res).round(4).to_csv(
        TABLES / "03_PLSDA_permutaciones_resumen.csv", index=False)
    pd.DataFrame({"Q2_permutado": perm_q2}).to_csv(
        TABLES / "03_PLSDA_permutaciones_Q2.csv", index=False)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    w = 0.27
    for j, (col, colr) in enumerate([("exactitud_LOOCV", "#2a78d6"),
                                     ("R2", "#eb6834"), ("Q2", "#1baf7a")]):
        ax[0].bar(cv.componentes + (j - 1) * w, cv[col], width=w - 0.02,
                  color=colr, label=col)
    ax[0].set_xlabel("Número de componentes")
    ax[0].set_ylabel("Desempeño")
    ax[0].set_ylim(min(0, cv.Q2.min() - 0.05), 1.05)
    ax[0].axhline(0, c="#52514e", lw=0.8)
    ax[0].set_title("PLS-DA - validación cruzada (LOOCV)")
    ax[0].legend(frameon=False, fontsize=8)
    ax[1].hist(perm_q2, bins=30, color="#b4b2a9", edgecolor="white")
    ax[1].axvline(q2_obs, c="#d03b3b", lw=2)
    ax[1].annotate(f"Q2 observado = {q2_obs:.2f}\np = {p_perm:.3f}",
                   (q2_obs, ax[1].get_ylim()[1] * 0.85), xytext=(-110, 0),
                   textcoords="offset points", fontsize=9)
    ax[1].set_xlabel("Q2 con etiquetas permutadas")
    ax[1].set_ylabel("Frecuencia")
    ax[1].set_title(f"Prueba de permutaciones ({N_PERM}x, {best_k} comp.)")
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    savefig(fig, "03_PLSDA_CV_permutacion.png")
    plt.close(fig)

    # ------------------------------ VIP -------------------------------------- #
    pls = PLSRegression(n_components=max(best_k, 2), scale=False).fit(Xv, yv)
    vip = pd.DataFrame({
        "VIP_comp1": vip_scores(pls, 1),
        f"VIP_{max(best_k, 2)}comp": vip_scores(pls),
        "mayor_en": np.where(X[y == HEAT].mean() > X[y == CONTROL].mean(),
                             HEAT, CONTROL),
    }, index=Xn.columns).sort_values("VIP_comp1", ascending=False)
    vip.round(3).to_csv(TABLES / "03_PLSDA_VIP.csv")
    print("  VIP >= 1 (comp. 1): " + ", ".join(vip.index[vip.VIP_comp1 >= 1]))

    fig, ax = plt.subplots(figsize=(6.5, 6))
    top = vip.head(15)[::-1]
    cols = [GROUP_COLORS[g] for g in top.mayor_en]
    ax.hlines(top.index, 0, top.VIP_comp1, color="#d3d1c7", lw=1)
    ax.scatter(top.VIP_comp1, top.index, s=70, c=cols, edgecolor="white", lw=1.5,
               zorder=3)
    ax.axvline(1, ls="--", c="#52514e", lw=1)
    ax.set_xlabel("VIP (componente 1)")
    ax.set_title("PLS-DA - metabolitos más importantes")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=GROUP_COLORS[g],
                              label=f"más alto en {g}") for g in (CONTROL, HEAT)],
              frameon=False, fontsize=8, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    savefig(fig, "03_PLSDA_VIP.png")
    plt.close(fig)

    # ----------------------------- Heatmap ----------------------------------- #
    z = (Xn - Xn.mean()) / Xn.std(ddof=1)
    g = sns.clustermap(z.T, method="ward", metric="euclidean", cmap="RdBu_r",
                       center=0, vmin=-3, vmax=3, figsize=(8, 8),
                       col_colors=[GROUP_COLORS[l] for l in y],
                       cbar_kws={"label": "z-score"}, dendrogram_ratio=0.12,
                       cbar_pos=(0.0, 0.9, 0.03, 0.08))
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    savefig(g.fig, "03_heatmap_clustering.png")
    plt.close(g.fig)

    # Subestructura del grupo 34 °C: ¿qué réplicas se separan del control?
    pc1 = pd.Series(scores[:, 0], index=Xn.index)
    print("  PC1 por muestra 34 °C: " +
          ", ".join(f"{s}={v:.2f}" for s, v in pc1[y == HEAT].items()))


if __name__ == "__main__":
    main()
