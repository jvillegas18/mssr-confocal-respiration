# -*- coding: utf-8 -*-
"""
stage4_statistics.py
====================
Etapa 4 — Análisis estadístico consolidado.

Reúne en un solo módulo las tres rutinas previas (EstadisticaPareadaMSSR,
EstadisticaMultiple_PareadaMSSR y ESTADISTICA_INT_BOD_MSSR), que eran variantes
del mismo análisis sobre peaks_results.csv:

  - Carga de todas las condiciones con promedio por célula (image_name).
  - Métricas Delta_FWHM y mejora porcentual de resolución.
  - ANOVA factorial sobre Delta_FWHM (factores definidos en config.ANOVA_FACTORS).
  - Tests pareados RAW vs MSSR (t pareada y Wilcoxon) por tratamiento.
  - Análisis sub-difractivo frente al límite de Rayleigh por fluorocromo.
  - Figuras clave (scatter con línea de identidad, boxplots, Bland-Altman,
    histograma de vecino más cercano) guardadas en BASE_DIR/STATS.

Las salidas gráficas se guardan en disco para uso reproducible en lugar de
depender de la visualización interactiva de Colab.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import ttest_rel, wilcoxon

import config as cfg

STATS_DIR = os.path.join(cfg.BASE_DIR, "STATS")


# ------------------------------------------------------------
# Carga de datos
# ------------------------------------------------------------
def load_peaks(per_cell=True):
    """Combina los peaks_results.csv de todas las condiciones en un DataFrame.

    per_cell=True promedia por célula (image_name); per_cell=False conserva
    todos los picos (necesario para el análisis de vecino más cercano)."""
    frames = []
    for cond in cfg.CONDITIONS:
        path = cfg.peaks_csv(cond)
        if not os.path.exists(path):
            print(f"  ⚠ No existe {path}")
            continue
        fd = pd.read_csv(path)

        if per_cell:
            fd = fd.dropna(subset=['fwhm_raw_nm', 'fwhm_mssr_nm'])
            fd = (fd.groupby("image_name")
                    .agg(fwhm_raw_nm=("fwhm_raw_nm", "mean"),
                         fwhm_mssr_nm=("fwhm_mssr_nm", "mean"))
                    .reset_index())

        fd["treatment"] = cond["treatment"]
        fd["dye"] = cond["dye"]
        fd["Tratamiento"] = cfg.label(cond)
        frames.append(fd)

    if not frames:
        raise FileNotFoundError("No se encontró ningún peaks_results.csv.")
    return pd.concat(frames, ignore_index=True)


# ------------------------------------------------------------
# Métricas
# ------------------------------------------------------------
def add_metrics(df):
    df = df.copy()
    df["Delta_FWHM"] = df["fwhm_raw_nm"] - df["fwhm_mssr_nm"]
    df["Improvement_%"] = 100 * df["Delta_FWHM"] / df["fwhm_raw_nm"]
    return df


# ------------------------------------------------------------
# Resumen e inferencia
# ------------------------------------------------------------
def descriptive_summary(df):
    print("\n=== ΔFWHM por condición (nm) ===")
    print(df.groupby(["treatment", "dye"])["Delta_FWHM"].describe())
    print("\n=== Mejora porcentual de resolución (media ± SD) ===")
    print(df.groupby(["treatment", "dye"])["Improvement_%"].agg(['mean', 'std']))


def factorial_anova(df):
    """ANOVA tipo II sobre Delta_FWHM con los factores de config.ANOVA_FACTORS.
    Soporta diseños de uno o dos factores."""
    factors = [f for f in cfg.ANOVA_FACTORS if df[f].nunique() > 1]
    if not factors:
        print("\nANOVA omitido: no hay factores con más de un nivel.")
        return None

    terms = [f"C({f})" for f in factors]
    if len(factors) == 2:
        terms.append(f"C({factors[0]}):C({factors[1]})")
    formula = "Delta_FWHM ~ " + " + ".join(terms)

    model = smf.ols(formula, data=df).fit()
    table = sm.stats.anova_lm(model, typ=2)
    print(f"\n=== ANOVA ({formula}) ===")
    print(table)
    return table


def paired_tests(df):
    """t pareada y Wilcoxon RAW vs MSSR por tratamiento. Requiere el detalle
    por pico, por lo que recarga con per_cell=False."""
    raw = load_peaks(per_cell=False).dropna(subset=['fwhm_raw_nm', 'fwhm_mssr_nm'])
    print("\n=== Tests pareados RAW vs MSSR por tratamiento ===")
    rows = []
    for trt in raw["Tratamiento"].unique():
        sub = raw[raw["Tratamiento"] == trt]
        t_stat, p_t = ttest_rel(sub["fwhm_raw_nm"], sub["fwhm_mssr_nm"])
        try:
            w_stat, p_w = wilcoxon(sub["fwhm_raw_nm"], sub["fwhm_mssr_nm"])
        except ValueError:
            w_stat, p_w = np.nan, np.nan
        delta = sub["fwhm_raw_nm"] - sub["fwhm_mssr_nm"]
        rows.append({"Tratamiento": trt, "n": len(sub),
                     "reduccion_media_nm": delta.mean(),
                     "t": t_stat, "p_t": p_t, "W": w_stat, "p_wilcoxon": p_w})
        print(f"\n{trt}  (n={len(sub)})")
        print(f"  t pareada : t={t_stat:.3f}, p={p_t:.3e}")
        if not np.isnan(w_stat):
            print(f"  Wilcoxon  : W={w_stat:.3f}, p={p_w:.3e}")
    return pd.DataFrame(rows)


def subdiffraction_analysis():
    """Fracción de separaciones de vecino más cercano por debajo del límite de
    Rayleigh, por fluorocromo."""
    df = load_peaks(per_cell=False).dropna(subset=["nearest_neighbor_nm"])
    print("\n=== Eventos sub-difractivos (vecino más cercano < Rayleigh) ===")
    rows = []
    for cond in cfg.CONDITIONS:
        limit = cfg.rayleigh_limit(cond)
        sub = df[(df["dye"] == cond["dye"]) & (df["treatment"] == cond["treatment"])]
        if len(sub) == 0:
            continue
        n_below = int((sub["nearest_neighbor_nm"] < limit).sum())
        pct = 100 * n_below / len(sub)
        rows.append({"Tratamiento": cfg.label(cond), "rayleigh_nm": round(limit, 1),
                     "n_total": len(sub), "n_sub": n_below, "pct_sub": round(pct, 1)})
        print(f"  {cfg.label(cond)}: Rayleigh≈{limit:.0f} nm | "
              f"{n_below}/{len(sub)} ({pct:.1f}%) bajo el límite")
    return pd.DataFrame(rows)


# ------------------------------------------------------------
# Figuras
# ------------------------------------------------------------
def _save(fig_name):
    os.makedirs(STATS_DIR, exist_ok=True)
    path = os.path.join(STATS_DIR, fig_name)
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  figura → {path}")


def plot_identity_scatter(df):
    plt.figure(figsize=(7, 7))
    sns.scatterplot(data=df, x="fwhm_raw_nm", y="fwhm_mssr_nm",
                    hue="Tratamiento", palette="Set1", s=70)
    lo = min(df["fwhm_raw_nm"].min(), df["fwhm_mssr_nm"].min())
    hi = max(df["fwhm_raw_nm"].max(), df["fwhm_mssr_nm"].max())
    plt.plot([lo, hi], [lo, hi], 'k--', linewidth=2)
    plt.xlabel("FWHM Raw (nm)"); plt.ylabel("FWHM MSSR (nm)")
    plt.title("FWHM RAW vs MSSR (línea de identidad y = x)")
    plt.legend(title="Tratamiento")
    _save("scatter_identidad_fwhm.png")


def plot_delta_box(df):
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=df, x="treatment", y="Delta_FWHM", hue="dye")
    sns.stripplot(data=df, x="treatment", y="Delta_FWHM", hue="dye",
                  dodge=True, palette="dark:black", alpha=0.3, legend=False)
    plt.ylabel("Δ FWHM (nm)")
    plt.title("Reducción de FWHM tras MSSR")
    _save("boxplot_delta_fwhm.png")


def plot_bland_altman(per_cell_df):
    raw = per_cell_df
    mean_vals = (raw["fwhm_raw_nm"] + raw["fwhm_mssr_nm"]) / 2
    diff = raw["fwhm_raw_nm"] - raw["fwhm_mssr_nm"]
    md, sd = diff.mean(), diff.std(ddof=1)
    plt.figure(figsize=(8, 6))
    plt.scatter(mean_vals, diff, alpha=0.7)
    plt.axhline(md, color='red', linestyle='--', label=f"Media = {md:.2f}")
    plt.axhline(md + 1.96 * sd, color='blue', linestyle='--', label=f"+1.96 SD = {md + 1.96 * sd:.2f}")
    plt.axhline(md - 1.96 * sd, color='blue', linestyle='--', label=f"-1.96 SD = {md - 1.96 * sd:.2f}")
    plt.xlabel("Media de FWHM (nm)"); plt.ylabel("Diferencia RAW - MSSR (nm)")
    plt.title("Bland–Altman FWHM"); plt.legend()
    _save("bland_altman_fwhm.png")


def plot_nearest_neighbor():
    df = load_peaks(per_cell=False).dropna(subset=["nearest_neighbor_nm"])
    plt.figure(figsize=(8, 6))
    sns.histplot(data=df, x="nearest_neighbor_nm", hue="dye", bins=60,
                 element="step", stat="density", common_norm=False)
    for cond in {c["dye"]: c for c in cfg.CONDITIONS}.values():
        limit = cfg.rayleigh_limit(cond)
        plt.axvline(limit, linestyle="--", linewidth=2,
                    label=f"Rayleigh {cond['dye']} ≈ {limit:.0f} nm")
    plt.xlabel("Distancia al vecino más cercano (nm)"); plt.ylabel("Densidad")
    plt.title("Separaciones vs límite de difracción")
    plt.legend()
    _save("histograma_vecino_cercano.png")


# ------------------------------------------------------------
# Orquestador de la etapa
# ------------------------------------------------------------
def run_all():
    df = add_metrics(load_peaks(per_cell=True))
    descriptive_summary(df)
    factorial_anova(df)
    paired_tests(df)
    subdiffraction_analysis()
    plot_identity_scatter(df)
    plot_delta_box(df)
    plot_bland_altman(df)
    plot_nearest_neighbor()
    print(f"\nEstadística completa. Figuras en {STATS_DIR}")
    return df


if __name__ == "__main__":
    run_all()
