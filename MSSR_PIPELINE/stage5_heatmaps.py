# -*- coding: utf-8 -*-
"""
stage5_heatmaps.py
==================
Etapa 5 — Heatmaps de intensidad MSSR enmascarados por la señal RAW.

Calcula un máximo global (percentil 99 sobre todos los píxeles enmascarados de
todas las condiciones) para que la escala de color sea comparable entre
condiciones, y genera por cada recorte:
  - un heatmap con colorbar y contorno de máscara,
  - un PNG limpio 96x96 sin ejes (para montajes).

Equivale a HEATMAPS, parametrizado por condición y con el nombre del archivo
MSSR derivado de config.MSSR_SUFFIX.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from skimage import io
from skimage.filters import threshold_otsu

import config as cfg


def _mssr_name(raw_file):
    """Nombre del TIFF MSSR a partir del recorte _raw.tif."""
    return raw_file.replace("_raw.tif", "_raw" + cfg.MSSR_SUFFIX + ".tif")


def _iter_pairs(condition):
    """Genera pares (ruta_raw, ruta_mssr) existentes para una condición."""
    raw_folder = cfg.crops_dir(condition)
    mssr_folder = cfg.mssr_dir(condition)
    for f in os.listdir(raw_folder):
        if not f.endswith("_raw.tif"):
            continue
        mssr_path = os.path.join(mssr_folder, _mssr_name(f))
        if os.path.exists(mssr_path):
            yield f, os.path.join(raw_folder, f), mssr_path


def compute_global_vmax(percentile=99):
    """Percentil global de intensidad MSSR dentro de la máscara RAW."""
    print("Calculando vmax global (P%d)..." % percentile)
    all_pixels = []
    for cond in cfg.CONDITIONS:
        for _, raw_path, mssr_path in _iter_pairs(cond):
            raw = io.imread(raw_path).astype(float)
            mssr = io.imread(mssr_path).astype(float)
            mask = raw > threshold_otsu(raw)
            pixels = mssr[mask]
            if len(pixels):
                all_pixels.extend(pixels)
    vmax = np.percentile(np.array(all_pixels), percentile) if all_pixels else 1.0
    print(f"vmax global P{percentile}: {vmax:.3f}")
    return vmax


def run_heatmaps(condition, vmax, vmin=0):
    heatmap_folder = cfg.heatmap_dir(condition)
    clean_folder = os.path.join(heatmap_folder, "clean_png")
    os.makedirs(heatmap_folder, exist_ok=True)
    os.makedirs(clean_folder, exist_ok=True)

    n = 0
    for fname, raw_path, mssr_path in _iter_pairs(condition):
        raw = io.imread(raw_path).astype(float)
        mssr = io.imread(mssr_path).astype(float)
        mask = raw > threshold_otsu(raw)
        masked = mssr * mask

        # Heatmap con colorbar y contorno
        plt.figure(figsize=(4, 4))
        im = plt.imshow(masked, cmap="turbo", vmin=vmin, vmax=vmax, interpolation="bicubic")
        plt.contour(mask, levels=[0.5], colors="cyan", linewidths=1)
        plt.colorbar(im, label="MSSR intensity")
        plt.axis("off")
        plt.savefig(os.path.join(heatmap_folder, fname.replace("_raw.tif", "_MSSR_heatmap.png")),
                    bbox_inches="tight", pad_inches=0)
        plt.close()

        # PNG limpio 96x96
        fig = plt.figure(figsize=(1, 1), dpi=cfg.CROP_SIZE)
        plt.imshow(masked, cmap="turbo", vmin=vmin, vmax=vmax, interpolation="nearest")
        plt.axis("off")
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        fig.savefig(os.path.join(clean_folder, fname.replace("_raw.tif", "_MSSR_heatmap_clean.png")),
                    dpi=cfg.CROP_SIZE)
        plt.close(fig)
        n += 1

    print(f"[{cfg.label(condition)}] {n} heatmaps → {heatmap_folder}")


def run_all():
    vmax = compute_global_vmax()
    for condition in cfg.CONDITIONS:
        run_heatmaps(condition, vmax=vmax)
    print("Heatmaps generados correctamente.")


if __name__ == "__main__":
    run_all()
