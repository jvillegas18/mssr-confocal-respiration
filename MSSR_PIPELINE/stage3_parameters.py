# -*- coding: utf-8 -*-
"""
stage3_parameters.py
====================
Etapa 3 — Extracción de parámetros cuantitativos por par RAW/MSSR.

Para cada recorte detecta picos en la imagen MSSR (peak_local_max), mide el
FWHM por ajuste gaussiano del perfil en RAW y en MSSR, y calcula la distancia
al vecino más cercano.

Entrada : recortes *_raw.tif (crops_dir) y sus MSSR (mssr_dir)
Salida  : peaks_results.csv y cell_summary.csv por condición.

Equivale a Parametros_before_statistics. Mejora frente al original: el
emparejamiento RAW→MSSR exige extensión .tif para no seleccionar por error los
PNG de salida.
"""

import os
import numpy as np
import pandas as pd
import tifffile as tiff
from scipy.optimize import curve_fit
from scipy.spatial.distance import cdist
from skimage.feature import peak_local_max
from skimage.measure import label, regionprops

import config as cfg


# ------------------------------------------------------------
# Tamaño celular por media altura (comparable RAW vs MSSR)
# ------------------------------------------------------------
def measure_cell_size(image):
    """Mide el tamaño del/los objeto(s) de un recorte usando el contorno a media
    altura (medio máximo), criterio autoescalante y por tanto comparable entre
    RAW y MSSR pese a su distinto rango de intensidad.

    Devuelve, promediando sobre los objetos válidos del recorte:
      - diam_eqd_nm  : diámetro de círculo equivalente, 2*sqrt(area/pi)
      - diam_feret_nm: diámetro de Feret máximo (eje mayor)
    o (nan, nan) si no hay objeto medible."""
    img = image.astype(np.float32)
    peak = float(img.max())
    if peak <= 0:
        return np.nan, np.nan

    # Contorno a media altura respecto al fondo del propio recorte.
    base = float(np.median(img))
    half_level = base + 0.5 * (peak - base)
    mask = img >= half_level

    lbl = label(mask)
    eqd, feret = [], []
    for r in regionprops(lbl):
        if r.area < cfg.CELL_SIZE_MIN_AREA_PX:
            continue
        eqd.append(2.0 * np.sqrt(r.area / np.pi) * cfg.PIXEL_SIZE_NM)
        # feret_diameter_max disponible en skimage >=0.18; fallback al eje mayor.
        fd = getattr(r, "feret_diameter_max", None)
        feret.append((fd if fd is not None else r.major_axis_length) * cfg.PIXEL_SIZE_NM)

    if not eqd:
        return np.nan, np.nan
    return float(np.mean(eqd)), float(np.mean(feret))


def measure_fwhm_precise(profile):
    """Devuelve el FWHM en píxeles a partir de un perfil 1D; NaN si el ajuste
    falla o cae fuera del rango físico configurado."""
    x = np.arange(len(profile))
    p0 = [np.max(profile) - np.min(profile), np.argmax(profile), 1.5, np.min(profile)]
    try:
        popt, _ = curve_fit(gaussian, x, profile, p0=p0, maxfev=cfg.FWHM_MAXFEV)
        sigma = abs(popt[2])
        fwhm_px = 2.3548 * sigma
        fwhm_nm = fwhm_px * cfg.PIXEL_SIZE_NM
        if fwhm_nm > cfg.FWHM_MAX_NM or fwhm_nm < cfg.FWHM_MIN_NM:
            return np.nan
        return fwhm_px
    except Exception:
        return np.nan


def fwhm_at_peak(image, y, x, window=None):
    """FWHM en nm del perfil horizontal centrado en (y, x)."""
    window = cfg.FWHM_WINDOW if window is None else window
    if y - window < 0 or y + window >= image.shape[0] or \
       x - window < 0 or x + window >= image.shape[1]:
        return np.nan
    profile = image[y, x - window:x + window]
    fwhm_px = measure_fwhm_precise(profile)
    if np.isnan(fwhm_px):
        return np.nan
    return fwhm_px * cfg.PIXEL_SIZE_NM


# ------------------------------------------------------------
# Análisis de un par RAW / MSSR
# ------------------------------------------------------------
def analyze_pair(raw_path, mssr_path):
    image_name = os.path.basename(raw_path)
    raw = tiff.imread(raw_path)
    mssr = tiff.imread(mssr_path)
    if raw.ndim > 2:
        raw = raw.squeeze()
    if mssr.ndim > 2:
        mssr = mssr.squeeze()
    raw = raw.astype(np.float32)
    mssr = mssr.astype(np.float32)

    coordinates = peak_local_max(
        mssr, min_distance=cfg.PEAK_MIN_DISTANCE, threshold_rel=cfg.PEAK_THRESHOLD_REL
    )

    # Distancia al vecino más cercano (nm)
    if len(coordinates) > 1:
        distances = cdist(coordinates, coordinates)
        np.fill_diagonal(distances, np.inf)
        nearest = np.min(distances, axis=1) * cfg.PIXEL_SIZE_NM
    else:
        nearest = [np.nan] * len(coordinates)

    results = []
    for i, (y, x) in enumerate(coordinates):
        results.append({
            "image_name": image_name,
            "peak_id": i,
            "x_px": x,
            "y_px": y,
            "fwhm_raw_nm": fwhm_at_peak(raw, y, x),
            "fwhm_mssr_nm": fwhm_at_peak(mssr, y, x),
            "nearest_neighbor_nm": nearest[i],
        })

    summary = {
        "image_name": image_name,
        "n_peaks": len(coordinates),
        "total_intensity_raw": raw.sum(),
        "total_intensity_mssr": mssr.sum(),
    }

    # Tamaño celular pre/post MSSR (media altura), comparable entre ambas.
    if getattr(cfg, "MEASURE_CELL_SIZE", False):
        eqd_raw, feret_raw = measure_cell_size(raw)
        eqd_mssr, feret_mssr = measure_cell_size(mssr)
        summary.update({
            "diam_eqd_raw_nm": eqd_raw,
            "diam_eqd_mssr_nm": eqd_mssr,
            "diam_feret_raw_nm": feret_raw,
            "diam_feret_mssr_nm": feret_mssr,
        })

    return pd.DataFrame(results), pd.DataFrame([summary])


def _find_mssr(base_name, mssr_path):
    """Localiza el TIFF MSSR correspondiente a un recorte (exige .tif)."""
    candidates = [
        f for f in os.listdir(mssr_path)
        if f.startswith(base_name) and "mssr" in f.lower() and f.lower().endswith(".tif")
    ]
    return candidates[0] if candidates else None


# ------------------------------------------------------------
# Procesamiento por condición
# ------------------------------------------------------------
def run_parameters(condition):
    raw_path = cfg.crops_dir(condition)
    mssr_path = cfg.mssr_dir(condition)
    cfg.ensure_dirs(condition)

    peaks_out = cfg.peaks_csv(condition)
    summary_out = cfg.summary_csv(condition)
    for p in (peaks_out, summary_out):
        if os.path.exists(p):
            os.remove(p)

    raw_files = [f for f in os.listdir(raw_path) if f.lower().endswith("_raw.tif")]
    print(f"[{cfg.label(condition)}] {len(raw_files)} recortes RAW para parametrizar")

    for raw_file in raw_files:
        base_name = raw_file.replace(".tif", "")
        mssr_name = _find_mssr(base_name, mssr_path)
        if mssr_name is None:
            print(f"  ⚠ Sin MSSR para {raw_file}")
            continue

        df_peaks, df_summary = analyze_pair(
            os.path.join(raw_path, raw_file),
            os.path.join(mssr_path, mssr_name)
        )

        df_peaks.to_csv(peaks_out, mode='a',
                        header=not os.path.exists(peaks_out), index=False)
        df_summary.to_csv(summary_out, mode='a',
                          header=not os.path.exists(summary_out), index=False)

    print(f"[{cfg.label(condition)}] Parámetros → {peaks_out}")


def run_all():
    for condition in cfg.CONDITIONS:
        run_parameters(condition)


if __name__ == "__main__":
    run_all()
