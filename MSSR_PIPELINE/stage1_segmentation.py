# -*- coding: utf-8 -*-
"""
stage1_segmentation.py
======================
Etapa 1 — Segmentación y extracción de recortes celulares.

Entrada : imágenes .tif crudas en raw_dir(condition)
Salida  : recortes 96x96 (_raw.tif para MSSR, _vis.png con barra de escala),
          overlays etiquetados y results_cells.csv por condición.

Equivale a SegmentacionTerminado140226, parametrizado por condición.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cv2
import tifffile as tiff
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops, regionprops_table
from skimage.color import label2rgb
from skimage.transform import resize
from skimage.morphology import white_tophat, disk

import config as cfg


# ------------------------------------------------------------
# Utilidades de imagen
# ------------------------------------------------------------
def load_tiff(path):
    img = tiff.imread(path)
    if img.ndim > 2:
        img = img[0]
    return img.astype(np.float32)


def segment(img):
    thresh = threshold_otsu(img)
    return img > thresh


def measure(mask, img_raw):
    labels = label(mask)
    props = regionprops_table(
        labels,
        intensity_image=img_raw,
        properties=('label', 'area', 'mean_intensity', 'centroid')
    )
    return labels, pd.DataFrame(props)


def normalize(img):
    img = img.astype(np.float32)
    return (img - img.min()) / (img.max() - img.min())


def add_scale_bar(img, pixel_size_um, bar_um=1):
    img_out = img.copy()
    bar_pixels = int(bar_um / pixel_size_um)
    h, w = img_out.shape
    margin = int(0.05 * w)
    bar_height = int(0.02 * h)
    x_start = w - margin - bar_pixels
    x_end = w - margin
    y_start = h - margin - bar_height
    y_end = h - margin
    img_out[y_start:y_end, x_start:x_end] = np.percentile(img_out, 99)
    return img_out


def to_display(img):
    img_disp = cv2.normalize(img, None, 0, 255, cv2.NORM_MINMAX)
    return img_disp.astype(np.uint8)


def passes_object_filter(region):
    """True si el objeto cumple los criterios de 'célula clara'. Devuelve también
    la razón de descarte (None si pasa) para diagnóstico agregado."""
    if region.area < cfg.OBJECT_MIN_AREA:
        return False, "area_min"
    if cfg.OBJECT_FILTER:
        if region.area > cfg.OBJECT_MAX_AREA:
            return False, "area_max"
        if region.solidity < cfg.OBJECT_MIN_SOLIDITY:
            return False, "solidez"
        if region.eccentricity > cfg.OBJECT_MAX_ECCENTRICITY:
            return False, "excentricidad"
    return True, None


def extract_fixed_size_crops_strict(img_raw, labels, crop_size=96):
    """Recortes cuadrados de tamaño fijo centrados en cada objeto que pasa el
    filtro de objetos. Descarta los que tocan el borde (para garantizar el
    tamaño exacto) y devuelve un conteo de descartes por criterio."""
    crops = []
    rejected = {"area_min": 0, "area_max": 0, "solidez": 0,
                "excentricidad": 0, "borde": 0}
    half = crop_size // 2
    H, W = img_raw.shape
    for region in regionprops(labels):
        ok, reason = passes_object_filter(region)
        if not ok:
            rejected[reason] += 1
            continue
        cy, cx = map(int, region.centroid)
        minr, maxr = cy - half, cy + half
        minc, maxc = cx - half, cx + half
        if minr < 0 or minc < 0 or maxr > H or maxc > W:
            rejected["borde"] += 1
            continue
        crop = img_raw[minr:maxr, minc:maxc]
        if crop.shape == (crop_size, crop_size):
            crops.append(crop.astype(np.float32))
    return crops, rejected


def analyze_image(path):
    img_raw = load_tiff(path)

    # La segmentación opera sobre una imagen con sustracción de fondo (white
    # top-hat), que elimina el fondo de fluorescencia suave y conserva los
    # objetos brillantes. El crudo (img_raw) NO se modifica: de él salen los
    # recortes para MSSR y las medidas de intensidad.
    if getattr(cfg, "BACKGROUND_SUBTRACTION", False):
        seg_input = white_tophat(img_raw, disk(cfg.TOPHAT_RADIUS))
    else:
        seg_input = img_raw

    img_norm = normalize(seg_input)
    binary_mask = segment(img_norm)
    labels, props = measure(binary_mask, img_raw)
    props["pixel_size_um"] = cfg.PIXEL_SIZE_UM
    props["area_um2"] = props["area"] * (cfg.PIXEL_SIZE_UM ** 2)
    return img_raw, binary_mask, labels, props


def estimate_snr(img_raw, binary_mask):
    """Relación señal/fondo aproximada de un campo: contraste de la señal
    (mediana dentro de la máscara) sobre el fondo (mediana fuera), normalizado
    por el ruido del fondo. Útil para QC e inclusión/exclusión de imágenes."""
    bg = img_raw[~binary_mask]
    sig = img_raw[binary_mask]
    if bg.size == 0 or sig.size == 0:
        return np.nan
    noise = np.std(bg)
    if noise <= 0:
        return np.nan
    return float((np.median(sig) - np.median(bg)) / noise)


# ------------------------------------------------------------
# Procesamiento por condición
# ------------------------------------------------------------
def run_segmentation(condition, show=False):
    folder_path = cfg.raw_dir(condition)
    out_path = cfg.salida_dir(condition)
    crop_path = cfg.crops_dir(condition)
    cfg.ensure_dirs(condition)

    tif_files = [f for f in os.listdir(folder_path)
                 if f.lower().endswith(".tif")
                 and os.path.isfile(os.path.join(folder_path, f))]

    print(f"[{cfg.label(condition)}] {len(tif_files)} imágenes RAW en {folder_path}")
    all_results = []

    for fname in tif_files:
        img_path = os.path.join(folder_path, fname)
        img_raw, binary, labels, props = analyze_image(img_path)

        snr = estimate_snr(img_raw, binary)

        cell_crops, rejected = extract_fixed_size_crops_strict(
            img_raw, labels, crop_size=cfg.CROP_SIZE
        )
        rej_str = ", ".join(f"{k}:{v}" for k, v in rejected.items() if v > 0)
        print(f"  {fname}: {len(cell_crops)} recortes | SNR≈{snr:.1f}"
              + (f" | descartados → {rej_str}" if rej_str else ""))

        # Control de calidad: descarta imágenes con relación señal/fondo baja
        # (típicamente fallos de adquisición con fondo de fluorescencia alto).
        if cfg.QC_MIN_SNR is not None and (np.isnan(snr) or snr < cfg.QC_MIN_SNR):
            print(f"  {fname} descartada por bajo SNR (< {cfg.QC_MIN_SNR})")
            continue

        # Filtro anti-ruido: descarta campos con demasiados objetos. En campos
        # densos (procarioplancton) este filtro no aplica; se omite si es None.
        if cfg.MAX_OBJECTS_PER_IMAGE is not None and len(cell_crops) > cfg.MAX_OBJECTS_PER_IMAGE:
            print(f"  {fname} descartada por exceso de objetos")
            continue

        props["image"] = fname
        all_results.append(props)

        for i, crop in enumerate(cell_crops):
            base_name = f"{fname.replace('.tif', '')}_cell_{i:03d}"
            # TIFF puro para MSSR (sin modificar)
            tiff.imwrite(os.path.join(crop_path, base_name + "_raw.tif"),
                         crop.astype(np.float32))
            # Versión visual con barra de escala
            crop_vis = add_scale_bar(crop, cfg.PIXEL_SIZE_UM, bar_um=cfg.SCALE_BAR_UM)
            plt.imsave(os.path.join(crop_path, base_name + "_vis.png"),
                       to_display(crop_vis), cmap='gray')

        # Overlay etiquetado
        overlay = label2rgb(labels, image=img_raw, alpha=0.4)
        overlay = overlay / overlay.max()
        plt.imsave(os.path.join(out_path, fname.replace(".tif", "_overlay.png")), overlay)
        tiff.imwrite(os.path.join(out_path, fname.replace(".tif", "_overlay.tif")),
                     (overlay * 65535).astype(np.uint16))

    if not all_results:
        print(f"[{cfg.label(condition)}] No se detectaron células / ninguna imagen pasó el QC.")
        return None

    final_df = pd.concat(all_results, ignore_index=True)
    csv_path = os.path.join(out_path, "results_cells.csv")
    final_df.to_csv(csv_path, index=False)
    kept = len(all_results)
    total = len(tif_files)
    print(f"[{cfg.label(condition)}] conservadas {kept}/{total} imágenes | "
          f"{len(final_df)} células | CSV: {csv_path}")
    return final_df


def run_all(show=False):
    for condition in cfg.CONDITIONS:
        run_segmentation(condition, show=show)


if __name__ == "__main__":
    run_all()
