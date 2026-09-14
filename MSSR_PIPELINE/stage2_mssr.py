# -*- coding: utf-8 -*-
"""
stage2_mssr.py
==============
Etapa 2 — Aplicación de MSSR a los recortes celulares.

Entrada : recortes *_raw.tif en crops_dir(condition)
Salida  : TIFF MSSR float32 (para análisis), PNG con barra de escala y PNG de
          presentación (upscale + 600 dpi) en mssr_dir(condition).

Equivale a MultiparametricoMSSRFernando, parametrizado por condición. El núcleo
algorítmico vive en mssr_core (ver nota sobre tMSSR/tMean en ese archivo).
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import cv2
from skimage.io import imread, imsave
from skimage.color import rgb2gray

import config as cfg
from mssr_core import apply_mssr, tMSSR, tMean


def _prepare_2d(image):
    """Normaliza la entrada a una imagen 2D lista para MSSR; devuelve None si
    la estructura no es compatible."""
    if image.ndim == 2:
        return image
    if image.ndim == 3:
        if image.shape[2] == 3:                       # RGB
            return rgb2gray(image)
        if image.shape[0] < 10:                        # pila (frames, H, W)
            stack = tMSSR(image, fwhm=cfg.MSSR_FWHM, amp=cfg.MSSR_AMP, order=cfg.MSSR_ORDER)
            return tMean(stack)
    return None


def _save_scalebar_pngs(image_mssr, out_path, base_filename):
    """Guarda un PNG con barra de escala de 1 µm y una versión de presentación."""
    img_norm = cv2.normalize(image_mssr, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    img_rgb = (cv2.cvtColor(img_norm, cv2.COLOR_GRAY2BGR)
               if img_norm.ndim == 2 else img_norm.copy())

    bar_length_px = int(cfg.SCALE_BAR_UM * cfg.PIXELS_PER_MICRON)
    h, w = img_rgb.shape[:2]
    margin, thickness = 10, 2
    x_start, x_end = w - margin - bar_length_px, w - margin
    y = h - margin
    cv2.line(img_rgb, (x_start, y), (x_end, y), (255, 255, 255), thickness)
    cv2.putText(img_rgb, f"{cfg.SCALE_BAR_UM} um", (x_start - 4, y - 5),
                cv2.FONT_HERSHEY_SIMPLEX, 0.25, (255, 255, 255), 1, cv2.LINE_AA)

    png_out = os.path.join(out_path, base_filename + ".png")
    cv2.imwrite(png_out, img_rgb)

    # Presentación: upscale bicúbico + 600 dpi
    img_big = cv2.resize(img_rgb, None, fx=2, fy=2, interpolation=cv2.INTER_CUBIC)
    pres_out = os.path.join(out_path, base_filename + "_presentation.png")
    plt.figure(figsize=(6, 6))
    plt.imshow(cv2.cvtColor(img_big, cv2.COLOR_BGR2RGB))
    plt.axis('off')
    plt.tight_layout(pad=0)
    plt.savefig(pres_out, dpi=600, bbox_inches='tight', pad_inches=0)
    plt.close()


def run_mssr(condition, show=False):
    folder_path = cfg.crops_dir(condition)
    out_path = cfg.mssr_dir(condition)
    cfg.ensure_dirs(condition)

    files = [f for f in os.listdir(folder_path) if f.lower().endswith('_raw.tif')]
    print(f"[{cfg.label(condition)}] {len(files)} recortes a procesar con MSSR")

    for filename in files:
        file_path = os.path.join(folder_path, filename)
        image = imread(file_path)
        img2d = _prepare_2d(image)
        if img2d is None:
            print(f"  ⚠ Estructura no compatible: {filename}. Saltando.")
            continue

        image_mssr = apply_mssr(img2d, amp=cfg.MSSR_AMP, fwhm=cfg.MSSR_FWHM, order=cfg.MSSR_ORDER)

        if show:
            fig, ax = plt.subplots(1, 2, figsize=(12, 6))
            ax[0].imshow(img2d, cmap='gray'); ax[0].set_title('Original'); ax[0].axis('off')
            ax[1].imshow(image_mssr, cmap='gray'); ax[1].set_title('MSSR'); ax[1].axis('off')
            plt.tight_layout(); plt.show()

        base = filename.replace('.tif', cfg.MSSR_SUFFIX)
        imsave(os.path.join(out_path, base + ".tif"), image_mssr.astype('float32'))
        _save_scalebar_pngs(image_mssr, out_path, base)

    print(f"[{cfg.label(condition)}] MSSR completado → {out_path}")


def run_all(show=False):
    for condition in cfg.CONDITIONS:
        run_mssr(condition, show=show)


if __name__ == "__main__":
    run_all()
