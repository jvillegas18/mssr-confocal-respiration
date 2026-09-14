# -*- coding: utf-8 -*-
"""
stage0_convert.py
=================
Etapa 0 — Conversión por lote de .oib (Olympus) a .tif y archivado.

Lee los .oib crudos de cfg.RAW_OIB_DIR (carpeta plana), deduce la condición
(tratamiento × fluorocromo) a partir del nombre, conserva SOLO el canal de
fluorescencia (los .oib traen además un canal de luz de transmisión) y lo
escribe como .tif 2D en cfg.raw_dir(condition). Tras esta etapa, las etapas 1–5
corren sin cambios sobre la jerarquía generada.

Codificación del nombre (confirmada para el experimento 260327):
    Muestras BOD (sin INT):  <muestra><letra>-BOD-<excitacion>-<idx>.oib
    Muestras INT:            <muestra><letra>-<excitacion>-<idx>.oib  (sin token)
        letra       : D = DAPI, S = SYBR Green   (último carácter del token 0)
        tratamiento : token 'BOD' presente = BOD ; ausente = INT
        excitacion  : 405 = DAPI ; 488 = SYBR     (validación cruzada)
    Ejemplos:
        13D-BOD-405-01.oib  → BOD + DAPI
        13S-BOD-488-01.oib  → BOD + SYBR
        13S-488-01.oib      → INT + SYBR
        13D-405-01.oib      → INT + DAPI

Canal de fluorescencia: cfg.OIB_FLUORESCENCE_CHANNEL (None = autodetección por
asimetría; o un índice fijo 0/1). El ensayo en seco reporta el canal elegido y
avisa si no es unánime entre archivos.

Recomendación: ejecutar primero convert_all(dry_run=True) para confirmar la
clasificación (4 condiciones, nada 'sin clasificar') y el canal elegido.
"""

import os
import numpy as np
import tifffile as tiff
import config as cfg

try:
    import oiffile
except ImportError:
    oiffile = None


DYE_FROM_LETTER = {"D": "DAPI", "S": "SYBR"}
DYE_FROM_WAVELENGTH = {"405": "DAPI", "488": "SYBR"}


# ------------------------------------------------------------
# Clasificación por nombre de archivo
# ------------------------------------------------------------
def parse_condition(fname):
    """Devuelve el dict de condición de cfg.CONDITIONS que corresponde al
    archivo, o None si el nombre no encaja. Valida letra ↔ longitud de onda.

    Regla de tratamiento de este experimento: las muestras BOD llevan el token
    'BOD' en el nombre y carecen de INT; las muestras INT NO llevan token de
    tratamiento (p. ej. 13S-488-01), de modo que INT = ausencia de 'BOD'."""
    stem = os.path.splitext(fname)[0]
    tokens = stem.split("-")
    if not tokens or not tokens[0]:
        return None

    # Fluorocromo: último carácter del primer token (D = DAPI, S = SYBR).
    dye = DYE_FROM_LETTER.get(tokens[0][-1].upper())
    if dye is None:
        return None

    # Tratamiento: BOD si aparece el token 'BOD'; INT si aparece 'INT'; en
    # ausencia de ambos se asume el valor configurado (en 260327, las INT no
    # llevan token de tratamiento, así que ausencia de 'BOD' = INT).
    toks_up = [tok.upper() for tok in tokens]
    if any("BOD" in tok for tok in toks_up):
        treatment = "BOD"
    elif any("INT" in tok for tok in toks_up):
        treatment = "INT"
    else:
        treatment = cfg.DEFAULT_TREATMENT_WHEN_UNLABELED

    # Validación cruzada con la excitación, si está presente.
    wl = next((tok for tok in tokens if tok in DYE_FROM_WAVELENGTH), None)
    if wl and DYE_FROM_WAVELENGTH[wl] != dye:
        print(f"  ⚠ {fname}: letra→{dye} pero {wl} nm→{DYE_FROM_WAVELENGTH[wl]}. "
              f"Se usa la letra ({dye}).")

    for c in cfg.CONDITIONS:
        if c["treatment"] == treatment and c["dye"] == dye:
            return c
    return None


# ------------------------------------------------------------
# Selección del canal de fluorescencia
# ------------------------------------------------------------
from scipy.stats import skew


def _channel_metrics(arr_cyx):
    """Métricas por canal para distinguir fluorescencia de luz de transmisión.
    La fluorescencia es señal dispersa sobre fondo oscuro: media baja y
    asimetría (skew) alta. La transmisión es brillante y uniforme: media alta
    y asimetría baja o negativa."""
    metrics = []
    for c in range(arr_cyx.shape[0]):
        flat = arr_cyx[c].ravel().astype(np.float64)
        metrics.append({"mean": float(flat.mean()), "skew": float(skew(flat))})
    return metrics


def select_fluorescence(arr, channel=None):
    """Devuelve (imagen_2d, idx_canal, metricas).

    - 2D puro: se devuelve tal cual (idx None).
    - (C, H, W) con C pequeño: se elige el canal de fluorescencia. Si channel es
      None se autodetecta por mayor asimetría; si es un entero, se usa ese índice.
    - Dimensiones mayores (p. ej. C, Z, H, W): se elige canal y se proyecta z por
      máxima intensidad (caso no esperado en 260327, incluido por robustez)."""
    arr = np.squeeze(np.asarray(arr)).astype(np.float32)

    if arr.ndim == 2:
        return arr, None, None

    if arr.ndim >= 3 and arr.shape[0] <= 4:
        chan_axis = arr[:]  # (C, ...)
        if arr.ndim > 3:    # colapsa z dentro de cada canal
            chan_axis = arr.reshape(arr.shape[0], -1, arr.shape[-2], arr.shape[-1]).max(axis=1)
        metrics = _channel_metrics(chan_axis)
        idx = int(channel) if channel is not None else int(
            np.argmax([m["skew"] for m in metrics]))
        return chan_axis[idx].astype(np.float32), idx, metrics

    # Sin eje de canal reconocible: proyección de máxima intensidad.
    out = arr
    while out.ndim > 2:
        out = out.max(axis=0)
    return out.astype(np.float32), None, None


# ------------------------------------------------------------
# Conversión por lote
# ------------------------------------------------------------
def convert_all(channel=None, dry_run=False):
    """Convierte todos los .oib conservando el canal de fluorescencia.

    channel: None usa cfg.OIB_FLUORESCENCE_CHANNEL (None=autodetección por
             archivo). Un entero (0/1) fija el índice del canal.
    dry_run: solo reporta clasificación, forma y canal elegido, sin escribir."""
    if oiffile is None:
        raise ImportError("Falta 'oiffile'. Instala con: pip install oiffile")

    chan = cfg.OIB_FLUORESCENCE_CHANNEL if channel is None else channel

    src = cfg.RAW_OIB_DIR
    files = sorted(f for f in os.listdir(src) if f.lower().endswith(".oib"))
    print(f"{len(files)} archivos .oib en {src}")
    print(f"Canal de fluorescencia: {'autodetección' if chan is None else f'fijo = {chan}'}\n")

    counts, skipped, chosen = {}, [], {}
    for f in files:
        cond = parse_condition(f)
        if cond is None:
            skipped.append(f)
            continue
        counts[cfg.label(cond)] = counts.get(cfg.label(cond), 0) + 1

        arr = oiffile.imread(os.path.join(src, f))
        img2d, idx, metrics = select_fluorescence(arr, channel=chan)
        if idx is not None:
            chosen[idx] = chosen.get(idx, 0) + 1

        if dry_run:
            det = ""
            if metrics is not None:
                det = " | skew=[" + ", ".join(f"{m['skew']:.2f}" for m in metrics) + \
                      f"] → canal {idx}"
            print(f"  {f} → {cfg.label(cond):8s} | {np.asarray(arr).shape}→{img2d.shape}{det}")
            continue

        cfg.ensure_dirs(cond)
        out = os.path.join(cfg.raw_dir(cond), os.path.splitext(f)[0] + ".tif")
        tiff.imwrite(out, img2d)

    # ---- Resumen ----
    print("\nResumen por condición:")
    for k in sorted(counts):
        print(f"  {k}: {counts[k]}")

    if chosen:
        print("\nCanal de fluorescencia elegido (conteo de archivos):")
        for idx in sorted(chosen):
            print(f"  canal {idx}: {chosen[idx]}")
        if len(chosen) > 1:
            print("  ⚠ La elección NO fue unánime entre archivos. Revisa el "
                  "ensayo en seco e idealmente fija OIB_FLUORESCENCE_CHANNEL.")
        elif chan is None:
            unico = next(iter(chosen))
            print(f"  ✓ Unánime: el canal {unico} es fluorescencia en todos. "
                  f"Puedes fijar OIB_FLUORESCENCE_CHANNEL = {unico} en config.")

    if skipped:
        print(f"\n⚠ {len(skipped)} sin clasificar (revisar codificación del nombre):")
        for f in skipped[:15]:
            print("   ", f)
    if dry_run:
        print("\n(dry_run: no se escribió ningún .tif)")
    return counts, skipped


def run_all():
    convert_all(dry_run=False)


if __name__ == "__main__":
    run_all()
