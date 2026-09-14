# -*- coding: utf-8 -*-
"""
config.py
=========
Única fuente de verdad del pipeline MSSR. Para correr el análisis sobre un
nuevo conjunto de imágenes, en la mayoría de los casos basta con editar
BASE_DIR y, si cambia el diseño experimental, la lista CONDITIONS.

Estructura de carpetas esperada (se construye automáticamente a partir de
BASE_DIR + cada condición):

    BASE_DIR/
      <TRATAMIENTO>/
        RAW/<CARPETA_FLUOROCROMO>/                imágenes .tif de entrada
        RAW/<CARPETA_FLUOROCROMO>/Salida/         overlays + results_cells.csv
        RAW/<CARPETA_FLUOROCROMO>/Salida/Recortes/ recortes _raw.tif y _vis.png
        MSSR/<CARPETA_FLUOROCROMO>/               salidas MSSR + peaks/summary.csv
        HEATMAPS_MSSR/<CARPETA_FLUOROCROMO>/       heatmaps

Esta es exactamente la jerarquía de tu proyecto NuevoAnalisis-BOD_INT.
"""

import os

# ============================================================
# 1. RAÍZ DEL PROYECTO  ←  EDITAR PARA UN NUEVO DATASET
# ============================================================

# En Colab apunta a Google Drive; en local apunta a una carpeta del disco.
# Raíz donde la etapa 0 construye la jerarquía de análisis (se crea sola).
BASE_DIR = "/content/drive/MyDrive/2025 Josue Villegas/LNMA 2026/260327/Analisis_MSSR"

# Carpeta plana con los .oib crudos del confocal (entrada de la etapa 0).
RAW_OIB_DIR = "/content/drive/MyDrive/2025 Josue Villegas/LNMA 2026/260327/260327"

# Si True, run_pipeline montará Google Drive antes de correr.
USE_COLAB = True

# ------------------------------------------------------------
# Conversión .oib (etapa 0)
# ------------------------------------------------------------
# Los .oib de 260327 traen 2 canales: fluorescencia + luz de transmisión.
# La etapa 0 conserva SOLO el de fluorescencia.
#   None → detección automática por archivo (canal con mayor asimetría:
#          señal dispersa y brillante sobre fondo oscuro).
#   0 / 1 → fija el índice del canal de fluorescencia (recomendado una vez
#          confirmado con el ensayo en seco, para evitar cualquier ambigüedad).
OIB_FLUORESCENCE_CHANNEL = None

# Tratamiento asumido cuando el nombre no lleva token 'BOD' ni 'INT'.
# En 260327 las muestras INT se nombran sin token de tratamiento (p. ej.
# 13S-488-01), por lo que ausencia de 'BOD' = INT.
DEFAULT_TREATMENT_WHEN_UNLABELED = "INT"

# ============================================================
# 2. CALIBRACIÓN ESPACIAL
# ============================================================

PIXELS_PER_MICRON = 18.8679                       # calibración del microscopio
PIXEL_SIZE_UM = 1.0 / PIXELS_PER_MICRON
PIXEL_SIZE_NM = 1000.0 / PIXELS_PER_MICRON

# ============================================================
# 3. PARÁMETROS MSSR (etapa 2)
# ============================================================

MSSR_AMP = 1        # factor de amplificación
MSSR_FWHM = 5       # FWHM del PSF en píxeles
MSSR_ORDER = 0      # orden del realce

# Sufijo de los archivos MSSR generados, derivado de los parámetros.
MSSR_SUFFIX = f"_MSSR_order{MSSR_ORDER}_amp{MSSR_AMP}"

# ============================================================
# 4. SEGMENTACIÓN (etapa 1)
# ============================================================

CROP_SIZE = 96               # tamaño del recorte cuadrado (px)
MIN_AREA = 60                # área mínima de un objeto para conservarlo (px)
# Filtro anti-ruido a nivel de imagen: descarta campos con más objetos que este
# valor. Era 20 para el experimento previo (pocas células por campo). Para
# campos densos de procarioplancton (cientos de células) se desactiva con None.
MAX_OBJECTS_PER_IMAGE = None
SCALE_BAR_UM = 1             # longitud de la barra de escala en visualizaciones

# Sustracción de fondo para SEGMENTAR (no altera los recortes crudos de MSSR).
# Corrige el fondo de fluorescencia alto de las muestras sin INT (BOD), que de
# otro modo Otsu sobre-segmenta. White top-hat con elemento estructurante de
# radio TOPHAT_RADIUS px: debe ser mayor que el radio celular y menor que la
# escala del fondo. Ajustar revisando los overlays.
BACKGROUND_SUBTRACTION = True
TOPHAT_RADIUS = 25

# Control de calidad por relación señal/fondo (SNR) por imagen.
#   None  → solo reporta el SNR, no excluye.
#   float → excluye imágenes con SNR por debajo del umbral (fallos de adquisición).
# En 260327 las imágenes buenas tienen SNR > 13 y las fallidas SNR < 5.3, con un
# hueco limpio entre ambos grupos; 8 separa sin ambigüedad. (Alternativa por
# conteo: poner QC_MIN_SNR = None y MAX_OBJECTS_PER_IMAGE = 3.)
QC_MIN_SNR = 8

# ------------------------------------------------------------
# Filtro de objetos a nivel individual (etapa 1, modo estricto)
# ------------------------------------------------------------
# Conserva solo objetos claramente celulares. Se aplica sobre la imagen CRUDA
# segmentada (el objeto físico que pasó el prefiltrado de 1 µm de poro), no
# sobre MSSR. Cada criterio mapea a una propiedad de skimage.regionprops:
#   área           → 'area'
#   solidez        → 'solidity'      (área / área del casco convexo)
#   excentricidad  → 'eccentricity'  (0 = círculo, 1 = línea)
OBJECT_FILTER = True

# Techo de tamaño: el prefiltrado físico es de 1 µm de poro. Se fija a 1.3 µm
# para no recortar células de ~1 µm que, por difracción, se ven algo mayores en
# RAW antes de MSSR (justamente el efecto que mide la etapa 3).
CELL_MAX_DIAMETER_UM = 1.3
# Área máxima derivada del diámetro máximo (círculo equivalente), en px².
OBJECT_MAX_AREA = 3.14159265 * (0.5 * CELL_MAX_DIAMETER_UM * PIXELS_PER_MICRON) ** 2
# Área mínima: se reutiliza MIN_AREA (≈0.46 µm de diámetro).
OBJECT_MIN_AREA = MIN_AREA
# Forma (modo estricto): compactas y convexas.
OBJECT_MIN_SOLIDITY = 0.90
OBJECT_MAX_ECCENTRICITY = 0.85

# ============================================================
# 5. DETECCIÓN DE PICOS Y FWHM (etapa 3)
# ============================================================

PEAK_MIN_DISTANCE = 1        # distancia mínima entre picos (px)
PEAK_THRESHOLD_REL = 0.3     # umbral relativo para peak_local_max
FWHM_WINDOW = 8              # semiventana del perfil para el ajuste gaussiano (px)
FWHM_MAXFEV = 800            # iteraciones máximas de curve_fit
FWHM_MIN_NM = 1              # rango físico aceptable del FWHM ajustado
FWHM_MAX_NM = 1000

# Medición de tamaño celular pre/post MSSR (etapa 3). El borde se define por
# media altura (medio máximo dentro de cada recorte), criterio autoescalante y
# por tanto comparable entre RAW y MSSR pese al reescalado de intensidad de MSSR.
# Se reportan dos diámetros por objeto y por imagen (RAW y MSSR):
#   - diámetro de círculo equivalente: 2*sqrt(area_half_max/pi)
#   - diámetro de Feret (máximo): eje mayor del objeto a media altura
MEASURE_CELL_SIZE = True
CELL_SIZE_MIN_AREA_PX = 20   # área mínima (px) de la región a media altura para medir

# ============================================================
# 6. ÓPTICA / LÍMITE DE DIFRACCIÓN (etapa 4)
# ============================================================

NUMERICAL_APERTURE = 1.3     # apertura numérica del objetivo
RAYLEIGH_FACTOR = 0.61       # criterio de Rayleigh: d = 0.61 * lambda / NA

# ============================================================
# 7. DISEÑO EXPERIMENTAL  ←  EDITAR SI CAMBIAN LAS CONDICIONES
# ============================================================
# Cada condición combina un tratamiento y un fluorocromo. 'dye_folder' es el
# nombre literal de la carpeta en disco (ojo: en tu proyecto es "SYBER"),
# mientras que 'dye' es la etiqueta usada en tablas y figuras ("SYBR").
# 'emission_nm' es la longitud de onda de emisión para el límite de Rayleigh.

# Diseño 2x2: tratamiento (INT/BOD) x fluorocromo (DAPI/SYBR).
# La etapa 0 deduce la condición del nombre del .oib y crea las carpetas, por lo
# que 'dye_folder' puede ser un nombre limpio. 'emission_nm' alimenta el límite
# de Rayleigh (etapa 4); son las emisiones aproximadas, no la excitación.
#
# Codificación del nombre confirmada para 260327:
#   <muestra><letra>-<tratamiento>-<excitacion>-<idx>.oib
#     letra        : D = DAPI, S = SYBR          (último carácter del token 0)
#     tratamiento  : BOD = sin INT, INT = con INT (token 1)
#     excitacion   : 405 = DAPI, 488 = SYBR       (validación cruzada)
#   Ej.: 13D-BOD-405-01 → BOD+DAPI ; 13S-BOD-488-01 → BOD+SYBR

CONDITIONS = [
    {"treatment": "INT", "dye": "DAPI", "dye_folder": "DAPI", "emission_nm": 460},
    {"treatment": "INT", "dye": "SYBR", "dye_folder": "SYBR", "emission_nm": 520},
    {"treatment": "BOD", "dye": "DAPI", "dye_folder": "DAPI", "emission_nm": 460},
    {"treatment": "BOD", "dye": "SYBR", "dye_folder": "SYBR", "emission_nm": 520},
]

# Factores del ANOVA. Si el diseño deja de ser 2x2, ajusta esta lista.
ANOVA_FACTORS = ["treatment", "dye"]

# ============================================================
# 8. CONSTRUCTORES DE RUTAS  (no suele requerir edición)
# ============================================================

def label(condition):
    """Etiqueta legible de un tratamiento, p. ej. 'INT-DAPI'."""
    return f"{condition['treatment']}-{condition['dye']}"

def raw_dir(condition):
    return os.path.join(BASE_DIR, condition["treatment"], "RAW", condition["dye_folder"])

def salida_dir(condition):
    return os.path.join(raw_dir(condition), "Salida")

def crops_dir(condition):
    return os.path.join(salida_dir(condition), "Recortes")

def mssr_dir(condition):
    return os.path.join(BASE_DIR, condition["treatment"], "MSSR", condition["dye_folder"])

def heatmap_dir(condition):
    return os.path.join(BASE_DIR, condition["treatment"], "HEATMAPS_MSSR", condition["dye_folder"])

def peaks_csv(condition):
    return os.path.join(mssr_dir(condition), "peaks_results.csv")

def summary_csv(condition):
    return os.path.join(mssr_dir(condition), "cell_summary.csv")

def rayleigh_limit(condition):
    """Límite de resolución de Rayleigh en nm para la condición."""
    return RAYLEIGH_FACTOR * condition["emission_nm"] / NUMERICAL_APERTURE

def ensure_dirs(condition):
    """Crea las carpetas de salida de una condición si no existen."""
    for d in (salida_dir(condition), crops_dir(condition),
              mssr_dir(condition), heatmap_dir(condition)):
        os.makedirs(d, exist_ok=True)
