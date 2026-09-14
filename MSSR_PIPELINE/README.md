# Pipeline MSSR — análisis cuantitativo de super-resolución

Reorganización de la rutina de análisis de imágenes (originalmente siete
*notebooks* de Colab) en un paquete modular con configuración centralizada.
El objetivo es correr el flujo completo sobre un nuevo conjunto de imágenes
editando un único archivo.

## Inicio rápido — experimento 260327 (.oib)

`config.py` ya viene apuntado a este experimento:
`RAW_OIB_DIR` = carpeta plana con los 100 `.oib`; `BASE_DIR` = `.../260327/Analisis_MSSR`
(la etapa 0 la crea). En Colab, tras copiar el paquete a Drive y añadir su ruta
con `sys.path.append(...)`:

```python
import stage0_convert as s0

# 1) Ensayo en seco: clasifica y muestra la forma de cada .oib, sin escribir.
s0.convert_all(dry_run=True)

# 2) Si el resumen por condición y las formas son correctos, convierte de verdad.
s0.convert_all(dry_run=False)
```

Confirma en el ensayo en seco que (a) los 100 archivos se reparten entre las
cuatro condiciones y ninguno cae en "sin clasificar", y (b) las formas de
entrada son monocanal (un solo canal por archivo); si algún `.oib` trae dos
canales, la proyección de máxima intensidad los mezclaría y habría que separar
canales antes. Hecho esto, corre el resto: `python run_pipeline.py 1 2 3 4 5`.

## Estructura del paquete

| Archivo | Etapa | Función |
|---|---|---|
| `config.py` | — | Rutas, calibración, parámetros y diseño experimental. **Único archivo a editar.** |
| `mssr_core.py` | — | Algoritmo MSSR (espacial + temporal). |
| `stage0_convert.py` | 0 | Conversión por lote `.oib`→`.tif` y archivado por condición. |
| `stage1_segmentation.py` | 1 | Segmentación Otsu y extracción de recortes 96×96. |
| `stage2_mssr.py` | 2 | Aplicación de MSSR a los recortes. |
| `stage3_parameters.py` | 3 | FWHM (ajuste gaussiano), detección de picos, vecino más cercano. |
| `stage4_statistics.py` | 4 | Descriptivos, ANOVA factorial, tests pareados, Bland–Altman, Rayleigh. |
| `stage5_heatmaps.py` | 5 | Heatmaps de intensidad con escala global común. |
| `run_pipeline.py` | — | Orquestador (monta Drive y ejecuta las etapas). |

El flujo es secuencialmente dependiente: **0 → 1 → 2 → 3 → {4, 5}**. La etapa 0
solo es necesaria si se parte de archivos crudos `.oib`; si ya hay `.tif`
organizados en la jerarquía, se empieza en la etapa 1.

## Cómo correrlo sobre un nuevo conjunto de imágenes

1. **Coloca las imágenes RAW** siguiendo la jerarquía de carpetas (ver abajo).
2. **Edita `config.py`**: como mínimo `BASE_DIR`; si cambia el diseño, también
   `CONDITIONS`.
3. **Ejecuta**:

   ```bash
   python run_pipeline.py        # pipeline completo
   python run_pipeline.py 1      # solo segmentación
   python run_pipeline.py 3 4    # parámetros + estadística
   ```

   En Colab, copia el paquete a Drive, añade su ruta con
   `sys.path.append(...)` y llama a `run_pipeline.main([...])`, o ejecuta cada
   `stageN.run_all()` por celda.

## Jerarquía de carpetas esperada

Construida automáticamente a partir de `BASE_DIR` y cada condición:

```
BASE_DIR/
  <TRATAMIENTO>/
    RAW/<CARPETA_FLUOROCROMO>/                ← imágenes .tif de entrada
    RAW/<CARPETA_FLUOROCROMO>/Salida/          overlays + results_cells.csv
    RAW/<CARPETA_FLUOROCROMO>/Salida/Recortes/ recortes _raw.tif y _vis.png
    MSSR/<CARPETA_FLUOROCROMO>/                salidas MSSR + peaks/summary.csv
    HEATMAPS_MSSR/<CARPETA_FLUOROCROMO>/       heatmaps
  STATS/                                       figuras de la etapa 4
```

Solo necesitas crear y poblar las carpetas `RAW/<fluorocromo>/`; el resto se
genera. Para un diseño distinto al 2×2 (INT/BOD × DAPI/SYBR), redefine la lista
`CONDITIONS`: admite cualquier número de condiciones. Para un único grupo de
imágenes, basta una sola entrada.

## Parámetros principales (en `config.py`)

| Parámetro | Valor por defecto | Significado |
|---|---|---|
| `PIXELS_PER_MICRON` | 18.8679 | Calibración espacial. |
| `MSSR_AMP / FWHM / ORDER` | 1 / 5 / 0 | Parámetros del realce MSSR. |
| `CROP_SIZE` | 96 | Lado del recorte (px). |
| `MIN_AREA` | 60 | Área mínima de objeto (px). |
| `MAX_OBJECTS_PER_IMAGE` | 20 | Filtro anti-ruido por campo. |
| `PEAK_THRESHOLD_REL` | 0.3 | Umbral relativo de detección de picos. |
| `FWHM_WINDOW` | 8 | Semiventana del perfil para el ajuste. |
| `NUMERICAL_APERTURE` | 1.3 | NA del objetivo (límite de Rayleigh). |

## Cambios respecto a los notebooks originales

- Rutas centralizadas: se elimina el *hardcoding* repetido en cada script.
- Las tres rutinas de estadística se unifican en `stage4_statistics.py`, con
  ANOVA factorial que se adapta al número de factores con varios niveles.
- **Corrección**: el emparejamiento RAW→MSSR en la etapa 3 ahora exige
  extensión `.tif`, evitando seleccionar por error los PNG de salida.
- Los nombres de archivo MSSR se derivan de `MSSR_SUFFIX` en vez de fijar
  `order0_amp1` a mano.
- Las figuras de la etapa 4 se guardan en disco (`BASE_DIR/STATS`) para
  reproducibilidad.

## Pendiente de verificación

- **`mssr_module` / MSSR temporal.** El módulo original de tu Drive no estaba
  entre los archivos compartidos. `mssr_core.py` conserva tu MSSR espacial sin
  cambios e incluye una reconstrucción estándar de `tMSSR`/`tMean` (MSSR por
  fotograma + media), usada solo para pilas multi-fotograma. Si tu módulo
  oficial difiere, reemplaza ese bloque o ajusta el import en `stage2_mssr.py`.
