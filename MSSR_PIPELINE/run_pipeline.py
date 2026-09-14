# -*- coding: utf-8 -*-
"""
run_pipeline.py
===============
Orquestador del pipeline MSSR completo.

Uso (Colab o local):

    # todo el pipeline sobre las condiciones de config.py
    python run_pipeline.py

    # etapas sueltas
    python run_pipeline.py 0        # solo conversión .oib→.tif
    python run_pipeline.py 3 4      # parámetros + estadística

Las etapas son secuencialmente dependientes:
    0 conversión → 1 segmentación → 2 MSSR → 3 parámetros → 4 estadística
                                                         └→ 5 heatmaps
"""

import sys
import config as cfg


def mount_drive():
    """Monta Google Drive si se ejecuta en Colab y USE_COLAB=True."""
    if not cfg.USE_COLAB:
        return
    try:
        from google.colab import drive
        drive.mount('/content/drive')
    except ImportError:
        print("No es un entorno Colab; se omite el montaje de Drive.")


STAGES = {
    "0": ("Conversión .oib→.tif", lambda: __import__("stage0_convert").run_all()),
    "1": ("Segmentación", lambda: __import__("stage1_segmentation").run_all()),
    "2": ("MSSR",         lambda: __import__("stage2_mssr").run_all()),
    "3": ("Parámetros",   lambda: __import__("stage3_parameters").run_all()),
    "4": ("Estadística",  lambda: __import__("stage4_statistics").run_all()),
    "5": ("Heatmaps",     lambda: __import__("stage5_heatmaps").run_all()),
}


def main(stage_args):
    mount_drive()
    stages = stage_args if stage_args else list(STAGES.keys())
    for s in stages:
        if s not in STAGES:
            print(f"Etapa desconocida: {s}")
            continue
        name, fn = STAGES[s]
        print(f"\n{'=' * 60}\nETAPA {s} — {name}\n{'=' * 60}")
        fn()
    print("\nPipeline finalizado.")


if __name__ == "__main__":
    main(sys.argv[1:])
