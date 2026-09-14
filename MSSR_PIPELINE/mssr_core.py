# -*- coding: utf-8 -*-
"""
mssr_core.py
============
Implementación de Mean-Shift Super-Resolution (MSSR).

Las funciones espaciales (bicInter, ftInterp, meshing, sfMSSR, apply_mssr)
son tu código original sin modificaciones algorítmicas: el comportamiento es
idéntico al de MultiparametricoMSSRFernando.

NOTA IMPORTANTE sobre tMSSR / tMean
-----------------------------------
Tu código original importaba tMSSR y tMean desde un 'mssr_module' que vive en
tu Drive y que NO está entre los archivos compartidos. Aquí se incluye una
reconstrucción estándar (MSSR temporal = sfMSSR por fotograma + promedio),
consistente con la formulación publicada del método. Solo se usa para pilas
multi-fotograma (image.ndim == 3 con pocos planos). Si tu mssr_module define
estas funciones de otra forma, reemplaza este bloque o ajusta el import en
stage2_mssr.py para usar tu módulo oficial.
"""

import math
import numpy as np
from scipy.interpolate import RectBivariateSpline


# ------------------------------------------------------------
# Interpolación bicúbica
# ------------------------------------------------------------
def bicInter(img, amp, mesh):
    width, height = img.shape
    y = np.linspace(1, width, width)
    x = np.linspace(1, height, height)
    imgInter = RectBivariateSpline(y, x, img, kx=3, ky=3)
    y2 = np.linspace(1, width, width * amp)
    x2 = np.linspace(1, height, height * amp)
    Z2 = imgInter(y2, x2)
    if mesh:
        Z2 = meshing(Z2, amp)
    return Z2


# ------------------------------------------------------------
# Interpolación de Fourier
# ------------------------------------------------------------
def ftInterp(img, amp, mesh):
    width, height = img.shape
    mdX = math.ceil(width / 2) + 1
    mdY = math.ceil(height / 2) + 1
    extraBorder = math.ceil(amp / 2)
    Nwidth = (width * amp) + extraBorder
    Nheight = (height * amp) + extraBorder
    lnX = len(np.arange(mdX, width))
    lnY = len(np.arange(mdY, height))
    imgFt = np.fft.fft2(img)
    imgFt = imgFt * (Nwidth / width) * (Nheight / height)
    fM = np.zeros((Nwidth, Nheight), dtype=complex)
    fM[0:mdX, 0:mdY] = imgFt[0:mdX, 0:mdY]
    fM[0:mdX, (Nheight - lnY):Nheight] = imgFt[0:mdX, mdY:height]
    fM[(Nwidth - lnX):Nwidth, 0:mdY] = imgFt[mdX:width, 0:mdY]
    fM[(Nwidth - lnX):Nwidth, (Nheight - lnY):Nheight] = imgFt[mdX:width, mdY:height]
    Z2 = (np.fft.ifft2(fM)).real
    Z2 = Z2[0:(width * amp), 0:(height * amp)]
    if mesh:
        Z2 = meshing(Z2, amp)
    return Z2


# ------------------------------------------------------------
# Compensación de malla
# ------------------------------------------------------------
def meshing(img, amp):
    width, height = img.shape
    desp = math.ceil(amp / 2)
    imgPad = np.pad(img, desp, 'symmetric')
    imgS1 = imgPad[0:width, desp:height + desp]
    imgS2 = imgPad[(desp * 2):width + (desp * 2), desp:height + desp]
    imgS3 = imgPad[desp:width + desp, 0:height]
    imgS4 = imgPad[desp:width + desp, (desp * 2):height + (desp * 2)]
    imgF = (img + imgS1 + imgS2 + imgS3 + imgS4) / 5
    return imgF


# ------------------------------------------------------------
# MSSR espacial
# ------------------------------------------------------------
def sfMSSR(img, fwhm, amp, order, mesh=True, ftI=False, intNorm=True):
    hs = round(0.5 * fwhm * amp)
    if hs < 1:
        hs = 1
    if amp > 1 and not ftI:
        img = bicInter(img, amp, mesh)
    elif amp > 1 and ftI:
        img = ftInterp(img, amp, mesh)
    width, height = img.shape
    xPad = np.pad(img, hs, 'symmetric')
    M = np.zeros((width, height))
    for i in range(-hs, hs + 1):
        for j in range(-hs, hs + 1):
            if i != 0 or j != 0:
                xThis = xPad[hs + i:width + hs + i, hs + j:height + hs + j]
                M = np.maximum(M, np.abs(img - xThis))
    weightAccum = np.zeros((width, height))
    yAccum = np.zeros((width, height))
    for i in range(-hs, hs + 1):
        for j in range(-hs, hs + 1):
            if i != 0 or j != 0:
                spatialkernel = np.exp(-(i ** 2 + j ** 2) / (hs ** 2))
                xThis = xPad[hs + i:width + hs + i, hs + j:height + hs + j]
                xDiffSq0 = ((img - xThis) / M) ** 2
                intensityKernel = np.exp(-xDiffSq0)
                weightThis = spatialkernel * intensityKernel
                weightAccum += weightThis
                yAccum += xThis * weightThis
    MS = img - (yAccum / weightAccum)
    MS[MS < 0] = 0
    MS[np.isnan(MS)] = 0
    I3 = MS / np.max(MS)
    x3 = img / np.max(img)
    for _ in range(order):
        I4 = x3 - I3
        I5 = np.max(I4) - I4
        I5 = I5 / np.max(I5)
        I6 = I5 * I3
        I7 = I6 / np.max(I6)
        x3 = I3
        I3 = I7
    I3[np.isnan(I3)] = 0
    IMSSR = I3 * img if intNorm else I3
    return IMSSR


def apply_mssr(image, amp=1, fwhm=5, order=0):
    """Wrapper de MSSR espacial sobre una imagen 2D."""
    return sfMSSR(image, fwhm=fwhm, amp=amp, order=order)


# ------------------------------------------------------------
# MSSR temporal  (RECONSTRUCCIÓN — ver nota de cabecera)
# ------------------------------------------------------------
def tMSSR(stack, fwhm, amp, order):
    """Aplica MSSR espacial a cada fotograma de una pila (frames, H, W)."""
    return np.stack([sfMSSR(frame, fwhm=fwhm, amp=amp, order=order)
                     for frame in stack], axis=0)


def tMean(stack):
    """Proyección de media temporal de una pila MSSR."""
    return np.mean(stack, axis=0)
