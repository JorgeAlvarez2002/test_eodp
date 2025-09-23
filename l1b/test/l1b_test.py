import sys, os
import numpy as np
import matplotlib.pyplot as plt

# si necesitas forzar el path para encontrar common.io.writeToa
sys.path.append(r"C:\Users\jorge\PycharmProjects\EODP\EODP_CODE\test_eodp")
from common.io.writeToa import readToa

def to_1d_spectrum(arr):
    """
    Convierte arr a un vector 1D de longitud de onda:
    - si es xarray.DataArray -> usa .values
    - si es 1D ya lo devuelve
    - si es 2D -> media sobre la primera dimensión (espacial)
    - si es 3D -> media sobre las dos primeras dimensiones (p. ej. ALT,CROSS)
    """
    # si es xarray
    try:
        import xarray as xr
        if isinstance(arr, xr.DataArray) or isinstance(arr, xr.Dataset):
            # si es Dataset intentar coger variable obvia, si no, convertir a numpy
            if isinstance(arr, xr.Dataset):
                # intentar escoger la primera variable del dataset
                var = list(arr.data_vars)[0]
                arr = arr[var].values
            else:
                arr = arr.values
    except Exception:
        pass

    arr = np.asarray(arr)
    if arr.ndim == 1:
        return arr
    elif arr.ndim == 2:
        # forma típica: [spatial, wavelength] -> promediamos spatial
        return np.mean(arr, axis=0)
    elif arr.ndim == 3:
        # forma típica: [alt, cross, wavelength] -> promediamos alt y cross
        return np.mean(arr, axis=(0, 1))
    else:
        raise ValueError(f"Array con ndim inesperado: {arr.ndim}")


bands = ["VNIR-0", "VNIR-1", "VNIR-2", "VNIR-3"]
fig, axes = plt.subplots(2, 2, figsize=(10, 6), constrained_layout=True)
axes = axes.flatten()

for i, band in enumerate(bands):
    # cargar
    toa_isrf = readToa(r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-ISM\output", f"ism_toa_isrf_{band}.nc")
    toa_true = readToa(r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-L1B\jorgeoutputs", f"l1b_toa_{band}.nc")
    toa_false = readToa(r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-L1B\jorgeoutputs_SIN_EQ", f"l1b_toa_{band}.nc")

    # convertir a espectros 1D
    spec_isrf = to_1d_spectrum(toa_isrf)
    spec_true = to_1d_spectrum(toa_true)
    spec_false = to_1d_spectrum(toa_false)

    ax = axes[i]
    ax.plot(spec_isrf, label="TOA after ISRF", linewidth=2)
    ax.plot(spec_true, label="Restored (EQ=True)", linewidth=2)
    ax.plot(spec_false, label="Restored (EQ=False)", linewidth=2)
    ax.set_title(band)
    ax.set_xlabel("Wavelength index")
    ax.set_ylabel("TOA Radiance")
    ax.grid(True)

fig.legend(
    ["TOA after ISRF", "Restored (EQ=True)", "Restored (EQ=False)"],
    loc="upper center", ncol=3
)
plt.show()
###############################################################################
# Ejercicio 1
###############################################################################
# LISTA DE BANDAS
bands = ["VNIR-0", "VNIR-1", "VNIR-2", "VNIR-3"]

# LEER TODOS LOS DATOS USANDO BUCLES
toa_true = {band: readToa(r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-L1B\jorgeoutputs",
                          f"l1b_toa_{band}.nc")
            for band in bands}

toa_original = {band: readToa(r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-L1B\output",
                              f"l1b_toa_{band}.nc")
                for band in bands}

toa_noeq = {band: readToa(r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-L1B\jorgeoutputs_SIN_EQ",
                          f"l1b_toa_{band}.nc")
            for band in bands}

# TEST 1: DIFERENCIAS (μ + 3σ)
def diff_band(l1b_toa, l1b_toa_original):
    tol = 0.01  # en %
    diff = np.abs(l1b_toa_original - l1b_toa) / np.maximum(l1b_toa_original, 1e-12) * 100
    mu = np.mean(diff)
    sigma = np.std(diff)
    threshold = mu + 3 * sigma  # valor hasta 3σ
    ok = threshold < tol
    return threshold, ok

print("=== COMPROBACIÓN μ+3σ < 0.01% ===")
for band in bands:
    threshold, ok = diff_band(toa_true[band], toa_original[band])
    print(f"{band}: μ+3σ = {threshold:.6f}%  -> {'✅ CUMPLE' if ok else '❌ NO CUMPLE'}")

