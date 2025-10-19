import sys
import numpy as np
from pathlib import Path

sys.path.append(r"C:\Users\jorge\PycharmProjects\EOPD\EODP_CODE\test_eodp\common\io")
from writeToa import readToa


def to_1d_spectrum(arr):
    """
    Convierte arr a un vector 1D de longitud de onda:
    - si es xarray.DataArray -> usa .values
    - si es Dataset -> coge la primera variable
    - si es 2D -> media sobre la primera dimensión
    - si es 3D -> media sobre las dos primeras dimensiones
    """
    try:
        import xarray as xr
        if isinstance(arr, xr.DataArray) or isinstance(arr, xr.Dataset):
            if isinstance(arr, xr.Dataset):
                var = list(arr.data_vars)[0]
                arr = arr[var].values
            else:
                arr = arr.values
    except ImportError:
        pass

    arr = np.asarray(arr)
    if arr.ndim == 1:
        return arr
    elif arr.ndim == 2:
        return np.nanmean(arr, axis=0)
    elif arr.ndim == 3:
        return np.nanmean(arr, axis=(0, 1))
    else:
        raise ValueError(f"Array con ndim inesperado: {arr.ndim}")


def diff_band(l1b_toa, toa_ref):
    """Calcula μ + 3σ del error relativo en %"""
    diff = np.abs(l1b_toa - toa_ref) / np.maximum(toa_ref, 1e-12) * 100
    mu = np.nanmean(diff)
    sigma = np.nanstd(diff)
    threshold = mu + 3 * sigma
    return threshold, threshold < 0.01  # retorna si cumple <0.01%


def check_bands_against_isrf(bands, base_paths):
    """Comprueba μ + 3σ < 0.01% con respecto a TOA"""
    print("=== COMPROBACIÓN μ+3σ < 0.01% vs OPTICAL ===")
    for band in bands:
        toa_isrf = to_1d_spectrum(readToa(base_paths["optical"], f"ism_toa_optical_{band}.nc"))
        toa_referencia = to_1d_spectrum(readToa(base_paths["referencia"], f"ism_toa_optical_{band}.nc"))
        threshold, ok = diff_band(toa_referencia, toa_isrf)
        print(f"{band}: μ+3σ = {threshold:.6f}%  -> {'CUMPLE' if ok else 'NO CUMPLE'}")


def main():
    bands = ["VNIR-0", "VNIR-1", "VNIR-2", "VNIR-3"]
    base_paths = {
        "optical": Path(r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-ISM\jorgeoutputs"),
        "referencia": Path(r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-ISM\output")
    }

    check_bands_against_isrf(bands, base_paths)


if __name__ == "__main__":
    main()
