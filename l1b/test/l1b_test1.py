import sys
import numpy as np
import matplotlib.pyplot as plt
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


def plot_bands(bands, base_paths):
    """Dibuja las 4 bandas en subplots, comparando EQ=True, EQ=False e ISRF"""
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), constrained_layout=True)
    axes = axes.flatten()

    for i, band in enumerate(bands):
        toa_isrf = readToa(base_paths["isrf"], f"ism_toa_isrf_{band}.nc")
        toa_true = readToa(base_paths["eq_true"], f"l1b_toa_{band}.nc")
        toa_false = readToa(base_paths["eq_false"], f"l1b_toa_{band}.nc")

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

    fig.legend(["TOA after ISRF", "Restored (EQ=True)", "Restored (EQ=False)"],
               loc="upper center", ncol=3)
    plt.show()


def diff_band(l1b_toa, l1b_toa_original):
    """Calcula μ + 3σ del error relativo en % y lo compara con 0.01%"""
    tol = 0.01  # en %
    diff = np.abs(l1b_toa_original - l1b_toa) / np.maximum(l1b_toa_original, 1e-12) * 100
    mu = np.nanmean(diff)
    sigma = np.nanstd(diff)
    threshold = mu + 3 * sigma
    return threshold, threshold < tol


def check_bands(bands, base_paths):
    """Ejecuta el test μ+3σ < 0.01% para todas las bandas"""
    print("=== COMPROBACIÓN μ+3σ < 0.01% ===")
    for band in bands:
        l1b = readToa(base_paths["eq_true"], f"l1b_toa_{band}.nc")
        l1b_ref = readToa(base_paths["original"], f"l1b_toa_{band}.nc")
        threshold, ok = diff_band(l1b, l1b_ref)
        print(f"{band}: μ+3σ = {threshold:.6f}%  -> {'CUMPLE' if ok else 'NO CUMPLE'}")


def main():
    bands = ["VNIR-0", "VNIR-1", "VNIR-2", "VNIR-3"]
    base_paths = {
        "isrf": Path(r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-ISM\output"),
        "eq_true": Path(r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-L1B\jorgeoutputs"),
        "eq_false": Path(r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-L1B\jorgeoutputs_SIN_EQ"),
        "original": Path(r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-L1B\output")
    }

    plot_bands(bands, base_paths)
    check_bands(bands, base_paths)


if __name__ == "__main__":
    main()


