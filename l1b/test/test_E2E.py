import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# reader
sys.path.append(r"C:\Users\jorge\PycharmProjects\EOPD\EODP_CODE\test_eodp\common\io")
from writeToa import readToa

BANDS = ["VNIR-0", "VNIR-1", "VNIR-2", "VNIR-3"]

PATHS = {
    "isrf": Path(r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-E2E\jorgeoutputs_END2END"),
    "l1b":  Path(r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-E2E\l1b_jorgeoutputs"),
}

# --- helpers ---
def central_alt_line(arr):
    """Devuelve la línea ALT central como vector 1D."""
    a = np.asarray(arr)
    # asegura 2D: (ALT, ACT)
    if a.ndim == 1:
        return a
    if a.ndim == 3:
        a = a.squeeze()
    if a.shape[0] < a.shape[1]:
        # típico: ALT x ACT (ok)
        line = a.shape[0] // 2
        v = a[line, :]
    else:
        # si viniera traspuesto, invierte
        line = a.shape[1] // 2
        v = a[:, line]
    return np.asarray(v, dtype=float)

def clean_and_crop(v, crop=2):
    """Quita no finitos y recorta bordes para evitar efectos del kernel."""
    v = np.asarray(v, dtype=float)
    v[~np.isfinite(v)] = np.nan
    if crop > 0 and v.size > 2*crop:
        v = v[crop:-crop]
    return v

def sort_pair(a, b):
    """Ordena ambas series con la misma permutación (por la serie de referencia)."""
    idx = np.argsort(a)
    return a[idx], b[idx]

def rel_stats(ref, est):
    """μ, σ y μ+3σ del error relativo en % (enmascara ref≈0)."""
    ref = np.asarray(ref, dtype=float)
    est = np.asarray(est, dtype=float)
    mask = np.isfinite(ref) & np.isfinite(est) & (np.abs(ref) > 1e-12)
    if not np.any(mask):
        return np.nan, np.nan, np.nan
    diff = np.abs(est[mask] - ref[mask]) / np.abs(ref[mask]) * 100.0
    mu = float(np.nanmean(diff))
    sigma = float(np.nanstd(diff))
    return mu, sigma, mu + 3*sigma

def compare_band(band, axs):
    # lee
    toa_isrf = readToa(PATHS["isrf"], f"ism_toa_isrf_{band}.nc")
    toa_l1b  = readToa(PATHS["l1b"],  f"l1b_toa_{band}.nc")

    # línea central
    y_isrf = central_alt_line(toa_isrf)
    y_l1b  = central_alt_line(toa_l1b)

    # limpia y recorta bordes
    y_isrf = clean_and_crop(y_isrf, crop=3)
    y_l1b  = clean_and_crop(y_l1b, crop=3)

    # iguala longitudes si difieren por 1–2 px (bordes)
    n = min(y_isrf.size, y_l1b.size)
    y_isrf, y_l1b = y_isrf[:n], y_l1b[:n]

    # ordena ambas (el test permite distinto orden de ejecución)
    x = np.arange(n)
    # si quieres estrictamente ordenar por amplitud:
    y_isrf_sorted, y_l1b_sorted = sort_pair(y_isrf, y_l1b)

    # stats
    mu, sigma, thr = rel_stats(y_isrf_sorted, y_l1b_sorted)

    # plot
    ax = axs
    ax.plot(y_isrf_sorted, label="TOA after ISRF", linewidth=2)
    ax.plot(y_l1b_sorted,  label="L1B TOA radiance", linewidth=2, linestyle="--")
    ax.plot(y_l1b_sorted - y_isrf_sorted, label="Difference (L1B−ISRF)", linewidth=1, alpha=0.8)
    ax.set_title(f"{band}  |  μ={mu:.4f}%  σ={sigma:.4f}%  μ+3σ={thr:.4f}%")
    ax.set_xlabel("Sample index (central ALT line, sorted)")
    ax.set_ylabel("TOA Radiance [mW·m⁻²·sr⁻¹]")
    ax.grid(True)
    return mu, sigma, thr

def main():
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)
    axes = axes.ravel()
    summary = []
    for i, b in enumerate(BANDS):
        mu, sigma, thr = compare_band(b, axes[i])
        summary.append((b, mu, sigma, thr))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3)
    plt.suptitle("Central ALT line – L1B vs ISRF (sorted, cropped edges)")
    plt.show()

    print("\n=== Relative error stats (L1B vs ISRF) ===")
    print(f"{'Band':<8}{'mu %':>10}{'sigma %':>12}{'mu+3σ %':>12}")
    for b, mu, s, t in summary:
        print(f"{b:<8}{mu:10.6f}{s:12.6f}{t:12.6f}")

if __name__ == "__main__":
    main()
