import sys
from pathlib import Path
import numpy as np

# === Ajusta estas rutas a tu entorno ===
# Módulo lector de TOA (el mismo que ya usas)
sys.path.append(r"C:\Users\jorge\PycharmProjects\EOPD\EODP_CODE\test_eodp\common\io")
from writeToa import readToa

# Directorio con tus salidas tras VCU (ficheros ism_toa_<BAND>.nc)
OUT_DIR = Path(r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-ISM\jorgeoutputs")

# Bandas a analizar
BANDS = ["VNIR-0", "VNIR-1", "VNIR-2", "VNIR-3"]

# Parámetros de cuantificación
BIT_DEPTH = 12
MAX_DN = (1 << BIT_DEPTH) - 1  # 2^bit - 1

def analyze_band(band: str):
    fname = f"ism_toa_{band}.nc"
    arr = readToa(OUT_DIR, fname)

    # Convertimos a array de enteros (DN deberían ser enteros tras la digitización)
    dn = np.asarray(arr)
    # Si por algún motivo está en float, lo llevamos a entero con seguridad
    if np.issubdtype(dn.dtype, np.floating):
        dn = np.rint(dn).astype(np.int64)
    else:
        dn = dn.astype(np.int64)

    # Enmascaramos NaN/Inf si existieran
    valid = np.isfinite(dn)
    total = dn.size
    n_valid = int(valid.sum())
    if n_valid == 0:
        raise RuntimeError(f"[{band}] No hay píxeles válidos en {fname}")

    # Saturaciones
    sat_hi = (dn == MAX_DN)
    sat_lo = (dn == 0)

    # Ojo: contamos solo sobre válidos
    sat_hi_cnt = int((sat_hi & valid).sum())
    sat_lo_cnt = int((sat_lo & valid).sum())
    sat_any_cnt = int(((sat_hi | sat_lo) & valid).sum())

    pct_hi = 100.0 * sat_hi_cnt / n_valid
    pct_lo = 100.0 * sat_lo_cnt / n_valid
    pct_any = 100.0 * sat_any_cnt / n_valid

    return {
        "band": band,
        "n_total": total,
        "n_valid": n_valid,
        "sat_hi_cnt": sat_hi_cnt,
        "sat_lo_cnt": sat_lo_cnt,
        "sat_any_cnt": sat_any_cnt,
        "pct_hi": pct_hi,
        "pct_lo": pct_lo,
        "pct_any": pct_any,
    }

def main():
    print(f"=== Saturation check (DN stage) — bit_depth={BIT_DEPTH}, MAX_DN={MAX_DN} ===\n")
    header = f"{'Band':<8}{'%SatHi':>10}{'%SatLo':>10}{'%SatAny':>12}{'Valid/Total':>16}"
    print(header)
    print("-" * len(header))

    for band in BANDS:
        try:
            r = analyze_band(band)
            print(f"{band:<8}{r['pct_hi']:10.4f}{r['pct_lo']:10.4f}{r['pct_any']:12.4f}"
                  f"{(str(r['n_valid']) + '/' + str(r['n_total'])):>16}")
        except Exception as e:
            print(f"{band:<8}ERROR: {e}")

    print("\nNotas:")
    print("  • %SatHi = porcentaje DN == 2^bit − 1 (saturación superior).")
    print("  • %SatLo = porcentaje DN == 0 (saturación inferior).")
    print("  • %SatAny = porcentaje DN == 0 o DN == 2^bit − 1.")
    print("  • Se ignoran valores no finitos si los hubiera.")

if __name__ == "__main__":
    main()