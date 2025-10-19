#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
EODP-TS-ISM-0002: Factores de conversión y trazabilidad numérica (print only)

- Irradiancia (mW/m^2) -> Fotones [ph]
- Fotones -> Electrones [e-]
- Electrones -> Voltios [V]
- Voltios -> Digital Numbers [DN]

Imprime:
  1) Factores por banda (I->ph y el combinado I->DN)
  2) Factores comunes (ph->e, e->V, V->DN)
  3) Trazabilidad con I=1 mW/m^2: ph, e-, V y DN resultantes por banda
"""

from dataclasses import dataclass

@dataclass(frozen=True)
class Config:
    # Geometría / tiempo
    pix_size_m: float = 30e-6        # [m]
    t_int_s: float = 0.00672         # [s]
    # Bandas: longitudes de onda centrales [m]
    lambda_m: tuple = (0.49e-6, 0.665e-6, 0.865e-6, 0.945e-6)
    band_names: tuple = ("VNIR-0", "VNIR-1", "VNIR-2", "VNIR-3")
    # Detección
    QE: float = 0.8                  # [e-/ph]
    # Cadena de vídeo
    OCF_V_per_e: float = 5.4e-6      # [V/e-]
    ADC_gain: float = 0.56           # [-]
    bit_depth: int = 12              # [-]
    Vmin: float = 0.0                # [V]
    Vmax: float = 0.86               # [V]
    # Constantes físicas
    h_planck: float = 6.62607015e-34 # [J·s]
    c_light: float = 299792458.0     # [m/s]

def main():
    cfg = Config()

    # Derivados (comunes)
    A_pix = cfg.pix_size_m ** 2                      # [m^2]
    K_ph_to_e = cfg.QE                               # [e-/ph]
    K_e_to_V = cfg.OCF_V_per_e * cfg.ADC_gain        # [V/e-]
    K_V_to_DN = (2**cfg.bit_depth - 1) / (cfg.Vmax - cfg.Vmin)  # [DN/V]

    print("\n=== PARÁMETROS ===")
    print(f"A_pix = {A_pix:.6e} m^2,  t_int = {cfg.t_int_s} s")
    print(f"QE = {cfg.QE} e-/ph")
    print(f"OCF = {cfg.OCF_V_per_e} V/e-,  ADC_gain = {cfg.ADC_gain}  -> K_e->V = {K_e_to_V:.6e} V/e-")
    print(f"bit_depth = {cfg.bit_depth}, Vmin = {cfg.Vmin} V, Vmax = {cfg.Vmax} V  -> K_V->DN = {K_V_to_DN:.10f} DN/V")
    print(f"h = {cfg.h_planck:.8e} J·s,  c = {cfg.c_light:.3f} m/s\n")

    # Cabecera compacta
    head = (
        f"{'Banda':<8}{'λ [µm]':>10}"
        f"{'K_I->ph [ph/(mW·m^-2)]':>28}"
        f"{'K_ph->e [e-/ph]':>16}"
        f"{'K_e->V [V/e-]':>16}"
        f"{'K_V->DN [DN/V]':>18}"
        f"{'K_I->DN [DN/(mW·m^-2)]':>30}"
    )
    print("=== FACTORES (entrada en mW·m^-2) ===")
    print(head)
    print("-" * len(head))

    # Por banda: sólo I->ph depende de λ; los demás factores son comunes
    results = []  # guardamos para la trazabilidad
    for name, lam in zip(cfg.band_names, cfg.lambda_m):
        # Irradiancia -> Fotones (mW/m^2 -> ph)
        # Nph = (I_mWm2 * 1e-3 [W/m^2]) * A_pix * t_int * lam / (h*c)
        K_I_to_ph = (1e-3) * A_pix * cfg.t_int_s * lam / (cfg.h_planck * cfg.c_light)

        # Factores combinados
        K_I_to_e  = K_I_to_ph * K_ph_to_e            # [e-/(mW·m^-2)]
        K_I_to_V  = K_I_to_e  * K_e_to_V             # [V/(mW·m^-2)]
        K_I_to_DN = K_I_to_V  * K_V_to_DN            # [DN/(mW·m^-2)]

        print(f"{name:<8}{lam*1e6:10.3f}{K_I_to_ph:28.6e}{K_ph_to_e:16.6f}"
              f"{K_e_to_V:16.6e}{K_V_to_DN:18.6f}{K_I_to_DN:30.6f}")

        results.append((name, lam, K_I_to_ph, K_I_to_e, K_I_to_V, K_I_to_DN))

    print("\nUnidades del TOA por etapa:")
    print("  1) Irradiancia: [mW·m⁻²]")
    print("  2) Fotones: [ph]")
    print("  3) Electrones: [e⁻]")
    print("  4) Voltaje: [V]")
    print("  5) Digital Number: [DN]\n")

if __name__ == "__main__":
    main()
