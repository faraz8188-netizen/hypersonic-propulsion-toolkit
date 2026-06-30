# Hypersonic Propulsion Toolkit

> **Open-source Python toolkit for ramjet and scramjet thermodynamic performance analysis (Mach 2-15)**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21054986.svg)](https://doi.org/10.5281/zenodo.21054986)
[![Validated vs EUCASS 2017-031](https://img.shields.io/badge/validated-EUCASS%202017--031-orange)]()
[![Chemistry vs Marzouk 2023](https://img.shields.io/badge/chemistry-Marzouk%202023-blueviolet)]()

---

## v1.1 Update - Real Equilibrium Combustion Chemistry

The original v1.0 combustor model used a fixed H2 lower heating value (LHV)
assumption. v1.1 adds genuine NASA-polynomial equilibrium combustion
thermodynamics, replacing that assumption with real chemistry:

- **nasa7_thermo.py**: NASA 7-coefficient polynomial thermodynamic database
  (H, H2, O, O2, OH, H2O, N2), exact coefficients from GRI-Mech 3.0 /
  NASA TM-4513. Validated: computed H2O formation enthalpy at 298.15K
  matches the NIST WebBook reference value exactly (-241.826 kJ/mol).
- **equilibrium_combustion.py**: H2/air adiabatic flame temperature solver
  with dissociation effects, validated to **0.06% error** against published
  CEARUN equilibrium-chemistry data (Marzouk 2023, DOI:10.48084/etasr.6132).
- **scramjet_combustor_v2.py**: integrates the above using the ACTUAL
  combustor static temperature (not freestream stagnation temperature) as
  the reactant initial condition.

**Key finding from this integration:** at realistic stoichiometric
combustion, the validated chemistry range (NASA polynomials valid to 3500K)
is exceeded starting at Mach 9, not Mach 15 as a naive check would suggest -
found by checking both inlet AND output flame temperature, after an earlier
version of the code missed this by only checking inlet temperature. At
realistic lean operation (phi=0.4), the validated range extends to M=10-11.
This finding, including the bug that revealed it, is documented in the code.

---

## Companion Thesis

This toolkit is the computational companion to:

**A Technical Study of Ramjet-Scramjet Propulsion Systems and Their Applications in Hypersonic Flight**
Faraz Shaikh | B.Tech Aerospace Engineering | Sandip University, Nashik, India | 2025

All simulation results in Chapter 11 of the thesis are implemented and validated by this toolkit. The toolkit is fully documented in Chapter 16.

---

## Scientific Basis - No Fabricated Values

Every constant is sourced from peer-reviewed literature:

| Parameter | Value | Source |
|---|---|---|
| H2 Lower Heating Value | 120.0 MJ/kg | Heiser & Pratt (1994) |
| Kerosene JP-7 LHV | 43.2 MJ/kg | Standard aviation reference |
| Ramjet inlet recovery | eta=0.92 | NATO AVT-185 / DTIC ADA596056 |
| Scramjet combustion efficiency | 0.88-0.92 | EUCASS 2017-031 |
| Scramjet f at M=10 | 0.039 | EUCASS 2017-031 (validated) |
| Scramjet Isp at M=10 | ~2607 s | EUCASS 2017-031 |
| Scramjet Fs at M=10 | ~1000 N.s/kg | EUCASS 2017-031 |
| Inlet recovery M=5-15 | 0.38-0.80 | Heiser & Pratt (1994) Fig 4.3 |
| Altitude trajectory | Constant dynamic-q | Heiser & Pratt (1994) |
| H2O formation enthalpy | -241.826 kJ/mol | NIST WebBook (exact match) |
| Stoichiometric air-H2 AFT | 2520.33 K (complete), 2378.62 K (equilibrium) | Marzouk (2023), DOI:10.48084/etasr.6132 |

---

## Repository Structure

```
hypersonic-propulsion-toolkit/
|- atmosphere.py                             # US Standard Atmosphere 1976 (NASA TM-X-74335)
|- ramjet_performance.py                     # Brayton cycle model, Mach 2-5
|- scramjet_performance.py                   # Enthalpy nozzle model, Mach 5-15
|- comparative_analysis.py                   # Full Mach 2-15 comparison
|- nasa7_thermo.py                           # v1.1: NASA polynomial thermo database
|- equilibrium_combustion.py                 # v1.1: equilibrium AFT solver
|- scramjet_combustor_v2.py                  # v1.1: real-chemistry combustor integration
|- run_all.py                                # Main runner
|- requirements.txt                          # Python dependencies
|- data/
|  |- stagnation_temperature_validated.csv   # Thesis Ch.11 Table vs toolkit: 0.0% error
|  |- ramjet_performance_results.csv         # Ramjet H2 and Kerosene results
|  |- scramjet_performance_results.csv       # Scramjet H2 results (EUCASS validated)
|- figures/                                  # Output plots (generated on run)
```

---

## Research Paper and DOI

This toolkit and accompanying research are archived on Zenodo, with both
v1.0.0 and v1.1.0 published as separate citable versions:

**Latest (v1.1.0): [10.5281/zenodo.21054986](https://doi.org/10.5281/zenodo.21054986)**
**v1.0.0: [10.5281/zenodo.21053028](https://doi.org/10.5281/zenodo.21053028)**

---

## Validated Results

### Stagnation Temperature: Thesis Ch.11 vs Toolkit

| Mach | Thesis (K) | Toolkit (K) | Error |
|------|-----------|-------------|-------|
| 2    | 396       | 396         | 0.0%  |
| 6    | 1,804     | 1,804       | 0.0%  |
| 10   | 4,620     | 4,620       | 0.0%  |
| 15   | 10,120    | 10,120      | 0.0%  |

### Scramjet vs EUCASS 2017-031 (M=10)

| Parameter | EUCASS Ref | Toolkit | Difference |
|---|---|---|---|
| Isp | 2,607 s | 2,921 s | +12% (altitude schedule) |
| Fs | ~1,000 N.s/kg | 1,117 N.s/kg | +12% |
| f | 0.039 | 0.039 | 0% |

### Equilibrium Combustion vs Marzouk (2023) - v1.1

| Parameter | Marzouk 2023 | Toolkit | Error |
|---|---|---|---|
| Complete combustion AFT | 2,520.33 K | 2,518.83 K | 0.06% |
| Equilibrium AFT (with dissociation) | 2,378.62 K | 2,377.12 K | 0.06% |

---

## Installation

```bash
git clone https://github.com/faraz8188-netizen/hypersonic-propulsion-toolkit.git
cd hypersonic-propulsion-toolkit
pip install -r requirements.txt
python run_all.py
```

Run v1.1 chemistry modules individually:
```bash
python nasa7_thermo.py            # Thermo database self-test
python equilibrium_combustion.py  # Equilibrium solver + Marzouk validation
python scramjet_combustor_v2.py   # Full combustor survey across Mach 5-15
```

---

## Key References

- W. H. Heiser and D. T. Pratt, Hypersonic Airbreathing Propulsion. AIAA, 1994.
- EUCASS 2017-031: Thermodynamic performance of scramjet at wide Mach numbers.
- NATO AVT-185 / DTIC ADA596056: Ramjet Intakes.
- NASA TM-X-74335: US Standard Atmosphere 1976.
- NASA TM-4513: McBride, Gordon and Reno, Coefficients for Calculating Thermodynamic and Transport Properties of Individual Species, 1993.
- Marzouk, O.A., Adiabatic Flame Temperatures for Oxy-Methane, Oxy-Hydrogen, Air-Methane, and Air-Hydrogen Stoichiometric Combustion, Eng. Technol. Appl. Sci. Res., 13(4), 2023. DOI:10.48084/etasr.6132
- G. P. Sutton and O. Biblarz, Rocket Propulsion Elements, 9th ed. Wiley, 2017.
- J. D. Anderson Jr., Hypersonic and High Temperature Gas Dynamics, 2nd ed. AIAA, 2006.

---

## Author

**Faraz Shaikh**
B.Tech Aerospace Engineering | Sandip University, Nashik, India (2025) | CGPA: 7.84/10

Research interests: Hypersonic propulsion, spacecraft propulsion, detonation-based engines.

---

## Citation

```bibtex
@software{shaikh2025hypersonic,
  author = {Shaikh, Faraz},
  title  = {Hypersonic Propulsion Toolkit},
  year   = {2026},
  version = {v1.1.0},
  doi    = {10.5281/zenodo.21054986},
  url    = {https://doi.org/10.5281/zenodo.21054986}
}
```

## License

MIT License
