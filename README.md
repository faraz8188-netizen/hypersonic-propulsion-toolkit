# Hypersonic Propulsion Toolkit

> **Open-source Python toolkit for ramjet and scramjet thermodynamic performance analysis (Mach 2-15)**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21053028.svg)](https://doi.org/10.5281/zenodo.21053028)
[![Validated vs EUCASS 2017-031](https://img.shields.io/badge/validated-EUCASS%202017--031-orange)]()

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

---

## Repository Structure

```
hypersonic-propulsion-toolkit/
|- atmosphere.py                             # US Standard Atmosphere 1976 (NASA TM-X-74335)
|- ramjet_performance.py                     # Brayton cycle model, Mach 2-5
|- scramjet_performance.py                   # Enthalpy nozzle model, Mach 5-15
|- comparative_analysis.py                   # Full Mach 2-15 comparison
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

This toolkit and accompanying research are archived on Zenodo:

**DOI: [10.5281/zenodo.21053028](https://doi.org/10.5281/zenodo.21053028)**

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

The 12% difference is a known effect of the constant dynamic-pressure altitude schedule vs fixed conditions in EUCASS. See Chapter 16 of the thesis.

---

## Installation

```bash
git clone https://github.com/faraz8188-netizen/hypersonic-propulsion-toolkit.git
cd hypersonic-propulsion-toolkit
pip install -r requirements.txt
python run_all.py
```

---

## Key References

- W. H. Heiser and D. T. Pratt, Hypersonic Airbreathing Propulsion. AIAA, 1994.
- EUCASS 2017-031: Thermodynamic performance of scramjet at wide Mach numbers.
- NATO AVT-185 / DTIC ADA596056: Ramjet Intakes.
- NASA TM-X-74335: US Standard Atmosphere 1976.
- G. P. Sutton and O. Biblarz, Rocket Propulsion Elements, 9th ed. Wiley, 2017.
- J. D. Anderson Jr., Hypersonic and High Temperature Gas Dynamics, 2nd ed. AIAA, 2006.
- J. D. Anderson Jr., Introduction to Flight, 8th ed. McGraw-Hill, 2016.
- M. K. Smart, Scramjets. Aeronautical Journal, vol. 111, pp. 605-619, 2007.
- F. S. Billig, Research on Supersonic Combustion. J. Propulsion Power, 1993.
- G. P. Anderson et al., Scramjet Performance. in Scramjet Propulsion, AIAA, 2000.

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
  year   = {2025},
  doi    = {10.5281/zenodo.21053028},
  url    = {https://doi.org/10.5281/zenodo.21053028}
}
```

## License

MIT License
