"""
Hypersonic Propulsion Toolkit - Main Runner
Runs all analysis modules and generates publication-quality figures.

Usage:
    python run_all.py

Outputs saved to figures/
"""
import os
os.makedirs('figures', exist_ok=True)

print('=' * 60)
print('HYPERSONIC PROPULSION TOOLKIT')
print('Ramjet & Scramjet Analysis | Mach 2-15')
print('Author: Faraz Shaikh, Sandip University (2025)')
print('=' * 60)

print('\n[1/3] Ramjet performance (Mach 2-5, Brayton cycle)...')
import ramjet_performance
ramjet_performance.plot_ramjet_performance()

print('\n[2/3] Scramjet performance (Mach 5-15, H2 fuel)...')
import scramjet_performance
scramjet_performance.plot_scramjet_performance()

print('\n[3/3] Comparative analysis (Mach 2-15)...')
import comparative_analysis
comparative_analysis.plot_combined()

print('\nAll figures saved to figures/')
