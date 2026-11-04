"""Small example datasets used by the tests and the gallery.

Ternary points are ``(solute, carrier, solvent)`` in wt%. Row *i* of a
raffinate list pairs with row *i* of the extract list as one tie line.

The tie-line tables were transcribed from coursework (unit operations
tutorials); verify against the original handbook source before citing them.
"""

#: Acetone (A) – water (B) – 1,1,2-trichloroethane (C)
SYSTEM_1 = {
    "labels": ("Acetone", "Water", "Trichloroethane"),
    "raffinate": [(5.96, 93.52, 0.52), (17.04, 82.23, 0.73), (26.92, 72.06, 1.02),
                  (30.88, 67.95, 1.17), (35.73, 62.67, 1.60), (40.90, 57.00, 2.10),
                  (46.05, 50.20, 3.75), (51.78, 41.70, 6.52)],
    "extract": [(8.75, 0.32, 90.93), (25.14, 1.10, 73.76), (38.52, 2.27, 59.21),
                (42.97, 3.11, 53.92), (48.21, 4.26, 47.53), (53.95, 6.05, 40.00),
                (57.40, 8.90, 33.70), (60.34, 13.40, 26.26)],
}

#: Solute (A) – water (B) – isopropyl ether (C); distribution coefficient K < 1
SYSTEM_2 = {
    "labels": ("Solute", "Water", "Isopropyl ether"),
    "raffinate": [(0.69, 98.1, 1.2), (2.89, 95.5, 1.6), (6.42, 91.7, 1.9),
                  (13.30, 84.4, 2.3), (25.50, 71.7, 3.4), (36.70, 58.9, 4.4),
                  (44.30, 45.1, 10.6), (46.40, 37.1, 16.5)],
    "extract": [(0.18, 0.5, 99.3), (0.79, 0.8, 98.4), (1.93, 1.0, 97.1),
                (4.82, 1.9, 93.3), (11.40, 3.9, 84.7), (21.60, 6.9, 71.5),
                (31.10, 10.8, 58.1), (36.20, 15.1, 48.7)],
}

#: Countercurrent extraction, solute-free basis (X, Y in kg/kg)
COUNTERCURRENT = {
    "eq_x": [0, 1.011e-3, 2.460e-3, 5.02e-3, 7.51e-3, 9.98e-3, 12.4e-3],
    "eq_y": [0, 1.907e-3, 3.961e-3, 5.96e-3, 7.86e-3, 9.31e-3, 9.7e-3],
    "raffinate_end": (1.058e-3, 0.0),
    "feed_end": (10.101e-3, 7.461e-3),
}

#: Batch drying: mean free moisture X (kg/kg) and rate R (kg/(m²·h))
DRYING = {
    "X": [0.2548, 0.2368, 0.2121, 0.1784, 0.1392, 0.0976, 0.0639, 0.0344, 0.0116, 0.0031],
    "R": [0.79, 1.03, 0.98, 0.97, 1.01, 0.73, 0.61, 0.35, 0.11, 0.04],
    "X_c": 0.14,
    "R_c": 1.0,
}
