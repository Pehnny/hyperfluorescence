from math import pi

BOLTZMANN : float = 8.617333262 * 10**(-5)                                              # [eV/K]
VACUUM_PERMITTIVITY : float = 55263.49406                                               # [e²/eV.nm]
RELATIVE_PERMITTIVITY : float = 3.
ELECTROSTATIC : float = 1. / (4. * pi * VACUUM_PERMITTIVITY * RELATIVE_PERMITTIVITY)    # [eV.nm]