NEURON {
    POINT_PROCESS GABAB
    :POINTER v_pre                 : Presynaptic voltage
    RANGE gGABAB, EGABAB, rGABAB, sGABAB, v_pre
    RANGE K1, K2, K3, K4
    RANGE normalizing_factor, area_post
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (S) = (siemens)
}

PARAMETER {
    gGABAB = 0.001e-3 (S/cm2)  : Conductance (0.001 mS/cm^2 -> S/cm^2)
    EGABAB = -95 (mV)          : Reversal potential
    K1 = 0.5 (/ms)             : Rate constant for rGABAB activation
    K2 = 0.0012 (/ms)          : Rate constant for rGABAB decay
    K3 = 0.18 (/ms)            : Rate constant for sGABAB activation
    K4 = 0.034 (/ms)           : Rate constant for sGABAB decay
    normalizing_factor = 1 (1) : Normalization factor
    area_post = 0.00029 (cm2)  : Postsynaptic neuron area
}

ASSIGNED {
    v (mV)                     : Postsynaptic voltage
    v_pre (mV)                 : Presynaptic voltage via pointer
    i (nA)                     : Synaptic current
}

STATE {
    rGABAB (1)                 : Receptor binding state (dimensionless)
    sGABAB (1)                 : Synaptic activation state (dimensionless)
}

INITIAL {
    rGABAB = rGABAB  : Random initial value between 0 and 0.1 (set from Python)
    sGABAB = sGABAB  : Random initial value between 0 and 0.1 (set from Python)
    i = 0
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    i = (gGABAB / normalizing_factor) * (sGABAB^4 / (sGABAB^4 + 100)) * (v - EGABAB) * 1e6 * area_post
}

DERIVATIVE states {
    rGABAB' = K1 * (2 * (1 + tanh(v_pre / 4))) * (1 - rGABAB) - K2 * rGABAB
    sGABAB' = K3 * rGABAB - K4 * sGABAB
}