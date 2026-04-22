NEURON {
    POINT_PROCESS GABAA
    :POINTER v_pre                 : Presynaptic voltage
    RANGE gGABAA, EGABAA, sGABAA, v_pre
    RANGE alpha, tauGABAA
    RANGE normalizing_factor, area_post
    RANGE propoCondMult, propoTauMult
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (S) = (siemens)
}

PARAMETER {
    gGABAA = 0.1e-3 (S/cm2)    : Conductance (0.1 mS/cm^2 -> S/cm^2)
    EGABAA = -80 (mV)          : Reversal potential
    alpha = 2 (/ms)            : Activation rate constant
    tauGABAA = 5 (ms)          : Decay time constant (base value)
    propoCondMult = 1          : Conductance multiplier
    propoTauMult = 1           : Tau multiplier
    normalizing_factor = 1 (1) : Normalization factor
    area_post = 0.000143 (cm2) : Postsynaptic neuron area
}

ASSIGNED {
    v (mV)                     : Postsynaptic voltage
    v_pre (mV)                 : Presynaptic voltage via pointer
    i (nA)                     : Synaptic current
}

STATE {
    sGABAA (1)                 : Synaptic activation (dimensionless)
}

INITIAL {
    sGABAA = sGABAA            : Random initial value between 0 and 0.1 (set from Python)
    i = 0
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    i = propoCondMult * (gGABAA / normalizing_factor) * sGABAA * (v - EGABAA) * 1e6 * area_post
}

DERIVATIVE states {
    sGABAA' = alpha * (1 + tanh(v_pre / 4)) * (1 - sGABAA) - sGABAA / (tauGABAA * propoTauMult)
}
