NEURON {
    POINT_PROCESS AMPA
    :POINTER v_pre                 : Presynaptic voltage
    RANGE gAMPA, EAMPA, sAMPA, alpha, tauAMPA, v_pre
    RANGE normalizing_factor, area_post
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (S) = (siemens)
}

PARAMETER {
    gAMPA = 0.4e-3 (S/cm2)        : Conductance (0.4 mS/cm^2 -> S/cm^2)
    EAMPA = 1 (mV)                : Reversal potential
    alpha = 5 (/ms)               : Activation rate constant
    tauAMPA = 2 (ms)        : Decay time constant
    normalizing_factor = 1 (1)        : Normalization factor
    area_post = 0.000143 (cm2)    : Postsynaptic neuron area
}

ASSIGNED {
    v (mV)                        : Postsynaptic voltage
    v_pre (mV)                    : Presynaptic voltage via pointer
    i (nA)                        : Synaptic current

}

STATE {
    sAMPA (1)                     : Synaptic activation (dimensionless)
}

INITIAL {
    sAMPA = sAMPA  : Random initial value between 0 and 0.1 (set from Python)
    i = 0
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    i = (gAMPA / normalizing_factor) * sAMPA * (v - EAMPA) * 1e6 * area_post
}

DERIVATIVE states {
    sAMPA' = alpha * (1 + tanh(v_pre / 4)) * (1 - sAMPA) - sAMPA / tauAMPA
}