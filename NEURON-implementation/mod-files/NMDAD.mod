NEURON {
    POINT_PROCESS NMDAD
    :POINTER v_pre
    RANGE gNMDA, ENMDA, sNMDA, xNMDA, res_NMDA, v_pre, spike_activity
    RANGE alphaS, tauS, alphaX, tauX, deprFactor, tauRes, normalizing_factor
    RANGE Npre, Npost, radius, remove_recurrent_bool, area_post
    RANGE release_time   : Delay for "down" event
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (S) = (siemens)
}

PARAMETER {
    gNMDA = 0.00257e-3 (S/cm2)
    ENMDA = 0 (mV)
    alphaS = 0.5 (/ms)
    tauS = 100 (ms)
    alphaX = 3.48 (/ms)
    tauX = 2 (ms)
    tauRes = 400 (ms)
    
    deprFactor = 0.9
    release_time = 0.02 (ms)  : Delay between "up" and "down" events (e.g., 2 * dt)

    Npre = 1
    Npost = 1
    radius = 10
    remove_recurrent_bool = 1
    area_post = 0.00035 (cm2)
}

ASSIGNED {
    v (mV)
    v_pre (mV)
    i (nA)
    normalizing_factor (1)
    spike_activity (1)        : Tracks spike state (0 or 1)
}

STATE {
    sNMDA (1)
    xNMDA (1)
    res_NMDA (1)
}

INITIAL {
    sNMDA = sNMDA
    xNMDA = xNMDA
    res_NMDA = res_NMDA
    normalizing_factor = clip((2 * radius + (1 - remove_recurrent_bool)) / (Npost / Npre), 0, Npre)
    i = 0
    spike_activity = 0
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    i = (gNMDA / normalizing_factor) * (res_NMDA * sNMDA) * (v - ENMDA) * 1e6 * area_post
}

DERIVATIVE states {
    : Continuous Activation
    sNMDA' = alphaS * xNMDA * (1 - sNMDA) - sNMDA / tauS
    xNMDA' = alphaX / (1 + exp(-(v_pre - 20) / 2)) - xNMDA / tauX

    : Continuous Recovery
    :res_NMDA' = (1 - res_NMDA) / tauRes
    res_NMDA' = (1 - res_NMDA) / tauRes + spike_activity * (-(1 - res_NMDA) / tauRes + (deprFactor * res_NMDA - res_NMDA) / dt)
}

NET_RECEIVE(w) {
    : Apply instantaneous depression on spike arrival
    :res_NMDA = res_NMDA * deprFactor
    if (flag == 0) {  : External event (presynaptic spike)
        spike_activity = 1  : "up" event
        net_send(release_time, 1)            : Schedule "down" event
    } else if (flag == 1) {  : Self-event for "down"
        spike_activity = 0  : "down" event
    }
}

FUNCTION clip(x, lo, hi) {
    if (x < lo) { clip = lo }
    else if (x > hi) { clip = hi }
    else { clip = x }
}