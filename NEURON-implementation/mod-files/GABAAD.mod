NEURON {
    POINT_PROCESS GABAAD
    :POINTER v_pre
    RANGE gGABAA, EGABAA, sGABAA, res_GABAA, v_pre, spike_activity
    RANGE alpha, tauGABAA, tauRes, deprFactor, normalizing_factor
    RANGE Npre, Npost, radius, remove_recurrent_bool, area_post
    RANGE propoCondMult, propoTauMult
    RANGE release_time   : Delay for "down" event
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (S) = (siemens)
}

PARAMETER {
    gGABAA = 0.1e-3 (S/cm2)
    EGABAA = -70 (mV)
    alpha = 1 (/ms)
    tauGABAA = 5 (ms)
    tauRes = 400 (ms)
    
    deprFactor = 0.9
    release_time = 0.02 (ms)  : Delay between "up" and "down" events (e.g., 2 * dt)
    
    Npre = 1
    Npost = 1
    radius = 10
    remove_recurrent_bool = 1
    propoCondMult = 1
    propoTauMult = 1
    area_post = 0.00015 (cm2)
}

ASSIGNED {
    v (mV)
    v_pre (mV)
    i (nA)
    normalizing_factor (1)
    spike_activity (1)        : Tracks spike state (0 or 1)
}

STATE {
    sGABAA (1)
    res_GABAA (1)
}

INITIAL {
    sGABAA = sGABAA
    res_GABAA = res_GABAA
    normalizing_factor = clip((2 * radius + (1 - remove_recurrent_bool)) / (Npost / Npre), 0, Npre)
    i = 0
    spike_activity = 0
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    i = propoCondMult * (gGABAA / normalizing_factor) * (res_GABAA * sGABAA) * (v - EGABAA) * 1e6 * area_post
}

DERIVATIVE states {
    : Continuous Activation
    sGABAA' = alpha / (1 + exp(-(v_pre - 20) / 2)) - sGABAA / (tauGABAA * propoTauMult)

    : Continuous Recovery
    :res_GABAA' = (1 - res_GABAA) / tauRes
    res_GABAA' = (1 - res_GABAA) / tauRes + spike_activity * (-(1 - res_GABAA) / tauRes + (deprFactor * res_GABAA - res_GABAA) / dt)
}

NET_RECEIVE(w) {
    : Apply instantaneous depression on spike arrival
    :res_GABAA = res_GABAA * deprFactor
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