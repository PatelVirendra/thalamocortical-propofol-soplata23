NEURON {
    POINT_PROCESS AMPAD
    :POINTER v_pre
    RANGE gAMPA, EAMPA, sAMPA, res_AMPA, v_pre, spike_activity
    RANGE alpha, tau, tauRes, deprFactor, normalizing_factor
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
    gAMPA = 0.005e-3 (S/cm2)
    EAMPA = 0 (mV)
    
    : Kinetic Parameters
    alpha = 3.48 (/ms)
    tau = 2 (ms)
    tauRes = 400 (ms)
    
    : Depression
    deprFactor = 0.9     : Instantaneous depression jump
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
    sAMPA (1)
    res_AMPA (1)
}

INITIAL {
    sAMPA = sAMPA
    res_AMPA = res_AMPA
    normalizing_factor = clip((2 * radius + (1 - remove_recurrent_bool)) / (Npost / Npre), 0, Npre)
    i = 0
    spike_activity = 0
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    i = (gAMPA / normalizing_factor) * (res_AMPA * sAMPA) * (v - EAMPA) * 1e6 * area_post
}

DERIVATIVE states {
    : 1. Continuous Activation (Driven by v_pre pointer)
    sAMPA' = alpha / (1 + exp(-(v_pre - 20) / 2)) - sAMPA / tau

    : 2. Continuous Recovery
    :res_AMPA' = (1 - res_AMPA) / tauRes
    res_AMPA' = (1 - res_AMPA) / tauRes + spike_activity * (-(1 - res_AMPA) / tauRes + (deprFactor * res_AMPA - res_AMPA) / dt)
}

: --- NET_RECEIVE BLOCK ---
: This block is triggered ONLY when a NetCon sends a spike event.
: We ignore the weight 'w' because our conductance is defined by gAMPA param.
NET_RECEIVE(w) {
    : Apply instantaneous depression
    :res_AMPA = res_AMPA * deprFactor
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