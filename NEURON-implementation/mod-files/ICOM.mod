NEURON {
    POINT_PROCESS ICOM
    :POINTER v_pre
    RANGE gCOM, area_post, v_pre
    NONSPECIFIC_CURRENT iCOM
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (nA) = (nanoamp)
    (S) = (siemens)
}

PARAMETER {
    gCOM = 0 (S/cm2)
    area_post = 0.00035 (cm2) : Postsynaptic neuron area
}

ASSIGNED {
    v (mV)
    v_pre (mV)
    iCOM (nA)
}

BREAKPOINT {
    iCOM = gCOM * (v - v_pre) * 1e6 * area_post
}