NEURON {
    SUFFIX IN
    :USEION na WRITE ina 
    :USEION k WRITE ik 
    NONSPECIFIC_CURRENT iNa, iK, ileak
    RANGE gNa, ENa, hNa, gk, Ek, nK, gl, El, minf, ina, ik, ileak
}

UNITS {
    (mA) = (milliamp)
    (uA) = (microamp)
    (mV) = (millivolt)
    (mS) = (millisiemens)
}

PARAMETER {
    gNa = 0.035 (S/cm2)
    ENa = 55 (mV)
    gk = 0.009 (S/cm2)
    Ek = -90 (mV)
    gl = 0.0001025 (S/cm2)
    El = -63.8 (mV)
}

STATE {
    hNa 
    nK
}

ASSIGNED {
    v (mV)
    iNa (mA/cm2)
    iK (mA/cm2)
    ileak (mA/cm2)
    alpha_n (1/ms)
    beta_n (1/ms)
    alpha_m (1/ms)
    beta_m (1/ms)
    alpha_h (1/ms)
    beta_h (1/ms)
    minf 
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    iNa = gNa * minf^3 * hNa * (v- ENa)
    iK = gk * nK^4 * (v - Ek)
    ileak = gl * (v - El)
}

DERIVATIVE states {
    evaluate_fct(v)
    hNa' = alpha_h * (1 - hNa) - beta_h * hNa
    nK' = alpha_n * (1 - nK) - beta_n * nK
}

UNITSOFF

INITIAL {
    hNa = hNa
    nK = nK
}

PROCEDURE evaluate_fct(v(mV)) {
    alpha_m = 0.5 * (v + 35) / (1 - exp(-(v + 35) / 10))
    beta_m = 20 * exp(-(v + 60) / 18)
    alpha_h = 0.35 * exp(-(v + 58) / 20)
    beta_h = 5 / (1 + exp(-(v + 28) / 10))
    alpha_n = 0.05 * (v + 34) / (1 - exp(-(v + 34) / 10))
    beta_n = 0.625 * exp(-(v + 44) / 80)
    minf = alpha_m / (alpha_m + beta_m)
}

UNITSON

:VERBATIM
:extern double scop_random();
:extern void set_seed(double);
:ENDVERBATIM
