NEURON {
    SUFFIX RE
    :USEION na WRITE ina 
    :USEION k WRITE ik 
    :USEION ca WRITE ica
    NONSPECIFIC_CURRENT iNa, iK, iCa, iLeak, iKLeak 
    RANGE gNa, ENa, gK, EK, gLeak, ELeak, gKLeak, EKLeak, gT, ET
    RANGE iNa, iK, iLeak, iKLeak, iT
    RANGE mNa, hNa, nK, mT, hT
    RANGE Minf, Hinf, tauM, tauH
    RANGE phiM, phiH, vShiftNa, vShiftK, vShiftT
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}

PARAMETER {
    gNa = 0.200 (S/cm2)
    ENa = 50 (mV)
    vShiftNa = 55 (mV)
    gK = 0.020 (S/cm2)
    EK = -100 (mV)
    vShiftK = 55 (mV)
    gLeak = 0.00005 (S/cm2)
    ELeak = -90 (mV)
    gKLeak = 0 (S/cm2)
    EKLeak = -95 (mV)
    gT = 0.003 (S/cm2)
    ET = 120 (mV)
    vShiftT = 4 (mV)
    phiM = 6.81
    phiH = 3.73
}

STATE {
    mNa 
    hNa 
    nK 
    mT 
    hT
}

ASSIGNED {
    v (mV)
    iNa (mA/cm2)
    iK (mA/cm2)
    iCa (mA/cm2)
    iLeak (mA/cm2)
    iKLeak (mA/cm2)
    iPoisson_RE (mA/cm2)
    Minf 
    Hinf
    tauM (ms) 
    tauH (ms)
    alphaM (1/ms)
    betaM (1/ms)
    alphaH (1/ms)
    betaH (1/ms)
    alphaN (1/ms)
    betaN (1/ms)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    iNa = gNa * mNa^3 * hNa * (v - ENa)
    iK = gK * nK^4 * (v - EK)
    iLeak = gLeak * (v - ELeak)
    iKLeak = gKLeak * (v - EKLeak)
    iCa = gT * mT^2 * hT * (v - ET)
}

DERIVATIVE states {
    evaluate_fct(v)
    mNa' = alphaM * (1 - mNa) - betaM * mNa
    hNa' = alphaH * (1 - hNa) - betaH * hNa
    nK' = alphaN * (1 - nK) - betaN * nK
    mT' = (Minf - mT) / tauM
    hT' = (Hinf - hT) / tauH
}

UNITSOFF

INITIAL {
    mNa = mNa
    hNa = hNa
    nK = nK
    mT = mT
    hT = hT
}

PROCEDURE evaluate_fct(v(mV)) {
    alphaM = 0.32 * (13 - (v + vShiftNa)) / (exp((13 - (v + vShiftNa)) / 4) - 1)
    betaM = 0.28 * ((v + vShiftNa) - 40) / (exp(((v + vShiftNa) - 40) / 5) - 1)
    alphaH = 0.128 * exp((17 - (v + vShiftNa)) / 18)
    betaH = 4 / (exp((40 - (v + vShiftNa)) / 5) + 1)
    alphaN = 0.032 * (15 - (v + vShiftK)) / (exp((15 - (v + vShiftK)) / 5) - 1)
    betaN = 0.5 * exp((10 - (v + vShiftK)) / 40)
    Minf = 1 / (1 + exp((-(v + vShiftT + 50)) / 7.4))
    tauM = ((3.0 + 1.0/(exp((v + vShiftT + 25) / 10) + exp(-(v + vShiftT + 100) / 15))) / phiM)
    Hinf = 1 / (1 + exp((v + vShiftT + 78) / 5))
    tauH = ((85.0 + 1.0/(exp((v + vShiftT + 46) / 4)  + exp(-(v + vShiftT + 405) / 50))) / phiH)
}

UNITSON

:VERBATIM
:#include <math.h>
:extern double scop_random();
:extern void set_seed(double);
:ENDVERBATIM

