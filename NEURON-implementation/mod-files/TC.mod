NEURON {
    SUFFIX TC
    :USEION na WRITE ina
    :USEION k WRITE ik
    :USEION ca READ cai WRITE ica
    NONSPECIFIC_CURRENT iNa, iK, iCa, iH, ileak, ikleak, iPoisson_TC
    RANGE gNa, ENa, vShiftNa, gK, EK, vShiftK, gLeak, ELeak, gKLeak, EKLeak, gT, vShiftT
    RANGE phiH, gH, gInc, EH, Cac, nca, k2, k4, PC, nexp
    RANGE Ca_inf, tauR, CaIC, CaNoiseIC
    RANGE Minf, tauH, Hinf, tauS
    RANGE mNa, hNa, nK, hT, Open, Pone, OpenLocked, cai
}

UNITS {
    (mA) = (milliamp)
    (uA) = (microamp)
    (mV) = (millivolt)
    (molar) = (1/liter)
    (mM) = (millimolar)
}

PARAMETER {
    gNa = 0.090 (S/cm2)
    ENa = 50 (mV)
    vShiftNa = -40 (mV)
    gK = 0.010 (S/cm2)
    EK = -100 (mV)
    vShiftK = -25 (mV)
    gLeak = 0.00001 (S/cm2)
    ELeak = -70 (mV)
    gKLeak = 0.0000172 (S/cm2)
    EKLeak = -100 (mV)
    gT = 0.002 (S/cm2)
    vShiftT = 2 (mV)
    phiH = 3.73
    Ca_inf = 0.00024 (mM)
    tauR = 5 (ms)
    CaIC = 0.0003 (mM)
    CaNoiseIC = 0.00001 (mM)
    gH = 0.000005 (S/cm2)
    gInc = 2
    EH = -40 (mV)
    Cac = 0.002 (mM)
    nca = 4
    k2 = 0.0004 (1/ms)
    k4 = 0.001 (1/ms)
    PC = 0.007
    nexp = 1
}

STATE {
    mNa 
    hNa 
    nK 
    hT 
    Open 
    Pone 
    OpenLocked 
    cai (mM)
}

ASSIGNED {
    v (mV)
    iNa (mA/cm2)
    iK (mA/cm2)
    ileak (mA/cm2)
    ikleak (mA/cm2)
    iCa (mA/cm2)
    iH (mA/cm2)
    iPoisson_TC (mA/cm2)
    alphaM (1/ms)
    betaM (1/ms)
    iNa_alphaH (1/ms)
    iNa_betaH (1/ms)
    alphaN (1/ms)
    betaN (1/ms)
    ET (mV)
    Minf (1)
    T_Hinf (1)
    tauH (ms)
    tauS (ms)
    Hinf (1)
    alphaH (1/ms)
    betaH (1/ms)
    A (mM*cm2/mA/ms)
    drive_channel (mM/ms)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ET = (1000 * 8.31441 * (273.15 + 36) / (2 * 96846) * log(2 / cai))
    iNa = gNa * mNa^3 * hNa * (v - ENa)
    iK = gK * nK^4 * (v - EK)
    ileak = gLeak * (v - ELeak)
    ikleak = gKLeak * (v - EKLeak)
    iCa = gT * Minf^2 * hT * (v - ET)
    iH = gH * (Open + gInc * OpenLocked) * (v - EH)
}

DERIVATIVE states { 
    evaluate_fct(v)
    drive_channel = -A*iCa
    if (drive_channel <= 0.) { drive_channel = 0. } : cannot pump inward
    mNa' = alphaM * (1 - mNa) - betaM * mNa
    hNa' = iNa_alphaH * (1 - hNa) - iNa_betaH * hNa
    nK' = alphaN * (1 - nK) - betaN * nK
    hT' = (T_Hinf - hT) / tauH
    cai' = drive_channel + (Ca_inf - cai) / tauR
    Open' = alphaH * (1 - Open - OpenLocked) - betaH * Open
    Pone' = k2 * pow((cai / Cac), nca) * (1 - Pone) - k2 * Pone
    OpenLocked' = k4 * pow((Pone / PC), nexp) * Open - k4 * OpenLocked
}

UNITSOFF

INITIAL {
    A = 10000 / (2 * 96489)
    cai = cai
    mNa = mNa
    hNa = hNa
    nK = nK
    hT = hT
    Open = Open
    Pone = Pone
    OpenLocked = OpenLocked
}

PROCEDURE evaluate_fct(v (mV)) {
    LOCAL vNa, vK, vT
    vNa = v - vShiftNa
    vK = v - vShiftK
    vT = v + vShiftT

    alphaM = 0.32 * (13 - vNa) / (exp((13 - vNa) / 4) - 1)
    betaM = 0.28 * (vNa - 40) / (exp((vNa - 40) / 5) - 1)
    iNa_alphaH = 0.128 * exp((17 - vNa) / 18)
    iNa_betaH = 4 / (exp((40 - vNa) / 5) + 1)

    alphaN = 0.032 * (15 - vK) / (exp((15 - vK) / 5) - 1)
    betaN = 0.5 * exp((10 - vK) / 40)

    Minf = 1 / (1 + exp((-(vT + 57)) / 6.2))
    T_Hinf = 1 / (1 + exp((vT + 81) / 4))
    tauH = ((30.8 + (211.4 + exp((vT + 113.2) / 5)) / (1 + exp((vT + 84) / 3.2))) / phiH)

    tauS = (20 + 1000 / (exp((v + 71.5) / 14.2) + exp(-(v + 89) / 11.6)))
    Hinf = 1 / (1 + exp((v + 75) / 5.5))
    alphaH = Hinf / tauS
    betaH = (1 - Hinf) / tauS
}

UNITSON

FUNCTION max(x, y) {
    if (x > y) {
        max = x
    } else {
        max = y
    }
}

:VERBATIM
:#include <math.h>
:extern double scop_random();
:extern void set_seed(double);
:ENDVERBATIM
