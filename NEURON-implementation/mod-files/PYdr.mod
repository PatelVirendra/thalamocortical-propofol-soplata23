TITLE PYdr compartment with Poisson input (Safe Version)

NEURON {
    SUFFIX PYdr
    NONSPECIFIC_CURRENT iLeak, iNaP, iAR, iKCa, iHVA
    RANGE gLeak, ELeak, gNaP, ENaP, gAR, EAR, gKCa, EKCa, KD, gHVA, EHVA, tauCa, alphaCa, areaDR
    RANGE CaBuffer
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
    (molar) = (1/liter)
    (mM) = (millimolar)
    (uM) = (micromolar)
}

PARAMETER {
    gLeak = 0.000005 (S/cm2)
    ELeak = -60.95   (mV)
    gNaP  = 0.0000686 (S/cm2)
    ENaP  = 55       (mV)
    gAR   = 0.0000257 (S/cm2)
    EAR   = -100     (mV)
    gKCa  = 0.00057  (S/cm2)
    EKCa  = -100     (mV)
    KD    = 30       (uM)
    gHVA  = 0.00043  (S/cm2)
    EHVA  = 120      (mV)
    tauCa = 150      (ms)
    alphaCa = 5000   (uM/mA/ms)
    areaDR = 0.00035 (cm2)
}

ASSIGNED {
    v (mV)
    iLeak      (mA/cm2)
    iNaP       (mA/cm2)
    iAR        (mA/cm2)
    iKCa       (mA/cm2)
    iHVA       (mA/cm2)
    Minf_NaP
    Hinf_AR
    Minf_HVA
}

STATE {
    CaBuffer   (uM)
}

INITIAL {
    CaBuffer = CaBuffer
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    iLeak = gLeak * (v - ELeak)
    iNaP  = gNaP * Minf_NaP^3 * (v - ENaP)
    iAR   = gAR * Hinf_AR * (v - EAR)
    iKCa  = gKCa * (CaBuffer / (CaBuffer + KD)) * (v - EKCa)
    iHVA  = gHVA * Minf_HVA^2 * (v - EHVA)
}

DERIVATIVE states { 
    evaluate_fct(v)
    CaBuffer' = -alphaCa * areaDR * iHVA - CaBuffer / tauCa
}

PROCEDURE evaluate_fct(v (mV)) {
    Minf_NaP = 1 / (1 + exp(-(v + 55.7) / 7.7))
    Hinf_AR  = 1 / (1 + exp((v + 75) / 4))
    Minf_HVA = 1 / (1 + exp(-(v + 20) / 9))
}

:VERBATIM
:extern double scop_random();
:extern void set_seed(double);
:ENDVERBATIM