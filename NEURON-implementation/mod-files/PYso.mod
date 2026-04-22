NEURON {
    SUFFIX PYso
    NONSPECIFIC_CURRENT iLeak, iNa, iK, iA, iKS
    RANGE gLeak, ELeak, gNa, ENa, gK, EK, gA, EA, tauH, gKS, EKS
    RANGE phi_Na, hNaIC, hNaNoiseIC, phi_K, nKIC, nKNoiseIC, phi_A, hAIC, hANoiseIC, phi_KS, mKSIC, mKSNoiseIC
    RANGE hNa, nK, hA, mKS 
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (uA) = (microamp)
    (molar) = (1/liter)
    (mM) = (millimolar)
    (uM) = (micromolar)
    (S) = (siemens)
}

PARAMETER {
    gLeak = 0.0000667 (S/cm2)
    ELeak = -60.95 (mV)

    gNa = 0.050 (S/cm2)
    ENa = 55 (mV)
    phi_Na = 4
    hNaIC = 0.7
    hNaNoiseIC = 0.1

    gK = 0.0105 (S/cm2)
    EK = -100 (mV)
    phi_K = 4
    nKIC = 0.05
    nKNoiseIC = 0.05

    gA = 0.001 (S/cm2)
    EA = -100 (mV)
    tauH = 15 (ms)
    phi_A = 1
    hAIC = 0.1
    hANoiseIC = 0.1

    gKS = 0.000576 (S/cm2)
    EKS = -100 (mV)
    phi_KS = 1
    mKSIC = 0.005
    mKSNoiseIC = 0.001
}

ASSIGNED {
    v (mV)
    iLeak (mA/cm2)
    iNa (mA/cm2)
    iK (mA/cm2)
    iA (mA/cm2)
    iKS (mA/cm2)
    Minf_Na
    alphaM_Na (1/ms)
    betaM_Na (1/ms)
    alphaH_Na (1/ms)
    betaH_Na (1/ms)
    Minf_K
    alphaM_K (1/ms)
    betaM_K (1/ms)
    Minf_A
    Hinf_A
    Minf_KS
    tauM_KS (ms)
}

STATE {
    hNa
    nK
    hA
    mKS
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    iLeak = gLeak * (v - ELeak)
    iNa = gNa * Minf_Na^3 * hNa * (v - ENa) 
    iK = gK * nK^4 * (v - EK)
    iA = gA * Minf_A^3 * hA * (v - EA)
    iKS = gKS * mKS^3 * (v - EKS)
}

DERIVATIVE states {
    evaluate_fct(v) 
    hNa' = phi_Na * (alphaH_Na * (1 - hNa) - betaH_Na * hNa)
    nK' = phi_K * (alphaM_K * (1 - nK) - betaM_K * nK)
    hA' = phi_A * (Hinf_A - hA) / tauH
    mKS' = phi_KS * (Minf_KS - mKS) / tauM_KS
}

UNITSOFF

INITIAL {
    hNa = hNa
    nK = nK
    hA = hA
    mKS = mKS
}

PROCEDURE evaluate_fct(v (mV)) {
    alphaM_Na = 0.1 * (v + 33) / (1 - exp(-(v + 33) / 10))
    betaM_Na = 4 * exp(-(v + 53.7) / 12)
    Minf_Na = alphaM_Na / (alphaM_Na + betaM_Na)
    alphaH_Na = 0.07 * exp(-(v + 50) / 10)
    betaH_Na = 1 / (1 + exp(-(v + 20) / 10))

    alphaM_K = 0.01 * (v + 34) / (1 - exp(-(v + 34) / 10))
    betaM_K = 0.125 * exp(-(v + 44) / 25)
    Minf_K = alphaM_K / (alphaM_K + betaM_K)

    Minf_A = 1 / (1 + exp(-(v + 50) / 20))
    Hinf_A = 1 / (1 + exp((v + 80) / 6))

    Minf_KS = 1 / (1 + exp(-(v + 34) / 6.5))
    tauM_KS = 8 / (exp(-(v + 55) / 30) + exp((v + 55) / 30))
}

UNITSON

:VERBATIM
:#include <math.h>
:extern double scop_random();
:extern void set_seed(double);
:ENDVERBATIM
