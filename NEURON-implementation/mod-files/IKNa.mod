NEURON {
    POINT_PROCESS IKNa
    :POINTER v_pre             : Access presynaptic voltage from PYdr (from target_var)
    NONSPECIFIC_CURRENT iKNa  : Only iKNa affects postsynaptic v
    RANGE gKNa, EKNa, concNaIC, concNaNoiseIC, alphaNa, RPump, eqNa, v_pre
    RANGE gNa, ENa, gNaP, ENaP, phi_Na, areaPYdr, areaPYso, area_post
    RANGE concNa, hNalocal, iNa_local, iNaP_local  : State and local currents
}

UNITS {
    (nA) = (nanoamp)
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (molar) = (1/liter)
    (mM) = (millimolar)
    (uM) = (micromolar)
}

PARAMETER {
    gKNa = 0.00133 (S/cm2)
    EKNa = -100 (mV)
    concNaIC = 12 (mM)
    concNaNoiseIC = 3 (mM)
    alphaNa = 10000 (mM/mA/ms)  : Adjusted for NEURON units
    RPump = 0.018 (mM/ms)
    eqNa = 9.5 (mM)
    : Sodium channel parameters (for iNa_local)
    gNa = 0.05 (S/cm2)
    ENa = 55 (mV)
    : Persistent sodium channel parameters (for iNaP_local)
    gNaP = 0.0000686 (S/cm2)
    ENaP = 55 (mV)
    phi_Na = 4  : Temperature scaling factor
    areaPYdr = 0.00035 (cm2)
    areaPYso = 0.00015 (cm2)
    area_post = 0.00015 (cm2) : Postsynaptic neuron area
    mm = 1 (mM)
}

ASSIGNED {
    v (mV)  : Postsynaptic voltage (PYso)
    v_pre (mV)  : Presynaptic voltage (PYdr)
    iKNa (nA)  : Current affecting v
    iNa_local (mA/cm2)  : Local Na+ current for concNa calculation
    iNaP_local (mA/cm2)  : Local persistent Na+ current for concNa calculation
    : Local Na+ channel gating (for iNa_local)
    Minf_Na_local (1)
    alphaM_Na_local (/ms)
    betaM_Na_local (/ms)
    alphaH_Na_local (/ms)
    betaH_Na_local (/ms)
    : Persistent Na+ channel gating (for iNaP_local)
    Minf_NaP_local (1)
    eqNaPumpTerm (1)
}

STATE {
    concNa (mM)  : Sodium concentration
    hNalocal (1)  : Inactivation gate for local Na+ current
}

INITIAL {
    : Initialize with default values (can be overridden by Python)
    concNa = concNa
    hNalocal = hNalocal
    
    : Compute initial values for assigned variables
    evaluate_fct(v)
    eqNaPumpTerm = pow(eqNa, 3) / (pow(eqNa, 3) + pow(15 * mm, 3))
}

BREAKPOINT {
    : Evaluate gating variables first
    evaluate_fct(v)
    
    : Solve ODEs using derivimplicit (CVode compatible for complex equations)
    SOLVE states METHOD derivimplicit
    
    : Calculate local currents
    iNa_local = gNa * Minf_Na_local * Minf_Na_local * Minf_Na_local * hNalocal * (v - ENa)
    
    : Persistent Na+ channel (presynaptic, PYdr)
    iNaP_local = gNaP * Minf_NaP_local * Minf_NaP_local * Minf_NaP_local * (v_pre - ENaP)
    
    : K(Na) current (postsynaptic, PYso)
    iKNa = gKNa * (0.37 / (1 + pow(38.7 * mm / concNa, 3.5))) * (v - EKNa) * 1e6 * area_post
}

DERIVATIVE states {
    LOCAL pumpTerm, eqPumpTerm
    
    : hNalocal dynamics
    hNalocal' = phi_Na * (alphaH_Na_local * (1 - hNalocal) - betaH_Na_local * hNalocal)
    
    : Sodium concentration dynamics
    eqPumpTerm = pow(eqNa, 3) / (pow(eqNa, 3) + pow(15 * mm, 3))
    pumpTerm = pow(concNa, 3) / (pow(concNa, 3) + pow(15 * mm, 3))
    
    concNa' = -alphaNa * (areaPYso * iNa_local + areaPYdr * iNaP_local) - RPump * (pumpTerm - eqPumpTerm)
}

PROCEDURE evaluate_fct(v (mV)) {
    : Local Na+ channel (postsynaptic, PYso)
    alphaM_Na_local = 0.1 * (v + 33) / (1 - exp(-(v + 33) / 10))
    betaM_Na_local = 4 * exp(-(v + 53.7) / 12)
    Minf_Na_local = alphaM_Na_local / (alphaM_Na_local + betaM_Na_local)
    alphaH_Na_local = 0.07 * exp(-(v + 50) / 10)
    betaH_Na_local = 1 / (1 + exp(-(v + 20) / 10))
    
    : Persistent Na+ channel (presynaptic, PYdr)
    Minf_NaP_local = 1 / (1 + exp(-(v_pre + 55.7) / 7.7))
}

COMMENT
================================================================================
IKNa - Sodium-dependent Potassium Current

This mechanism computes the K(Na) current based on intracellular sodium
concentration, which is modulated by Na+ influx through voltage-gated
and persistent sodium channels.

SOLVE METHOD:
Uses 'derivimplicit' instead of 'cnexp' because the concNa' equation
contains cubic terms that are not compatible with cnexp. The derivimplicit
method works well with CVode for complex, potentially stiff equations.

POINTER v_pre:
Receives presynaptic voltage from PYdr cell via target_var for computing
the persistent sodium current contribution to sodium concentration.

USAGE IN PYTHON:
    syn = h.IKNa(pyso_cell.soma(0.5))
    syn.gKNa = config.gKNa
    syn.concNa = 12 + random * 3  # Initial sodium concentration
    syn.hNalocal = 0.5 + random * 0.1  # Initial inactivation
    # ... other parameters ...
    
    # Voltage transfer from PYdr
    config.pc.target_var(syn, syn._ref_v_pre, pydr_src_gid + gid_offset)

================================================================================
ENDCOMMENT
