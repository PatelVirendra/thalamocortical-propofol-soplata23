"""
cells.py — Brian2 NeuronGroup Builders
=======================================
Defines one factory function per cell type.  Each function returns a
Brian2 NeuronGroup together with its associated PoissonGroup background
drive and the synapse that couples them.

Cell types
----------
PYdr  — Pyramidal dendrite (two-compartment dendritic compartment)
PYso  — Pyramidal soma     (two-compartment somatic compartment)
IN    — Cortical inhibitory interneuron
TC    — Thalamocortical relay neuron
RE    — Reticular nucleus neuron  (TRN)

All parameters are read directly from config.py so this file never
needs to be edited for condition sweeps.
"""

import numpy as np
import brian2 as b2
import config
from utils import get_independent_draws

# =============================================================================
# HELPER
# =============================================================================

def _make_poisson_drive(target_group, n_cells, baseline_rate, gext, tau, Eext,
                        var_name):
    """
    Create an independent per-cell Poisson background drive that mimics
    the NEURON NetStim → ExpSyn pattern.

    Returns (PoissonGroup, Synapses) — both must be added to the Network.
    """
    stim = b2.PoissonGroup(n_cells, rates=baseline_rate)
    syn = b2.Synapses(
        stim, target_group,
        on_pre=f'{var_name} += {gext!r}'
    )
    syn.connect(j='i')   # 1-to-1 mapping
    return stim, syn


# =============================================================================
# PYdr — Pyramidal Dendritic Compartment
# =============================================================================

def make_PYdr_group(params):
    N   = config.Npydr
    uM  = config.uM

    # --- Biophysical constants ---
    Cm_PYdr      = 1     * b2.uF  / b2.cm**2
    gLeak_PYdr   = 0.005 * b2.msiemens / b2.cm**2
    ELeak_PYdr   = -60.95 * b2.mV

    gNaP_PYdr    = 0.0686 * b2.msiemens / b2.cm**2
    ENaP_PYdr    = 55     * b2.mV

    gAR_PYdr     = 0.0257 * b2.msiemens / b2.cm**2
    EAR_PYdr     = -100   * b2.mV

    tauCa_PYdr   = 150    * b2.ms
    alphaCa_PYdr = 0.005 * 1000 * uM / (b2.uA * b2.ms)
    areaDR       = 0.00035 * b2.cm**2

    gKCa_PYdr    = 0.57   * b2.msiemens / b2.cm**2
    EKCa_PYdr    = -100   * b2.mV
    KD_PYdr      = 30     * uM

    gHVA_PYdr    = 0.43   * b2.msiemens / b2.cm**2
    EHVA_PYdr    = 120    * b2.mV

    # Background noise conductance
    gext_PYdr = params['gext_PYdr'] * b2.siemens / b2.cm**2

    # --- Neuron equations ---
    eqs_PYdr = '''
    dv/dt = (iLeak_PYdr + iNaP + iAR + iKCa + iHVA + iCOM
             + isyn_PYdr + iPoisson_PYdr) / Cm_PYdr : volt

    isyn_PYdr         = iAMPA_PYso_PYdr + iNMDA_PYso_PYdr + iAMPA_TC_PYdr : amp/meter**2
    iAMPA_PYso_PYdr   : amp/meter**2
    iNMDA_PYso_PYdr   : amp/meter**2
    iAMPA_TC_PYdr     : amp/meter**2

    iCOM              : amp/meter**2

    iLeak_PYdr = -gLeak_PYdr * (v - ELeak_PYdr) : amp/meter**2

    Minf_NaP = 1 / (1 + exp(-(v/mV + 55.7) / 7.7)) : 1
    iNaP = -gNaP_PYdr * Minf_NaP**3 * (v - ENaP_PYdr)  : amp/meter**2

    Hinf_AR = 1 / (1 + exp((v/mV + 75) / 4)) : 1
    iAR = -gAR_PYdr * Hinf_AR * (v - EAR_PYdr)          : amp/meter**2

    dCaBuffer/dt = alphaCa_PYdr * areaDR * iHVA - CaBuffer / tauCa_PYdr : mM

    iKCa = -gKCa_PYdr * (CaBuffer / (CaBuffer + KD_PYdr)) * (v - EKCa_PYdr) : amp/meter**2

    Minf_HVA = 1 / (1 + exp(-(v/mV + 20) / 9)) : 1
    iHVA = -gHVA_PYdr * Minf_HVA**2 * (v - EHVA_PYdr)   : amp/meter**2

    dgPoisson_PYdr/dt = -gPoisson_PYdr / tau_PYdr : siemens/meter**2
    iPoisson_PYdr = -gPoisson_PYdr * (v - Eext_PYdr)     : amp/meter**2
    '''

    group = b2.NeuronGroup(
        N, eqs_PYdr,
        method='euler',
        threshold='v > -25*mV and not_refractory',
        reset='',
        refractory='v > -25*mV',
        namespace={
            'Cm_PYdr': Cm_PYdr,
            'gLeak_PYdr': gLeak_PYdr, 'ELeak_PYdr': ELeak_PYdr,
            'gNaP_PYdr': gNaP_PYdr,   'ENaP_PYdr': ENaP_PYdr,
            'gAR_PYdr': gAR_PYdr,     'EAR_PYdr': EAR_PYdr,
            'tauCa_PYdr': tauCa_PYdr, 'alphaCa_PYdr': alphaCa_PYdr,
            'areaDR': areaDR,
            'gKCa_PYdr': gKCa_PYdr,   'EKCa_PYdr': EKCa_PYdr,
            'KD_PYdr': KD_PYdr,
            'gHVA_PYdr': gHVA_PYdr,   'EHVA_PYdr': EHVA_PYdr,
            'tau_PYdr': config.tau_PYdr,
            'Eext_PYdr': config.Eext_PYdr,
        },
        name='PYdr_group'
    )

    # --- Exact Random123 Initialization (Type ID: 1) ---
    gids = np.arange(config.Npyso, config.Npyso + config.Npydr)
    draws = get_independent_draws(gids, type_id=1, global_seed=config.randSeed, num_draws=2)
    
    group.CaBuffer = (0.001 + 0.01 * draws[0]) * uM
    group.v        = (-68 + 20 * draws[1]) * b2.mV
    group.gPoisson_PYdr = 0 * b2.msiemens / b2.cm**2

    # --- Poisson background drive ---
    poisson_stim, poisson_syn = _make_poisson_drive(
        group, N,
        config.baseline_PYdr, gext_PYdr, config.tau_PYdr, config.Eext_PYdr,
        var_name='gPoisson_PYdr'
    )

    return group, poisson_stim, poisson_syn


# =============================================================================
# PYso — Pyramidal Somatic Compartment
# =============================================================================

def make_PYso_group(params):
    N = config.Npyso

    # --- Biophysical constants ---
    Cm_PYso      = 1      * b2.uF / b2.cm**2
    gLeak_PYso   = 0.0667 * b2.msiemens / b2.cm**2
    ELeak_PYso   = -60.95 * b2.mV

    gNa_PYso     = 50    * b2.msiemens / b2.cm**2
    ENa_PYso     = 55    * b2.mV
    phi_Na_PYso  = 4

    gK_PYso      = 10.5  * b2.msiemens / b2.cm**2
    EK_PYso      = -100  * b2.mV
    phi_K_PYso   = 4

    gA_PYso      = 1     * b2.msiemens / b2.cm**2
    EA_PYso      = -100  * b2.mV
    tauH_A_PYso  = 15    * b2.ms
    phi_A_PYso   = 1

    gKS_PYso     = 0.576 * b2.msiemens / b2.cm**2
    EKS_PYso     = -100  * b2.mV
    phi_KS_PYso  = 1

    eqs_PYso = '''
    dv/dt = (iLeak_PYso + iNa_PYso + iK_PYso + iA_PYso + iKS_PYso
             + iKNa + iCOM + isyn_PYso) / Cm_PYso : volt

    isyn_PYso       = iGABAA_IN_PYso : amp/meter**2
    iGABAA_IN_PYso  : amp/meter**2

    iCOM            : amp/meter**2
    iKNa            : amp/meter**2

    iLeak_PYso = -gLeak_PYso * (v - ELeak_PYso) : amp/meter**2

    alphaM_Na = (0.1  * (v/mV + 33)  / (1 - exp(-(v/mV + 33)  / 10))) / ms : Hz
    betaM_Na  = (4    * exp(-(v/mV + 53.7) / 12)) / ms                      : Hz
    alphaH_Na = (0.07 * exp(-(v/mV + 50)   / 10)) / ms                      : Hz
    betaH_Na  = (1 / (1 + exp(-(v/mV + 20) / 10))) / ms                     : Hz
    Minf_Na   = alphaM_Na / (alphaM_Na + betaM_Na)                           : 1
    dhNa/dt   = phi_Na_PYso * (alphaH_Na * (1 - hNa) - betaH_Na * hNa)      : 1
    iNa_PYso  = -gNa_PYso * Minf_Na**3 * hNa * (v - ENa_PYso)               : amp/meter**2

    alphaM_K  = (0.01 * (v/mV + 34) / (1 - exp(-(v/mV + 34) / 10))) / ms   : Hz
    betaM_K   = (0.125 * exp(-(v/mV + 44) / 25)) / ms                       : Hz
    Minf_K    = alphaM_K / (alphaM_K + betaM_K)                             : 1
    dnK/dt    = phi_K_PYso * (alphaM_K * (1 - nK) - betaM_K * nK)           : 1
    iK_PYso   = -gK_PYso * nK**4 * (v - EK_PYso)                            : amp/meter**2

    Minf_A    = 1 / (1 + exp(-(v/mV + 50) / 20))                            : 1
    Hinf_A    = 1 / (1 + exp( (v/mV + 80) / 6))                             : 1
    dhA/dt    = phi_A_PYso * (Hinf_A - hA) / tauH_A_PYso                    : 1
    iA_PYso   = -gA_PYso * Minf_A**3 * hA * (v - EA_PYso)                   : amp/meter**2

    Minf_KS   = 1 / (1 + exp(-(v/mV + 34) / 6.5))                          : 1
    tauM_KS   = (8 / (exp(-(v/mV + 55) / 30) + exp((v/mV + 55) / 30))) * ms: second
    dmKS/dt   = phi_KS_PYso * (Minf_KS - mKS) / tauM_KS                    : 1
    iKS_PYso  = -gKS_PYso * mKS**3 * (v - EKS_PYso)                         : amp/meter**2
    '''

    group = b2.NeuronGroup(
        N, eqs_PYso,
        method='euler',
        threshold='v > -25*mV and not_refractory',
        reset='',
        refractory='v > -25*mV',
        namespace={
            'Cm_PYso': Cm_PYso,
            'gLeak_PYso': gLeak_PYso, 'ELeak_PYso': ELeak_PYso,
            'gNa_PYso': gNa_PYso,     'ENa_PYso': ENa_PYso,
            'phi_Na_PYso': phi_Na_PYso,
            'gK_PYso': gK_PYso,       'EK_PYso': EK_PYso,
            'phi_K_PYso': phi_K_PYso,
            'gA_PYso': gA_PYso,       'EA_PYso': EA_PYso,
            'tauH_A_PYso': tauH_A_PYso, 'phi_A_PYso': phi_A_PYso,
            'gKS_PYso': gKS_PYso,     'EKS_PYso': EKS_PYso,
            'phi_KS_PYso': phi_KS_PYso,
        },
        name='PYso_group'
    )

    # --- Exact Random123 Initialization (Type ID: 0) ---
    gids = np.arange(0, config.Npyso)
    draws = get_independent_draws(gids, type_id=0, global_seed=config.randSeed, num_draws=5)
    
    group.hNa = 0.7 + 0.1 * draws[0]
    group.nK  = 0.05 + 0.05 * draws[1]
    group.hA  = 0.1 + 0.1 * draws[2]
    group.mKS = 0.005 + 0.001 * draws[3]
    group.v   = (-68 + 20 * draws[4]) * b2.mV

    return group


# =============================================================================
# IN — Cortical Inhibitory Interneuron
# =============================================================================

def make_IN_group(params):
    N = config.Ninh

    Cm_IN      = 1      * b2.uF / b2.cm**2
    gNa_IN     = 35     * b2.msiemens / b2.cm**2
    ENa_IN     = 55     * b2.mV
    gK_IN      = 9      * b2.msiemens / b2.cm**2
    EK_IN      = -90    * b2.mV
    gLeak_IN   = 0.1025 * b2.msiemens / b2.cm**2
    ELeak_IN   = -63.8  * b2.mV

    eqs_IN = '''
    dv/dt = (iNa_IN + iK_IN + iLeak_IN + isyn_IN) / Cm_IN : volt

    isyn_IN         = iAMPA_PYso_IN + iNMDA_PYso_IN + iGABAA_IN_IN + iAMPA_TC_IN : amp/meter**2
    iAMPA_PYso_IN   : amp/meter**2
    iNMDA_PYso_IN   : amp/meter**2
    iGABAA_IN_IN    : amp/meter**2
    iAMPA_TC_IN     : amp/meter**2

    iLeak_IN = -gLeak_IN * (v - ELeak_IN) : amp/meter**2

    alphaM_IN = (0.5  * (v/mV + 35) / (1 - exp(-(v/mV + 35) / 10))) / ms  : Hz
    betaM_IN  = (20   * exp(-(v/mV + 60) / 18))                       / ms  : Hz
    alphaH_IN = (0.35 * exp(-(v/mV + 58) / 20))                       / ms  : Hz
    betaH_IN  = (5 / (1 + exp(-(v/mV + 28) / 10)))                    / ms  : Hz
    Minf_IN   = alphaM_IN / (alphaM_IN + betaM_IN)                          : 1
    alphaN_IN = (0.05  * (v/mV + 34) / (1 - exp(-(v/mV + 34) / 10))) / ms  : Hz
    betaN_IN  = (0.625 * exp(-(v/mV + 44) / 80))                      / ms  : Hz
    dhNa_IN/dt = alphaH_IN * (1 - hNa_IN) - betaH_IN * hNa_IN               : 1
    dnK_IN/dt  = alphaN_IN * (1 - nK_IN)  - betaN_IN * nK_IN                : 1
    iNa_IN = -gNa_IN * Minf_IN**3 * hNa_IN * (v - ENa_IN)                   : amp/meter**2

    iK_IN  = -gK_IN * nK_IN**4 * (v - EK_IN)                                : amp/meter**2
    '''

    group = b2.NeuronGroup(
        N, eqs_IN,
        method='euler',
        threshold='v > -25*mV and not_refractory',
        reset='',
        refractory='v > -25*mV',
        namespace={
            'Cm_IN': Cm_IN,
            'gNa_IN': gNa_IN, 'ENa_IN': ENa_IN,
            'gK_IN':  gK_IN,  'EK_IN':  EK_IN,
            'gLeak_IN': gLeak_IN, 'ELeak_IN': ELeak_IN,
        },
        name='IN_group'
    )

    # --- Exact Random123 Initialization (Type ID: 2) ---
    gids = np.arange(config.Npyso + config.Npydr, config.Npyso + config.Npydr + config.Ninh)
    draws = get_independent_draws(gids, type_id=2, global_seed=config.randSeed, num_draws=3)
    
    group.hNa_IN = 0.7 + 0.1 * draws[0]
    group.nK_IN  = 0.13 + 0.01 * draws[1]
    group.v      = (-68 + 20 * draws[2]) * b2.mV

    return group


# =============================================================================
# TC — Thalamocortical Relay Neuron
# =============================================================================

def make_TC_group(params):
    N   = config.Ntc
    uM  = config.uM

    Cm_TC        = 1     * b2.uF  / b2.cm**2
    gNa_TC       = 90    * b2.msiemens / b2.cm**2
    ENa_TC       = 50    * b2.mV
    vShiftNa_TC  = -40

    gK_TC        = 10    * b2.msiemens / b2.cm**2
    EK_TC        = -100  * b2.mV
    vShiftK_TC   = -25

    gLeak_TC     = 0.01  * b2.msiemens / b2.cm**2
    ELeak_TC     = -70   * b2.mV
    EKLeak_TC    = -100  * b2.mV

    gT_TC        = 2     * b2.msiemens / b2.cm**2
    vShiftT_TC   = 2
    phiH_TC      = 3.73
    Ca_inf_TC    = 0.00024 * b2.mM
    A_TC         = (10 / (2 * 96489)) * b2.mM * b2.cm**2 / (b2.uA * b2.ms)
    tauR_TC      = 5  * b2.ms

    gInc_TC      = 2
    EH_TC        = -40   * b2.mV
    Cac_TC       = 0.002 * b2.mM
    nca_TC       = 4
    k2_TC        = 0.0004 / b2.ms
    k4_TC        = 0.001  / b2.ms
    pc_TC        = 0.007
    nexp_TC      = 1

    gH_TC     = params['gH']     * b2.siemens / b2.cm**2
    gKLeak_TC = params['gKLeak'] * b2.siemens / b2.cm**2
    gext_TC   = params['gext_TC'] * b2.siemens / b2.cm**2

    eqs_TC = '''
    dv/dt = (iNa_TC + iK_TC + iLeak_TC + iKLeak_TC
             - iT_TC + iH_TC + isyn_TC + iPoisson_TC) / Cm_TC : volt

    isyn_TC          = iGABAA_RE_TC + iGABAB_RE_TC + iAMPA_PYso_TC : amp/meter**2
    iGABAA_RE_TC     : amp/meter**2
    iGABAB_RE_TC     : amp/meter**2
    iAMPA_PYso_TC    : amp/meter**2

    iLeak_TC  = -gLeak_TC  * (v - ELeak_TC)  : amp/meter**2
    iKLeak_TC = -gKLeak_TC * (v - EKLeak_TC) : amp/meter**2

    alphaM_TC = 0.32 * (13 - (v/mV - vShiftNa_TC)) / (exp((13 - (v/mV - vShiftNa_TC)) / 4) - 1) / ms : Hz
    betaM_TC  = 0.28 * ((v/mV - vShiftNa_TC) - 40) / (exp(((v/mV - vShiftNa_TC) - 40) / 5) - 1) / ms : Hz
    alphaHNa_TC = 0.128 * exp((17 - (v/mV - vShiftNa_TC)) / 18) / ms : Hz
    betaHNa_TC  = 4 / (exp((40 - (v/mV - vShiftNa_TC)) / 5) + 1) / ms : Hz
    dmNa_TC/dt  = alphaM_TC * (1 - mNa_TC) - betaM_TC * mNa_TC : 1
    dhNa_TC/dt  = alphaHNa_TC * (1 - hNa_TC) - betaHNa_TC * hNa_TC : 1
    iNa_TC = -gNa_TC * mNa_TC**3 * hNa_TC * (v - ENa_TC)           : amp/meter**2

    alphaN_TC = 0.032 * (15 - (v/mV - vShiftK_TC)) / (exp((15 - (v/mV - vShiftK_TC)) / 5) - 1) / ms : Hz
    betaN_TC  = 0.5   * exp((10 - (v/mV - vShiftK_TC)) / 40) / ms : Hz
    dnK_TC/dt = alphaN_TC * (1 - nK_TC) - betaN_TC * nK_TC : 1
    iK_TC = -gK_TC * nK_TC**4 * (v - EK_TC) : amp/meter**2

    ET_TC     = (1000 * 8.31441 * (273.15 + 36) / (2 * 96846) * log((2 * mM)/ Ca_TC)) * mV : volt
    Minf_T_TC = 1 / (1 + exp((-(v/mV + vShiftT_TC + 57)) / 6.2))  : 1
    Ca_Hinf   = 1 / (1 + exp( (v/mV + vShiftT_TC + 81) / 4))      : 1
    tauHT_TC  = ((30.8 + (211.4 + exp((v/mV + vShiftT_TC + 113.2) / 5)) /
                  (1 + exp((v/mV + vShiftT_TC + 84) / 3.2))) / phiH_TC) * ms : second
    dhT_TC/dt = (Ca_Hinf - hT_TC) / tauHT_TC : 1
    iT_TC     = gT_TC * Minf_T_TC**2 * hT_TC * (v - ET_TC)        : amp/meter**2
    dCa_extra = clip(-A_TC * iT_TC, 0*mM/second, inf*mM/second)   : mM/second (constant over dt)
    dCa_TC/dt = dCa_extra + (Ca_inf_TC - Ca_TC) / tauR_TC         : mM

    tauS_TC = (20 + 1000 / (exp((v/mV + 71.5) / 14.2) + exp(-(v/mV + 89) / 11.6))) * ms : second
    Hinf_TC = 1 / (1 + exp((v/mV + 75) / 5.5))  : 1
    alphaH_TC = Hinf_TC / tauS_TC                : Hz
    betaH_TC  = (1 - Hinf_TC) / tauS_TC          : Hz
    dOpen/dt  = alphaH_TC * (1 - Open - OpenLocked) - betaH_TC * Open : 1
    dPone/dt  = k2_TC * (Ca_TC / Cac_TC)**nca_TC * (1 - Pone) - k2_TC * Pone : 1
    dOpenLocked/dt = k4_TC * (Pone / pc_TC)**nexp_TC * Open - k4_TC * OpenLocked : 1
    iH_TC = -gH_TC * (Open + gInc_TC * OpenLocked) * (v - EH_TC) : amp/meter**2

    dgPoisson_TC/dt = -gPoisson_TC / tau_TC : siemens/meter**2
    iPoisson_TC = -gPoisson_TC * (v - Eext_TC) : amp/meter**2
    '''

    group = b2.NeuronGroup(
        N, eqs_TC,
        method='euler',
        threshold='v > -25*mV and not_refractory',
        reset='',
        refractory='v > -25*mV',
        namespace={
            'Cm_TC': Cm_TC,
            'gNa_TC': gNa_TC, 'ENa_TC': ENa_TC, 'vShiftNa_TC': vShiftNa_TC,
            'gK_TC': gK_TC,   'EK_TC': EK_TC,   'vShiftK_TC': vShiftK_TC,
            'gLeak_TC': gLeak_TC,   'ELeak_TC': ELeak_TC,
            'gKLeak_TC': gKLeak_TC, 'EKLeak_TC': EKLeak_TC,
            'gT_TC': gT_TC, 'vShiftT_TC': vShiftT_TC, 'phiH_TC': phiH_TC,
            'Ca_inf_TC': Ca_inf_TC, 'A_TC': A_TC, 'tauR_TC': tauR_TC,
            'gH_TC': gH_TC, 'gInc_TC': gInc_TC, 'EH_TC': EH_TC,
            'Cac_TC': Cac_TC, 'nca_TC': nca_TC,
            'k2_TC': k2_TC, 'k4_TC': k4_TC, 'pc_TC': pc_TC, 'nexp_TC': nexp_TC,
            'tau_TC': config.tau_TC, 'Eext_TC': config.Eext_TC,
        },
        name='TC_group'
    )

    # --- Exact Random123 Initialization (Type ID: 3) ---
    gids = np.arange(config.Npyso + config.Npydr + config.Ninh, config.Npyso + config.Npydr + config.Ninh + config.Ntc)
    draws = get_independent_draws(gids, type_id=3, global_seed=config.randSeed, num_draws=9)

    group.Ca_TC      = (0.0003 + 0.00001 * draws[0]) * b2.mM
    group.mNa_TC     = 0.00007 + 0.00001 * draws[1]
    group.hNa_TC     = 0.8 + 0.1 * draws[2]
    group.nK_TC      = 0.00025 + 0.00001 * draws[3]
    group.hT_TC      = 0.01 + 0.005 * draws[4]
    group.Open       = 0.05 + 0.01 * draws[5]
    group.Pone       = 0.06 + 0.01 * draws[6]
    group.OpenLocked = 0.55 + 0.01 * draws[7]
    group.v          = (-68 + 20 * draws[8]) * b2.mV
    group.gPoisson_TC = 0 * b2.msiemens / b2.cm**2

    poisson_stim, poisson_syn = _make_poisson_drive(
        group, N,
        config.baseline_TC, gext_TC, config.tau_TC, config.Eext_TC,
        var_name='gPoisson_TC'
    )

    return group, poisson_stim, poisson_syn


# =============================================================================
# RE — Reticular Nucleus Neuron  (TRN)
# =============================================================================

def make_RE_group(params):
    N = config.Nre

    Cm_RE        = 1     * b2.uF  / b2.cm**2
    gNa_RE       = 200   * b2.msiemens / b2.cm**2
    ENa_RE       = 50    * b2.mV
    vShiftNa_RE  = -55

    gK_RE        = 20    * b2.msiemens / b2.cm**2
    EK_RE        = -100  * b2.mV
    vShiftK_RE   = -55

    gLeak_RE     = 0.05  * b2.msiemens / b2.cm**2
    ELeak_RE     = -90   * b2.mV
    gKLeak_RE    = 0.0   * b2.msiemens / b2.cm**2
    EKLeak_RE    = -95   * b2.mV

    gT_RE        = 3     * b2.msiemens / b2.cm**2
    ET_RE        = 120   * b2.mV
    vShiftT_RE   = 4
    phiM_RE      = 6.81
    phiH_RE      = 3.73

    gext_RE = params['gext_RE'] * b2.siemens / b2.cm**2

    eqs_RE = '''
    dv/dt = (iNa_RE + iK_RE + iLeak_RE + iKLeak_RE
             + iT_RE + isyn_RE + iPoisson_RE) / Cm_RE : volt

    isyn_RE         = iAMPA_TC_RE + iGABAA_RE_RE + iAMPA_PYso_RE : amp/meter**2
    iAMPA_TC_RE     : amp/meter**2
    iGABAA_RE_RE    : amp/meter**2
    iAMPA_PYso_RE   : amp/meter**2

    iLeak_RE  = -gLeak_RE  * (v - ELeak_RE)  : amp/meter**2
    iKLeak_RE = -gKLeak_RE * (v - EKLeak_RE) : amp/meter**2

    alphaM_RE = 0.32 * (13 - (v/mV - vShiftNa_RE)) / (exp((13 - (v/mV - vShiftNa_RE)) / 4) - 1) / ms : Hz
    betaM_RE  = 0.28 * ((v/mV - vShiftNa_RE) - 40) / (exp(((v/mV - vShiftNa_RE) - 40) / 5) - 1) / ms : Hz
    alphaH_RE = 0.128 * exp((17 - (v/mV - vShiftNa_RE)) / 18) / ms : Hz
    betaH_RE  = 4 / (exp((40 - (v/mV - vShiftNa_RE)) / 5) + 1) / ms : Hz
    dmNa_RE/dt = alphaM_RE * (1 - mNa_RE) - betaM_RE * mNa_RE : 1
    dhNa_RE/dt = alphaH_RE * (1 - hNa_RE) - betaH_RE * hNa_RE : 1
    iNa_RE = -gNa_RE * mNa_RE**3 * hNa_RE * (v - ENa_RE) : amp/meter**2

    alphaN_RE = 0.032 * (15 - (v/mV - vShiftK_RE)) / (exp((15 - (v/mV - vShiftK_RE)) / 5) - 1) / ms : Hz
    betaN_RE  = 0.5   * exp((10 - (v/mV - vShiftK_RE)) / 40) / ms : Hz
    dnK_RE/dt = alphaN_RE * (1 - nK_RE) - betaN_RE * nK_RE : 1
    iK_RE = -gK_RE * nK_RE**4 * (v - EK_RE) : amp/meter**2

    Minf_T_RE  = 1 / (1 + exp((-(v/mV + vShiftT_RE + 50)) / 7.4)) : 1
    tauM_T_RE  = ((3.0 + 1.0 / (exp((v/mV + vShiftT_RE + 25) / 10)
                  + exp(-(v/mV + vShiftT_RE + 100) / 15))) / phiM_RE) * ms : second
    Hinf_T_RE  = 1 / (1 + exp((v/mV + vShiftT_RE + 78) / 5))      : 1
    tauH_T_RE  = ((85.0 + 1.0 / (exp((v/mV + vShiftT_RE + 46) / 4)
                  + exp(-(v/mV + vShiftT_RE + 405) / 50))) / phiH_RE) * ms : second
    dmT_RE/dt  = (Minf_T_RE - mT_RE) / tauM_T_RE : 1
    dhT_RE/dt  = (Hinf_T_RE - hT_RE) / tauH_T_RE : 1
    iT_RE = -gT_RE * mT_RE**2 * hT_RE * (v - ET_RE)               : amp/meter**2

    dgPoisson_RE/dt = -gPoisson_RE / tau_RE : siemens/meter**2
    iPoisson_RE = -gPoisson_RE * (v - Eext_RE) : amp/meter**2
    '''

    group = b2.NeuronGroup(
        N, eqs_RE,
        method='euler',
        threshold='v > -25*mV and not_refractory',
        reset='',
        refractory='v > -25*mV',
        namespace={
            'Cm_RE': Cm_RE,
            'gNa_RE': gNa_RE,   'ENa_RE': ENa_RE,   'vShiftNa_RE': vShiftNa_RE,
            'gK_RE': gK_RE,     'EK_RE': EK_RE,     'vShiftK_RE': vShiftK_RE,
            'gLeak_RE': gLeak_RE,   'ELeak_RE': ELeak_RE,
            'gKLeak_RE': gKLeak_RE, 'EKLeak_RE': EKLeak_RE,
            'gT_RE': gT_RE, 'ET_RE': ET_RE, 'vShiftT_RE': vShiftT_RE,
            'phiM_RE': phiM_RE, 'phiH_RE': phiH_RE,
            'tau_RE': config.tau_RE, 'Eext_RE': config.Eext_RE,
        },
        name='RE_group'
    )

    # --- Exact Random123 Initialization (Type ID: 4) ---
    gids = np.arange(config.Npyso + config.Npydr + config.Ninh + config.Ntc, config.Npyso + config.Npydr + config.Ninh + config.Ntc + config.Nre)
    draws = get_independent_draws(gids, type_id=4, global_seed=config.randSeed, num_draws=6)

    group.mNa_RE = 0.00002 + 0.00001 * draws[0]
    group.hNa_RE = 0.8 + 0.1 * draws[1]
    group.nK_RE  = 0.00015 + 0.00001 * draws[2]
    group.mT_RE  = 0.01 + 0.01 * draws[3]
    group.hT_RE  = 0.6 + 0.01 * draws[4]
    group.v      = (-68 + 20 * draws[5]) * b2.mV
    group.gPoisson_RE = 0 * b2.msiemens / b2.cm**2

    poisson_stim, poisson_syn = _make_poisson_drive(
        group, N,
        config.baseline_RE, gext_RE, config.tau_RE, config.Eext_RE,
        var_name='gPoisson_RE'
    )

    return group, poisson_stim, poisson_syn
