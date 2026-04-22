"""
synapses.py — Brian2 Synaptic Connections
"""

import numpy as np
import brian2 as b2
import config
from netcon import netcon_nearest_neighbors
from utils import get_independent_draws

# =============================================================================
# SHARED SYNAPSE EQUATION TEMPLATES 
# =============================================================================

_AMPA_STP_EQNS = '''
normalizing_factor = clip((2*radius + (1 - int(remove_recurrent_bool)))
                          / (Npost / Npre), 0, Npre)             : 1
iAMPA_post = -(gAMPA / normalizing_factor) * res_AMPA * sAMPA
             * (v_post - EAMPA)                                  : amp/meter**2 (summed)
dsAMPA/dt  = alpha_AMPA / (1 + exp(-(v_pre - 20*mV) / (2*mV)))
             - sAMPA / tauAMPA                                   : 1 (clock-driven)
dres_AMPA/dt = (1 - res_AMPA) / tauRes_AMPA
               + spike_activity * (-(1 - res_AMPA) / tauRes_AMPA
               + (deprFactor * res_AMPA - res_AMPA) / dt)        : 1 (clock-driven)
spike_activity : 1
Npre           : 1
Npost          : 1
alpha_AMPA     : 1/second
tauAMPA        : second
deprFactor     : 1
tauRes_AMPA    : second
gAMPA          : siemens/meter**2
EAMPA          : volt
radius         : 1
remove_recurrent_bool : boolean
'''

_NMDA_STP_EQNS = '''
normalizing_factor = clip((2*radius + (1 - int(remove_recurrent_bool)))
                          / (Npost / Npre), 0, Npre)             : 1
iNMDA_post = -(gNMDA / normalizing_factor) * res_NMDA * sNMDA
             * (v_post - ENMDA)                                  : amp/meter**2 (summed)
dsNMDA/dt  = alphaS_NMDA * xNMDA * (1 - sNMDA) - sNMDA / tauS_NMDA : 1 (clock-driven)
dxNMDA/dt  = alphaX_NMDA / (1 + exp(-(v_pre - 20*mV) / (2*mV)))
             - xNMDA / tauX_NMDA                                 : 1 (clock-driven)
dres_NMDA/dt = (1 - res_NMDA) / tauRes_NMDA
               + spike_activity * (-(1 - res_NMDA) / tauRes_NMDA
               + (deprFactor * res_NMDA - res_NMDA) / dt)        : 1 (clock-driven)
spike_activity : 1
Npre           : 1
Npost          : 1
alphaS_NMDA    : 1/second
tauS_NMDA      : second
alphaX_NMDA    : 1/second
tauX_NMDA      : second
deprFactor     : 1
tauRes_NMDA    : second
gNMDA          : siemens/meter**2
ENMDA          : volt
radius         : 1
remove_recurrent_bool : boolean
'''

_GABAA_STP_EQNS = '''
normalizing_factor = clip((2*radius + (1 - int(remove_recurrent_bool)))
                          / (Npost / Npre), 0, Npre)             : 1
iGABAA_post = -propoCondMult * gGABAA / normalizing_factor
              * res_GABAA * sGABAA * (v_post - EGABAA)           : amp/meter**2 (summed)
dsGABAA/dt  = alpha_GABAA / (1 + exp(-(v_pre - 20*mV) / (2*mV)))
              - sGABAA / (tauGABAA * propoTauMult)               : 1 (clock-driven)
dres_GABAA/dt = (1 - res_GABAA) / tauRes_GABAA
                + spike_activity * (-(1 - res_GABAA) / tauRes_GABAA
                + (deprFactor * res_GABAA - res_GABAA) / dt)     : 1 (clock-driven)
spike_activity : 1
Npre           : 1
Npost          : 1
alpha_GABAA    : 1/second
tauGABAA       : second
deprFactor     : 1
tauRes_GABAA   : second
gGABAA         : siemens/meter**2
EGABAA         : volt
radius         : 1
remove_recurrent_bool : boolean
propoCondMult  : 1
propoTauMult   : 1
'''

_STP_ON_PRE = {'up': 'spike_activity = 1', 'down': 'spike_activity = 0'}
_STP_DELAY = {'up': config.transmission_delay, 'down': config.transmission_delay + config.release_time}

def _nn_connect(syn, n_pre, n_post, n_neighbors, remove_recurrent):
    mat = netcon_nearest_neighbors(n_neighbors * 2, n_pre, n_post, remove_recurrent)
    src, tgt = mat.nonzero()
    syn.connect(i=src, j=tgt)
    return src, tgt

def _all_to_all_connect(syn):
    syn.connect()

# =============================================================================
# COMPARTMENTAL CONNECTIONS
# =============================================================================

def make_ICOM(PYdr_group, PYso_group):
    N = config.Npyso
    random_factors = get_independent_draws(np.arange(N), type_id=5, global_seed=config.randSeed, dist='normal')[0]

    ICOM_PYdr_PYso = b2.Synapses(PYdr_group, PYso_group, model='iCOM_post = -gCOM * (v_post - v_pre) : amp/meter**2 (summed)\ngCOM : siemens/meter**2', name='ICOM_PYdr_PYso')
    ICOM_PYdr_PYso.connect(j='i')
    ICOM_PYdr_PYso.gCOM = (11.667 + 0.667 * random_factors) * b2.msiemens / b2.cm**2

    ICOM_PYso_PYdr = b2.Synapses(PYso_group, PYdr_group, model='iCOM_post = -gCOM * (v_post - v_pre) : amp/meter**2 (summed)\ngCOM : siemens/meter**2', name='ICOM_PYso_PYdr')
    ICOM_PYso_PYdr.connect(j='i')
    ICOM_PYso_PYdr.gCOM = (5.0 + 0.286 * random_factors) * b2.msiemens / b2.cm**2

    return ICOM_PYdr_PYso, ICOM_PYso_PYdr

def make_IKNa(PYdr_group, PYso_group):
    N = config.Npyso
    gKNa_val = config.current_state_params['gKNa'] * b2.siemens / b2.cm**2
    EKNa = -100 * b2.mV
    alphaNa = 0.01 * 1000 * b2.mM / (b2.uA * b2.ms)
    RPump = 0.018 * b2.mM / b2.ms
    eqNa = 9.5 * b2.mM
    eqNaPumpTerm = eqNa**3 / (eqNa**3 + (15 * b2.mM)**3)
    gNa_local, ENa_local, phi_Na = 50 * b2.msiemens / b2.cm**2, 55 * b2.mV, 4
    gNaP_local, ENaP_local = 0.0686 * b2.msiemens / b2.cm**2, 55 * b2.mV
    areaSO, areaDR = 0.00015 * b2.cm**2, 0.00035 * b2.cm**2

    IKNa_PYso_PYdr = b2.Synapses(
        PYdr_group, PYso_group,
        model='''
        alphaM_Na_local  = (0.1  * (v_post/mV + 33)  / (1 - exp(-(v_post/mV + 33)  / 10))) / ms : Hz
        betaM_Na_local   = (4    * exp(-(v_post/mV + 53.7) / 12)) / ms                            : Hz
        alphaH_Na_local  = (0.07 * exp(-(v_post/mV + 50)   / 10)) / ms                            : Hz
        betaH_Na_local   = (1 / (1 + exp(-(v_post/mV + 20) / 10))) / ms                           : Hz
        Minf_Na_local    = alphaM_Na_local / (alphaM_Na_local + betaM_Na_local)                   : 1
        dhNa_local/dt    = phi_Na * (alphaH_Na_local * (1 - hNa_local) - betaH_Na_local * hNa_local) : 1 (clock-driven)
        iNa_local        = gNa_local * Minf_Na_local**3 * hNa_local * (v_post - ENa_local)        : amp/meter**2
        Minf_NaP_local   = 1 / (1 + exp(-(v_pre/mV + 55.7) / 7.7))                                : 1
        iNaP_local       = gNaP_local * Minf_NaP_local**3 * (v_pre - ENaP_local)                  : amp/meter**2
        dconcNa/dt = -alphaNa * (areaSO * iNa_local + areaDR * iNaP_local) - RPump * (concNa**3 / (concNa**3 + (15*mM)**3) - eqNaPumpTerm) : mM (clock-driven)
        iKNa_post = -gKNa_val * (0.37 / (1 + (38.7*mM / concNa)**3.5)) * (v_post - EKNa) : amp/meter**2 (summed)
        ''',
        method='euler', namespace={'gKNa_val': gKNa_val, 'EKNa': EKNa, 'alphaNa': alphaNa, 'RPump': RPump, 'eqNaPumpTerm': eqNaPumpTerm, 'gNa_local': gNa_local, 'ENa_local': ENa_local, 'phi_Na': phi_Na, 'gNaP_local': gNaP_local, 'ENaP_local': ENaP_local, 'areaSO': areaSO, 'areaDR': areaDR}, name='IKNa_PYso_PYdr'
    )
    IKNa_PYso_PYdr.connect(j='i')
    
    draws = get_independent_draws(np.arange(N), type_id=6, global_seed=config.randSeed, num_draws=2)
    IKNa_PYso_PYdr.hNa_local = 0.5 + 0.1 * draws[0]
    IKNa_PYso_PYdr.concNa = (12.0 + 3.0 * draws[1]) * b2.mM
    return IKNa_PYso_PYdr

# =============================================================================
# INTRACORTICAL SYNAPSES
# =============================================================================

def make_IAMPA_PYso_PYdr(PYso_group, PYdr_group, params, n_neighbors=10, remove_recurrent=True):
    syn = b2.Synapses(PYso_group, PYdr_group, model=_AMPA_STP_EQNS.replace('iAMPA_post', 'iAMPA_PYso_PYdr_post'), on_pre=_STP_ON_PRE, delay=_STP_DELAY, method='euler', name='IAMPA_PYso_PYdr')
    src, _ = _nn_connect(syn, config.Npyso, config.Npydr, n_neighbors, remove_recurrent)
    syn.Npre = len(syn.source); syn.Npost = len(syn.target)
    syn.gAMPA, syn.EAMPA, syn.alpha_AMPA, syn.tauAMPA = params['gAMPA_PYso_PYdr'] * b2.siemens / b2.cm**2, 0 * b2.mV, 3.48 / b2.ms, 2 * b2.ms
    syn.deprFactor, syn.tauRes_AMPA, syn.radius, syn.remove_recurrent_bool = 0.9, 400 * b2.ms, n_neighbors, remove_recurrent
    
    draws = get_independent_draws(np.arange(len(src)), type_id=7, global_seed=config.randSeed, num_draws=2)
    syn.sAMPA = 0.0 + 0.1 * draws[0]
    syn.res_AMPA = 0.8 - 0.1 * draws[1]
    syn.spike_activity = 0
    return syn

def make_INMDA_PYso_PYdr(PYso_group, PYdr_group, params, n_neighbors=10, remove_recurrent=True):
    syn = b2.Synapses(PYso_group, PYdr_group, model=_NMDA_STP_EQNS.replace('iNMDA_post', 'iNMDA_PYso_PYdr_post'), on_pre=_STP_ON_PRE, delay=_STP_DELAY, method='euler', name='INMDA_PYso_PYdr')
    src, _ = _nn_connect(syn, config.Npyso, config.Npydr, n_neighbors, remove_recurrent)
    syn.Npre = len(syn.source); syn.Npost = len(syn.target)
    syn.gNMDA, syn.ENMDA, syn.alphaS_NMDA, syn.tauS_NMDA, syn.alphaX_NMDA, syn.tauX_NMDA = params['gNMDA_PYso_PYdr'] * b2.siemens / b2.cm**2, 0 * b2.mV, 0.5 / b2.ms, 100 * b2.ms, 3.48 / b2.ms, 2 * b2.ms
    syn.deprFactor, syn.tauRes_NMDA, syn.radius, syn.remove_recurrent_bool = 0.9, 400 * b2.ms, n_neighbors, remove_recurrent
    
    draws = get_independent_draws(np.arange(len(src)), type_id=8, global_seed=config.randSeed, num_draws=3)
    syn.sNMDA = 0.0 + 0.1 * draws[0]
    syn.xNMDA = 0.0 + 0.1 * draws[1]
    syn.res_NMDA = 0.8 - 0.1 * draws[2]
    syn.spike_activity = 0
    return syn

def make_IAMPA_PYso_IN(PYso_group, IN_group, params, n_neighbors=10, remove_recurrent=False):
    syn = b2.Synapses(PYso_group, IN_group, model=_AMPA_STP_EQNS.replace('iAMPA_post', 'iAMPA_PYso_IN_post'), on_pre=_STP_ON_PRE, delay=_STP_DELAY, method='euler', name='IAMPA_PYso_IN')
    src, _ = _nn_connect(syn, config.Npyso, config.Ninh, n_neighbors, remove_recurrent)
    syn.Npre = len(syn.source); syn.Npost = len(syn.target)
    syn.gAMPA, syn.EAMPA, syn.alpha_AMPA, syn.tauAMPA = params['gAMPA_PYso_IN'] * b2.siemens / b2.cm**2, 0 * b2.mV, 3.48 / b2.ms, 2 * b2.ms
    syn.deprFactor, syn.tauRes_AMPA, syn.radius, syn.remove_recurrent_bool = 0.9, 400 * b2.ms, n_neighbors, remove_recurrent
    
    draws = get_independent_draws(np.arange(len(src)), type_id=9, global_seed=config.randSeed, num_draws=2)
    syn.sAMPA = 0.0 + 0.1 * draws[0]
    syn.res_AMPA = 0.8 - 0.1 * draws[1]
    syn.spike_activity = 0
    return syn

def make_INMDA_PYso_IN(PYso_group, IN_group, params, n_neighbors=10, remove_recurrent=False):
    syn = b2.Synapses(PYso_group, IN_group, model=_NMDA_STP_EQNS.replace('iNMDA_post', 'iNMDA_PYso_IN_post'), on_pre=_STP_ON_PRE, delay=_STP_DELAY, method='euler', name='INMDA_PYso_IN')
    src, _ = _nn_connect(syn, config.Npyso, config.Ninh, n_neighbors, remove_recurrent)
    syn.Npre = len(syn.source); syn.Npost = len(syn.target)
    syn.gNMDA, syn.ENMDA, syn.alphaS_NMDA, syn.tauS_NMDA, syn.alphaX_NMDA, syn.tauX_NMDA = params['gNMDA_PYso_IN'] * b2.siemens / b2.cm**2, 0 * b2.mV, 0.5 / b2.ms, 100 * b2.ms, 3.48 / b2.ms, 2 * b2.ms
    syn.deprFactor, syn.tauRes_NMDA, syn.radius, syn.remove_recurrent_bool = 0.9, 400 * b2.ms, n_neighbors, remove_recurrent
    
    draws = get_independent_draws(np.arange(len(src)), type_id=10, global_seed=config.randSeed, num_draws=3)
    syn.sNMDA = 0.0 + 0.1 * draws[0]
    syn.xNMDA = 0.0 + 0.1 * draws[1]
    syn.res_NMDA = 0.8 - 0.1 * draws[2]
    syn.spike_activity = 0
    return syn

def make_IGABAA_IN_PYso(IN_group, PYso_group, params, n_neighbors=10, remove_recurrent=False):
    syn = b2.Synapses(IN_group, PYso_group, model=_GABAA_STP_EQNS.replace('iGABAA_post', 'iGABAA_IN_PYso_post'), on_pre=_STP_ON_PRE, delay=_STP_DELAY, method='euler', name='IGABAA_IN_PYso')
    src, _ = _nn_connect(syn, config.Ninh, config.Npyso, n_neighbors, remove_recurrent)
    syn.Npre = len(syn.source); syn.Npost = len(syn.target)
    syn.gGABAA, syn.EGABAA, syn.alpha_GABAA, syn.tauGABAA = params['gGABAA_IN_PYso'] * b2.siemens / b2.cm**2, -70 * b2.mV, 1 / b2.ms, 5 * b2.ms
    syn.deprFactor, syn.tauRes_GABAA, syn.radius, syn.remove_recurrent_bool = 0.9, 400 * b2.ms, n_neighbors, remove_recurrent
    syn.propoCondMult, syn.propoTauMult = params['propoCondMult'], params['propoTauMult']
    
    draws = get_independent_draws(np.arange(len(src)), type_id=11, global_seed=config.randSeed, num_draws=2)
    syn.sGABAA = 0.0 + 0.1 * draws[0]
    syn.res_GABAA = 0.8 - 0.1 * draws[1]
    syn.spike_activity = 0
    return syn

def make_IGABAA_IN_IN(IN_group, params, n_neighbors=10, remove_recurrent=True):
    syn = b2.Synapses(IN_group, IN_group, model=_GABAA_STP_EQNS.replace('iGABAA_post', 'iGABAA_IN_IN_post'), on_pre=_STP_ON_PRE, delay=_STP_DELAY, method='euler', name='IGABAA_IN_IN')
    src, _ = _nn_connect(syn, config.Ninh, config.Ninh, n_neighbors, remove_recurrent)
    syn.Npre = len(syn.source); syn.Npost = len(syn.target)
    syn.gGABAA, syn.EGABAA, syn.alpha_GABAA, syn.tauGABAA = 0.000825 * b2.msiemens / b2.cm**2, -70 * b2.mV, 1 / b2.ms, 5 * b2.ms
    syn.deprFactor, syn.tauRes_GABAA, syn.radius, syn.remove_recurrent_bool = 0.9, 400 * b2.ms, n_neighbors, remove_recurrent
    syn.propoCondMult, syn.propoTauMult = params['propoCondMult'], params['propoTauMult']
    
    draws = get_independent_draws(np.arange(len(src)), type_id=12, global_seed=config.randSeed, num_draws=2)
    syn.sGABAA = 0.0 + 0.1 * draws[0]
    syn.res_GABAA = 0.8 - 0.1 * draws[1]
    syn.spike_activity = 0
    return syn

# =============================================================================
# INTRATHALAMIC SYNAPSES
# =============================================================================

def make_IAMPA_TC_RE(TC_group, RE_group, params):
    eqns = '''
    iAMPA_TC_RE_post = -(gAMPA / Npre) * sAMPA * (v_post - EAMPA) : amp/meter**2 (summed)
    dsAMPA/dt = P_AMPA * (1 + tanh(v_pre / (4*mV))) * (1 - sAMPA) - sAMPA / tauAMPA : 1 (clock-driven)
    Npre   : 1
    Npost  : 1
    gAMPA  : siemens/meter**2
    EAMPA  : volt
    tauAMPA : second
    P_AMPA  : 1/second
    '''
    syn = b2.Synapses(TC_group, RE_group, model=eqns, method='euler', name='IAMPA_TC_RE')
    _all_to_all_connect(syn)
    syn.Npre = len(syn.source); syn.Npost = len(syn.target)
    syn.gAMPA, syn.EAMPA, syn.tauAMPA, syn.P_AMPA = 0.4 * b2.msiemens / b2.cm**2, 1 * b2.mV, 2 * b2.ms, 5 / b2.ms
    
    num_synapses = config.Ntc * config.Nre
    draws = get_independent_draws(np.arange(num_synapses), type_id=13, global_seed=config.randSeed, num_draws=1)
    syn.sAMPA = 0.0 + 0.1 * draws[0]
    return syn

def make_IGABAA_RE_RE(RE_group, params):
    eqns = '''
    iGABAA_RE_RE_post = -propoCondMult * gGABAA / Npre * sGABAA * (v_post - EGABAA) : amp/meter**2 (summed)
    dsGABAA/dt = P_GABAA * (1 + tanh(v_pre / (4*mV))) * (1 - sGABAA) - sGABAA / (tauGABAA * propoTauMult) : 1 (clock-driven)
    Npre  : 1
    Npost : 1
    gGABAA : siemens/meter**2
    EGABAA : volt
    tauGABAA : second
    P_GABAA  : 1/second
    propoCondMult : 1
    propoTauMult  : 1
    '''
    syn = b2.Synapses(RE_group, RE_group, model=eqns, method='euler', name='IGABAA_RE_RE')
    _all_to_all_connect(syn)
    syn.Npre = len(syn.source); syn.Npost = len(syn.target)
    syn.gGABAA, syn.EGABAA, syn.tauGABAA, syn.P_GABAA = 0.1 * b2.msiemens / b2.cm**2, -80 * b2.mV, 5 * b2.ms, 2 / b2.ms
    syn.propoCondMult, syn.propoTauMult = params['propoCondMult'], params['propoTauMult']
    
    num_synapses = config.Nre * config.Nre
    draws = get_independent_draws(np.arange(num_synapses), type_id=14, global_seed=config.randSeed, num_draws=1)
    syn.sGABAA = 0.0 + 0.1 * draws[0]
    return syn

def make_IGABAA_RE_TC(RE_group, TC_group, params):
    eqns = '''
    iGABAA_RE_TC_post = -propoCondMult * gGABAA / Npre * sGABAA * (v_post - EGABAA) : amp/meter**2 (summed)
    dsGABAA/dt = P_GABAA * (1 + tanh(v_pre / (4*mV))) * (1 - sGABAA) - sGABAA / (tauGABAA * propoTauMult) : 1 (clock-driven)
    Npre  : 1
    Npost : 1
    gGABAA : siemens/meter**2
    EGABAA : volt
    tauGABAA : second
    P_GABAA  : 1/second
    propoCondMult : 1
    propoTauMult  : 1
    '''
    syn = b2.Synapses(RE_group, TC_group, model=eqns, method='euler', name='IGABAA_RE_TC')
    _all_to_all_connect(syn)
    syn.Npre = len(syn.source); syn.Npost = len(syn.target)
    syn.gGABAA, syn.EGABAA, syn.tauGABAA, syn.P_GABAA = 0.1 * b2.msiemens / b2.cm**2, -80 * b2.mV, 5 * b2.ms, 2 / b2.ms
    syn.propoCondMult, syn.propoTauMult = params['propoCondMult'], params['propoTauMult']
    
    num_synapses = config.Nre * config.Ntc
    draws = get_independent_draws(np.arange(num_synapses), type_id=15, global_seed=config.randSeed, num_draws=1)
    syn.sGABAA = 0.0 + 0.1 * draws[0]
    return syn

def make_IGABAB_RE_TC(RE_group, TC_group, params):
    eqns = '''
    iGABAB_RE_TC_post = -gGABAB / Npre * (sGABAB**4 / (sGABAB**4 + 100)) * (v_post - EGABAB) : amp/meter**2 (summed)
    drGABAB/dt = K1 * (2 * (1 + tanh(v_pre / (4*mV)))) * (1 - rGABAB) - K2 * rGABAB : 1 (clock-driven)
    dsGABAB/dt = K3 * rGABAB - K4 * sGABAB : 1 (clock-driven)
    Npre   : 1
    Npost  : 1
    gGABAB : siemens/meter**2
    EGABAB : volt
    K1     : 1/second
    K2     : 1/second
    K3     : 1/second
    K4     : 1/second
    '''
    syn = b2.Synapses(RE_group, TC_group, model=eqns, method='euler', name='IGABAB_RE_TC')
    _all_to_all_connect(syn)
    syn.Npre = len(syn.source); syn.Npost = len(syn.target)
    syn.gGABAB, syn.EGABAB, syn.K1, syn.K2, syn.K3, syn.K4 = 0.001 * b2.msiemens / b2.cm**2, -95 * b2.mV, 0.5 / b2.ms, 0.0012 / b2.ms, 0.18 / b2.ms, 0.034 / b2.ms
    
    num_synapses = config.Nre * config.Ntc
    draws = get_independent_draws(np.arange(num_synapses), type_id=16, global_seed=config.randSeed, num_draws=2)
    syn.rGABAB = 0.0 + 0.1 * draws[0]
    syn.sGABAB = 0.0 + 0.1 * draws[1]
    return syn

# =============================================================================
# THALAMOCORTICAL & CORTICOTHALAMIC SYNAPSES
# =============================================================================

def make_IAMPA_TC_PYdr(TC_group, PYdr_group, params, n_neighbors=10, remove_recurrent=False):
    eqns = '''
    normalizing_factor = clip((2*radius + (1 - int(remove_recurrent_bool))) / (Npost / Npre), 0, Npre) : 1
    iAMPA_TC_PYdr_post = -(gAMPA / normalizing_factor) * sAMPA * (v_post - EAMPA) : amp/meter**2 (summed)
    dsAMPA/dt = P_AMPA * (1 + tanh(v_pre / (4*mV))) * (1 - sAMPA) - sAMPA / tauAMPA : 1 (clock-driven)
    Npre   : 1
    Npost  : 1
    gAMPA  : siemens/meter**2
    EAMPA  : volt
    tauAMPA : second
    P_AMPA : 1/second
    radius : 1
    remove_recurrent_bool : boolean
    '''
    syn = b2.Synapses(TC_group, PYdr_group, model=eqns, method='euler', name='IAMPA_TC_PYdr')
    src, _ = _nn_connect(syn, config.Ntc, config.Npydr, n_neighbors, remove_recurrent)
    syn.Npre = len(syn.source); syn.Npost = len(syn.target)
    syn.gAMPA, syn.EAMPA, syn.tauAMPA, syn.P_AMPA, syn.radius, syn.remove_recurrent_bool = params['gAMPA_TC_PYdr'] * b2.siemens / b2.cm**2, 1 * b2.mV, 2 * b2.ms, 5 / b2.ms, n_neighbors, remove_recurrent
    
    draws = get_independent_draws(np.arange(len(src)), type_id=17, global_seed=config.randSeed, num_draws=1)
    syn.sAMPA = 0.0 + 0.1 * draws[0]
    return syn

def make_IAMPA_TC_IN(TC_group, IN_group, params, n_neighbors=10, remove_recurrent=False):
    eqns = '''
    normalizing_factor = clip((2*radius + (1 - int(remove_recurrent_bool))) / (Npost / Npre), 0, Npre) : 1
    iAMPA_TC_IN_post = -(gAMPA / normalizing_factor) * sAMPA * (v_post - EAMPA) : amp/meter**2 (summed)
    dsAMPA/dt = P_AMPA * (1 + tanh(v_pre / (4*mV))) * (1 - sAMPA) - sAMPA / tauAMPA : 1 (clock-driven)
    Npre   : 1
    Npost  : 1
    gAMPA  : siemens/meter**2
    EAMPA  : volt
    tauAMPA : second
    P_AMPA : 1/second
    radius : 1
    remove_recurrent_bool : boolean
    '''
    syn = b2.Synapses(TC_group, IN_group, model=eqns, method='euler', name='IAMPA_TC_IN')
    src, _ = _nn_connect(syn, config.Ntc, config.Ninh, n_neighbors, remove_recurrent)
    syn.Npre = len(syn.source); syn.Npost = len(syn.target)
    syn.gAMPA, syn.EAMPA, syn.tauAMPA, syn.P_AMPA, syn.radius, syn.remove_recurrent_bool = params['gAMPA_TC_IN'] * b2.siemens / b2.cm**2, 1 * b2.mV, 2 * b2.ms, 5 / b2.ms, n_neighbors, remove_recurrent
    
    draws = get_independent_draws(np.arange(len(src)), type_id=18, global_seed=config.randSeed, num_draws=1)
    syn.sAMPA = 0.0 + 0.1 * draws[0]
    return syn

def make_IAMPA_PYso_TC(PYso_group, TC_group, params, n_neighbors=10, remove_recurrent=False):
    eqns = '''
    normalizing_factor = clip((2*radius + (1 - int(remove_recurrent_bool))) / (Npost / Npre), 0, Npre) : 1
    iAMPA_PYso_TC_post = -(gAMPA / normalizing_factor) * sAMPA * (v_post - EAMPA) : amp/meter**2 (summed)
    dsAMPA/dt = P_AMPA * (1 + tanh(v_pre / (4*mV))) * (1 - sAMPA) - sAMPA / tauAMPA : 1 (clock-driven)
    Npre   : 1
    Npost  : 1
    gAMPA  : siemens/meter**2
    EAMPA  : volt
    tauAMPA : second
    P_AMPA : 1/second
    radius : 1
    remove_recurrent_bool : boolean
    '''
    syn = b2.Synapses(PYso_group, TC_group, model=eqns, method='euler', name='IAMPA_PYso_TC')
    src, _ = _nn_connect(syn, config.Npyso, config.Ntc, n_neighbors, remove_recurrent)
    syn.Npre = len(syn.source); syn.Npost = len(syn.target)
    syn.gAMPA, syn.EAMPA, syn.tauAMPA, syn.P_AMPA, syn.radius, syn.remove_recurrent_bool = params['gAMPA_PYso_TC'] * b2.siemens / b2.cm**2, 1 * b2.mV, 2 * b2.ms, 5 / b2.ms, n_neighbors, remove_recurrent
    
    draws = get_independent_draws(np.arange(len(src)), type_id=19, global_seed=config.randSeed, num_draws=1)
    syn.sAMPA = 0.0 + 0.1 * draws[0]
    return syn

def make_IAMPA_PYso_RE(PYso_group, RE_group, params, n_neighbors=10, remove_recurrent=False):
    eqns = '''
    normalizing_factor = clip((2*radius + (1 - int(remove_recurrent_bool))) / (Npost / Npre), 0, Npre) : 1
    iAMPA_PYso_RE_post = -(gAMPA_PYso_RE / normalizing_factor) * sAMPA_PYso_RE * (v_post - EAMPA_PYso_RE) : amp/meter**2 (summed)
    dsAMPA_PYso_RE/dt = P_AMPA * (1 + tanh(v_pre / (4*mV))) * (1 - sAMPA_PYso_RE) - sAMPA_PYso_RE / tauAMPA : 1 (clock-driven)
    Npre   : 1
    Npost  : 1
    gAMPA_PYso_RE   : siemens/meter**2
    EAMPA_PYso_RE   : volt
    tauAMPA : second
    P_AMPA : 1/second
    radius : 1
    remove_recurrent_bool : boolean
    '''
    syn = b2.Synapses(PYso_group, RE_group, model=eqns, method='euler', name='IAMPA_PYso_RE')
    src, _ = _nn_connect(syn, config.Npyso, config.Nre, n_neighbors, remove_recurrent)
    syn.Npre = len(syn.source); syn.Npost = len(syn.target)
    syn.gAMPA_PYso_RE, syn.EAMPA_PYso_RE, syn.tauAMPA, syn.P_AMPA, syn.radius, syn.remove_recurrent_bool = params['gAMPA_PYso_RE'] * b2.siemens / b2.cm**2, 1 * b2.mV, 2 * b2.ms, 5 / b2.ms, n_neighbors, remove_recurrent
    
    draws = get_independent_draws(np.arange(len(src)), type_id=20, global_seed=config.randSeed, num_draws=1)
    syn.sAMPA_PYso_RE = 0.0 + 0.1 * draws[0]
    return syn
