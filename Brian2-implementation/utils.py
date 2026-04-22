# =============================================================================
# Robust Randomization for Reproducibility
# =============================================================================

import numpy as np

def get_independent_draws(ids, type_id, global_seed, num_draws=1, dist='uniform', **kwargs):
    """
    Replicates NEURON's Random123(id1, id2, seed) logic for Brian 2 arrays.
    
    Args:
        ids: An iterable or 1D NumPy array of specific integer IDs (e.g., GIDs or indices).
        type_id: Unique integer for the cell/synapse population (NEURON's id2).
        global_seed: The global trial seed.
        num_draws: How many sequential values to draw per item.
        dist: 'uniform' or 'normal'.
    
    Returns:
        A 2D NumPy array of shape (num_draws, len(ids)).
    """
    ids = np.asarray(ids, dtype=int)
    results = np.zeros((num_draws, len(ids)))
    
    for idx, unique_id in enumerate(ids):
        # Hash the exact 3-tuple used in NEURON
        seq = np.random.SeedSequence((int(unique_id), int(type_id), int(global_seed)))
        
        # Philox is the exact core engine powering NEURON's Random123
        rng = np.random.Generator(np.random.Philox(seq))
        
        if dist == 'uniform':
            draws = rng.uniform(kwargs.get('low', 0.0), kwargs.get('high', 1.0), size=num_draws)
        elif dist == 'normal':
            draws = rng.normal(kwargs.get('loc', 0.0), kwargs.get('scale', 1.0), size=num_draws)
            
        results[:, idx] = draws
        
    return results
