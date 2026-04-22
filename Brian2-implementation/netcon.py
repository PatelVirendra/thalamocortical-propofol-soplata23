"""
netcon.py — Nearest-Neighbour Connectivity Matrix
==================================================
Generates a binary (n_pre × n_post) connectivity matrix for
nearest-neighbour (ring-topology) synaptic projections.

Identical to the NEURON netcon.py so that connectivity patterns match
exactly between the two simulator implementations.

Author  : Erik Roberts (converted from MATLAB)
Modified: Austin Soplata (AES)
"""

import numpy as np

def netcon_nearest_neighbors(n_neighbors, n_pre, n_post, remove_recurrent_bool):
    """
    Build a nearest-neighbour connectivity matrix for one- or two-population
    projections.

    Each pre-synaptic neuron i connects to the n_neighbors post-synaptic
    neurons that are topographically closest to it.  The population is
    arranged on a 1-D ring (periodic boundary conditions).

    Parameters
    ----------
    n_neighbors : int
        Connective diameter — number of post-synaptic targets per pre cell.
        Always rounded down to the nearest even number so that the fan-out
        is symmetric (n_neighbors/2 on each side).
    n_pre : int
        Number of pre-synaptic neurons.
    n_post : int
        Number of post-synaptic neurons.
    remove_recurrent_bool : bool
        If True, remove the diagonal entries (autapses).  Only meaningful
        when n_pre == n_post.

    Returns
    -------
    netcon : np.ndarray, shape (n_pre, n_post), dtype float
        Binary connectivity matrix.  Entry [i, j] == 1 means neuron i
        synapses onto neuron j.
    """
    netcon = np.zeros((n_pre, n_post))

    # Ensure the diameter is even for symmetric connectivity
    n_neighbors = round(n_neighbors - (n_neighbors % 2))
    n_half      = n_neighbors // 2

    if n_pre > n_neighbors or n_post > n_neighbors:

        if n_pre == n_post:
            # --- Same-size populations: direct ring topology ---
            for i in range(n_pre):
                j = list(range(i - n_half, i + n_half + 1))
                j = [x + n_post if x < 0       else x for x in j]
                j = [x - n_post if x >= n_post  else x for x in j]
                netcon[i, j] = 1

        elif n_pre > n_post:
            # --- More pre cells: compress the pre index ---
            spacing = round(n_pre / n_post)
            for i in range(n_pre):
                centre = round(i / spacing)
                j = list(range(centre - n_half, centre + n_half + 1))
                j = [x + n_post if x < 0       else x for x in j]
                j = [x - n_post if x >= n_post  else x for x in j]
                netcon[i, j] = 1

        else:
            # --- More post cells: expand the pre index ---
            spacing = round(n_post / n_pre)
            for i in range(n_pre):
                centre = i * spacing
                j = list(range(centre - n_half, centre + n_half + 1))
                j = [x + n_post if x < 0       else x for x in j]
                j = [x - n_post if x >= n_post  else x for x in j]
                netcon[i, j] = 1

    else:
        # --- Populations are smaller than the diameter: all-to-all ---
        netcon = np.ones((n_pre, n_post))

    # Remove autapses / recurrent connections if requested
    if remove_recurrent_bool:
        netcon = netcon - np.eye(n_pre, n_post)

    return netcon
