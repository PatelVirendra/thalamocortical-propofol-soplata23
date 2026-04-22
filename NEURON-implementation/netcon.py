import numpy as np

def netcon_nearest_neighbors(n_neighbors, n_pre, n_post, remove_recurrent_bool):
    """
    Calculate netcon for radius-type connections.

    Version 2
    Author: Erik Roberts (converted from MATLAB)
    Some modifications by Austin Soplata (AES)

    Purpose: Makes a connectivity matrix for 'nearest neighbors' connections,
             either within a single population or between two populations.
             This works for connections where the two populations are of
             either the same or different sizes.

    Args:
        n_neighbors: number of nearest neighbors to connect to, aka connective 'diameter'
        n_pre: number of presynaptic neurons
        n_post: number of postsynaptic neurons
        remove_recurrent_bool: Remove recurrent connections. Only meant for when
                              making a connection between a population and itself.

    Returns:
        netcon: the connection matrix
    """
    netcon = np.zeros((n_pre, n_post))

    # Make even
    n_neighbors = round(n_neighbors - (n_neighbors % 2))

    n_half = n_neighbors // 2

    if n_pre > n_neighbors or n_post > n_neighbors:
        if n_pre == n_post:
            # If height < width, then i needs to wrap around
            for i in range(n_pre):
                j = list(range(i - n_half, i + n_half + 1))
                # Handle wrapping for negative indices
                j = [x + n_post if x < 0 else x for x in j]
                # Handle wrapping for indices exceeding n_post
                j = [x - n_post if x >= n_post else x for x in j]
                netcon[i, j] = 1
        elif n_pre > n_post:
            spacing = round(n_pre / n_post)
            for i in range(n_pre):
                j = list(range(round(i / spacing) - n_half, round(i / spacing) + n_half + 1))
                # Handle wrapping for negative indices
                j = [x + n_post if x < 0 else x for x in j]
                # Handle wrapping for indices exceeding n_post
                j = [x - n_post if x >= n_post else x for x in j]
                netcon[i, j] = 1
        elif n_pre < n_post:
            spacing = round(n_post / n_pre)
            for i in range(n_pre):
                j = list(range(i * spacing - n_half, i * spacing + n_half + 1))
                # Handle wrapping for negative indices
                j = [x + n_post if x < 0 else x for x in j]
                # Handle wrapping for indices exceeding n_post
                j = [x - n_post if x >= n_post else x for x in j]
                netcon[i, j] = 1
    else:
        netcon = np.ones((n_pre, n_post))

    # Remove recurrent connections
    if remove_recurrent_bool:
        netcon = netcon - np.eye(n_pre, n_post)

    return netcon