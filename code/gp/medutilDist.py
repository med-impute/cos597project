import numpy as np
# testing
from sklearn.utils import array2d, check_random_state
from sklearn.metrics.pairwise import manhattan_distances

def med_cross_distances_test(X):
	"""
	Computes the nonzero componentwise L1 cross-distances between the vectors
	in X.

	Parameters
	----------

	X: array_like
	    An array with shape (n_samples, n_features)

	Returns
	-------

	D: array with shape (n_samples * (n_samples - 1) / 2, n_features)
	    The array of componentwise L1 cross-distances.

	ij: arrays with shape (n_samples * (n_samples - 1) / 2, 2)
	    The indices i and j of the vectors in X associated to the cross-
	    distances in D: D[k] = np.abs(X[ij[k, 0]] - Y[ij[k, 1]]).
	"""
	X = array2d(X)
	n_samples, n_features = X.shape
	n_nonzero_cross_dist = n_samples * (n_samples - 1) / 2
	ij = np.zeros((n_nonzero_cross_dist, 2), dtype=np.int)
	D = np.zeros((n_nonzero_cross_dist, n_features))
	ll_1 = 0
	for k in range(n_samples - 1):
	    ll_0 = ll_1
	    ll_1 = ll_0 + n_samples - k - 1
	    ij[ll_0:ll_1, 0] = k
	    ij[ll_0:ll_1, 1] = np.arange(k + 1, n_samples)
	    D[ll_0:ll_1] = np.abs(X[k] - X[(k + 1):n_samples])

	return D, ij.astype(np.int)

def med_distance_test(X, Y):
	print 'medutilDist: testing cross distances function (equivalent to manhattan_distances)'
	D = manhattan_distances(X, Y, sum_over_features=False)

	return D

def med_cross_distnaces(X):
    """
    Computes the nonzero componentwise distances for heterogenous medical data X.
    Parameters
    ----------
    X: array_like
        An array with shape (n_samples, n_features)
        n_features: (V + 3), where V is the intended # of features for cross-individual distance
        0: PAN, ID of the observation
        1: Feature type/index of the observation
        2: Time stamp of the observation
        3 - (V + 2): Cross-individual features
    Returns
    -------
    D: array with shape (n_samples * (n_samples - 1) / 2, n_features)
        The array of componentwise L1 cross-distances.
    ij: arrays with shape (n_samples * (n_samples - 1) / 2, 2)
        The indices i and j of the vectors in X associated to the cross-
        distances in D: D[k] = np.abs(X[ij[k, 0]] - Y[ij[k, 1]]).
    """
    # Now the distance is handled in the kernel, so we only wrap the data into correct format here
    X = check_array(X)
    n_samples, n_features = X.shape
    n_nonzero_cross_dist = n_samples * (n_samples - 1) // 2
    ij = np.zeros((n_nonzero_cross_dist, 2), dtype=np.int)
    D = np.zeros((n_nonzero_cross_dist, n_features))
    ll_1 = 0
    for k in range(n_samples - 1):
        ll_0 = ll_1
        ll_1 = ll_0 + n_samples - k - 1
        ij[ll_0:ll_1, 0] = k
        ij[ll_0:ll_1, 1] = np.arange(k + 1, n_samples)
        D[ll_0:ll_1] = np.concatenate([X[k], X[(k + 1):n_samples]])

    return D, ij

def med_distnaces(X):
	"""
	"""
	print 'medutilDist: custome distance fuction'
	D = manhattan_distances(X, Y, sum_over_features=False)
	
	return D

