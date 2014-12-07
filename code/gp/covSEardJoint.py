import numpy as np
from scipy.spatial.distance import squareform

def covSEardJoint(theta, d):
    """
    by-pass test: Squared exponential correlation model
    """
    print 'covSEardJoint: by-pass mode of custom kernel (Squared Exponential)'
    theta = np.asarray(theta, dtype=np.float)
    d = np.asarray(d, dtype=np.float)

    if d.ndim > 1:
        n_features = d.shape[1]
    else:
        n_features = 1

    if theta.size == 1:
        return np.exp(-theta[0] * np.sum(d ** 2, axis=1))
    elif theta.size != n_features:
        raise ValueError("Length of theta must be 1 or %s" % n_features)
    else:
        return np.exp(-np.sum(theta.reshape(1, n_features) * d ** 2, axis=1))
    """
    Self-defined Squared exponential correlation model (Radial Basis Function).
    
    Parameters
    ----------
    theta : array_like
        An array with shape 1 (isotropic) or n (anisotropic) giving the
        autocorrelation parameter(s).
        theta = [ 
        log(lambda_ii)  // 1 : class-wise individual length scale (1 unique coef. so far)
        log(theta_1)    // 1 + (1 to D): signal variance for each covariate
        log(theta_2)  
        ...
        log(theta_D)    
        log(lambda_1)   // 1 + D + (1 to D): temporal length scale for each covariate
        log(lambda_2)  
        ...
        log(lambda_D)    
        log(lambda_12)   // 1 + 2*D + (1 to 1/2*D*(D-1)):
        log(lambda_13)   // pairwise length scale for each pair of covariate
        ...
        log(lambda_D-1D) // Total # of hyperparameter now: 1 + 2*D + 1/2*D*(D-1)
        ]
    dx : array_like
        An array with shape (n_eval, n_features) giving the componentwise
        distances between locations x and x' at which the correlation model
        should be evaluated.
        
        *** 
        NOTE: Original GP library uses manhattan_distances to calculate the distance,
        and then call the correlation function to calculate the values.
        
        We need a custom distance function also, so will need to overwrite the GP library also.
        ***
        
        n_eval: total number of points (across individual, covariate, and time)
        n_features = 2*(V + 3), where
        0: PAN (ID) of training observation
        1: Feature index of training observation
        2: Elapsed timestamp of training observation (unit: accumulated hours)
        3 - (V+2): V-dim demographic data of training observation (*NOT YET CODED*)
        (V+3) - 2*(V+3)-1: for testing observation
        
    Returns
    -------
    r : array_like
        An array with shape (n_eval, ) containing the values of the
        autocorrelation model.

    #     if d.ndim > 1:
    #         n_features = d.shape[1]
    #     else:
    #         n_features = 1
    
    #     if theta.size == 1:
    #         return np.exp(-theta[0] * np.sum(d ** 2, axis=1))
    #     elif theta.size != n_features:
    #         raise ValueError("Length of theta must be 1 or %s" % n_features)
    #     else:
    #         return np.exp(-np.sum(theta.reshape(1, n_features) * d ** 2, axis=1))
    ## hard-coded variables (optimized later)
    train_fnum = 3
    test_fnum = train_fnum
    
    ## Parse Data
    print('covSEardJoint: parsing data...')
    
    theta = np.asarray(theta, dtype = np.float)
    d = np.asarray(d, dtype = np.float)
    
    N, V = d.shape
    V = V/2
    V = V - 3
    print(['covSEardJoint: current # of dimension for inter-individual kernel: V = ', str(V)])
    
    train_pid = d[:, 0]
    train_fid = d[:, 2]
    train_time = d[:, 4]
    train_demogr = d[:, 6:(V + 5)]
        
    test_pid = d[:, 1]
    test_fid = d[:, 3]
    test_time = d[:, 5]
    test_demogr = d[:, (V + 6):(2*V + 5)]
    
    ## Parse hyperparameter
    print('covSEardJoint: parsing hyperparam...')
    lambda_ii = np.exp(hyp[0]) # length scale for inter-individual term (not yet finished)
    theta = np.exp(2*hyp[1:(train_fnum)]) # scale factor (variance) for each covariate
    lambda_j = np.exp(hyp[(train_fnum + 1):(2*train_fnum)]) 
    lambda_jj = np.exp(squareform(hyp[(2*train_fnum + 1):(train_fnum*(train_fnum - 1)/2 + 2*train_fnum)])
    
    # debugging print
    print lambda_ii
    print theta
    print lambda_j
    print lambda_jj
    
    print('covSEardJoint: computing autocorrelations...')
    Kc = 0.5*np.ones(N)
    Kj = np.zeros(N)
    Kjj = np.zeros(N)
    Kii = np.zeros(N)
    for obs in range(0, N):
        # if (train_pid == test_pid) && (train_fid == test_fid) && (train_time == test_time):
        # diag. case, log(kernel val.) = 0
        if(train_pid[obs] == test_pid[obs]):
            # intra-individual case
            if(train_fid[obs] == test_fid[obs]):
                # intra-covariate case
                if(train_time[obs] == test_time[obs]):
                    print(['dulplicate time stamps for patient: ', str(train_pid[obs])])
                else:
                    # TODO: add periodic kernel or causal kernel here
                    Kj[obs] = ((train_time[obs] - test_time[obs])/lambda_j[train_fid[obs]]) ** 2
            else:
                # inter-covariate case
                if(train_time[obs] == test_time[obs]):
                    # need to find out correlation for inter-covariate cases (hard-coded as 0.4 now)
                    Kjj[obs] = ((0.4)/lambda_jj[train_fid[obs], test_fid[obs]]) ** 2
                else:
                    # need to find out correlation for inter-covariate cases (hard-coded as 0.4 now)
                    # may need to re-define joint time length scale
                    # (use product of two lambda_j instead now)
                    Kj[obs] = ((train_time[obs] - test_time[obs]) ** 2)/lambda_j[train_fid[obs]]/lambda_j[test_fid[obs]]
                    Kjj[obs] = ((0.4)/lambda_jj[train_fid[obs], test_fid[j]]) ** 2
        else:
            # inter-individual case
            Kj[obs] = 10 ** 20
            Kjj[obs] = 10 ** 20
            if(train_fid[obs] == test_fid[obs]):
                # intra-covariate case: need to re-define distance
                Kii[obs] = 10 # (demogr_dist[obs]/lambda_ii) ** 2
            else:
                # inter-covariate case: need to re-define distance
                Kii[obs] = 10 ** 2 # make it a near-zero number now
    
    Kvj = np.exp(-Kj/2);
    Kvjj = np.exp(-Kjj/2);
    Kvii = np.exp(-Kii/2);
    K = Kc * (Kvj * Kvjj + Kvii)
    return K
    """