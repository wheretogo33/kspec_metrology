import numpy as np
from scipy.spatial import cKDTree
from kspec_metrology.analysis.utils import transform, camera2focal_coeff
from kspec_metrology.logging.log import get_logger

def matchfiber(x, y
               , xobs_raw, yobs_raw
               , nbuffer=10):
    log = get_logger()

    coeff_temp = np.copy(camera2focal_coeff)
    xobs, yobs = transform(xobs_raw, yobs_raw, coeff_temp)

    nhunt = 720
    theta_grid = np.linspace(0., 2.*np.pi, nhunt)
    dd_sum = np.zeros(nhunt)

    log.info("Finding field rotation angle")
    for ihunt, theta_temp in enumerate(theta_grid):
        xobs_rot = np.cos(theta_temp)*xobs - np.sin(theta_temp)*yobs
        yobs_rot = np.sin(theta_temp)*xobs + np.cos(theta_temp)*yobs
        obs_flag = np.concatenate( (np.full(xobs.size, 0), np.full(xobs.size, 1)) )
        pos_tot = np.concatenate( (np.vstack((x, y)).T
                                 , np.vstack((xobs_rot, yobs_rot)).T) )

        tree = cKDTree(pos_tot)
        dd, ii = tree.query(pos_tot, k=nbuffer)

        for ipeak in range(xobs.size):
            dd_sum[ihunt] += dd[ipeak][obs_flag[ii[ipeak]] == 1].min()

    theta_guess = theta_grid[dd_sum.argmin()]

    xobs_rot = np.cos(theta_guess)*xobs - np.sin(theta_guess)*yobs
    yobs_rot = np.sin(theta_guess)*xobs + np.cos(theta_guess)*yobs

    ngrid = 41
    offset_grid = np.linspace(-10., 10., ngrid)
    dsum_temp = np.zeros((ngrid, ngrid))
    for i in range(ngrid):
        for j in range(ngrid):
            tree = cKDTree(np.c_[xobs_rot+offset_grid[i], yobs_rot+offset_grid[j]])
            d, _ = tree.query(np.c_[x, y], k=1)

            dsum_temp[i,j] = d.sum()

    imin, jmin = np.unravel_index(dsum_temp.argmin(), dsum_temp.shape)

    xobs_rot += offset_grid[imin]
    yobs_rot += offset_grid[jmin]
    #==========================================================================
                 
    pos_tot = np.concatenate( (np.vstack((x, y)).T
                                 , np.vstack((xobs_rot, yobs_rot)).T) )
    tree = cKDTree(pos_tot)
    dd, ii = tree.query(pos_tot, k=nbuffer)   

    imatch = np.zeros(xobs.size, dtype=np.int32)
    for ipeak in range(xobs.size):
        itemp = ii[ipeak]
        imatch[ipeak] = itemp[obs_flag[itemp] == 1][dd[ipeak][obs_flag[itemp] == 1].argmin()] - xobs.size

    if np.unique(imatch).size != xobs.size:
        log.warning("Fiber matching is not unique")
        # log about which number of fibers are not matched

    return imatch, theta_guess, (offset_grid[imin], offset_grid[jmin])
