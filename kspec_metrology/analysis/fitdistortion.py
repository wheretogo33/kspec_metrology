import numpy as np
from scipy.optimize import curve_fit
from kspec_metrology.analysis.utils import transform_polynomial, transform, camera2focal_coeff

def fitdistortion(x, y, fid_flag
                  , xobs, yobs
                  , xorigin, yorigin
                  , imatch, theta_guess):

    xobs_match = xobs[imatch]
    yobs_match = yobs[imatch]

    xfit = np.concatenate((xobs_match[fid_flag], xobs_match[fid_flag]))
    yfit = np.concatenate((yobs_match[fid_flag], yobs_match[fid_flag]))
    ccd_fit = np.concatenate((x[fid_flag], y[fid_flag]))
    
    coeff_temp = np.copy(camera2focal_coeff)
    coeff_temp[1] = theta_guess
    inv_popt_obs, _ = curve_fit(transform_polynomial, (xfit, yfit), ccd_fit
                       , maxfev=10000
                       , p0=coeff_temp
                       )

    xfocal_obs, yfocal_obs = transform(xobs_match, yobs_match, inv_popt_obs)

    dx = xfocal_obs-x
    dy = yfocal_obs-y

    theta_true, phi_true = np.full(x.size, 1e99), np.full(x.size, 1e99)
    theta_obs, phi_obs = np.full(x.size, 1e99), np.full(x.size, 1e99)
    for i in range(x.size):
      if fid_flag[i]:
        continue
      theta_true[i], phi_true[i] = find_angle_double_method2(xorigin[i], yorigin[i],
                                                             x[i], y[i])
      theta_obs[i], phi_obs[i] = find_angle_double_method2(xorigin[i], yorigin[i],
                                                           xfocal_obs[i], yfocal_obs[i] )
    theta_new, phi_new = (theta_true + theta_true - theta_obs), (phi_true + phi_true, phi_obs)
    print(theta_new[~fid_flag], phi_new[~fid_flag])
                    
    print(dx*1e3, dy*1e3)
                    
    return xfocal_obs, yfocal_obs, dx, dy, theta_true, phi_true, theta_new, phi_new, inv_popt_obs
