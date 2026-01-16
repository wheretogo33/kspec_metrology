import numpy as np
from kspec_metrology.analysis.findpeak import findpeak
from kspec_metrology.analysis.matchfiber import matchfiber
from kspec_metrology.analysis.fitdistortion import fitdistortion

def mtlcal(data_dir='./MTL/data/'):
    # Load Fiber positions and Fiducial flag------------------------------------------------
    #x, y, fid_flag = np.load(data_dir+"pos.npy")
    #fid_flag = fid_flag.astype(bool)
    
    
    
    npeaks = x.size

    _, xobs, yobs = findpeak(npeaks)

    imatch, theta_guess = matchfiber(x, y, xobs, yobs)

    _, _, dx, dy, _, _, theta_new, phi_new, _ = fitdistortion(x, y, fid_flag, xobs, yobs, xorigin, yorigin, imatch, theta_guess)

    return dx, dy, theta_new, phi_new 
