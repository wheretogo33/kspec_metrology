import numpy as np
from kspec_metrology.analysis.findpeak import findpeak
from kspec_metrology.analysis.matchfiber import matchfiber
from kspec_metrology.analysis.fitdistortion import fitdistortion
from astropy.io import ascii
import json

def mtlcal(data_dir='./MTL/data/'):
    tab = ascii.read('Fiber_Configuration_250415.txt')
    tab.rename_columns(tab.colnames[:4], ["ID", "X", "Y", "FiducialFlag"])
    
    ids = [
           'B5' , 'C2', 'D0',
           'E5' , 'E8', 'G4',
           'G11', 'H7', 'I2',
           'I9' , 'K3', 'K6',
           'L10', 'M7'
          ]
    
    
    fid_ids = tab["ID"][tab["FiducialFlag"] == 1]
    pick_ids = list(dict.fromkeys(list(ids)) + list(fid_ids))
    
    #pick_ids.remove('Z0')
    
    xy_map = {r["ID"]: (r["X"], r["Y"]) for r in tab}
    xy = np.array([xy_map[i] for i in pick_ids], dtype=float)
    xorigin, yorigin = xy[:, 0], xy[:, 1]
    
    x, y = np.copy(xorigin), np.copy(yorigin)
    fid_flag = np.zeros(x.size, dtype=bool)
    fid_flag[14:] = True
    
    #---------------------------------------------------------------
    x[~fid_flag] += 16.8
    
    #---------------------------------------------------------------
    with open('./target/object.info', 'r') as ff:
        target = json.load(ff)
        xall, yall = target['x'], target['y']
    x[:14], y[:14] = xall[:14], yall[:14]


    #---------------------------------------------------------------
    npeaks = x.size
    _, xobs, yobs = findpeak(npeaks)

    imatch, theta_guess = matchfiber(x, y, xobs, yobs)

    _, _, dx, dy, _, _, theta_new, phi_new, _ = fitdistortion(x, y, fid_flag, xobs, yobs, xorigin, yorigin, imatch, theta_guess)

    return dx, dy, theta_new, phi_new 
