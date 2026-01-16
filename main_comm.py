from kspec_metrology.exposure.mtlexp import mtlexp
from kspec_metrology.analysis.mtlcal import mtlcal
from astropy.io import ascii
import numpy as np

tab = ascii.read('Fiber_Configuration_250415.txt')
tab.rename_columns(tab.colnames[:4], ["ID", "X", "Y", "FiducialFlag"])

ids = ['B5' , 'C2', 'D0',
       'E5' , 'E8', 'G4',
       'G11', 'H7', 'I2',
       'I9' , 'K3', 'K6',
       'L10', 'M7']


fid_ids = tab["ID"][tab["FiducialFlag"] == 1]
pick_ids = list(dict.fromkeys(list(ids)) + list(fid_ids))

#pick_ids.remove('Z0')

xy_map = {r["ID"]: (r["X"], r["Y"]) for r in tab}
xy = np.array([xy_map[i] for i in pick_ids], dtype=float)
xorigin, yorigin = xy[:, 0], xy[:, 1]


mtlexp(exptime=3000.,
        readmode=1,
        nexposure=1,
        data_dir='./MTL/data/')


dx, dy, theta_new, phi_new = mtlcal() 