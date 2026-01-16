from kspec_metrology.exposure.mtlexp import mtlexp
from kspec_metrology.analysis.mtlcal import mtlcal

mtlexp(exptime=3000.,
        readmode=1,
        nexposure=1,
        data_dir='./MTL/data/')


dx, dy, theta_new, phi_new = mtlcal() 
