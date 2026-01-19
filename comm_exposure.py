from kspec_metrology.exposure.mtlexp import mtlexp


mtlexp(
        exptime=2e3, readmode=2, nexposure=1,
        data_dir='./data/', head='Fiducial_short_'
        )

#mtlexp(
#        exptime=2e4, readmode=2, nexposure=1,
#        data_dir='./data/', head='Fiducial_long_'
#        )

mtlexp(
        exptime=5e4, readmode=2, nexposure=1, 
        data_dir='./data/', head='Fiber_trial1_'
        )


