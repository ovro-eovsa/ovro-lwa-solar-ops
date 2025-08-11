import sys
sys.path.append('/data07/msurajit/ovro-lwa-solar-ops')
import caltable_generator
import glob,os
'''
This code is an example on how to solve the crosshand phase using OVRO-LWA data.
We shall use CygA transit for this purpose. We will use data in the time_interval
(t_0,t_1). t_0 is the time after which the elevation of CygA exceeds 60 degrees
t_1 is the time after which the elevation of CygA decreases below 60 degrees. The
60 degree is a heurestic which I have seen works well. It balances the region of
the sky where the beam model is fairly well known and also the beam induced
polarization is sufficient to detect CygA in all Stokes.

starttime: t_0 This should be in Astropy isot format
endtime: t_1   This should also be in Astropy isot format
tdt: 900  This is the cadence within which data will be used
            in the interval t_0 and t_1. This is the time in 
            seconds. I advise using 900s to facilitate some
            averaging in time. 3600 is good for testing purposes
Note that by we use the beam_caltables for the initial bandpass calibration.
This is to ensure that the Stokes V bilobe is not there. However, even if 
regular caltables are provided, the crosshand phase solver explicitly will
flag the outrigger antennas. 

do_apply: False (by default) This ensures that if the user is not confident
                             then the phase solutions are not applied. However
                             it should be noted that while the solutions are plotted,
                             the solutions themselves are not saved, if do_apply=False
                             So, if the user is not sure, it is better to keep a copy
                             of the caltables, and then set do_apply=True. If the
                             user does not like the solutions, they should delete the
                             corresponding key from the database, and also delete the
                             corrupted caltables.
bands: all bands (by default). 
database: Default is None. Not required if do_apply is False. However, if do_apply=True,
            the user must specify a database. If not present, the database is created. 
            Database is a hdf5 file. The solutions will be written in a new group, where
            the group key is the YYYYMMDD of the starttime
By default the solutions will be plotted. The plotting always happens in a non-GUI environment
The figures will have the following 2 names by default : DI_leakage_variation.png
and crosshand_theta_variation.png . 
The database will not be overwritten by default. 
'''
starttime='2025-08-11T04:00:00'
endtime='2025-08-11T09:30:00'
tdt=3600 
cygA_coord='19h59m28.35663s +40d44m02.0970s'
caltable_folder='/fast/msurajit/slow_calib/realtime_code_output/investigate_caltable_issue/ovrolwasolar_code/caltables/'
beam_caltable_folder='/fast/msurajit/slow_calib/realtime_code_output/investigate_caltable_issue/ovrolwasolar_code/caltables_beam/'

caltable_prefix='20250804_071500'
doapply=True
bands=['32MHz','50MHz','69MHz']
database='test_crosshand_db.hdf5'

caltable_generator.crosshand_phase_solver(starttime,endtime,tdt,cygA_coord,\
                                        caltable_folder=caltable_folder,\
                                        beam_caltable_folder=beam_caltable_folder,\
                                        caltable_prefix=caltable_prefix,\
                                        doapply=True,\
                                        bands=bands,\
                                        database=database)

