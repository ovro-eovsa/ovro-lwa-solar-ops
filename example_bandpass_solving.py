import sys
sys.path.append('/data07/msurajit/ovro-lwa-solar-ops')
import caltable_generator
from astropy.time import Time
'''
This code is an example on how to solve the crosshand phase using OVRO-LWA data.
We shall use CygA transit for this purpose. We will use data in the time_interval
(t_0,t_1). t_0 is the time after which the elevation of CygA exceeds 80 degrees
t_1 is the time after which the elevation of CygA decreases below 80 degrees.

The high elevation is chosen to ensure that the beam model is very accurate, and
also the A-team only model which is used in the bandpass task is accurate. 

There is a function in caltable_generator, which can be used to find these times

time=Time('2025-08-11T10:00:00',format='isot')
rise, set1=caltable_generator.source_riseset(cygA_coord,time,altitude_limit=80)

starttime: t_0 This should be in Astropy isot format
endtime: t_1   This should also be in Astropy isot format
workdir: Path where all the data processing will be done. By default, this is current
         directory

By default caltables will be stored in workdir/caltables and workdir/caltables_beam. 
However these locations can also be changed. The code will also create a waterfall plot
which can be used to check the quality of the solutions.
'''



starttime='2025-08-11T07:45:00'
endtime='2025-08-11T07:50:00'

workdir='/fast/msurajit/slow_calib/realtime_code_output/investigate_caltable_issue/ovrolwasolar_code/'

bands=['55MHz']

bcaltbs,bcaltbs_bm=caltable_generator.gen_multitime_caltable(starttime,endtime, workdir=workdir, bands=bands)
