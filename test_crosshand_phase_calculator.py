import numpy as np
import matplotlib.pyplot as plt
from ovrolwasolar.beam_polcalib import image_polcal_astronomical_source as img_polcal
from ovrolwasolar import beam_polcalib,utils
import glob,os

dynamic_spectrum=np.load("avged_dynamic_spectrum.npy")
mjds=np.load("time_mjds.npy")
freqs=np.load("freqs_MHz.npy")/1e6

img_pol=img_polcal(dynamic_spectrum=dynamic_spectrum,\
                    freqs=freqs,\
                    tim_mjds=mjds,\
                    sky_coord='19h59m28.35663s +40d44m02.0970s')
                    
img_pol.crosshand_phase_solver()
'''
fig,ax=plt.subplots(nrows=1,ncols=3,figsize=[12,4],constrained_layout=True)
for j,(ax1,stokes) in enumerate(zip(ax,['Q','U','V'])):
    ax1.plot(freqs,img_pol.leakage[j+1,:],'o-')
    ax1.set_ylabel(stokes+" leakage fraction")
    ax1.set_ylim([-0.2,0.2])
fig.supxlabel("Frequency (MHz)")
fig.suptitle("April 2, 2025, CygA transit")
plt.show()


plt.plot(img_pol.freqs,img_pol.crosshand_theta,'o-')
plt.show()
'''
inplace=False
crosshand_phase=img_pol.crosshand_theta
crosshand_freqs=img_pol.freqs

caltables=glob.glob("/fast/msurajit/slow_calib/realtime_code_output/caltables/20250325_161009_*.bcal")
    
crosshand_applied=[None]*len(caltables)

for j,caltable in enumerate(caltables):
    if not inplace:
        copied_caltable=caltable+".crosshand"
        crosshand_applied[j]=copied_caltable
        os.system("cp -r "+caltable+" "+copied_caltable)
        beam_polcalib.combine_crosshand_theta_on_caltable(copied_caltable,crosshand_phase,crosshand_freqs)
