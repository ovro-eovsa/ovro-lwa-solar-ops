import numpy as np
import matplotlib.pyplot as plt
from ovrolwasolar.beam_polcalib import image_polcal_astronomical_source as img_polcal
from ovrolwasolar import beam_polcalib,utils
import glob,os
import sys
sys.path.append('/data07/msurajit/ovro-lwa-solar-ops2/')
import caltable_generator

cygA_coord='19h59m28.35663s +40d44m02.0970s'
starttime='2025-04-17T13:00:00'
endtime='2025-04-17T13:10:00'
tdt=900  ### in seconds
caltable_folder='/lustre/solarpipe/realtime_pipeline/caltables_latest' ### ensure that there is one caltable for frequency band in this folder

dynamic_spectrum,freqs,tim_mjds=caltable_generator.get_source_DS(starttime,endtime,tdt,cygA_coord,\
                                    caltable_folder=caltable_folder)

img_pol=img_polcal(dynamic_spectrum=dynamic_spectrum,\
                    freqs=freqs,\
                    tim_mjds=tim_mjds,\
                    sky_coord=cygA_coord)

img_pol.fit_UV=False

img_pol.crosshand_phase_solver()

crosshand_V=img_pol.crosshand_theta ### obtained by matching V

np.savetxt("crosshand_phase.txt",np.array([img_pol.freqs,img_pol.crosshand_theta]).T)

plt.plot(img_pol.freqs,img_pol.crosshand_theta,'bs-',label='match V')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Crosshand phase (radian)")
plt.savefig("crosshand_theta_variation.png")
plt.close()
#### Depending on the crosshand phase variation obtained we can choose which is better. My tests show that these are practically same in
#### both cases. If very similar, I will prefer the phases obtained by matching both U and V


### This is a bonus plot, which shows the direction independent leakage

np.savetxt("DI_leakage.txt",np.array([img_pol.freqs,img_pol.leakage[0,:],\
                                    img_pol.leakage[1,:],img_pol.leakage[2,:],\
                                    img_pol.leakage[3,:]]).T)
fig,ax=plt.subplots(nrows=1,ncols=3,figsize=[12,4],constrained_layout=True)
for j,(ax1,stokes) in enumerate(zip(ax,['Q','U','V'])):
    ax1.plot(freqs,img_pol.leakage[j+1,:],'o-')
    ax1.set_ylabel(stokes+" leakage fraction")
    ax1.set_ylim([-0.2,0.2])
fig.supxlabel("Frequency (MHz)")
plt.savefig("DI_leakage_variation.png")
plt.close()
