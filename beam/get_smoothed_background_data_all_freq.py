import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.signal import savgol_filter
import matplotlib.dates as mdates
from astropy.time import Time

background_file='OVRO_LWA_background_calibrated.hdf5'
outfile='OVRO_LWA_smoothed_beam_background.hdf5'

mjd_all=[]
background_all=[]

freq_ind=700

with h5py.File(background_file,'r') as hf:
    group_keys=hf.keys()
    for key in group_keys:
        if key=='freq_MHz':
            freqs=np.array(hf[key])
        else:
            hf_group=hf[key]
            background=np.array(hf_group['background'])
            mjd=np.array(hf_group['mjd'])
            
        
            background_all+=background.tolist()
            mjd_all+=mjd.tolist()
            

background_all=np.array(background_all)
mjd_all=np.array(mjd_all)

datetimes=Time(mjd_all,format='mjd').datetime

datetime_plot=mdates.date2num(datetimes)

background_smooth=np.zeros_like(background_all)*np.nan
shape=background_all.shape

num_chan=shape[1]
for i in range(num_chan):
    background_freq=background_all[:,i]
    #print (background_freq)
    pos=np.where((background_freq<0.9) & (background_freq>0.4))
    try:
        filtered_background=savgol_filter(background_freq[pos],window_length=20,polyorder=2)
       
        background_smooth[:,i]=np.interp(mjd_all,mjd_all[pos],filtered_background)

    except Exception as e:
        print (i,e)
        raise e


with h5py.File(outfile,'w') as hf:
    hf.create_dataset("background_SFU",data=background_smooth.T)
    hf.create_dataset("freq_MHz",data=freqs)
    hf.create_dataset("mjd",data=mjd_all)


fig=plt.figure()
ax=fig.add_subplot(111)
plt.imshow(background_smooth.T,vmax=1,vmin=0.4,cmap='viridis',origin='lower',extent=(datetime_plot[0],datetime_plot[-1],freqs[0],freqs[-1]),aspect='auto')
plt.xlabel("Date")
plt.ylabel("Frequency (MHz)")
cbar=plt.colorbar()
cbar.ax.set_ylabel("Flux density (SFU)",rotation=270,labelpad=15)
ax.xaxis_date()
date_format = mdates.DateFormatter('%Y-%m-%d')
ax.xaxis.set_major_formatter(date_format)


fig.autofmt_xdate()

plt.savefig("background_variation_with_time.pdf")
plt.show()

