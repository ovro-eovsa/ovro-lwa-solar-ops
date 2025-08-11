#from ovrolwasolar import beam_polcalib
from astropy.io import fits
import numpy as np
from astropy.time import Time
import astropy.table
import matplotlib.pyplot as plt
import os,h5py
from scipy.interpolate import RegularGridInterpolator
import datetime as dt
import sys
sys.path.insert(0,'/data1/smondal/softwares/lib/python3.8/site-packages/')
from ovrolwasolar import beam_polcalib,utils

def do_beam_correction(filename):
    with fits.open(filename,mode='update') as hdul:
        if 'PBCOR' not in hdul[0].header:
            flux=hdul[0].data
            freqs=hdul[1].data
            tim = astropy.table.Table(hdul[2].data)
            tmjd= np.array(tim['mjd']) + np.array(tim['time']) /86400000
            
            
            freq_MHz=np.array(freqs.tolist())*1e3
            tim_mjd=Time(tmjd,format='mjd')
            
            az,alt=utils.get_solar_altaz_multiple_times(tim_mjd)
            
            good_pos=np.where(alt>0.001)[0]
            bad_pos=np.where(alt<=0.001)[0]
            
            beam=beam_polcalib.beam_polcal(filename=filename)
            beam.model_beam_file='/common/lwa/support_files/OVRO-LWA_MROsoil_updatedheight.h5'
            
            beam.stokes_data=np.squeeze(flux,axis=0)
            beam.freqs=freq_MHz
            beam.times=tim_mjd[good_pos]
            
            beam.get_primary_beam(freq_sep=1,tim_sep=300,overwrite=True,normalise_wrt_I=False)
            
            
            num_freqs=beam.freqs.size
            freq_arr=np.arange(0,num_freqs,1)


            #beam.stokes_data[0][freq_arr][good_pos]/=beam.primary_beam[0,...]
            #beam.stokes_data[0][freq_arr][bad_pos]=np.nan
           
            for i in range(num_freqs):
                beam.stokes_data[0,i,good_pos]/=beam.primary_beam[0,i,:]
                beam.stokes_data[0,i,bad_pos]=-1000

            hdul[0].data=np.expand_dims(beam.stokes_data,axis=0)
            hdul[0].header['PBCOR']=True
            hdul[0].header['BEAMFILE']=os.path.basename(beam.model_beam_file)
            hdul.flush()
        else:
            print ("Primary beam correction already done. Not repeating")

    return
    
def plot_beamfile(filename,vmax=2,vmin=0.001,figfile=None):
    with fits.open(filename) as hdul:
        flux=hdul[0].data
    
    plt.imshow(np.squeeze(flux),origin='lower',aspect='auto',vmax=vmax,vmin=vmin)
    plt.colorbar()
    if figfile:
        plt.savefig(figfile)
    plt.show()
    
def read_background_file(filename):
    with h5py.File(filename,'r') as hf:
        background_flux=np.array(hf['background_SFU'])
        freqs=np.array(hf['freq_MHz'])
        mjd=np.array(hf['mjd'])
    print (background_flux.shape,freqs.shape,mjd.shape)
    interp = RegularGridInterpolator((freqs,mjd), background_flux)
    return  interp

def do_background_subtraction(filename,interpolating_func):
    
    with fits.open(filename,mode='update') as hdul:
        if 'BKG_SUB' not in hdul[0].header:
            
            flux=np.squeeze(hdul[0].data)
            
            freqs=hdul[1].data
            tim = astropy.table.Table(hdul[2].data)
            tmjd= np.array(tim['mjd']) + np.array(tim['time']) /86400000
            
            
            freq_MHz=np.array(freqs.tolist())*1e3
            
            

            ### using the periodic behaviour
            try:
                background_flux=interpolating_func((freq_MHz,tmjd))
            except:
                try:
                    background_flux=interpolating_func((freq_MHz,tmjd+365))
                except:
                    background_flux=interpolating_func((freq_MHz,tmjd-365))
            
            

            flux-=background_flux
            
            hdul[0].data=np.expand_dims(flux,axis=(0,1))
            hdul[0].header['BKG_SUB']=True
            hdul[0].header['CTYPE1']='stokes'
            hdul[0].header['CRVAL1']=(1,'I') 
            hdul[0].header['CTYPE2']=('dummy','kept for historical reasons')
            hdul[0].header['CTYPE3']='freq'
            hdul[0].header['CTYPE4']='time'
            hdul[0].header['CUNIT3']='GHz'
            hdul[0].header['BUNIT']='SFU'
            
            hdul[1].header['TUNIT1']='GHz'
            
            hdul.flush()
        else:
            print ("Background is already subtracted. Not repeating.")
    return
    



basedir='/common/lwa/spec/fits'
background_filename='/common/lwa/support_files/smoothed_background.hdf5'

starttime=dt.datetime(2024,2,16,1,0,0) ##dt.datetime(2023,10,12,1,0,0)
tdt=dt.timedelta(days=1)

endtime=dt.datetime(2025,7,8,10,0,0)

st=starttime

os.chdir(basedir)




interp_func=read_background_file(background_filename)

current_year=st.year

datedir=os.path.join(basedir,st.strftime("%Y"))
print (datedir)
if os.path.isdir(datedir):
    os.chdir(datedir)


j=0

with open("/common/lwa/bad_times.txt","w") as f1:
    while st<endtime:
       
        if st.year>current_year:
            datedir=os.path.join(basedir,st.strftime("%Y"))
            if os.path.isdir(datedir):
                os.chdir(datedir)
            current_year=st.year        
        
        try:    
            file1=f'ovro-lwa.lev1_bmf_256ms_96kHz.{st.strftime("%Y-%m-%d")}.dspec_I.fits'
            
    
            if os.path.isfile(file1):
                 
                do_background_subtraction(file1,interp_func)
                
                do_beam_correction(file1)
                
                print (f'{file1} finalized')
        except Exception as e:
            f1.write(f"{file1}\n")
            pass
        
    
        st+=tdt
        j+=1
#        if j==2:
#            break
        




    
    

    
    
    
    
    
    
        
