import glob,os
import datetime as dt
from suncasa import dspec
import numpy as np
import h5py
from astropy.time import Time

basedir='/common/lwa/spec/fits_bcgrd/'
starttime=dt.datetime(2023,10,24,1,0,0)
tdt=dt.timedelta(days=1)

endtime=dt.datetime(2025,7,1,1,0,0)

st=starttime

home=os.getcwd()

hf=h5py.File("OVRO_LWA_background_calibrated.hdf5",'w')

os.chdir(basedir)



current_month=st.month
current_year=st.year

background_data=[]
times=[]

datedir=os.path.join(basedir,st.strftime("%Y"))
print (datedir)
if os.path.isdir(datedir):
    os.chdir(datedir)


try:
    while st<endtime:
        
        if st.year>current_year:
            datedir=os.path.join(basedir,st.strftime("%Y"))
            print (datedir)
            if os.path.isdir(datedir):
                os.chdir(datedir)
            current_year=st.year                
            
        if st.month!=current_month:
            background_data=np.array(background_data)
            times=np.array(times)
            
            prev_time=st-tdt
            group=hf.create_group(prev_time.strftime("%Y%m"))
            group.create_dataset("background",data=background_data)
            group.create_dataset("mjd",data=times)
            del background_data,times
            background_data=[]
            times=[]
            current_month=st.month
            
        
        file1=f'ovro-lwa.lev1_bmf_256ms_96kHz.{st.strftime("%Y-%m-%d")}.dspec_I.fits'
         
 #       print (file1)
        if os.path.isfile(file1):
            try:
                ds=dspec.Dspec()
                ds.read(file1,source='general',timebin=10,freqbin=4,stokes='I')
            
                data=np.squeeze(ds.data)
                shape=data.shape
                print (shape)
                num_times=shape[1]
                num_freqs=shape[0]
                freq_avged_data=np.mean(data,axis=1)
                background_data.append(freq_avged_data)
                
                mjd=Time(st,format='datetime').mjd
                freqs=ds.freq_axis
                #print (freqs[0],freqs[-1])
                times.append(mjd)
#                print (times)
            except Exception as e:
                print (st)
                print (os.getcwd())
                print (file1)
                raise e
                #pass
        else:
            print (f"{file1} is not present")
        st+=tdt
        #print (st)

    background_data=np.array(background_data)
    times=np.array(times)
    prev_time=st-tdt
    group=hf.create_group(prev_time.strftime("%Y%m"))
    group.create_dataset("background",data=background_data)
    group.create_dataset("mjd",data=times)
    hf.create_dataset('freq_MHz',data=np.array(freqs)*1e-6)
finally:
    hf.close()



os.chdir(home)



