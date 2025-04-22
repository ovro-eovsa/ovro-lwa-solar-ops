from astropy.coordinates import SkyCoord, EarthLocation, get_body, AltAz
import numpy as np
from astropy.time import Time, TimeDelta
import astropy.units as u
import sys,os,glob,logging
import data_downloader
from casatools import table
from ovrolwasolar import utils,flagging,calibration
from casatools import msmetadata
from casatasks import applycal
from ovrolwasolar import beam_polcalib,utils
from ovrolwasolar.beam_polcalib import image_polcal_astronomical_source as img_polcal


def source_riseset(skycoord, date_time,observatory='ovro', altitude_limit=15):
    '''
    :param date_time: input time in astropy.time.Time format
    :param observatory: name of the observatory recognized by astropy
    :param altitude_limit: lower limit of altitude to consider. Default to 15 degrees.
    :param skycoord: Source can be provided by using a astropy skycoordinate. 
                     
    :return trise, tset  ## has 30 min resolution
    '''
    try:
        date_mjd = Time(date_time).mjd
    except Exception as e:
        logging.error(e)

    obs = EarthLocation.of_site(observatory)
    t0 = Time(int(date_mjd), format='mjd')
    source_rise=t0
    

    alt = skycoord.transform_to(AltAz(obstime=t0, location=obs)).alt.degree

    source_risen=False
    if alt>altitude_limit:
        source_risen=True

    if not source_risen:    
        while alt < altitude_limit:
            source_rise += TimeDelta(1800., format='sec')
            alt = skycoord.transform_to(AltAz(obstime=source_rise, location=obs)).alt.degree
            
        source_set=source_rise+TimeDelta(1800., format='sec')
        while alt > altitude_limit:
            source_set += TimeDelta(1800., format='sec')
            alt = skycoord.transform_to(AltAz(obstime=source_set, location=obs)).alt.degree
        return source_rise,source_set
    else:
        while alt > altitude_limit:
            source_rise -= TimeDelta(1800., format='sec')
            alt = skycoord.transform_to(AltAz(obstime=source_rise, location=obs)).alt.degree
        
        source_set=t0
        alt = skycoord.transform_to(AltAz(obstime=source_set, location=obs)).alt.degree
        while alt>altitude_limit:
            source_set += TimeDelta(1800., format='sec')
            alt = skycoord.transform_to(AltAz(obstime=source_set, location=obs)).alt.degree

        return source_rise,source_set

def flag_ms(msfiles):
    for ms_file in msfiles:
        try:
           antflagfile = flagging.flag_bad_ants(ms_file)
        except:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            logging.error(exc_type, fname, exc_tb.tb_lineno)
    return

def gen_caltables(calib_in, bcaltb=None, uvrange='>10lambda', refant='202', flag_outrigger_antenna=True, 
        proc_dir='./',doplot=True):
    
    """
    Function to generate calibration tables for a list of calibration ms files
    :param calib_in: input used for calibration. This can be either a) a string of time stamp recognized by astropy.time.Time 
            or b) a specific list of ms files used for calibration
    :param bcaltb: name of the calibration tables. Use the default if None
    :param uvrange: uv range to be used, default to '>10lambda'
    :param refant: reference antenna, default to '202'
    :param flag_outrigger_antenna: if True, flag all outrigger antennas. These would be used for beamforming.
    :proc_dir: directory to process the data and hold output files.
    
    In the early days we were generating a normalising factor for each channel. Now the normalization is being
    done in the beamformer itself. I am removing that part of the code on April 4, 2025.
    """
    import pandas as pd
    
    msmd=msmetadata()
    
    download_fold = proc_dir + '/ms_calib/'
    caltable_fold = proc_dir + '/caltables/'
    beam_caltable_fold = proc_dir + '/caltables_beam/'

    if not os.path.exists(download_fold):
        os.makedirs(download_fold)

    if not os.path.exists(caltable_fold):
        os.makedirs(caltable_fold)

    if not os.path.exists(beam_caltable_fold):
        os.makedirs(beam_caltable_fold)

    bcaltbs = []
    chan_freqs = []
    if type(calib_in) == list:
        calib_in.sort()
        ms_calib = calib_in
    elif type(calib_in) == str:
        try:
            calib_time = Time(calib_in)
        except:
            print('<<',Time.now().isot,'>>','The input time needs to be astropy.time.Time format. Abort...')
            return -1
        ms_calib = data_downloader.download_calibms(calib_time, download_fold=download_fold)
    else:
        print('<<',Time.now().isot,'>>','Input not recognized. Abort...')
        return -1

    
    if len(ms_calib) > 0:
        # TODO: somehow the parallel processing failed if flagging has run. I have no idea why. Returning to the slow serial processing.
        #pool = multiprocessing.pool.Pool(processes=len(ms_calib))
        #gen_calib_partial = partial(calibration.gen_calibration, uvrange=uvrange, caltable_fold=caltable_fold,
        #                refant=refant)
        #result = pool.map_async(gen_calib_partial, ms_calib)
        #timeout = 2000.
        #result.wait(timeout=timeout)
        #bcaltbs = result.get()
        #pool.close()
        #pool.join()
        flag_ms(ms_calib)
        
        for ms_calib_ in ms_calib:
            try:
                bcaltb = calibration.gen_calibration(ms_calib_, uvrange=uvrange, caltable_fold=caltable_fold, refant=refant)
                msmd.open(ms_calib_)
                chan_freqs.append(msmd.chanfreqs(0))
                msmd.done()
                bcaltbs.append(bcaltb)
            except Exception as e:
                print('<<',Time.now().isot,'>>','Something is wrong when making calibrations for ', ms_calib_)
                print(e)
            
        chan_freqs = np.concatenate(chan_freqs)
        if doplot:
            create_waterfall_plot(bcaltbs,ms_calib)
    else:
        print('<<',Time.now().isot,'>>','The list of calibration ms files seems to be empty. Abort...')
        return -1
        

    if flag_outrigger_antenna:
        bcaltbs_bm = []
        for bcaltb in bcaltbs:
            bcaltb_bm = beam_caltable_fold + '/' + os.path.basename(bcaltb)
            os.system('cp -r ' + bcaltb + ' ' + bcaltb_bm)
            flag_outrigger(bcaltb_bm,ms_calib[0])
            
            bcaltbs_bm.append(bcaltb_bm)

        
        return bcaltbs, bcaltbs_bm
    else:
        return bcaltbs

def get_gain_amplitude_phase(caltable):
    tb=table()
    tb.open(caltable)
    try:
        data=tb.getcol('CPARAM')
        flag=tb.getcol('FLAG')
    finally:
        tb.close()

    pos=np.where(flag==True)
    data[pos]=np.nan
    return np.abs(data),np.angle(data)



def robust_linear_fit(x,y,num_trial=5,thresh=3):
    for i in range(num_trial):
        pos=np.where(np.isnan(y)==False)
        poly=np.polyfit(x[pos],y[pos],deg=1)
        predicted_y=np.poly1d(poly)(x)
        residual=predicted_y-y
        std=np.nanstd(residual)
        pos=np.where(residual>thresh*std)
        y[pos]=np.nan
    return poly
    
def find_delay(freqs,phase):
    pos=np.where(phase<0)[0]
    if pos.size>0:
        phase+=2*np.pi ### assume that phase is from -pi to pi

    pos=np.where((np.isnan(phase)==False) & (freqs>40*1e6))[0]
    if pos.size==0:
        return np.nan
    phase1=phase[pos]
    freqs1=freqs[pos]/1e6 ### converting to MHz
    unwrapped_phase=np.unwrap(phase1)
    delay=robust_linear_fit(freqs1,unwrapped_phase)[0]/(2*np.pi) ### delay is in microseconds
    return delay
    
def find_delay_all_ant_corr(freqs,phase):
    phase_shape=phase.shape
    num_corr=phase_shape[0]
    num_ant=phase_shape[2]
    
    delay=np.zeros((num_corr,num_ant))*np.nan
    
    for i in range(num_corr):
        for j in range(num_ant):
            wrapped_phase=phase[i,:,j]
            delay[i,j]=find_delay(freqs,wrapped_phase)
    return delay
    
def create_waterfall_plot(caltables,msnames,figname=None,num_chan=192,num_ant=352):
    import matplotlib.pyplot as plt
    
    bands=['13MHz', '18MHz', '23MHz', '27MHz', '32MHz', '36MHz', '41MHz', \
            '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', \
            '78MHz', '82MHz']
            
    msnames.sort()
    
    assert len(caltables)<=len(msnames),"The number of MS is smaller than the number of caltables. Please check"

    num_msnames=len(msnames)
    num_bands=len(bands)
    
    ms_freqs_str=[]
    caltable_freqs=np.zeros(len(caltables))
    
    amp=np.zeros((2,num_bands*num_chan,num_ant))
    phase=np.zeros((2,num_bands*num_chan,num_ant))
    freqs=np.zeros((num_bands*num_chan))

    for j,msname in enumerate(msnames):
        ms_freqs_str.append(utils.get_freqstr_from_name(msname))
    
    for j,caltable in enumerate(caltables):
        caltable_freqs[j]=utils.get_caltable_freq(caltable)[0]/1e6 ### converting to MHz
    

    j=0   
    for k,band in enumerate(bands):
        if band!=ms_freqs_str[j]:
            print ("MS file for "+band+" is missing")
            amp[:,k*num_chan:(k+1)*num_chan,:]=np.nan
            phase[:,k*num_chan:(k+1)*num_chan,:]=np.nan
            freqs[k*num_chan:(k+1)*num_chan]=np.nan
            continue    
        msname=msnames[j]
        ms_freq_MHz=int(ms_freqs_str[j].split('MHz')[0])
        ind=np.argmin(abs(ms_freq_MHz-caltable_freqs))
        #ind=caltable_freqs.index(ms_freqs[j])
        if abs(ms_freq_MHz-caltable_freqs[ind])>2:
            print ("Freq index not found. Some caltables have not been produced")
            amp[:,k*num_chan:(k+1)*num_chan,:]=np.nan
            phase[:,k*num_chan:(k+1)*num_chan,:]=np.nan
            freqs[k*num_chan:(k+1)*num_chan]=np.nan
        else:    
            amp[:,k*num_chan:(k+1)*num_chan,:], phase[:,k*num_chan:(k+1)*num_chan,:]=\
                                                            get_gain_amplitude_phase(caltables[ind])
            freqs[k*num_chan:(k+1)*num_chan]=utils.get_caltable_freq(caltables[ind])
        j+=1
    
    delays=find_delay_all_ant_corr(freqs,phase)
        
    fig,ax=plt.subplots(nrows=2,ncols=1,sharex=True,figsize=[15,10])

    for ind,pol in zip([0,1],['XX','YY']):
        im=ax[ind].imshow(amp[ind,:,:],origin='lower',interpolation='none',vmax=0.005,\
                            vmin=0.0002,cmap='viridis',extent=(0,351,13,87),aspect='auto')
        ax[ind].set_title(pol)
        plt.colorbar(im,ax=ax[ind])

    fig.suptitle(utils.get_timestr_from_name(caltables[0]))
    fig.supxlabel("Antenna")
    fig.supylabel("Gain amplitude")

    if not isinstance(figname,list):
        figname_amp='_'.join(caltables[0].split('_')[:-1])+"_amp.png"
    else:
        figname_amp=figname[0]
    plt.savefig(figname_amp)
    plt.close()
    
    fig,ax=plt.subplots(nrows=2,ncols=1,sharex=True,figsize=[12,10])
    
    for ind,pol in zip([0,1],['XX','YY']):
        im=ax[ind].plot(delays[ind],'ro')
        ax[ind].set_title(pol)

    fig.suptitle(utils.get_timestr_from_name(caltables[0]))
    fig.supxlabel("Antenna")
    fig.supylabel("Delay (microseconds)")
    if not isinstance(figname,list):
        figname_delay='_'.join(caltables[0].split('_')[:-1])+"_delay.png"
    else:
        figname_delay=figname[1]
    plt.savefig(figname_delay)
    plt.close()
    
    return

def flag_outrigger(dataset, ref_ms):
    '''
    dataset: can be a MS or caltable
    ref_MS: MS from which the antenna names will be read
    
    Note that this function assumes that there is a single timeslice
    in the dataset. If there are more than one timelsice, this
    function will give wrong results.
    '''
    tb=table()
    core_ant_ids, exp_ant_ids = flagging.get_antids(ref_ms)
    tb.open(dataset,nomodify=False)
    try:
        flag=tb.getcol('FLAG')
        flag[:,:,exp_ant_ids]=True
        tb.putcol('FLAG',flag)
        tb.flush()
    finally:
        tb.close()
    return

def crosshand_phase_solver(starttime,endtime,tdt,sky_coord,freq_avg=16,proc_dir='./',\
                            model_beam_file='/lustre/msurajit/beam_model_nivedita/OVRO-LWA_soil_pt.h5',\
                            caltable_folder='/lustre/solarpipe/realtime_pipeline/caltables_latest',\
                            ):
    '''
    This function is the wrapper function for determining the
    crosshand phase using a A-team source.
    
    :param starttime: Startime of the time period which will be
                      used for the solution
    :type starttime: Astropy isot format
    :param endtime: Endtime of the time period which will be
                      used for the solution
    :type endtime: Astropy isot format
    :param tdt: Time difference in seconds between consecutive datasets chosen for 
                solutions
    :type tdt: int
    :param sky_coord: sky coordinates of the object which will be used for solution
                     RA-Dec coordinates in the format "J2000 HHhMMmSS.Ss DDdMMmSS.Ss"
                     Note the separation by hms and dms in RA and DEC respectively.
    :type sky_coord:  str
    :param freq_avg: Number of channels to average
    :type freq_avg: int
    '''
    dynamic_spectrum,freqs,tim_mjds=get_source_DS(starttime,endtime,tdt,sky_coord,proc_dir=proc_dir,\
                                    caltable_folder=caltable_folder,model_beam_file=model_beam_file, freq_avg=freq_avg)
    
    img_pol=img_polcal(dynamic_spectrum=dynamic_spectrum,\
                    freqs=freqs,\
                    tim_mjds=mjds,\
                    sky_coord=sky_coord)
                    
    img_pol.crosshand_phase_solver()
    
    fig,ax=plt.subplots(nrows=1,ncols=3,figsize=[12,4],constrained_layout=True)
    for j,(ax1,stokes) in enumerate(zip(ax,['Q','U','V'])):
        ax1.plot(freqs,img_pol.leakage[j+1,:],'o-')
        ax1.set_ylabel(stokes+" leakage fraction")
        ax1.set_ylim([-0.2,0.2])
    fig.supxlabel("Frequency (MHz)")
    plt.savefig("DI_leakage_variation.png")
    plt.close()
    
    plt.plot(img_pol.freqs,img_pol.crosshand_theta,'o-')
    plt.xlabel("Frequency (MHz)")
    plt.ylabel("Crosshand phase (radian)")
    plt.savefig("crosshand_theta_variation.png")
    plt.close()
    
    return img_pol.freqs,img_pol.crosshand_theta
    
def apply_crosshand_phase_on_caltables(caltables,crosshand_phase,crosshand_freqs,inplace=False):
    
    crosshand_applied=[None]*len(caltables)
    
    for j,caltable in enumerate(caltables):
        if not inplace:
            copied_caltable=caltable+".crosshand"
            crosshand_applied[j]=copied_caltable
            os.system("cp "+caltable+" "+copied_caltable)
        else:
            crosshand_applied[j]=caltable
        beam_polcalib.combine_crosshand_theta_on_caltable(crosshand_applied[j],crosshand_phase,crosshand_freqs)
    
    return  crosshand_applied         
        
def get_source_DS(starttime,endtime,tdt,sky_coord,proc_dir='./',\
                     bands=['13MHz', '18MHz', '23MHz', '27MHz', '32MHz', '36MHz', '41MHz', \
                    '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', \
                    '78MHz', '82MHz'],caltable_folder='/lustre/solarpipe/realtime_pipeline/caltables_latest',\
                    freq_avg=16):
    '''
    This function will produce the dynamic spectrum of the source. It will also do this
    after flagging outrigger. Then it will phase up to the sky_coord supplied. Then it
    will add the visibilities to obtain the pseudo-stokes I,Q,U and V. The outrigger
    flagging is done to take care of the bi-lobe structure seen in Stokes V. It is 
    not seen if the angular resolution is poor. But the crosshand phase in independent
    of angular resolution.
    
    :param starttime: Startime of the time period which will be
                      used for the solution
    :type starttime: Astropy isot format
    :param endtime: Endtime of the time period which will be
                      used for the solution
    :type endtime: Astropy isot object
    :param tdt: Time difference in seconds between consecutive datasets chosen for 
                solutions
    :type tdt: int
    :param sky_coord: sky coordinates of the object which will be used for solution
                     RA-Dec coordinates in the format "HH:MM:SS.S DD.MM.SS.S"
    :type sky_coord:  str
    :param proc_dir: Path of directory where all processing will happen
    :type proc_dir: str
    '''
    
            
    
    num_bands=len(bands)
    
    
    tdt=TimeDelta(tdt,format='sec')
    starttime=Time(starttime,format='isot',scale='utc')
    endtime=Time(endtime,format='isot',scale='utc')
    
    st=starttime
    end_mjd=endtime.mjd
    
    dynamic_spectrum=[]
    mjds=[]
    
    download_fold = proc_dir + '/ms_calib/'
    if not os.path.exists(download_fold):
        os.makedirs(download_fold)
    
    while st.mjd<=end_mjd:
        msfiles=data_downloader.download_calibms(st,download_fold=download_fold)#glob.glob(download_fold+"/*.ms")#
        flag_ms(msfiles)
        
        if len(msfiles)==num_bands:
            stokes_data,freqs=get_source_spectrum_single_time(msfiles,sky_coord)
        else:
            stokes_data,_=get_source_spectrum_single_time(msfiles,sky_coord)
        dynamic_spectrum.append(stokes_data)
        mjds.append(st.mjd)
        st+=tdt
    
    concat_ds=np.array(dynamic_spectrum)
    concat_ds=np.swapaxes(np.swapaxes(concat_ds,0,1),1,2)
    mjds=np.array(mjds)
    
    DS=np.swapaxes(concat_ds,1,2) ### now frequency in axis=2
    shape=DS.shape
    num_chan=shape[2]
    num_times=shape[1]
    num_pol=shape[0]
    if num_chan%freq_avg!=0:
        num_chan=int((num_chan//freq_avg)*freq_avg)
    DS_reshaped=DS[:,:,:num_chan].reshape((num_pol,num_times,num_chan//freq_avg,freq_avg))
    DS_freq_avged=np.nanmean(DS_reshaped,axis=3)
    freqs_reshaped=freqs[:num_chan].reshape((num_chan//freq_avg,freq_avg))
    freqs_avged=np.nanmean(freqs_reshaped,axis=1)
    
    return np.swapaxes(DS_freq_avged,1,2),freqs_avged,mjds ## freqs in MHz
            
            
def get_source_spectrum_single_time(msfiles,sky_coord,caltable_folder='/lustre/solarpipe/realtime_pipeline/caltables_latest'):
    
    bands=['13MHz', '18MHz', '23MHz', '27MHz', '32MHz', '36MHz', '41MHz', \
            '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', \
            '78MHz', '82MHz']
            
    
    num_bands=len(bands)
    
    ms_freqs_str=[]
    
    
    msfiles.sort()
    
    
    num_msfiles=len(msfiles)
    
    for j,msname in enumerate(msfiles):
        ms_freqs_str.append(utils.get_freqstr_from_name(msname))

    ms_freqs=np.zeros(num_msfiles)
    for j,msname in enumerate(msfiles):
        frequencies=utils.get_caltable_freq(msname)[0]/1e6 ### converting to MHz
        ms_freqs[j]=frequencies
        
    msmd=msmetadata()
    msmd.open(msfiles[0])
    try:
        num_ant=msmd.nantennas()
        num_chan=msmd.nchan(0)
    finally:
        msmd.done()
        
    num_correlation=int(num_ant*(num_ant-1)/2+num_ant)
    stokes_data=np.zeros((4,num_bands*num_chan))
    freqs=np.zeros((num_bands*num_chan))
    
        
    j=0   
    for k,band in enumerate(bands):
        if j>len(ms_freqs_str)-1:
            logging.info("MS file for "+band+" is missing or bandpass caltable was not present")
            print ("MS file for "+band+" is missing or bandpass caltable was not present")
            stokes_data[:,k*num_chan:(k+1)*num_chan]=np.nan
            freqs[k*num_chan:(k+1)*num_chan]=np.nan
            continue    
            
        caltables=glob.glob(os.path.join(caltable_folder,"*"+ms_freqs_str[j]+"*.bcal"))
        
        if band!=ms_freqs_str[j] or len(caltables)==0:
            logging.info("MS file for "+band+" is missing or bandpass caltable was not present")
            print ("MS file for "+band+" is missing or bandpass caltable was not present")
            stokes_data[:,k*num_chan:(k+1)*num_chan]=np.nan
            freqs[k*num_chan:(k+1)*num_chan]=np.nan
            continue    
        msname=msfiles[j]
        ms_freq_MHz=int(ms_freqs_str[j].split('MHz')[0])
        ind=np.argmin(abs(ms_freq_MHz-ms_freqs))
        if abs(ms_freq_MHz-ms_freqs[ind])>2:
            logging.info("Freq index not found. Some MSfiles are missing")
            print ("Freq index not found. Some MSfiles are missing, ",ms_freq_MHz)
            stokes_data[:,k*num_chan:(k+1)*num_chan]=np.nan
            freqs[k*num_chan:(k+1)*num_chan]=np.nan
            continue
        else:    
            applycal(msname,gaintable=caltables[0])
            stokes_data[:,k*num_chan:(k+1)*num_chan]=get_IQUV(msname,sky_coord)                                                            
            freqs[k*num_chan:(k+1)*num_chan]=utils.get_caltable_freq(msname)*1e-6 #### converting to MHz
        j+=1
    return stokes_data,freqs
    
def get_baseline_index(ant_id1, ant_id2,num_ant=352,auto_corr=True):
    '''
    ant_id1 <=ant_id2
    '''
    if auto_corr==True:
        return num_ant*ant_id1+ant_id2-ant_id1*(ant_id1+1)/2
    else:
        return num_ant*ant_id1+ant_id2-(ant_id1+1)*(ant_id1+2)/2 
    
def get_IQUV(msname,sky_coord):
    '''
    Returns the pseudo-stokes IQUV flux of a source. First it flags the out_rigger
    antennas to ensure that the array has a rather poor spatial resolution. This
    also helps in ensuring that the source is within the PSF even when shifts due to
    ionospheric refraction is high.
    
    '''
    flag_outrigger(msname,msname)
    os.system("chgcentre " + msname + " " + sky_coord)
    
    datacolumn='DATA'
    corrected_present=utils.check_corrected_data_present(msname)
    if corrected_present:
        datacolumn='CORRECTED_DATA'
        
    tb=table()
    tb.open(msname)
    try:
        data=tb.getcol(datacolumn)
        flag=tb.getcol('FLAG')
    finally:
        tb.close()
    
    pos=np.where(flag==True)
    data[pos]=np.nan
    
    auto_corr_inds=np.array(get_baseline_index(np.arange(0,352),np.arange(0,352)),dtype=int)
    data[:,:,auto_corr_inds]=np.nan
    
    I=np.real(0.5*np.nanmean(data[0,:,:]+data[3,:,:],axis=1))
    Q=np.real(0.5*np.nanmean(data[0,:,:]-data[3,:,:],axis=1))
    U=np.real(0.5*np.nanmean(data[1,:,:]+data[2,:,:],axis=1))
    V=np.imag(0.5*np.nanmean(data[1,:,:]-data[2,:,:],axis=1))
    
    shape=data.shape
    stokes_data=np.zeros((shape[0],shape[1]))
    stokes_data[0,...]=I
    stokes_data[1,...]=Q
    stokes_data[2,...]=U
    stokes_data[3,...]=V
    
    return stokes_data
    
    
            
        
        
    
    
    
    
