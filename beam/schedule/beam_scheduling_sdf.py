'''
This code is primarily written by Prof. Dale Gary. Surajit Mondal
used his template and added a new function called make_sdf, which 
can be used to track an arbitrary ra-dec location.
'''
from astropy.coordinates import SkyCoord, EarthLocation, get_body, AltAz
from astropy.time import Time
import astropy.units as u
from time import sleep
import numpy as np
import os, subprocess
from copy import copy
from lwa_antpos.station import Station
station=Station('ovro')
obs = EarthLocation.from_geocentric(*station.ecef, unit=u.m)  # NB: The code below uses this global variable!
scandur = 60.   # Duration of a normal scan (minutes)
scangap = 1.    # Time gap between scans (minutes)
def approx_rise(t=None):
    '''  Given a time in MJD, determine the approximate sunrise time for that date.
         Sunrise is defined as rising above 10 degrees altitude.
         obs is the 'ovro' observatory object calculated above.
    '''
    if t is None:
        t = Time.now().mjd
    t1 = int(t) + 13*3600/86400.
    sc = get_body('sun', Time(t1, format='mjd'), location=obs)
    for dt in range(100):
        ts = t1 + dt*300./86400.
        tstart = Time(ts,format='mjd')
        alt = sc.transform_to(AltAz(obstime=tstart, location=obs)).alt.degree
        if alt > 10.:
            return Time(ts,format='mjd').iso

def approx_set(t=None):
    '''  Given a time in MJD, determine the approximate sunset time for that date.
         Sunrise is defined as setting below 11 degrees altitude.
         obs is the 'ovro' observatory object calculated above.
    '''
    if t is None:
        t = Time.now().mjd
    t1 = int(t) + 23.5*3600/86400.
    sc = get_body('sun', Time(t1, format='mjd'), location=obs)
    for dt in range(100):
        ts = t1 + dt*300./86400.
        tstart = Time(ts,format='mjd')
        alt = sc.transform_to(AltAz(obstime=tstart, location=obs)).alt.degree
        if alt < 11.:
            return Time(ts,format='mjd').iso

def noon_radec(t=None):
    if t is None:
        t = Time.now().mjd
    t1 = int(t) + 20*3600/86400.   # This is 20 UT (near noon) on requested date
    sc = get_body('sun', Time(t1, format='mjd'), location=obs)
    ra = sc.ra.hour
    dec = sc.dec.deg+20
    return ra, dec

def make_solar_sdf(trange=None, beam_to_use=2):
    ''' Make a standard solar sdf file for scheduling the solar beam observations.

        trange     An astropy Time() object representing a timerange (start and stop time).
                     If only one time is given, the rise and set times for that date are
                     used to calculate a timerange.  If None, the rise and set times for
                     the current day are used.
    '''
    from numpy import ceil
    if trange is None:
        trange = Time([approx_rise(), approx_set()])  # No time given => use today's date
    else:
        try:
            if len(trange) >= 2:
                pass   # Will use the first two times in an array of times as start/stop times
        except:
            # Case of trange existing but having no len(), which is hopefully due to it having only one (start) time.
            trange = Time([trange.iso, approx_set(t=trange.mjd)])
        # Ensure that start and stop times are bounded by rise and set times:
        mjd0 = trange[0].mjd
        mjd1 = trange[1].mjd
        mjdrise = Time(approx_rise(mjd0)).mjd
        mjdset = Time(approx_set(mjd0)).mjd
        trange = Time([max(mjd0,mjdrise),min(mjd1,mjdset)],format='mjd')

    # Create a template file using "lwaobserving" create command
    blah = subprocess.run(["lwaobserving","create-sdf","--n-obs","1","--sess-mode","POWER","--beam-num",str(beam_to_use),
       "--obs-mode","TRK_SOL", "--obs-start", trange[0].isot[:19],"--obs-dur","1800000","--obj-name","sun",
       "--int-time","64","--do-cal","/lustre/solarpipe/solar_beam_sdfs/template.sdf"],stdout = subprocess.PIPE)
    output = blah.stdout.decode('utf-8').split('\n')
    f = open("/lustre/solarpipe/solar_beam_sdfs/template.sdf",'r')
    lines = f.readlines()
    f.close()
    hdr_block = lines[:11]
    obs_block = lines[11:]
    obs_block.append('\n')   # Add a blank link to separate sections
    alt_block = copy(obs_block)  # Off Sun pointing block
    alt_block = alt_block[1:]    # Shift lines up one line to eliminate OBS_TARGET for offsun
    sdf_name = '/lustre/solarpipe/solar_beam_sdfs/solar_'+trange[0].iso[:19].replace('-','').replace(':','').replace(' ','_')+'.sdf'
    f = open(sdf_name,'w')
    fmt_changed = False  # Flag to indicate that the schedule file format has changed (error condition)
    if 'CONFIG_FILE' in hdr_block[7]:
        #hdr_block[7] = 'CONFIG_FILE      /opt/devel/dgary/lwa_config_calim_std.yaml\n'
        hdr_block[7] = 'CONFIG_FILE      /opt/devel/solarpipe/lwa_config_calim_std.yaml\n'
    else: fmt_changed = True
    if 'CAL_DIR' in hdr_block[8]:
        hdr_block[8] = 'CAL_DIR          /opt/devel/solarpipe/operation/caltab/caltables_beam_latest\n'
    else: fmt_changed = True
    for line in hdr_block:
        f.write(line)
    durmin = (trange[1].mjd - trange[0].mjd)*1440
    nblocks = int(ceil(durmin/(scandur+scangap)))
    offsun = False
    kblock = 0
    for i in range(nblocks):
        kblock += 1
        tstart = trange[0].mjd + i*(scandur+scangap)/1440.
        if tstart - int(trange[0].mjd) > 20/24. and not offsun:
            ra, dec = noon_radec(tstart)   # Get RA Dec of a point 10 degrees north of the Sun
            tend = min(tstart + 2/1440.,trange[1].mjd)
            alt_block[0] = 'OBS_ID          '+str(kblock)+'\n'
            alt_block[1] = 'OBS_START_MJD   '+str(int(tstart))+'\n'
            alt_block[2] = 'OBS_START_MPM   '+str(int(ceil((tstart % 1)*86400000.)))+'\n'
            alt_block[3] = 'OBS_START       UTC '+Time(tstart,format='mjd').iso[:19].replace('-',' ')+'\n'
            alt_block[4] = 'OBS_DUR         '+str(int((tend-tstart)*86400000.))+'\n'
            alt_block[5] = 'OBS_INT_TIME    64\n'
            alt_block[6] = 'OBS_DUR+        00:02:00\n'
            alt_block[7] = 'OBS_MODE        TRK_RADEC\n'
            alt_block[8] = 'OBS_RA          '+str(ra)+'\n'
            alt_block[9] = 'OBS_DEC         '+str(dec)+'\n'
            for line in alt_block:
                f.write(line)
            kblock += 1
            tstart +=  (2. + scangap)/1440.
            tend = min(tstart + (scandur - 2. - scangap)/1440.,trange[1].mjd)
            durminutes = (tend - tstart)*1440.
            hh = int(durminutes / 60)
            mm = int(durminutes % 60)
            ss = int(((durminutes % 60) - mm)*60)
            durstr = '{:02d}:{:02d}:{:02d}'.format(hh,mm,ss)
            obs_block[0] = 'OBS_ID          '+str(kblock)+'\n'
            obs_block[2] = 'OBS_START_MJD   '+str(int(tstart))+'\n'
            obs_block[3] = 'OBS_START_MPM   '+str(int(ceil((tstart % 1)*86400000.)))+'\n'
            obs_block[4] = 'OBS_START       UTC '+Time(tstart,format='mjd').iso[:19].replace('-',' ')+'\n'
            obs_block[5] = 'OBS_DUR         '+str(int((tend-tstart)*86400000.))+'\n'
            obs_block[7] = 'OBS_DUR+        '+durstr+'\n'
            for line in obs_block:
                f.write(line)
            offsun = True
        else:
            tend = min(tstart + scandur/1440.,trange[1].mjd)
            durminutes = (tend - tstart)*1440.
            if durminutes > 0:
                hh = int(durminutes / 60)
                mm = int(durminutes % 60)
                ss = int(((durminutes % 60) - mm)*60)
                durstr = '{:02d}:{:02d}:{:02d}'.format(hh,mm,ss)
                if 'OBS_ID' in obs_block[0]:
                    obs_block[0] = 'OBS_ID          '+str(kblock)+'\n'
                else: fmt_changed = True
                if 'OBS_START_MJD' in obs_block[2]:
                    obs_block[2] = 'OBS_START_MJD   '+str(int(tstart))+'\n'
                else: fmt_changed = True
                if 'OBS_START_MPM' in obs_block[3]:
                    obs_block[3] = 'OBS_START_MPM   '+str(int(ceil((tstart % 1)*86400000.)))+'\n'
                else: fmt_changed = True
                if 'OBS_START' in obs_block[4]:
                    obs_block[4] = 'OBS_START       UTC '+Time(tstart,format='mjd').iso[:19].replace('-',' ')+'\n'
                else: fmt_changed = True
                if 'OBS_DUR' in obs_block[5]:
                    obs_block[5] = 'OBS_DUR         '+str(int((tend-tstart)*86400000.))+'\n'
                else: fmt_changed = True
                if 'OBS_DUR+' in obs_block[7]:
                    obs_block[7] = 'OBS_DUR+        '+durstr+'\n'
                if fmt_changed:
                    print('The format of the sdf file has changed! Code must be modified to continue. Exiting.')
                    return
                else:
                    for line in obs_block:
                        f.write(line)
    f.close()
    return sdf_name

def multiday_obs(ndays=7, startday=0, send=True, beam_to_use=2):
    ''' Make and submit consecutive schedule (.sdf) files.

        Inputs:
          ndays         Integer number of consecutive days for which to create the schedule.  Default is 7.
          startday      Integer number of days ahead to start.  Default is 0 (i.e. start with today's
                          schedule). 
          send          Boolean stating whether to submit the schedule.  Default is True.  If False,
                          sdfs are created in /tmp, but not submitted.
    '''
    for i in range(startday, startday+ndays):
        if i == 0:
            # First time through, make sure initial time is at least 15 min in the future
            mjd = Time.now().mjd + 15.*60./86400.
        else:
            # On subsequent days, set start time to 12 UT (actual time will be sunrise)
            mjd = int(Time.now().mjd) + 0.5 + i
        sdf_name = make_solar_sdf(Time(mjd,format='mjd'), beam_to_use=beam_to_use)
        if send:
            print('Submitting',sdf_name)
            os.system('lwaobserving submit-sdf '+sdf_name)
            sleep(2)
            
def make_sdf(source_coord,trange,beam_to_use=5, scan_duration=1800,source_name='source',\
            integration_time=64, sdf_filename='template.sdf', config_file="/opt/devel/solarpipe/lwa_config_calim_std.yaml",\
            caltable_folder='/opt/devel/solarpipe/operation/caltab/caltables_beam_latest', send=False,\
            scandur = 60.,scangap = 1.):
    ''' Make a standard sdf file for scheduling the beam observations, where a source of given ra-dec is tracked

        :param trange: An astropy Time() object representing a timerange (start and stop time).
                     If only one time is given, the rise and set times for that date are
                     used to calculate a timerange.  If None, the rise and set times for
                     the current day are used.
        :param source_coord : Astropy skycoord of source to track
        :param beam_to_use: Specifies which beam to use for the observation. Default is 5
        :param scan_duration: Specifies the duration of a scan. Takes value in seconds. 
                                Default is 1800, which is 30min
        :param source_name: The source name is appended to the filename. Default is "source"
        :param integration_time: Time integration to do in ms. Default is 64 .
        :param sdf_filename: This is a file template which is created. Defaullt is template.sdf
        :param config_file: The configuration file which will be used. Default is 
                            /opt/devel/dgary/lwa_config_calim_std.yaml  . Remember to verify that the
                            beam you plan to use is listed at the bottom of the file
        :param caltable_folder: This is the folder from which caltables will be chose to calibrate 
                                the beam if required. Ensure that the folder contains one bcal file 
                                for each frequency band. Default is /lustre/solarpipe/realtime_pipeline/caltables_beam_latest
        :param send : If True, submits the sdf. Default is False
        :param scandur: Duration of a normal scan (minutes)
        :param scangap: Time gap between scans (minutes)
    '''
    from numpy import ceil
    
    
    source_ra=source_coord.ra.hour
    source_dec=source_coord.dec.deg

    # Create a template file using "lwaobserving" create command
    blah = subprocess.run(["lwaobserving","create-sdf","--n-obs","1","--sess-mode","POWER","--beam-num",str(beam_to_use),\
        "--obs-mode","TRK_RADEC", "--ra", str(source_ra), "--dec", str(source_dec),"--obs-start", trange[0].isot[:19],\
        "--obs-dur",str(int(scan_duration*1000)),\
        "--obj-name",source_name,"--do-cal","--int-time",str(integration_time),sdf_filename],stdout = subprocess.PIPE)
    output = blah.stdout.decode('utf-8').split('\n')
    f = open(sdf_filename,'r')
    lines = f.readlines()
    f.close()
    hdr_block = lines[:11]
    obs_block = lines[11:]
    obs_block.append('\n')   # Add a blank link to separate sections
    alt_block = copy(obs_block)  # Off Sun pointing block
    alt_block = alt_block[1:]    # Shift lines up one line to eliminate OBS_TARGET for offsun
    
    sdf_name = os.path.join(os.path.dirname(sdf_filename),'source_'+\
                            trange[0].datetime.strftime("%Y%m%d%H%M%S")+\
                            "_"+trange[1].datetime.strftime("%Y%m%d%H%M%S")+".sdf")
    with open(sdf_name,'w') as f:
        fmt_changed = False  # Flag to indicate that the schedule file format has changed (error condition)
        if 'CONFIG_FILE' in hdr_block[7]:
            hdr_block[7] = f'CONFIG_FILE      {config_file}\n'
        else: fmt_changed = True
        if 'CAL_DIR' in hdr_block[8]:
            hdr_block[8] = f'CAL_DIR          {caltable_folder}\n'
        else: fmt_changed = True
        for line in hdr_block:
            f.write(line)
        durmin = (trange[1].mjd - trange[0].mjd)*1440  ### what is 1440?
        scandur=min(scan_duration,scandur) 
        nblocks = int(ceil(durmin/(scandur+scangap)))
        offsun = False
        kblock = 0
        for i in range(nblocks):
            kblock += 1
            tstart = trange[0].mjd + i*(scandur+scangap)/1440.
            
            
            tend = min(tstart + scandur/1440.,trange[1].mjd)
            durminutes = (tend - tstart)*1440.
            if durminutes > 0:
                hh = int(durminutes / 60)
                mm = int(durminutes % 60)
                ss = int(((durminutes % 60) - mm)*60)
                durstr = '{:02d}:{:02d}:{:02d}'.format(hh,mm,ss)
                if 'OBS_ID' in obs_block[0]:
                    obs_block[0] = 'OBS_ID          '+str(kblock)+'\n'
                else: fmt_changed = True
                if 'OBS_START_MJD' in obs_block[2]:
                    obs_block[2] = 'OBS_START_MJD   '+str(int(tstart))+'\n'
                else: fmt_changed = True
                if 'OBS_START_MPM' in obs_block[3]:
                    obs_block[3] = 'OBS_START_MPM   '+str(int(ceil((tstart % 1)*86400000.)))+'\n'
                else: fmt_changed = True
                if 'OBS_START' in obs_block[4]:
                    obs_block[4] = 'OBS_START       UTC '+Time(tstart,format='mjd').iso[:19].replace('-',' ')+'\n'
                else: fmt_changed = True
                if 'OBS_DUR' in obs_block[5]:
                    obs_block[5] = 'OBS_DUR         '+str(int((tend-tstart)*86400000.))+'\n'
                else: fmt_changed = True
                if 'OBS_DUR+' in obs_block[7]:
                    obs_block[7] = 'OBS_DUR+        '+durstr+'\n'
                if fmt_changed:
                    print('The format of the sdf file has changed! Code must be modified to continue. Exiting.')
                    return
                else:
                    for line in obs_block:
                        f.write(line)
        if send:
            print('Submitting',sdf_name)
            os.system('lwaobserving submit-sdf '+sdf_name)
            sleep(2)
    return sdf_name
            
