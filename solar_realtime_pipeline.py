#!/usr/bin/env python
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
import os, sys, glob, getopt
from ovrolwasolar import solar_pipeline as sp
from ovrolwasolar.primary_beam import analytic_beam as beam
from ovrolwasolar import utils, calibration, flagging
from casatasks import clearcal, applycal, flagdata, tclean, exportfits, imsubimage
from casatools import msmetadata, quanta, measures, table
from suncasa.utils import helioimage2fits as hf
from suncasa.io import ndfits
from ovrolwasolar import file_handler
import logging
import timeit
import multiprocessing
from multiprocessing import TimeoutError
from astropy.time import Time, TimeDelta
import astropy.units as u
from sunpy.coordinates import frames
from astropy.coordinates import SkyCoord, EarthLocation, get_body, AltAz
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
from suncasa.utils import plot_mapX as pmX
import sunpy.map as smap
import shlex, subprocess
from functools import partial
from time import sleep
import socket
from matplotlib.patches import Ellipse
import argparse
import pandas as pd

from ovrolwasolar import visualization as ovis
from ovrolwasolar import refraction_correction as orefr

matplotlib.use('agg')

msmd = msmetadata()
qa = quanta()
me = measures()
tb = table()

def sun_riseset(date=Time.now(), observatory='ovro', altitude_limit=15.):
    '''
    Given a date in Time object, determine the sun rise and set time as viewed from OVRO
    :param date: input time in astropy.time.Time format
    :param observatory: name of the observatory recognized by astropy
    :param altitude_limit: lower limit of altitude to consider. Default to 15 degrees.
    '''
    try:
        date_mjd = Time(date).mjd
    except Exception as e:
        logging.error(e)

    obs = EarthLocation.of_site(observatory)
    t0 = Time(int(date_mjd) + 13. / 24., format='mjd')
    sun_loc = get_body('sun', t0, location=obs)
    alt = sun_loc.transform_to(AltAz(obstime=t0, location=obs)).alt.degree
    while alt < altitude_limit:
        t0 += TimeDelta(60., format='sec')
        alt = sun_loc.transform_to(AltAz(obstime=t0, location=obs)).alt.degree

    t1 = Time(int(date_mjd) + 22. / 24., format='mjd')
    sun_loc = get_body('sun', t1, location=obs)
    alt = sun_loc.transform_to(AltAz(obstime=t1, location=obs)).alt.degree
    while alt > altitude_limit:
        t1 += TimeDelta(60., format='sec')
        alt = sun_loc.transform_to(AltAz(obstime=t1, location=obs)).alt.degree

    return t0, t1


#### define data files #####
def list_msfiles_old(intime, server='lwacalim', distributed=True, file_path='slow',
            nodes = [1, 2, 3, 4, 5, 6, 7, 8], time_interval='10s'):
    """
    Return a list of visibilities to be copied for pipeline processing for a given time
    :param intime: astropy Time object
    :param time_interval: Options are '10s', '1min', '10min', '1hr', '1day'
    :param nband_min: minimum number of available subbands acceptable
    """
    intimestr = intime.isot[:-4].replace('-','').replace(':','').replace('T','_')
    if time_interval == '10s':
        tstr = intimestr[:-1]
    if time_interval == '1min':
        tstr = intimestr[:-2]
    if time_interval == '10min':
        tstr = intimestr[:-3]
    if time_interval == '1hr':
        tstr = intimestr[:-4]
    if time_interval == '1day':
        tstr = intimestr[:9]

    msfiles = []
    if not distributed:
        args = ['ssh', '{}'.format(server), 'ls', '{}'.format(file_path), '|', 'grep', '{}'.format(tstr)]
        p = subprocess.run(args, capture_output=True)
        filenames = p.stdout.decode('utf-8').split('\n')[:-1]
        for filename in filenames:
            if filename[-6:] == 'MHz.ms':
                pathstr = '{0:s}:{1:s}/{2:s}'.format(server, file_path, filename)
                tmpstr = filename[:15].replace('_', 'T')
                timestr = tmpstr[:4] + '-' + tmpstr[4:6] + '-' + tmpstr[6:11] + ':' + tmpstr[11:13] + ':' + tmpstr[13:]
                freqstr = filename[16:21]
                msfiles.append({'path': pathstr, 'name': filename, 'time': timestr, 'freq': freqstr})
    else:
        processes=[]
        for i in nodes:
            args = ['ssh', '{0:s}0{1:d}'.format(server, i), 'ls', '/data0{0:d}/{1:s}'.format(i, file_path), '|', 'grep', '{}'.format(tstr)]
            p = subprocess.Popen(args, stdout=subprocess.PIPE)
            processes.append(p)
            filenames = p.communicate()[0].decode('utf-8').split('\n')[:-1]
            for filename in filenames:
                if filename[-6:] == 'MHz.ms':
                    pathstr = '{0:s}0{1:d}:/data0{2:d}/{3:s}/{4:s}'.format(server, i, i, file_path, filename)
                    tmpstr = filename[:15].replace('_', 'T')
                    timestr = tmpstr[:4] + '-' + tmpstr[4:6] + '-' + tmpstr[6:11] + ':' + tmpstr[11:13] + ':' + tmpstr[13:]
                    freqstr = filename[16:21]
                    msfiles.append({'path': pathstr, 'name': filename, 'time': timestr, 'freq': freqstr})
    return msfiles

def list_msfiles(intime, lustre=True, file_path='slow', server=None, time_interval='10s', 
                 bands=['32MHz', '36MHz', '41MHz', '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', '78MHz', '82MHz']):
    """
    Return a list of visibilities to be copied for pipeline processing for a given time
    :param intime: astropy Time object
    :param server: name of the server to list available data
    :param lustre: if True, specific to lustre system on lwacalim nodes. If not, try your luck in combination with file_path
    :param file_path: file path to the data files. For lustre, it is either 'slow' or 'fast'. For other servers, provide full path to data.
    :param time_interval: Options are '10s', '1min', '10min'
    :param bands: bands to list/download. Default to 12 bands above 30 MHz. Full list of available bands is
            ['13MHz', '18MHz', '23MHz', '27MHz', '32MHz', '36MHz', '41MHz', '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', '78MHz', '82MHz']
    """
    intimestr = intime.isot[:-4].replace('-','').replace(':','').replace('T','_')
    datestr = intime.isot[:10]
    hourstr = intime.isot[11:13]
    if time_interval == '10s':
        tstr = intimestr[:-1]
    if time_interval == '1min':
        tstr = intimestr[:-2]
    if time_interval == '10min':
        tstr = intimestr[:-3]

    msfiles = []
    if lustre:
        processes=[]
        for b in bands:
            pathstr = '/lustre/pipeline/{0:s}/{1:s}/{2:s}/{3:s}/'.format(file_path, b, datestr, hourstr)
            if server:
                cmd = 'ssh ' + server + ' ls ' + pathstr + ' | grep ' + tstr
            else:
                cmd = 'ls ' + pathstr + ' | grep ' + tstr
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            processes.append(p)
            filenames = p.communicate()[0].decode('utf-8').split('\n')[:-1]
            for filename in filenames:
                if filename[-6:] == 'MHz.ms':
                    filestr = pathstr + filename
                    tmpstr = filename[:15].replace('_', 'T')
                    timestr = tmpstr[:4] + '-' + tmpstr[4:6] + '-' + tmpstr[6:11] + ':' + tmpstr[11:13] + ':' + tmpstr[13:]
                    freqstr = filename[16:21]
                    msfiles.append({'path': filestr, 'name': filename, 'time': timestr, 'freq': freqstr})
    else:
        if server:
            cmd = 'ssh ' + server + ' ls ' + file_path + ' | grep ' + tstr
        else:
            cmd = 'ls ' + file_path + ' | grep ' + tstr
        p = subprocess.run(cmd, capture_output=True, shell=True)
        filenames = p.stdout.decode('utf-8').split('\n')[:-1]
        for filename in filenames:
            if filename[-6:] == 'MHz.ms':
                if server:
                    pathstr = '{0:s}:{1:s}/{2:s}'.format(server, file_path, filename)
                else:
                    pathstr = '{0:s}/{1:s}'.format(file_path, filename)
                tmpstr = filename[:15].replace('_', 'T')
                timestr = tmpstr[:4] + '-' + tmpstr[4:6] + '-' + tmpstr[6:11] + ':' + tmpstr[11:13] + ':' + tmpstr[13:]
                freqstr = filename[16:21]
                msfiles.append({'path': pathstr, 'name': filename, 'time': timestr, 'freq': freqstr})
    return msfiles

def download_msfiles_cmd(msfile_path, server, destination):
    if server:
        p = subprocess.Popen(shlex.split('rsync -az --numeric-ids --info=progress2 --no-perms --no-owner --no-group {0:s}:{1:s} {2:s}'.format(server, msfile_path, destination)))
    else:
        p = subprocess.Popen(shlex.split('rsync -az --numeric-ids --info=progress2 --no-perms --no-owner --no-group {0:s} {1:s}'.format(msfile_path, destination)))
    std_out, std_err = p.communicate()
    if std_err:
        print(std_err)

def download_msfiles(msfiles, destination='/fast/bin.chen/realtime_pipeline/slow_working/', bands=None, verbose=True, server=None, maxthread=6):
    from multiprocessing.pool import ThreadPool
    """
    Parallelized downloading for msfiles returned from list_msfiles() to a destination.
    """
    inmsfiles_path = [f['path'] for f in msfiles]
    inmsfiles_name = [f['name'] for f in msfiles]
    inmsfiles_band = [f['freq'] for f in msfiles]
    omsfiles_path = []
    omsfiles_name = []
    if bands is None:
        omsfiles_path = inmsfiles_path
        omsfiles_name = inmsfiles_name
    else:
        for bd in bands:
            if bd in inmsfiles_band:
                idx = inmsfiles_band.index(bd)
                #omsfiles_server.append(inmsfiles_server[idx])
                omsfiles_path.append(inmsfiles_path[idx])
                omsfiles_name.append(inmsfiles_name[idx])

    nfile = len(omsfiles_path)
    if nfile == 0:
        print('No files to download. Abort...')
        return -1
    time_bg = timeit.default_timer() 
    if verbose:
        print('I am going to download {0:d} files'.format(nfile))

    tp = ThreadPool(maxthread)
    for omsfile_path in omsfiles_path:
        tp.apply_async(download_msfiles_cmd, args=(omsfile_path, server, destination))

    tp.close()
    tp.join()

    time_completed = timeit.default_timer() 
    if verbose:
        print('Downloading {0:d} files took in {1:.1f} s'.format(nfile, time_completed-time_bg))
    omsfiles = [destination + n for n in omsfiles_name]
    return omsfiles


def download_timerange(starttime, endtime, download_interval='1min', destination='/fast/bin.chen/20231027/slow/', 
                server=None, lustre=True, file_path='slow', bands=None, verbose=True, maxthread=5):
    '''
    :param download_interval: If str should either be 10s, 1min or 10min. If integer should be in seconds.
    '''
    time_bg = timeit.default_timer() 
    t_start = Time(starttime)
    t_end = Time(endtime)
    print('Start time: ', t_start.isot)
    print('End time: ', t_end.isot)
    if not os.path.exists(destination):
        os.makedirs(destination)

    if download_interval == '10s':
        dt = TimeDelta(10., format='sec')
    elif download_interval == '1min':
        dt = TimeDelta(60., format='sec')
    elif download_interval == '10min':
        dt = TimeDelta(600., format='sec')
    else:
        try:
            download_interval=int(download_interval)
            if download_interval%10!=0:
                logging.warning("Data is recorded with 10s cadence. Separation should be a multiple of 10."+\
                                "Setting download interval to nearest multiple of 10")
            dt=TimeDelta(download_interval,format='sec')
        except Exception as e:
            logging.error("download interval should either be 10s, 1min, 10min, or an integer in seconds")
            print ("download interval should either be 10s, 1min, 10min, or an integer in seconds")
            raise e
    nt = int(np.ceil((t_end - t_start) / dt))
    if isinstance(download_interval, int):
        download_interval=str(download_interval)+"s"
    print('====Will download {0:d} times at an interval of {1:s}===='.format(nt, download_interval))
    for i in range(nt):
        intime = t_start + i * dt
        msfiles = list_msfiles(intime, server=server, lustre=lustre, file_path=file_path, time_interval='10s', bands=bands)
        if verbose:
            print('Downloading time ', intime.isot)
        download_msfiles(msfiles, destination=destination, bands=bands, verbose=verbose, server=server, maxthread=maxthread)
    time_completed = timeit.default_timer() 
    if verbose:
        print('====Downloading all {0:d} times took {1:.1f} s===='.format(nt, time_completed-time_bg))


def download_calibms(calib_time, download_fold = '/lustre/bin.chen/realtime_pipeline/ms_calib/', doflag=True,
        bands=['13MHz', '18MHz', '23MHz', '27MHz', '32MHz', '36MHz', '41MHz', '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', '78MHz', '82MHz']):
    """
    Function to download calibration ms files for all or selected bands based on a given time
    :param calib_time: time selected for generating the calibration tables. A string recognized by astropy Time format
    :param download_fold: directory to hold the downloaded ms files 
    :param bands: band selection. Default to use all 16 bands.
    """
    if type(calib_time) == str:
        try:
            calib_time = Time(calib_time)
        except:
            print('The input time needs to be astropy.time.Time format')
    print(socket.gethostname(), '=======Calibrating Time {0:s}======='.format(calib_time.isot))
    ms_calib0 = list_msfiles(calib_time, file_path='slow', bands=bands)
    ms_calib = download_msfiles(ms_calib0, destination=download_fold, bands=bands)
    ms_calib.sort()
    if doflag:
        for ms_calib_ in ms_calib:
           antflagfile = flagging.flag_bad_ants(ms_calib_)
    return ms_calib


def gen_caltables(calib_in, bcaltb=None, uvrange='>10lambda', refant='202', flag_outrigger=True, 
        proc_dir='./'):
    """
    Function to generate calibration tables for a list of calibration ms files
    :param calib_in: input used for calibration. This can be either a) a string of time stamp recognized by astropy.time.Time 
            or b) a specific list of ms files used for calibration
    :param bcaltb: name of the calibration tables. Use the default if None
    :param uvrange: uv range to be used, default to '>10lambda'
    :param refant: reference antenna, default to '202'
    :param flag_outrigger: if True, flag all outrigger antennas. These would be used for beamforming.
    :proc_dir: directory to process the data and hold output files.
    """
    import pandas as pd
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
            print('The input time needs to be astropy.time.Time format. Abort...')
            return -1
        ms_calib = download_calibms(calib_time, doflag=True, download_fold=download_fold)
    else:
        print('Input not recognized. Abort...')
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
        for ms_calib_ in ms_calib:
            try:
                bcaltb = calibration.gen_calibration(ms_calib_, uvrange=uvrange, caltable_fold=caltable_fold, refant=refant)
                msmd.open(ms_calib_)
                chan_freqs.append(msmd.chanfreqs(0))
                msmd.done()
                bcaltbs.append(bcaltb)
            except Exception as e:
                print('Something is wrong when making calibrations for ', ms_calib_)
                print(e)
        chan_freqs = np.concatenate(chan_freqs)
    else:
        print('The list of calibration ms files seems to be empty. Abort...')
        return -1
        

    if flag_outrigger:
        core_ant_ids, exp_ant_ids = flagging.get_antids(ms_calib[0])
        bcaltbs_bm = []
        bmcalfac = []
        for bcaltb in bcaltbs:
            bcaltb_bm = beam_caltable_fold + '/' + os.path.basename(bcaltb)
            os.system('cp -r ' + bcaltb + ' ' + bcaltb_bm)
            tb.open(bcaltb_bm, nomodify=False)
            flags = tb.getcol('FLAG')
            npol, nch, nant = flags.shape
            for exp_ant_id in exp_ant_ids:
                flags[:, :, exp_ant_id] = True
            num_ant_per_chan = nant - np.sum(flags, axis=2)
            bmcalfac_per_chan = num_ant_per_chan ** 2.
            bmcalfac.append(bmcalfac_per_chan)
            tb.putcol('FLAG', flags)
            tb.close()
            bcaltbs_bm.append(bcaltb_bm)

        bmcalfac = np.concatenate(bmcalfac, axis=1)
        # write channel frequencies and corresponding beam scaling factors into a csv file
        df = pd.DataFrame({"chan_freqs":chan_freqs, "calfac_x":bmcalfac[0], "calfac_y":bmcalfac[1]})
        bcalfac_file = beam_caltable_fold + '/' + os.path.basename(bcaltb)[:15] + '_bmcalfac.csv'
        df.to_csv(bcalfac_file, index=False)
        return bcaltbs, bcaltbs_bm, bcalfac_file
    else:
        return bcaltbs
        
def convert_caltables_for_fast_vis(solar_ms,calib_ms,caltables):
    fast_caltables=[]
    for caltb in caltables:
        fast_caltables.append(calibration.make_fast_caltb_from_slow(calib_ms, solar_ms, caltb))
    return fast_caltables
    
def get_selfcal_table_to_apply(msname,caltable_folder):
    mstime = utils.get_time_from_name(msname)
    mstime_str = utils.get_timestr_from_name(msname)
    msfreq_str = utils.get_freqstr_from_name(msname)
    print (mstime,mstime_str,msfreq_str)

    caltables = glob.glob(os.path.join(caltable_folder,"*" + msfreq_str + "*.gcal"))
    print (caltables)
    if len(caltables) == 0:
        return []
    selfcal_time = utils.get_selfcal_time_to_apply(msname, caltables) ### Real time pipeline does not do DD cal.
                                                                    ### Hence caltables will only contain DI caltables
    print (selfcal_time)
    caltables = glob.glob(caltable_folder + "/" + selfcal_time + "*" + msfreq_str + "*.gcal")
    return caltables

def check_fast_ms(msname):
    msmd.open(msname)
    try:
        antids = msmd.antennaids()
    finally:
        msmd.done()
    num_ants=len(antids)
    if num_ants>50:
        return False
    return True
    
def run_calib(msfile, msfiles_cal=None, bcal_tables=None, do_selfcal=True, num_phase_cal=0, \
                num_apcal=1, caltable_folder=None, logger_file=None, visdir_slfcaled=None, \
                flagdir=None):
    msmd.open(msfile)
    trange = msmd.timerangeforobs(0)
    btime = qa.time(trange['begin']['m0'],form='fits')[0]
    etime = qa.time(trange['end']['m0'],form='fits')[0]
    msmd.close()
    cfreqidx = os.path.basename(msfile).find('MHz') - 2
    cfreq = os.path.basename(msfile)[cfreqidx:cfreqidx+2]+'MHz'
    msfile_cal_ = [m for m in msfiles_cal if cfreq in m]
    bcal_tables_ = [m for m in bcal_tables if cfreq in m]
    #### Generate calibrations ####
    imagename = os.path.basename(msfile)[:-3]+'_sun_selfcal'
    
    fast_vis=check_fast_ms(msfile)
    if not do_selfcal:
        gaintables = get_selfcal_table_to_apply(msfile,caltable_folder)
    if len(bcal_tables_) > 0:
        bcal_table = [bcal_tables_[0]]
        print('Found calibration table {0:s}'.format(bcal_table[0]))
        if not do_selfcal:
            for cal in gaintables:
                bcal_table.append(cal)
        
        try:
            outms, tmp = sp.image_ms_quick(msfile, calib_ms=None, bcal=bcal_table, do_selfcal=do_selfcal,\
                                        imagename=imagename, logging_level='info', \
                                        num_phase_cal=num_phase_cal, num_apcal=num_apcal,
                                        logfile=logger_file, caltable_folder=caltable_folder, \
                                        do_final_imaging=False, do_fluxscaling=False, freqbin=1, \
                                        fast_vis=fast_vis)
            os.system('cp -r '+ outms + ' ' + visdir_slfcaled + '/')
            if os.path.exists(msfile.replace('.ms', '.badants')):
                os.system('cp '+ msfile.replace('.ms', '.badants') + ' ' + flagdir + '/')
            msfile_slfcaled = visdir_slfcaled + '/' + os.path.basename(outms)
            return msfile_slfcaled
        except Exception as e:
            logging.error(e)
            return -1
    elif len(msfile_cal_) > 0:
        msfile_cal = msfile_cal_[0]
        logging.warning("No selfcal tables will be applied here. You cannot have possibly"+\
                        " have a previous selfcal table with getting bandpass")
        try:
            outms, tmp = sp.image_ms_quick(msfile, calib_ms=msfile_cal, do_selfcal=do_selfcal, imagename=imagename, logging_level='info', 
                        num_phase_cal=num_phase_cal, num_apcal=num_apcal,
                        logfile=logger_file, caltable_folder=caltable_folder, do_final_imaging=False, do_fluxscaling=False, freqbin=1)
            os.system('cp -r '+ outms + ' ' + visdir_slfcaled + '/')
            msfile_slfcaled = visdir_slfcaled + '/' + os.path.basename(outms)
            return msfile_slfcaled
        except Exception as e:
            logging.error(e)
            return -1
    else:
        print('No night time ms or caltable available for {0:s}. Skip...'.format(msfile))
        return -1


def run_imager(msfile_slfcaled, imagedir_allch=None, ephem=None, nch_out=12, stokes='I', beam_fit_size=2):
    blc = int(512 - 128)
    trc = int(512 + 128 - 1)
    region='box [ [ {0:d}pix , {1:d}pix] , [{2:d}pix, {3:d}pix ] ]'.format(blc, blc, trc, trc)
    try:
        msmd.open(msfile_slfcaled)
        trange = msmd.timerangeforobs(0)
        msmd.close()
        btime = Time(trange['begin']['m0']['value'], format='mjd')
        etime = Time(trange['end']['m0']['value'], format='mjd')
        tref_mjd = (btime.mjd + etime.mjd) / 2. 
        tref = Time(tref_mjd, format='mjd')
        tref_str = btime.isot+'~'+etime.isot
        msinfo = hf.read_msinfo(msfile_slfcaled, verbose=True)
        timeutc = me.epoch('UTC', '%fd' % tref.mjd)
        ovro = me.observatory('OVRO_MMA')
        me.doframe(ovro)
        me.doframe(timeutc)
        d0 = me.direction('SUN')
        d0_j2000 = me.measure(d0, 'J2000')
        azel = me.measure(d0, 'AZEL')
        elev = np.degrees(azel['m1']['value'])
        az = np.degrees(azel['m0']['value'])
        pb = beam(msfile_slfcaled)
        pb.srcjones(az=[az],el=[elev])
        jones_matrices = pb.get_source_pol_factors(pb.jones_matrices[0,:,:])
        sclfactor = 1. / jones_matrices[0][0]
        helio_imagename = imagedir_allch + os.path.basename(msfile_slfcaled).replace('.ms','.sun') 
        if nch_out>1:
            default_wscleancmd = ("wsclean -j 1 -mem 2 -no-reorder -no-dirty -no-update-model-required -horizon-mask 5deg -size 1024 1024 -scale 1.5arcmin -weight briggs -0.5 -minuv-l 10 -auto-threshold 3 -name " + 
                helio_imagename + " -niter 10000 -mgain 0.8 -beam-fitting-size " + str(beam_fit_size) + " -pol " + stokes + " -join-channels -channels-out " + str(nch_out) + ' ' + msfile_slfcaled)
        else:
            default_wscleancmd = ("wsclean -j 1 -mem 2 -no-reorder -no-dirty -no-update-model-required -horizon-mask 5deg -size 1024 1024 -scale 1.5arcmin -weight briggs -0.5 -minuv-l 10 -auto-threshold 3 -name " + 
                helio_imagename + " -niter 10000 -mgain 0.8 -beam-fitting-size " + str(beam_fit_size) + " -pol " + stokes + ' ' + msfile_slfcaled)
 
        os.system(default_wscleancmd)

        outfits = glob.glob(helio_imagename + '*-image.fits')
        outfits.sort()
        outfits_helio = hf.imreg(msfile_slfcaled, outfits, ephem=ephem, msinfo=msinfo, timerange=[tref_str] * len(outfits), 
                usephacenter=True, verbose=True, toTb=True, subregion=region, sclfactor=sclfactor)
        return outfits_helio
    except Exception as e:
        logging.error(e)
        #logging.error('{0:s}: Imaging for {1:s} failed'.format(socket.gethostname(), msfile_slfcaled))
        return -1


def daily_refra_correction(date, save_dir='/lustre/bin.chen/realtime_pipeline/', overwrite=True, overbright=2e6,
        dointerp=False, interp_method='linear', max_dt=600.):
    """
    Function for doing daily refraction corrections based on level 1 fits files produced in a given solar day
    :param date: format 'yyyy-mm-dd' or an astropy.time.Time object or an astropy.time.Time compatible string
    :param save_dir: directory to save the data prodcuts. Need to have a substructure of lev1/yyyy/mm/dd for level 1 files, and lev15/yyyy/mm/dd for level 1.5 files
    :param overwrite: if True, overwrite the existing level 1.5 and refraction coefficient csv file
    :param overbright: peak brightness temperature exceeding this value (in Kelvin) will be excluded for fitting
    :param interp: interpolation method used by scipy.interpolation.interp1d. Default to 'linear'
    :param max_dt: maximum time difference to perform the interpolation in seconds
    """
    if isinstance(date, str):
        try:
            date0 = Time(date)
        except:
            print("Input date not recognizable. Must be 'yyyy-mm-dd' or astropy format.")
    elif isinstance(date, Time):
        date0 = date
    else:
        print("Input date not recognizable. Must be 'yyyy-mm-dd' or astropy format.")

    # define output directories
    fits_dir = save_dir + '/fits/'
    fits_dir_lv10 = fits_dir + '/lev10/' 
    fits_dir_lv15 = fits_dir + '/lev15/' 
    hdf_dir = save_dir + '/hdf/'
    hdf_dir_lv10 = hdf_dir + '/lev1/' 
    hdf_dir_lv15 = hdf_dir + '/lev15/' 
    fig_mfs_dir = save_dir + '/figs_mfs/'
    fig_mfs_dir_lv10 = fig_mfs_dir + '/lev1/' 
    fig_mfs_dir_lv15 = fig_mfs_dir + '/lev15/' 
    fig_mfs_dir_synop = fig_mfs_dir + '/synop/' 
    refra_dir = save_dir + '/refra/'

    # Note the following reorganization is only for synoptic plots and refraction csv files
    # if UT time is before 3 UT, assign it to the earlier date. 
    datestr_synop = date0.isot.split('T')[0].replace('-','')

    datedir0 = date0.isot.split('T')[0].replace('-','/')+'/'
    date1 = date0 + TimeDelta(1., format='jd')
    datedir1 = date1.isot.split('T')[0].replace('-','/')+'/'
    fits_fch0_lv10 = glob.glob(save_dir + '/fits/lev1/' + datedir0 + '*fch*.fits')
    fits_fch0_lv10.sort()
    fits_fch1_lv10 = glob.glob(save_dir + '/fits/lev1/' + datedir1 + '*fch*.fits')
    fits_fch1_lv10.sort()
    # narrow down to all files between 13 UT and 03 UT of the second day 
    fits_fch_lv10_all_ = fits_fch0_lv10 + fits_fch1_lv10
    fits_fch_lv10_all = []


    for f in fits_fch_lv10_all_:
        datestr = os.path.basename(f).split('.')[2].split('T')[0]
        timestr = os.path.basename(f).split('.')[2].split('T')[1][:-1]
        datetimestr = datestr + 'T' + timestr[:2] + ':' + timestr[2:4] + ':' + timestr[4:]
        if Time(datetimestr).mjd > Time(date0).mjd + 13./24.  and Time(datetimestr).mjd < Time(date1).mjd + 3./24.:
            fits_fch_lv10_all.append(f)

    fits_fch_lv10_all.sort()

    refrafile = refra_dir + '/refra_coeff_' + datestr_synop + '.csv'
    if os.path.exists(refrafile):
        df = pd.read_csv(refrafile)
    else:
        print('Cannot find a refraction correction file in the designated directory. Make a new one from data.')  

    for fits_fch_lv10 in fits_fch_lv10_all:
        fits_mfs_lv10 = fits_fch_lv10.replace('fch', 'mfs')
        print('processing fits '+fits_fch_lv10)
        datestr = os.path.basename(fits_fch_lv10).split('.')[2].split('T')[0]
        timestr = os.path.basename(fits_fch_lv10).split('.')[2].split('T')[1][:-1]
        datetimestr = datestr + 'T' + timestr[:2] + ':' + timestr[2:4] + ':' + timestr[4:]
        datedir = datestr.replace('-','/') + '/'
        if not os.path.exists(fits_dir_lv15 + datedir):
            os.makedirs(fits_dir_lv15 + datedir)
        if not os.path.exists(hdf_dir_lv15 + datedir):
            os.makedirs(hdf_dir_lv15 + datedir)
        if not os.path.exists(fig_mfs_dir_lv15 + datedir):
            os.makedirs(fig_mfs_dir_lv15 + datedir)
        if not os.path.exists(fig_mfs_dir_synop + datedir):
            os.makedirs(fig_mfs_dir_synop + datedir)

        try:
            meta, data = ndfits.read(fits_fch_lv10)
            if 'df' in locals() and meta['header']['date-obs'][:19] in df.Time.values and not overwrite:
                print('Refraction correction record for '+ meta['header']['date-obs'] + ' already exists. Continue')
                continue
            else:
                refra_rec = orefr.refraction_fit_param(fits_fch_lv10, return_record=True, overbright=overbright)
                fits_mfs_lv15 = fits_dir_lv15 + datedir + os.path.basename(fits_mfs_lv10.replace('.lev1_mfs', '.lev1.5_mfs'))
                hdf_mfs_lv15 = hdf_dir_lv15 + datedir + os.path.basename(fits_mfs_lv15).replace('.fits', '.hdf')
                fits_fch_lv15 = fits_dir_lv15 + datedir + os.path.basename(fits_fch_lv10.replace('.lev1_fch', '.lev1.5_fch'))
                hdf_fch_lv15 = hdf_dir_lv15 + datedir + os.path.basename(fits_fch_lv15).replace('.fits', '.hdf')
                figname_lv15 = os.path.basename(fits_mfs_lv15).replace('.fits', '.png')
                if (not np.isnan(refra_rec['px0'])) and (not np.isnan(refra_rec['py0'])): 
                    df_new = pd.DataFrame(refra_rec, index=[0])
                    if 'df' in locals():
                        if (refra_rec['Time'] in df['Time'].unique()):
                            cols = list(df.columns)
                            df.loc[df.Time.isin(df_new.Time), cols] = df_new[cols].values
                            print('Refraction correction record for '+ refra_rec['Time'] + ' updated in '+ refrafile)
                        else:
                            df = pd.concat([df_new, df], ignore_index=True)
                            df = df.sort_values(by='Time')
                            df.to_csv(refrafile, index=False)
                            print('Refraction correction record for '+ refra_rec['Time'] + ' added to '+ refrafile)
                    else:
                        df = df_new
                        df.to_csv(refrafile, index=False)
                        print('Refraction correction record for '+ refra_rec['Time'] + ' added to '+ refrafile)

                    fits_mfs_lv15 = orefr.apply_refra_record(fits_mfs_lv10, refra_rec, fname_out=fits_mfs_lv15)
                    if fits_mfs_lv15:
                        utils.compress_fits_to_h5(fits_mfs_lv15, hdf_mfs_lv15)
                        fig, axes = ovis.slow_pipeline_default_plot(fits_mfs_lv15)
                        fig.savefig(fig_mfs_dir_lv15 + datedir + figname_lv15)
                        figname_synop = figname_lv15.replace('.lev1.5_mfs_10s.', '.synop_mfs_10s.')
                        os.system('cp '+ fig_mfs_dir_lv15 + datedir + figname_lv15 + ' ' + fig_mfs_dir_synop + datedir + figname_synop)

                        fits_fch_lv15 = orefr.apply_refra_record(fits_fch_lv10, refra_rec, fname_out=fits_fch_lv15)
                        utils.compress_fits_to_h5(fits_fch_lv15, hdf_fch_lv15)

                else:
                    print('Refraction correction failed for '+ datetimestr)
                    if dointerp:
                        print('Trying to interpolate from nearby times')
                        fits_mfs_lv15 = orefr.apply_refra_record(fits_mfs_lv10, df, fname_out=fits_mfs_lv15, interp=interp_method, max_dt=max_dt)
                        if fits_mfs_lv15:
                            print('Succeeded and updating level 1.5 files, but be cautious!')
                            utils.compress_fits_to_h5(fits_mfs_lv15, hdf_mfs_lv15)
                            fig, axes = ovis.slow_pipeline_default_plot(fits_mfs_lv15)
                            fig.savefig(fig_mfs_dir_lv15 + datedir + figname_lv15)
                            figname_synop = figname_lv15.replace('.lev1.5_mfs_10s.', '.synop_mfs_10s.')
                            os.system('cp '+ fig_mfs_dir_lv15 + datedir + figname_lv15 + ' ' + fig_mfs_dir_synop + datedir + figname_synop)
                            fits_fch_lv15 = orefr.apply_refra_record(fits_fch_lv10, df, fname_out=fits_fch_lv15, interp=interp_method, max_dt=max_dt)
                            utils.compress_fits_to_h5(fits_fch_lv15, hdf_fch_lv15)
                        else:
                            print('Interpolation failed')
                            continue
        except Exception as e:
            logging.error(e)
            logging.error('====Processing {0:s} failed'.format(fits_fch_lv10))



def pipeline_quick(image_time=Time.now() - TimeDelta(20., format='sec'), server=None, lustre=True, file_path='slow', 
            min_nband=6, nch_out=12, beam_fit_size=2, stokes='I', do_selfcal=True, num_phase_cal=0, num_apcal=1, 
            overwrite_ms=False, delete_ms_slfcaled=False, logger_file=None, compress_fits=True,
            proc_dir = '/fast/bin.chen/realtime_pipeline/',
            save_dir = '/lustre/bin.chen/realtime_pipeline/',
            calib_dir = '/lustre/bin.chen/realtime_pipeline/caltables/',
            calib_file = '20240117_145752',
            delete_working_ms=True, delete_working_fits=True, do_refra=True, overbright=2e6,
            do_imaging=True,
            bands = ['32MHz', '36MHz', '41MHz', '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', '78MHz', '82MHz']):
    """
    Pipeline for processing and imaging slow visibility data
    :param time_start: start time of the visibility data to be processed
    :param time_end: end time of the visibility data to be processed
    :param image_interval: time interval between adjacent visibilities to be processed
    :param server: server name on which the data is stored
    :param lustre: if True, specific to lustre system on lwacalim nodes. If not, try your luck in combination with file_path
    :param file_path: path to the data w.r.t. the server
    :param min_nband: minimum number of bands to be processed. Will skip if less than that.
    :param calib_file: calibration file to be used. Format yyyymmdd_hhmmss
    :param nch_out: number of channels to be imaged
    :param beam_fit_size: size of the beam area used for fitting to be passed to wsclean. See https://wsclean.readthedocs.io/en/v3.0/restoring_beam_size.html
    :param do_selfcal: if True, do selfcalibration
    :param num_phase_cal: number of phase calibration iterations
    :param num_apcal: number of amplitude calibration iterations
    :param overwrite_ms: if True, overwrite the ms files in the working directory
    :param delete_ms_slfcaled: if True, delete the selfcalibrated ms files after imaging
    :param logger_file: name of the log file
    :param compress_fits: if True, compress the fits (not loseless, reduce bitlen, but it's OK to do this in most of the cases of solar imaging) files after imaging
    :param proc_dir: directory to hold the working files
    :param save_dir: directory to hold the final products
    :param calib_dir: directory to hold the initial bandpass calibration tables
    :param calib_file: calibration file to be used. Format yyyymmdd_hhmmss
    :param delete_working_ms: if True, delete the working ms files after imaging (set False for debugging purpose)
    :param overbright: peak brightness temperature exceeding this value (in Kelvin) will be excluded for refraction correction fitting
    """

    time_begin = timeit.default_timer() 
    
    if file_path[-1]=='/':
        file_path=file_path[:-1]

    # caltable_folder is where the initial bandpass calibration tables are located 
    caltable_folder = calib_dir
    # gaintable_folder is where the intermediate gain tables are located
    gaintable_folder = proc_dir + '/caltables/'
    visdir_calib = proc_dir + '/'+file_path+'_slow_calib/'
    visdir_work = proc_dir + '/'+file_path+'_working/'
    visdir_slfcaled = proc_dir + '/'+file_path+'_slfcaled/'
    imagedir_allch = proc_dir + '/'+file_path+'_images_allch/'
    flagdir = save_dir + '/flags/'
    refradir = save_dir + '/refra/'

    imagedir_allch_combined = save_dir + '/'+file_path+'_fits/'
    hdf_dir = save_dir + '/'+ file_path+'_hdf/'
    fig_mfs_dir = save_dir + '/'+ file_path+'_figs_mfs/'

    ## Night-time MS files used for calibration ##
    msfiles_cal = glob.glob(visdir_calib + calib_file + '_*MHz.ms')
    msfiles_cal.sort()

    bcal_tables = glob.glob(caltable_folder + calib_file + '_*MHz.bcal')
    bcal_tables.sort()

    if not os.path.exists(visdir_work):
        os.makedirs(visdir_work)

    if not os.path.exists(visdir_slfcaled):
        os.makedirs(visdir_slfcaled)

    if not os.path.exists(gaintable_folder):
        os.makedirs(gaintable_folder)

    if not os.path.exists(flagdir):
        os.makedirs(flagdir)

    if not os.path.exists(refradir):
        os.makedirs(refradir)

    if not os.path.exists(imagedir_allch):
        os.makedirs(imagedir_allch)

    if not os.path.exists(imagedir_allch_combined):
        os.makedirs(imagedir_allch_combined)

    if not os.path.exists(hdf_dir):
        os.makedirs(hdf_dir)

    if not os.path.exists(fig_mfs_dir):
        os.makedirs(fig_mfs_dir)

    try:
        print(socket.gethostname(), '=======Processing Time {0:s}======='.format(image_time.isot))
        #logging.info('=======Processing Time {0:s}======='.format(image_time.isot))
        msfiles0 = list_msfiles(image_time, lustre=lustre, server=server, file_path=file_path)
        if len(msfiles0) < min(min_nband, len(bands)):
            print('This time only has {0:d} subbands. Check nearby +-10s time.'.format(len(msfiles0)))
            image_time_before = image_time - TimeDelta(10., format='sec')
            msfiles0_before = list_msfiles(image_time_before, lustre=lustre, server=server, file_path=file_path)
            image_time_after = image_time + TimeDelta(10., format='sec')
            msfiles0_after = list_msfiles(image_time_after, lustre=lustre, server=server, file_path=file_path)
            if len(msfiles0_before) < min(min_nband, len(bands)) and len(msfiles0_after) < min(min_nband, len(bands)):
                print('I cannot find a nearby time with at least {0:d} available subbands. Abort and wait for next time interval.'.format(min_nband))
                return False
            else:
                if len(msfiles0_before) > len(msfiles0_after):
                    msfiles0 = msfiles0_before
                    image_time = image_time_before
                else:
                    msfiles0 = msfiles0_after
                    image_time = image_time_after
        
        msfiles0_freq = [f['freq'] for f in msfiles0]
        msfiles0_name = [f['name'] for f in msfiles0]
        timestr = msfiles0_name[0][:15]
        
        prev_calfiles=glob.glob(os.path.join(gaintable_folder,"*.gcal"))#### these files will be deleted in this cycle
        
        msfiles_slfcaled = glob.glob(visdir_slfcaled + '/' + timestr + '_*MHz*.ms')
        msfiles_slfcaled.sort()
        if len(msfiles_slfcaled) == 0 or overwrite_ms:
            #msfiles0 = glob.glob(datadir_orig + timestr + '_*MHz.ms')
            #msfiles0.sort()
            # skip the first two bands (18-32 MHz)
            # msfiles0 = msfiles0[2:]

            #### copy files over to the working directory ####
            print('==Copying file over to working directory==')
            logging.debug('====Copying file over to working directory====')
            time1 = timeit.default_timer()
            msfiles = download_msfiles(msfiles0, destination=visdir_work, bands=bands)
            time2 = timeit.default_timer()
            logging.debug('Time taken to copy files is {0:.1f} s'.format(time2-time1))

            fitsfiles=[]
            msfiles_slfcaled = []

            # parallelized calibration, selfcalibration, and source subtraction
            logging.debug('Starting to calibrate all {0:d} bands'.format(len(msfiles)))
            time_cal1 = timeit.default_timer()
            pool = multiprocessing.pool.Pool(processes=len(msfiles))
            #result = pool.map_async(run_calib, msfiles)
            run_calib_partial = partial(run_calib, msfiles_cal=msfiles_cal, bcal_tables=bcal_tables, do_selfcal=do_selfcal, num_phase_cal=num_phase_cal, num_apcal=num_apcal, 
                    logger_file=logger_file, caltable_folder=gaintable_folder, visdir_slfcaled=visdir_slfcaled, flagdir=flagdir)
            result = pool.map_async(run_calib_partial, msfiles)
            timeout = 1200.
            result.wait(timeout=timeout)
            #if result.ready():
            #    time_cal2 = timeit.default_timer()
            #    logging.debug('Calibration for all {0:d} bands is done in {1:.1f} s'.format(len(msfiles), time_cal2-time_cal1))
            #else:
            #    logging.debug('Calibration for certain bands is incomplete in {0:.1f} s'.format(timeout))
            #    logging.debug('Proceed anyway')
            try:
                msfiles_slfcaled = result.get(timeout)
                time_cal2 = timeit.default_timer()
                logging.debug('Calibration for all {0:d} bands is done in {1:.1f} s'.format(len(msfiles), time_cal2-time_cal1))
                pool.close()
                pool.join()
            except TimeoutError as e:
                logging.error(e)
                logging.debug('Calibration for certain bands is incomplete in {0:.1f} s'.format(timeout))
                logging.debug('Proceed anyway')
                pool.terminate()
                pool.close()
            if delete_working_ms:
                for file1 in msfiles0_name:
                    timestr1 = utils.get_timestr_from_name(file1)
                    freqstr=utils.get_freqstr_from_name(file1)
                    os.system('rm -rf '+ visdir_work + '/' + timestr1 + '_*'+freqstr+'*')
                
        else:
            logging.debug('=====Selfcalibrated ms already exist for {0:s}. Proceed with imaging.========'.format(timestr)) 
            if delete_working_ms:
                for file1 in msfiles0_name:
                    timestr1 = utils.get_timestr_from_name(file1)
                    freqstr=utils.get_freqstr_from_name(file1)
                    os.system('rm -rf '+ visdir_work + '/' + timestr1 + '_*'+freqstr+'*')
                


        # Do imaging
        print('======= processed selfcaled ms files =====')
        success = [type(m) is str for m in msfiles_slfcaled]
        msfiles_slfcaled_success = []
        for m in msfiles_slfcaled:
            if type(m) is str:
                msfiles_slfcaled_success.append(m)
        if sum(success) > min(4,len(bands)-1):
            logging.info('{0:s}: Successfuly selfcalibrated {1:d} out of {2:d} bands'.format(socket.gethostname(), len(msfiles_slfcaled_success), len(bands)))
            time_img1 = timeit.default_timer()
            for i, m in enumerate(msfiles_slfcaled_success):
                try:
                    msmd.open(m)
                    trange = msmd.timerangeforobs(0)
                    msmd.close()
                    break
                except Exception as e:
                    if i < len(msfiles_slfcaled_success): 
                        logging.error('Reading file {0:s} has error {1:s}. Will try the next one'.format(m, e))
                        continue
                    else:
                        logging.error('Nothing seems to work. I will abort and continue to the next time')
                        for file1 in msfiles0_name:
                            timestr1 = utils.get_timestr_from_name(file1)
                            freqstr=utils.get_freqstr_from_name(file1)
                            os.system('rm -rf '+ visdir_slfcaled + '/' + timestr1 + '_*'+freqstr+'*.ms')
                            if do_selfcal:
                                os.system('rm -rf '+ gaintable_folder + '/' + timestr1 + '_*'+freqstr+'*')
                        return False
            fast_vis=check_fast_ms(msfiles_slfcaled[0])

            if do_imaging:
                if fast_vis:
                    nch_out=1
                fitsfiles=image_times(msfiles_slfcaled,imagedir_allch, nch_out=nch_out, \
                                   stokes=stokes, beam_fit_size=beam_fit_size)
            btime = Time(trange['begin']['m0']['value'], format='mjd')

        else:
            logging.error('For time {0:s}, less than 4 bands out of {1:d} bands were calibrated successfully. Abort....'.format(timestr, len(bands)))
            for file1 in msfiles0_name:
                timestr1 = utils.get_timestr_from_name(file1)
                freqstr=utils.get_freqstr_from_name(file1)
                os.system('rm -rf '+ visdir_slfcaled + '/' + timestr1 + '_*'+freqstr+'*.ms')
                if do_selfcal:
                    os.system('rm -rf '+ gaintable_folder + '/' + timestr1 + '_*'+freqstr+'*')
                os.system('rm -rf '+ visdir_work + '/' + timestr1 + '_*'+freqstr+'*')
            return False

        if delete_ms_slfcaled:
            new_caltables={}
            for file1 in msfiles0_name:
                timestr1 = utils.get_timestr_from_name(file1)
                freqstr=utils.get_freqstr_from_name(file1)
                os.system('rm -rf '+ visdir_slfcaled + '/' + timestr1 + '_*'+freqstr+'*.ms')
                os.system('rm -rf '+ visdir_work + '/' + timestr1 + '_*'+freqstr+'*')
                new_caltables[freqstr]=len(glob.glob(os.path.join(gaintable_folder,timestr1+"_"+freqstr+"*.gcal"))) 
                               
            for calfile in prev_calfiles:
                freqstr=utils.get_freqstr_from_name(calfile)
                try:
                    if new_caltables[freqstr]!=0 and do_selfcal:
                        os.system('rm -rf '+calfile) ### I am only deleting the previous calfile and
                                                 ### and that also when this round has exited successfully.
                except KeyError:
                    pass
        
        if 'fitsfiles' in locals():
            if (len(fitsfiles[0]) > 1 and not fast_vis) or (fast_vis and len(fitsfiles[0]) >= 1):
                datedir = btime.isot[:10].replace('-','/')+'/'
                
                # Note the following reorganization is only for synoptic plots and refraction csv files
                # if UT time is before 4 UT, assign it to the earlier date. 
    
                date_mjd = int(btime.mjd)
                if btime.mjd - date_mjd < 4./24.:
                    datestr_synop = Time(btime.mjd - 1., format='mjd').isot[:10].replace('-','')
                    datedir_synop = Time(btime.mjd - 1., format='mjd').isot[:10].replace('-','/') + '/'
                else:
                    datestr_synop = Time(btime.mjd, format='mjd').isot[:10].replace('-','')
                    datedir_synop = Time(btime.mjd, format='mjd').isot[:10].replace('-','/') + '/'
                    
                fig_mfs_dir_sub_synop = fig_mfs_dir + '/synop/' + datedir_synop
                
                if not os.path.exists(fig_mfs_dir_sub_synop):
                    os.makedirs(fig_mfs_dir_sub_synop)
                    
                print ("fast_vis compression:",fast_vis)
                fits_images,plotted_image=compress_plot_images(fitsfiles, btime, datedir, imagedir_allch_combined, hdf_dir, \
                                        fig_mfs_dir, stokes, fast_vis=fast_vis)
                
                logging.info("Level 1 images plotted ok")
                figname_to_copy=None
                
                if not fast_vis:    
                    if len(fits_images)==2:  #### 0  is mfs, 1 in fine channel
                        if do_refra:
                            refrafile = refradir + '/refra_coeff_' + datestr_synop + '.csv'
                            logging.info("Trying to do refraction correction")

                            refra_image,success=do_refraction_correction(fits_images, overbright, \
                                                        refrafile, datedir, imagedir_allch_combined, hdf_dir, \
                                                        fig_mfs_dir,btime)
                            if not success:
                                logging.info('Refraction correction failed for '+ btime.isot)
                                figname_to_copy=plotted_image
                                figname_synop = os.path.basename(plotted_image).replace('.lev1_mfs_10s.', '.synop_mfs_10s.')    
                            else:
                                figname_to_copy=refra_image
                                figname_synop = os.path.basename(refra_image).replace('.lev1.5_mfs_10s.', '.synop_mfs_10s.')
                      
                    if not figname_to_copy:
                        figname_to_copy=plotted_image
                        figname_synop = os.path.basename(plotted_image).replace('.lev1', '.synop')             
                    
                    synoptic_image=os.path.join(fig_mfs_dir_sub_synop, figname_synop)    
                    
                    os.system('cp '+ figname_to_copy + ' ' + synoptic_image)   
                
                if delete_working_fits:
                    os.system('rm -rf '+imagedir_allch + '*')
                time_completed= timeit.default_timer() 
                logging.debug('====All processing for time {0:s} is done in {1:.1f} minutes'.format(timestr, (time_completed-time_begin)/60.))
                return True
            else:
                time_exit = timeit.default_timer()
                logging.error('====Processing for time {0:s} failed in {1:.1f} minutes'.format(timestr, (time_exit-time_begin)/60.))
                return False
        elif do_imaging:
            time_exit = timeit.default_timer()
            logging.error('====Processing for time {0:s} failed in {1:.1f} minutes'.format(timestr, (time_exit-time_begin)/60.))
            return False
        else:
            time_completed= timeit.default_timer() 
            logging.debug('====All processing for time {0:s} is done in {1:.1f} minutes'.format(timestr, (time_completed-time_begin)/60.))
            return True
    except Exception as e:
        logging.error(e)
        time_exit = timeit.default_timer()
        logging.error('====Processing for time {0:s} failed in {1:.1f} minutes'.format(timestr, (time_exit-time_begin)/60.))
        return False


def image_times(msfiles_slfcaled, imagedir_allch, nch_out=12, stokes='I', beam_fit_size=2):
    msfiles_slfcaled_success = []
    for m in msfiles_slfcaled:
        if type(m) is str:
            msfiles_slfcaled_success.append(m)
    
    time_img1 = timeit.default_timer()
    for i, m in enumerate(msfiles_slfcaled_success):
        msmd.open(m)
        trange = msmd.timerangeforobs(0)
        msmd.close()
        break
                
    btime = Time(trange['begin']['m0']['value'], format='mjd')
    etime = Time(trange['end']['m0']['value'], format='mjd')
    tref_mjd = (btime.mjd + etime.mjd) / 2. 
    tref = Time(tref_mjd, format='mjd')
    ephem = hf.read_horizons(tref, dur=1./60./24., observatory='OVRO_MMA')
    pool = multiprocessing.pool.Pool(processes=len(msfiles_slfcaled_success))
    run_imager_partial = partial(run_imager, imagedir_allch=imagedir_allch, ephem=ephem, \
                nch_out=nch_out, stokes=stokes, beam_fit_size=beam_fit_size)
    results = [pool.apply_async(run_imager_partial, (file1,)) for file1 in msfiles_slfcaled_success]
    timeout = 200.

    fitsfiles=[]
    for res in results:
        try:
            fitsfiles.append(res.get(timeout))
            time_img2 = timeit.default_timer()
            logging.debug('Imaging for all {0:d} bands is done in {1:.1f} s'.format(len(msfiles_slfcaled_success), time_img2-time_img1))
            pool.close()
            
            
        except TimeoutError as e:
            logging.error(e)
            logging.debug('Imaging for certain bands is incomplete in {0:.1f} s'.format(timeout))
            logging.debug('Proceed anyway')
            
    pool.join()
    pool.terminate()

    return fitsfiles

def compress_plot_images(fitsfiles, starttime, datedir, imagedir_allch_combined, \
                            hdf_dir, fig_mfs_dir, stokes, fast_vis=False):    
    ## define subdirectories for storing the fits and png files
    
    btime=starttime
    imagedir_allch_combined_sub_lv10 = imagedir_allch_combined + '/lev1/' + datedir
    hdf_dir_sub_lv10 = hdf_dir + '/lev1/' + datedir
    fig_mfs_dir_sub_lv10 = fig_mfs_dir + '/lev1/' + datedir

    if not os.path.exists(imagedir_allch_combined_sub_lv10):
       os.makedirs(imagedir_allch_combined_sub_lv10)
    if not os.path.exists(hdf_dir_sub_lv10):
       os.makedirs(hdf_dir_sub_lv10)
    
    if not os.path.exists(fig_mfs_dir_sub_lv10):
        os.makedirs(fig_mfs_dir_sub_lv10)
    

    ## Wrap images
    timestr_iso = btime.isot[:-4].replace(':','')+'Z'
    
    
    
    
    # multi-frequency synthesis images
    fits_mfs = imagedir_allch_combined_sub_lv10 + '/ovro-lwa.lev1_mfs_10s.' + timestr_iso + '.image_'+stokes+'.fits' 
    #fitsfiles_mfs = glob.glob(imagedir_allch + '/' + timestr+ '*MFS-image.fits')
    fitsfiles_mfs = []
    for f in fitsfiles:

        if type(f) is list:
            if 'MFS' in f[-1] and (not fast_vis):
                fitsfiles_mfs.append(f[-1])
            elif fast_vis:
                fitsfiles_mfs+=f
             
        else:
            continue
    
    fitsfiles_mfs.sort()
    
    ndfits.wrap(fitsfiles_mfs, outfitsfile=fits_mfs)
    
    hdf_mfs = hdf_dir_sub_lv10 + os.path.basename(fits_mfs).replace('.fits', '.hdf')
    #utils.compress_fits_to_h5(fits_mfs, hdf_mfs)
    
    if not fast_vis:
        # fine channel spectral images
        fits_fch = imagedir_allch_combined_sub_lv10 + '/ovro-lwa.lev1_fch_10s.' + timestr_iso + '.image_'+stokes+'.fits' 
        fitsfiles_fch = []
        for f in fitsfiles:
            if type(f) is list:
                fitsfiles_fch += f[:-1]
            else:
                continue
        fitsfiles_fch.sort()
        ndfits.wrap(fitsfiles_fch, outfitsfile=fits_fch)
        hdf_fch = hdf_dir_sub_lv10 + os.path.basename(fits_fch).replace('.fits', '.hdf')
        #utils.compress_fits_to_h5(fits_fch, hdf_fch)
        
        
    
    
    fig, axes = ovis.slow_pipeline_default_plot(fits_mfs)
    figname_lv10 = os.path.basename(fits_mfs).replace('.fits', '.png')
    fig.savefig(fig_mfs_dir_sub_lv10 + '/' + figname_lv10)
    
    if not fast_vis:
        return [fits_mfs,fits_fch], \
                    os.path.join(fig_mfs_dir_sub_lv10, figname_lv10)
    else:
        return [fits_mfs], os.path.join(fig_mfs_dir_sub_lv10, figname_lv10)
    
def do_refraction_correction(fitsfiles, overbright, refrafile, datedir, imagedir_allch_combined, hdf_dir, \
                            fig_mfs_dir, image_time):
    btime=image_time                        
    imagedir_allch_combined_sub_lv15 = imagedir_allch_combined + '/lev15/' + datedir
    hdf_dir_sub_lv15 = hdf_dir + '/lev15/' + datedir
    fig_mfs_dir_sub_lv15 = fig_mfs_dir + '/lev15/' + datedir
    
    if not os.path.exists(fig_mfs_dir_sub_lv15):
        os.makedirs(fig_mfs_dir_sub_lv15)
        
    if not os.path.exists(imagedir_allch_combined_sub_lv15):
        os.makedirs(imagedir_allch_combined_sub_lv15)
    
    if not os.path.exists(hdf_dir_sub_lv15):
       os.makedirs(hdf_dir_sub_lv15)
   
    print ("directories made")
    fits_mfs,fits_fch=fitsfiles

    px, py = orefr.refraction_fit_param(fits_fch, overbright=overbright)

    if (not np.isnan(px).any()) and (not np.isnan(py).any()): 
        df_new = pd.DataFrame({"Time":btime.isot, "px0":px[0], "px1":px[1], "py0":py[0], "py1":py[1]}, index=[0])
        print (refrafile)
        if os.path.exists(refrafile):
            df = pd.read_csv(refrafile)
            if not (btime.isot in df['Time'].unique()):
                df = pd.concat([df_new, df], ignore_index=True)
                df = df.sort_values(by='Time')
                df.to_csv(refrafile, index=False)
                logging.info('Refraction correction record for '+ btime.isot + ' added to '+ refrafile)
            else:
                logging.info('Refraction correction record for '+ btime.isot + ' already exists. Do nothing.')
        else:
            df = df_new
            df.to_csv(refrafile, index=False)
            logging.info('Refraction correction record for '+ btime.isot + ' added to '+ refrafile)

        fits_mfs_lv15 = imagedir_allch_combined_sub_lv15 + os.path.basename(fits_mfs.replace('.lev1_mfs', '.lev1.5_mfs'))
        fits_mfs_lv15 = orefr.apply_refra_coeff(fits_mfs, px, py, fname_out=fits_mfs_lv15)
        hdf_mfs_lv15 = hdf_dir_sub_lv15 + os.path.basename(fits_mfs_lv15).replace('.fits', '.hdf')
        utils.compress_fits_to_h5(fits_mfs_lv15, hdf_mfs_lv15)

        fits_fch_lv15 = imagedir_allch_combined_sub_lv15 + os.path.basename(fits_fch.replace('.lev1_fch', '.lev1.5_fch'))
        fits_fch_lv15 = orefr.apply_refra_coeff(fits_fch, px, py, fname_out=fits_fch_lv15)
        hdf_fch_lv15 = hdf_dir_sub_lv15 + os.path.basename(fits_fch_lv15).replace('.fits', '.hdf')
        utils.compress_fits_to_h5(fits_fch_lv15, hdf_fch_lv15)

        fig, axes = ovis.slow_pipeline_default_plot(fits_mfs_lv15)
        figname_lv15 = os.path.basename(fits_mfs_lv15).replace('.fits', '.png')
        fig.savefig(fig_mfs_dir_sub_lv15 + '/' + figname_lv15)
        
        return os.path.join(fig_mfs_dir_sub_lv15,figname_lv15), True
      
    else:
        return None, False
        
            
    
        

def run_pipeline(time_start=Time.now(), time_end=None, time_interval=600., delay_from_now=180., do_selfcal=True, num_phase_cal=0, num_apcal=1, 
        server=None, lustre=True, file_path='slow', multinode=True, nodes='0123456789', delete_ms_slfcaled=True, 
        logger_dir = '/lustre/bin.chen/realtime_pipeline/logs/', logger_prefix='solar_realtime_pipeline', logger_level=20,
        proc_dir = '/fast/bin.chen/realtime_pipeline/',
        save_dir = '/lustre/bin.chen/realtime_pipeline/',
        calib_dir = '/lustre/bin.chen/realtime_pipeline/caltables/',
        calib_file = '20240117_145752', altitude_limit=15., 
        beam_fit_size = 2, 
        delete_working_ms=True, do_refra=True, delete_working_fits=True,
        do_imaging=True,
        bands = ['32MHz', '36MHz', '41MHz', '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', '78MHz', '82MHz']):
    '''
    Main routine to run the pipeline. Note each time stamp takes about 8.5 minutes to complete.
    "time_interval" needs to be set to something greater than that. 600 is recommended.
    :param time_start: time for starting the pipeline. astropy.time.Time object.
    :param time_end: time for ending the pipeline. astropy.time.Time object. 
                If not specified (default), it will never end until being killed manually.
    :param time_interval: interval between adjacent processing times in seconds for each session
    :param server: server name on which the data is stored
    :param lustre: if True, specific to lustre system on lwacalim nodes. If not, try your luck in combination with file_path
    :param file_path: path to the data w.r.t. the server
    :param delay_from_now: delay of the newest time to process compared to now.
    :param delete_ms_slfcaled: whether or not to delete the self-calibrated measurement sets.
    :param multinode: if True, will delay the start time by the node
    :param nodes: list of lwacalim nodes to be used. Default '0123456789'
    :param calib_file: calibration file to be used. Format yyyymmdd_hhmmss
    :param altitude_limit: lowest altitude to start the pipeline in degrees. Default to 15 deg.
    :param beam_fit_size: size of the beam area used for fitting to be passed to wsclean. See https://wsclean.readthedocs.io/en/v3.0/restoring_beam_size.html
    :param do_imaging: If True, imaging and other related steps will be done. Default: True
    '''
    try:
        time_start = Time(time_start)
    except Exception as e:
        logging.error(e)
        raise e
    
    if file_path[-1]=='/':
        file_path=file_path[:-1]
        
    (t_rise, t_set) = sun_riseset(time_start, altitude_limit=altitude_limit)
    # set up logging file
    server = socket.gethostname()
    datestr = Time(time_start.mjd, format='mjd').isot[:10].replace('-','')
    datedir = Time(time_start.mjd, format='mjd').isot[:10].replace('-','/') + '/'
    
    logger_file = logger_dir + datedir + logger_prefix + '_' + file_path+ '_'+ datestr + '_' + server + '.log'  

    if not os.path.exists(os.path.dirname(logger_file)):
        print('Path to logger file {0:s} does not exist. Attempting to create the directory tree.'.format(logger_file))
        os.makedirs(os.path.dirname(logger_file))

    logging.basicConfig(filename=logger_file, filemode='at',
        format='%(asctime)s %(funcName)s %(lineno)d %(levelname)-8s %(message)s',
        level=logger_level,
        datefmt='%Y-%m-%d %H:%M:%S', force=True)

    logging.info('{0:s}: I am asked to start imaging for {1:s}'.format(socket.gethostname(), time_start.isot))
    if multinode:
        nodenum = int(socket.gethostname()[-2:])
        nodes_list=[int(n) for n in list(nodes)]
        nnode = len(nodes_list)
        delay_by_node = nodes_list.index(nodenum) * (time_interval/nnode) 
    else:
        logging.info('{0:s}: I am running on a single node'.format(socket.gethostname()))
        nodenum = int(socket.gethostname()[-2:])
        nodes_list = [nodenum]
        delay_by_node = 0. 
    #while time_start > t_rise and time_start < Time.now() - TimeDelta(15.,format='sec'): 
    # find out when the Sun is high enough in the sky
    if time_start < t_rise:
        twait = t_rise - time_start
        logging.info('{0:s}: Start time {1:s} is before sunrise. Wait for {2:.1f} hours to start.'.format(socket.gethostname(), time_start.isot, twait.value * 24.))
        time_start += TimeDelta(twait.sec + 60., format='sec')
        sleep(twait.sec + 60.)
    else:
        logging.info("{0:s}: Start time {1:s} is after today's sunrise at {2:s}. Will try to proceed.".format(socket.gethostname(), time_start.isot, t_rise.isot))
    time_start += TimeDelta(delay_by_node, format='sec')
    logging.info('{0:s}: Delay {1:.1f} min to {2:s}'.format(socket.gethostname(), delay_by_node / 60., time_start.isot))
    sleep(delay_by_node)
    while True:
        if time_end:
            if time_start > Time(time_end):
                logging.info('The new imaging time now passes the provided end time. Ending the pipeline.'.format(Time(time_start).isot, Time(time_end).isot))
                print('The new imaging time now passes the provided end time. Ending the pipeline.'.format(Time(time_start).isot, Time(time_end).isot))
                break
        time1 = timeit.default_timer()
        if time_start > Time.now() - TimeDelta(delay_from_now, format='sec'):
            twait = time_start - Time.now()
            logging.info('{0:s}: Start time {1:s} is too close to current time. Wait {2:.1f} m to start.'.format(socket.gethostname(), time_start.isot, (twait.sec + delay_from_now) / 60.))
            sleep(twait.sec + delay_from_now)
        logging.info('{0:s}: Start processing {1:s}'.format(socket.gethostname(), time_start.isot))
        res = pipeline_quick(time_start, do_selfcal=do_selfcal, num_phase_cal=num_phase_cal, num_apcal=num_apcal, \
                            server=server, lustre=lustre, file_path=file_path, delete_ms_slfcaled=delete_ms_slfcaled,\
                            logger_file=logger_file, proc_dir=proc_dir, save_dir=save_dir, calib_dir=calib_dir, \
                            calib_file=calib_file, delete_working_ms=delete_working_ms,\
                            delete_working_fits=delete_working_fits, do_refra=do_refra, \
                            beam_fit_size=beam_fit_size, do_imaging=do_imaging, bands=bands)
        time2 = timeit.default_timer()
        if res:
            logging.info('{0:s}: Processing {1:s} was successful within {2:.1f}m'.format(socket.gethostname(), time_start.isot, (time2-time1)/60.))
        else:
            logging.info('{0:s}: Processing {1:s} was unsuccessful!!!'.format(socket.gethostname(), time_start.isot))

        if (time_interval - (time2-time1)) < 0:
            logging.info('{0:s}: Warning!! Processing {1:s} took {2:.1f}m to complete. This node may be falling behind'.format(socket.gethostname(), time_start.isot, (time2-time1)/60.))

        time_start += TimeDelta(time_interval, format='sec')

        if time_start > t_set:
            (t_rise_next, t_set_next) = sun_riseset(t_set + TimeDelta(6./24., format='jd'))
            date_mjd = int(time_start.mjd)
            if time_start.mjd - date_mjd < 4./24.:
                date_synop = Time(time_start.mjd - 1., format='mjd').isot[:10]
            else:
                date_synop = Time(time_start.mjd, format='mjd').isot[:10]
            
            if int(socket.gethostname()[-2:]) == nodes_list[-1]:
                logging.info('{0:s}: Sun is setting. Done for the day. Doing refraction corrections for the full day.'.format(socket.gethostname())) 
                daily_refra_correction(date_synop, save_dir=save_dir, overwrite=False, dointerp=True, interp_method='linear', max_dt=600.)
                twait = t_rise_next - Time.now() 
                logging.info('{0:s}: Refraction corrections done. Wait for {1:.1f} hours to start.'.format(socket.gethostname(), twait.value * 24.)) 
            else:
                twait = t_rise_next - Time.now() 
                logging.info('{0:s}: Sun is setting. Done for the day. Wait for {1:.1f} hours to start.'.format(socket.gethostname(), twait.value * 24.)) 

            time_start += TimeDelta(twait.sec + 60. + delay_by_node, format='sec')
            t_rise = t_rise_next
            t_set = t_set_next
            # updating the logger file
            datestr = Time(t_rise.mjd, format='mjd').isot[:10].replace('-','')
            datedir = Time(t_rise.mjd, format='mjd').isot[:10].replace('-','/') + '/'
            logger_file = logger_dir + datedir + logger_prefix + '_' + datestr + '_' + server + '.log'  

            if not os.path.exists(os.path.dirname(logger_file)):
                print('Path to logger file {0:s} does not exist. Attempting to create the directory tree.'.format(logger_file))
                os.makedirs(os.path.dirname(logger_file))

            logging.basicConfig(filename=logger_file, filemode='at',
                format='%(asctime)s %(funcName)s %(lineno)d %(levelname)-8s %(message)s',
                level=logger_level,
                datefmt='%Y-%m-%d %H:%M:%S', force=True)
            sleep(twait.sec + 60. + delay_by_node)



    



