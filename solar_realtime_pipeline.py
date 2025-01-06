#!/usr/bin/env python
import os, sys, glob, getopt
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
from ovrolwasolar import solar_pipeline as sp
from ovrolwasolar.primary_beam import analytic_beam as beam
from ovrolwasolar import utils, calibration, flagging
from casatasks import clearcal, applycal, flagdata, tclean, exportfits, imsubimage, split
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
import socket,glob
from matplotlib.patches import Ellipse
import argparse
import pandas as pd
import resource
import platform

from ovrolwasolar import visualization as ovis
from ovrolwasolar import refraction_correction as orefr

#import gc # garbage collection
#gc.enable()

matplotlib.use('agg')

msmd = msmetadata()
qa = quanta()
me = measures()
tb = table()

def get_memory():
    with open('/proc/meminfo', 'r') as mem:
        free_memory = 0
        for i in mem:
            sline = i.split()
            if str(sline[0]) in ('MemFree:', 'Buffers:', 'Cached:'):
                free_memory += int(sline[1])
    return free_memory

def set_memory_limit(percentage=0.1):
    """
    Only works in Linux systems
    """
    if platform.system() != "Linux":
        print('Only works on linux!')
        return
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (int(get_memory() * 1024 * percentage), hard))


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
            if len(filenames) > 0:
                for filename in filenames:
                    if filename[-6:] == 'MHz.ms':
                        filestr = pathstr + filename
                        tmpstr = filename[:15].replace('_', 'T')
                        timestr = tmpstr[:4] + '-' + tmpstr[4:6] + '-' + tmpstr[6:11] + ':' + tmpstr[11:13] + ':' + tmpstr[13:]
                        freqstr = filename[16:21]
                        msfiles.append({'path': filestr, 'name': filename, 'time': timestr, 'freq': freqstr})
            else:
                logging.info('Did not find any files at the given time {0:s}.'.format(intime.isot)) 
    else:
        if server:
            cmd = 'ssh ' + server + ' ls ' + file_path + ' | grep ' + tstr
        else:
            cmd = 'ls ' + file_path + ' | grep ' + tstr
        p = subprocess.run(cmd, capture_output=True, shell=True)
        filenames = p.stdout.decode('utf-8').split('\n')[:-1]
        if len(filenames) > 0:
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
        else:
            logging.info('Did not find any files at the given time {0:s}.'.format(intime.isot)) 
            msfiles = []
    return msfiles

def download_msfiles_cmd(msfile_path, server, destination):
    if server:
        p = subprocess.Popen(shlex.split('rsync -az --numeric-ids --info=progress2 --no-perms --no-owner --no-group {0:s}:{1:s} {2:s}'.format(server, msfile_path, destination)))
    else:
        #p = subprocess.Popen(shlex.split('rsync -az --numeric-ids --info=progress2 --no-perms --no-owner --no-group {0:s} {1:s}'.format(msfile_path, destination)))
        p = subprocess.Popen(shlex.split('cp -r {0:s} {1:s}'.format(msfile_path, destination)))
    std_out, std_err = p.communicate()
    if std_err:
        print('<<',Time.now().isot,'>>',std_err)

def download_msfiles(msfiles, destination='/fast/solarpipe/realtime_pipeline/slow_working/', bands=None, verbose=True, server=None, maxthread=3):
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
        print('<<',Time.now().isot,'>>','No files to download. Abort...')
        return -1
    time_bg = timeit.default_timer() 
    if verbose:
        print('<<',Time.now().isot,'>>','I am going to download {0:d} files'.format(nfile))

    tp = ThreadPool(maxthread)
    for omsfile_path in omsfiles_path:
        tp.apply_async(download_msfiles_cmd, args=(omsfile_path, server, destination))

    tp.close()
    tp.join()

    time_completed = timeit.default_timer() 
    if verbose:
        print('<<',Time.now().isot,'>>','Downloading {0:d} files took in {1:.1f} s'.format(nfile, time_completed-time_bg))
    omsfiles = [destination + n for n in omsfiles_name]
    return omsfiles


def download_timerange(starttime, endtime, download_interval='1min', destination='/fast/solarpipe/20231027/slow/', 
                server=None, lustre=True, file_path='slow', bands=None, verbose=True, maxthread=5):
    '''
    :param download_interval: If str should either be 10s, 1min or 10min. If integer should be in seconds.
    '''
    time_bg = timeit.default_timer() 
    t_start = Time(starttime)
    t_end = Time(endtime)
    print('<<',Time.now().isot,'>>','Start time: ', t_start.isot)
    print('<<',Time.now().isot,'>>','End time: ', t_end.isot)
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


def download_calibms(calib_time, download_fold = '/lustre/solarpipe/realtime_pipeline/ms_calib/', doflag=True,
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
            try:
               antflagfile = flagging.flag_bad_ants(ms_calib_)
            except:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                logging.error(exc_type, fname, exc_tb.tb_lineno)
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
            print('<<',Time.now().isot,'>>','The input time needs to be astropy.time.Time format. Abort...')
            return -1
        ms_calib = download_calibms(calib_time, doflag=True, download_fold=download_fold)
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
    else:
        print('<<',Time.now().isot,'>>','The list of calibration ms files seems to be empty. Abort...')
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

    caltables = glob.glob(os.path.join(caltable_folder,"*" + msfreq_str + "*.gcal"))
    if len(caltables) == 0:
        return []
    selfcal_time = utils.get_selfcal_time_to_apply(msname, caltables) ### Real time pipeline does not do DD cal.
                                                                    ### Hence caltables will only contain DI caltables
    caltables = glob.glob(caltable_folder + "/" + selfcal_time + "*" + msfreq_str + "*.gcal")
    return caltables

def get_allsky_image_to_use(msname, img_folder):
    mstime = utils.get_time_from_name(msname)
    mstime_str = utils.get_timestr_from_name(msname)
    msfreq_str = utils.get_freqstr_from_name(msname)
    
    imgs = glob.glob(os.path.join(img_folder,"*" + msfreq_str + "*allsky-image.fits"))
    if len(imgs) == 0:
        return []
    img_time = utils.get_selfcal_time_to_apply(msname, imgs) ## function matches time str only. Hence can be used here
    imgs = glob.glob(img_folder + "/" + img_time + "*" + msfreq_str + "*allsky-image.fits")
    return imgs
    

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

    
def run_calib(msfile, msfiles_cal=None, bcal_tables=None, do_selfcal=True, num_phase_cal=0, 
                num_apcal=1, caltable_folder=None, logger_file=None, visdir_slfcaled=None, 
                flagdir=None, delete_allsky=False, actively_rm_ms=True):
    
    try:
        msmd.open(msfile)
    except Exception as e:
        logging.error(e)
        return -1

    # do time average if the input ms file is fast visibility
    if check_fast_ms(msfile):
        omsfile = os.path.dirname(msfile) + '/' + os.path.basename(msfile).replace('.ms', '.10s.ms')
        split(msfile, omsfile, datacolumn='data', timebin='10s')
        os.system('rm -rf ' + msfile)
        os.system('mv ' + omsfile + ' ' + msfile)

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
        print('<<',Time.now().isot,'>>','Found calibration table {0:s}'.format(bcal_table[0]))
        if not do_selfcal:
            for cal in gaintables:
                bcal_table.append(cal)
        
        msfile_cal = None
        
        fast_vis_image_model_subtraction=False
        sky_image=None
        prev_allsky_img=get_allsky_image_to_use(msfile,os.path.dirname(msfile).replace('fast_working','slow_working'))
        if len(prev_allsky_img)!=0 and fast_vis:
            fast_vis_image_model_subtraction=True
            sky_image=prev_allsky_img[0].split('-image')[0]
            logging.debug("will use "+sky_image)
        
        try:
            outms, tmp = sp.image_ms_quick(msfile, calib_ms=None, bcal=bcal_table, do_selfcal=do_selfcal,\
                                        imagename=imagename, logging_level='info', \
                                        num_phase_cal=num_phase_cal, num_apcal=num_apcal,
                                        logfile=logger_file, caltable_folder=caltable_folder, \
                                        do_final_imaging=False, do_fluxscaling=False, freqbin=1, \
                                        fast_vis=fast_vis, delete_allsky=delete_allsky,\
                                        fast_vis_image_model_subtraction=fast_vis_image_model_subtraction,\
                                        sky_image=sky_image)
            if actively_rm_ms:
                os.system('rm -rf ' + msfile)
            os.system('cp -r '+ outms + ' ' + visdir_slfcaled + '/')
            if actively_rm_ms:
                os.system('rm -rf ' + outms)
            if os.path.exists(msfile.replace('.ms', '.badants')):
                os.system('cp '+ msfile.replace('.ms', '.badants') + ' ' + flagdir + '/')
            msfile_slfcaled = visdir_slfcaled + '/' + os.path.basename(outms)
            return msfile_slfcaled
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            logging.error(exc_type, fname, exc_tb.tb_lineno)
            return -1
    elif len(msfile_cal_) > 0:
        msfile_cal = msfile_cal_[0]
        logging.warning("No selfcal tables will be applied here. You cannot have possibly"+\
                        " have a previous selfcal table with getting bandpass")
        try:
            outms, tmp = sp.image_ms_quick(msfile, calib_ms=msfile_cal, do_selfcal=do_selfcal, imagename=imagename, logging_level='info', 
                        num_phase_cal=num_phase_cal, num_apcal=num_apcal,
                        logfile=logger_file, caltable_folder=caltable_folder, do_final_imaging=False, 
                        do_fluxscaling=False, freqbin=1, delete_allsky=delete_allsky)
            os.system('cp -r '+ outms + ' ' + visdir_slfcaled + '/')
            msfile_slfcaled = visdir_slfcaled + '/' + os.path.basename(outms)
            return msfile_slfcaled
        except Exception as e:
            logging.error(e)
            return -1
    else:
        print('<<',Time.now().isot,'>>','No night time ms or caltable available for {0:s}. Skip...'.format(msfile))
        return -1


def run_imager(msfile_slfcaled, imagedir_allch=None, ephem=None, nch_out=12, stokes='I', beam_fit_size=2, briggs=-0.5):
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
        default_wscleancmd = "wsclean -j 2 -mem 4 -quiet -no-reorder -no-dirty -no-update-model-required -horizon-mask 5deg -size 1024 1024 -scale 1.5arcmin -weight briggs " + str(briggs) + " -minuv-l 10 -auto-threshold 3 -name " + helio_imagename + " -niter 10000 -mgain 0.8 -beam-fitting-size " + str(beam_fit_size) + " -pol " + stokes

        if nch_out>1:
            # default to be used for slow visibility imaging for fine channel imaging
            default_wscleancmd += " -join-channels -channels-out " + str(nch_out)
        
        cmd= shlex.split(default_wscleancmd + ' ' + msfile_slfcaled)
        wsclean_proc=subprocess.run(cmd)

        outfits = glob.glob(helio_imagename + '*-image.fits')
        outfits.sort()
        if len(outfits) > 0:
            outfits_helio = hf.imreg(msfile_slfcaled, outfits, ephem=ephem, msinfo=msinfo, timerange=[tref_str] * len(outfits), 
                    usephacenter=True, verbose=True, toTb=True, subregion=region, sclfactor=sclfactor)
            return outfits_helio
        else:
            logging.error('No fits images produced.')
            return -1
    except Exception as e:
        logging.error(e)
        if wsclean_proc.poll() is None:
            wsclean_proc.terminate()
        return -1


def daily_refra_correction(date, save_dir='/lustre/solarpipe/realtime_pipeline/', overwrite=True, overbright=2e6,
        dointerp=False, interp_method='linear', max_dt=600., slowfast='slow'):
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
            print('<<',Time.now().isot,'>>',"Input date not recognizable. Must be 'yyyy-mm-dd' or astropy format.")
    elif isinstance(date, Time):
        date0 = date
    else:
        print('<<',Time.now().isot,'>>',"Input date not recognizable. Must be 'yyyy-mm-dd' or astropy format.")

    # define output directories
    fits_dir = save_dir + '/fits/' + slowfast + '/'
    fits_dir_lv10 = fits_dir + '/lev1/' 
    fits_dir_lv15 = fits_dir + '/lev15/' 
    hdf_dir = save_dir + '/hdf/' + slowfast + '/'
    hdf_dir_lv10 = hdf_dir + '/lev1/' 
    hdf_dir_lv15 = hdf_dir + '/lev15/' 
    fig_mfs_dir = save_dir + '/figs_mfs/' + slowfast + '/'
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
    fits_fch0_lv10 = glob.glob(fits_dir_lv10 + datedir0 + '*fch*.fits')
    fits_fch0_lv10.sort()
    fits_fch1_lv10 = glob.glob(fits_dir_lv10 + datedir1 + '*fch*.fits')
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
        print('<<',Time.now().isot,'>>','Cannot find a refraction correction file in the designated directory. Make a new one from data.')  

    for fits_fch_lv10 in fits_fch_lv10_all:
        fits_mfs_lv10 = fits_fch_lv10.replace('fch', 'mfs')
        print('<<',Time.now().isot,'>>','processing fits '+fits_fch_lv10)
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
            if 'df' in locals() and meta['header']['date-obs'][:19] in df.Time.values and not overwrite and slowfast == 'slow':
                print('<<',Time.now().isot,'>>','Refraction correction record for '+ meta['header']['date-obs'] + ' already exists. Continue')
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
                            print('<<',Time.now().isot,'>>','Refraction correction record for '+ refra_rec['Time'] + ' updated in '+ refrafile)
                        else:
                            df = pd.concat([df_new, df], ignore_index=True)
                            df = df.sort_values(by='Time')
                            df.to_csv(refrafile, index=False)
                            print('<<',Time.now().isot,'>>','Refraction correction record for '+ refra_rec['Time'] + ' added to '+ refrafile)
                    else:
                        df = df_new
                        df.to_csv(refrafile, index=False)
                        print('<<',Time.now().isot,'>>','Refraction correction record for '+ refra_rec['Time'] + ' added to '+ refrafile)

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
                    print('<<',Time.now().isot,'>>','Refraction correction failed for '+ datetimestr)
                    if dointerp:
                        print('<<',Time.now().isot,'>>','Trying to interpolate from nearby times')
                        fits_mfs_lv15 = orefr.apply_refra_record(fits_mfs_lv10, df, fname_out=fits_mfs_lv15, interp=interp_method, max_dt=max_dt)
                        if fits_mfs_lv15:
                            print('<<',Time.now().isot,'>>','Succeeded and updating level 1.5 files, but be cautious!')
                            utils.compress_fits_to_h5(fits_mfs_lv15, hdf_mfs_lv15)
                            fig, axes = ovis.slow_pipeline_default_plot(fits_mfs_lv15)
                            fig.savefig(fig_mfs_dir_lv15 + datedir + figname_lv15)
                            figname_synop = figname_lv15.replace('.lev1.5_mfs_10s.', '.synop_mfs_10s.')
                            os.system('cp '+ fig_mfs_dir_lv15 + datedir + figname_lv15 + ' ' + fig_mfs_dir_synop + datedir + figname_synop)
                            fits_fch_lv15 = orefr.apply_refra_record(fits_fch_lv10, df, fname_out=fits_fch_lv15, interp=interp_method, max_dt=max_dt)
                            utils.compress_fits_to_h5(fits_fch_lv15, hdf_fch_lv15)
                        else:
                            print('<<',Time.now().isot,'>>','Interpolation failed')
                            continue
        except Exception as e:
            logging.error(e)
            logging.error('====Processing {0:s} failed'.format(fits_fch_lv10))



def pipeline_quick(image_time=Time.now() - TimeDelta(20., format='sec'), server=None, lustre=True, file_path='slow', 
            min_nband=6, nch_out=12, beam_fit_size=2, briggs=-0.5, stokes='I', do_selfcal=True, num_phase_cal=0, num_apcal=1, 
            overwrite_ms=False, delete_ms_slfcaled=False, logger_file=None, compress_fits=True,
            proc_dir_mem = '/dev/shm/srtmp/', proc_dir = '/fast/solarpipe/realtime_pipeline/',
            save_dir = '/lustre/solarpipe/realtime_pipeline/',
            calib_dir = '/lustre/solarpipe/realtime_pipeline/caltables/',
            calib_file = '20240117_145752',
            delete_working_ms=True, delete_working_fits=True, do_refra=True, overbright=2e6, save_selfcaltab=False,
            slowfast='slow', do_imaging=True, delete_allsky=False, save_allsky=False,
            bands = ['32MHz', '36MHz', '41MHz', '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', '78MHz', '82MHz'],
            clear_old_files=True, clear_older_than=45, actively_rm_ms=True):
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
    :param briggs: briggs weighting parameter to be passed to wsclean. -1 is close to uniform and 1 is close to natural. See https://wsclean.readthedocs.io/en/latest/image_weighting.html
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
    :param slowfast: specify whether slow or fast visibilities are being processed
    :param delete_allsky: if True, delete the allsky image after each run. Otherwise keep the latest frame. 
    :param save_allsky: if True, save the allsky image FITS file into save_dir + 'allsky/'
    """

    time_begin = timeit.default_timer() 

    set_memory_limit()

    if slowfast.lower() != 'slow' and slowfast.lower() != 'fast':
        print("slowfast needs to be either 'slow' or 'fast'. Abort")
        return False

    if lustre:
        print('I am working with the default data archive on the lustre server.')
        if slowfast.lower() == 'slow':
            file_path = 'slow'
        elif slowfast.lower() == 'fast':
            file_path = 'fast'
        else:
            print("slowfast needs to be either 'slow' or 'fast'. Abort") 
            return False
    else:
        print('<<',Time.now().isot,'>>','I am not working with the default data archive on the lustre server. This is not fully tested. Good luck!')
    
    # caltable_folder is where the initial bandpass calibration tables are located 
    caltable_folder = calib_dir
    # gaintable_folder is where the intermediate gain tables are located
    gaintable_folder = proc_dir_mem + '/caltables/'
    visdir_calib = proc_dir_mem + '/slow_calib/'
    visdir_work = proc_dir_mem + '/' + slowfast + '_working/'
    visdir_slfcaled = proc_dir_mem + '/' + slowfast + '_slfcaled/'
    imagedir_allch = proc_dir + '/'+ slowfast +'_images_allch/'
    flagdir = save_dir + '/flags/'
    refradir = save_dir + '/refra/'
    caltabarchive = save_dir + '/caltablearchive/'

    imagedir_allch_combined = save_dir + '/fits/' + slowfast + '/'
    hdf_dir = save_dir + '/hdf/' + slowfast + '/'
    fig_mfs_dir = save_dir + '/figs_mfs/' + slowfast + '/'
    allsky_dir = save_dir + '/allsky/' 
    allsky_dir_fits = allsky_dir + 'fits/'
    allsky_dir_figs = allsky_dir + 'figs/'

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

    if not os.path.exists(allsky_dir_fits):
        os.makedirs(allsky_dir_fits)

    if not os.path.exists(allsky_dir_figs):
        os.makedirs(allsky_dir_figs)

    if not os.path.exists(caltabarchive):
        os.makedirs(caltabarchive)

    if clear_old_files:
        remove_old_items(visdir_work, clear_older_than)
        remove_old_items(visdir_slfcaled, clear_older_than)
        #remove_old_items(visdir_calib, clear_older_than)
        remove_old_items(imagedir_allch, clear_older_than)

    try:
        print('<<',Time.now().isot,'>>')
        print(socket.gethostname(), '=======Processing Time {0:s}======='.format(image_time.isot))
        #logging.info('=======Processing Time {0:s}======='.format(image_time.isot))
        msfiles0 = list_msfiles(image_time, lustre=lustre, server=server, file_path=file_path, time_interval='10s')
        
        if len(msfiles0) < len(bands):
            # try to find missing times from nearby times that are +-4 s within the reference time
            msfiles0_before = list_msfiles(image_time - TimeDelta(10., format='sec'), lustre=lustre, server=server, file_path=file_path, time_interval='10s')
            if len(msfiles0_before) > 0:
                msfiles0_before_ts = [f['time'] for f in msfiles0_before]
            else:
                msfiles0_before_ts = []
            msfiles0_after = list_msfiles(image_time + TimeDelta(10., format='sec'), lustre=lustre, server=server, file_path=file_path, time_interval='10s')
            if len(msfiles0_after) > 0:
                msfiles0_after_ts = [f['time'] for f in msfiles0_after]
            else:
                msfiles0_after_ts = []
            msfiles0_all = msfiles0_before + msfiles0 + msfiles0_after
            
            if len(msfiles0_all) > 0:
                msfiles0_ts = [f['time'] for f in msfiles0]
                msfiles0_tref = Time(np.median(Time(msfiles0_ts).mjd), format='mjd')
                msfiles0_all_ts = msfiles0_before_ts + msfiles0_ts + msfiles0_after_ts
                idxs, = np.where(np.abs((Time(msfiles0_all_ts).mjd - msfiles0_tref.mjd)) < 4./86400.)
                msfiles0 = [msfiles0_all[idx] for idx in idxs]
            else:
                logging.info('I cannot find any available subbands. Abort and wait for next time interval.')
                return False


        # check if the currently requested time has enough number of bands
        if len(msfiles0) < min(min_nband, len(bands)):
            #print('This time only has {0:d} subbands. Check nearby +-10s time.'.format(len(msfiles0)))
            #if slowfast.lower()=='slow':
            #    image_time_before = image_time - TimeDelta(10., format='sec')
            #    msfiles0_before = list_msfiles(image_time_before, lustre=lustre, server=server, file_path=file_path)
            #    image_time_after = image_time + TimeDelta(10., format='sec')
            #    msfiles0_after = list_msfiles(image_time_after, lustre=lustre, server=server, file_path=file_path)
            #    if len(msfiles0_before) < min(min_nband, len(bands)) and len(msfiles0_after) < min(min_nband, len(bands)):
            #        print('I cannot find a nearby time with at least {0:d} available subbands. Abort and wait for next time interval.'.format(min_nband))
            #        return False
            #    else:
            #        if len(msfiles0_before) > len(msfiles0_after):
            #            msfiles0 = msfiles0_before
            #            image_time = image_time_before
            #        else:
            #            msfiles0 = msfiles0_after
            #            image_time = image_time_after
            #else:
            #    print('I cannot find a nearby time with at least {0:d} available subbands. Abort and wait for next time interval.'.format(min_nband))
            #    return False
            logging.info('I cannot find at least {0:d} available subbands. Abort and wait for next time interval.'.format(min_nband))
            return False

        
        msfiles0_freq = [f['freq'] for f in msfiles0]
        msfiles0_name = [f['name'] for f in msfiles0]
        timestr = msfiles0_name[0][:15]
        logging.info('====Processing {0:s}===='.format(timestr))
        
        prev_calfiles=glob.glob(os.path.join(gaintable_folder,"*.gcal"))#### these files will be deleted in this cycle
        prev_allsky_imgs = glob.glob(os.path.join(visdir_work,'*_allsky-image.fits'))
        
                
        if save_selfcaltab and slowfast.lower()=='slow':

            if len(prev_calfiles) > 0:
                logging.info('Saving gain tables to '+caltabarchive)
                logging.info('prev_files : '+str(prev_calfiles))
                for g in prev_calfiles:
                    
                    date_caltable = g.split('/')[-1].split('_')[0]
                    time_caltable = g.split('/')[-1].split('_')[1]
                    gaintables_sub = caltabarchive + '/' + date_caltable + '/' + date_caltable + "_" + time_caltable +'/'
                    if not os.path.exists(gaintables_sub):
                        os.makedirs(gaintables_sub)
                
                    os.system('cp -r ' + g + ' ' + gaintables_sub)
            else:
                logging.info('No gain tables found for {0:s}'.format(timestr))     
        else:
            logging.info('No gain tables saved for {0:s}'.format(timestr))
            logging.info(str(save_selfcaltab)+ ' '+str(slowfast.lower()=='slow'))
            print('<<',Time.now().isot,'>>','No gain tables saved for {0:s}'.format(timestr))


        msfiles_slfcaled = glob.glob(visdir_slfcaled + '/' + timestr + '_*MHz*.ms')
        msfiles_slfcaled.sort()
        if len(msfiles_slfcaled) == 0 or overwrite_ms:
            #msfiles0 = glob.glob(datadir_orig + timestr + '_*MHz.ms')
            #msfiles0.sort()
            # skip the first two bands (18-32 MHz)
            # msfiles0 = msfiles0[2:]

            #### copy files over to the working directory ####
            print('<<',Time.now().isot,'>>','==Copying file over to working directory==')
            logging.debug('====Copying file over to working directory====')
            time1 = timeit.default_timer()
            msfiles = download_msfiles(msfiles0, destination=visdir_work, bands=bands)
            time2 = timeit.default_timer()
            logging.debug('Time taken to copy files is {0:.1f} s'.format(time2-time1))


            
            # parallelized calibration, selfcalibration, and source subtraction
            logging.debug('Starting to calibrate all {0:d} bands'.format(len(msfiles)))
            time_cal1 = timeit.default_timer()

            #result = pool.map_async(run_calib, msfiles)
            run_calib_partial = partial(run_calib, msfiles_cal=msfiles_cal, bcal_tables=bcal_tables, do_selfcal=do_selfcal, 
                    num_phase_cal=num_phase_cal, num_apcal=num_apcal, logger_file=logger_file, caltable_folder=gaintable_folder, 
                    visdir_slfcaled=visdir_slfcaled, flagdir=flagdir, delete_allsky=delete_allsky, actively_rm_ms=actively_rm_ms)

            if slowfast.lower()=='slow':
                timeout = 800.
            else:
                timeout = 300.
            
            msfiles_slfcaled=parallel_task_runner(run_calib_partial,msfiles,timeout=timeout)
            if len(msfiles_slfcaled)==len(msfiles):
                time_cal2 = timeit.default_timer()
                logging.debug('Calibration for all {0:d} bands is done in {1:.1f} s'.format(len(msfiles), time_cal2-time_cal1))
            else:
                logging.debug('Calibration for certain bands is incomplete in {0:.1f} s'.format(timeout))
                logging.debug('Proceed anyway')
                
            
            allsky_fitsfiles=[]
            for file1 in msfiles0_name:
                timestr1 = utils.get_timestr_from_name(file1)
                freqstr=utils.get_freqstr_from_name(file1)
                if save_allsky and slowfast.lower()=='slow':
                    allsky_dir_fits_sub = allsky_dir_fits + '/'+ timestr1[0:4]+'/'+timestr1[4:6]+'/'+timestr1[6:8]+'/'
                    if not os.path.exists(allsky_dir_fits_sub):
                        os.makedirs(allsky_dir_fits_sub)
                    timestr_out = timestr1[:4]+'-'+timestr1[4:6]+'-'+timestr1[6:8]+'T'+timestr1[9:15]
                    allsky_fitsfile_search = glob.glob(visdir_work + '/' + timestr1 + '_' + freqstr + '*_allsky-image.fits')
                    if len(allsky_fitsfile_search) == 1:
                        logging.debug('All sky image exists. Move it to '+allsky_dir_fits_sub)
                        allsky_fitsfile0 = allsky_fitsfile_search[0]
                        allsky_fitsfile = 'ovro-lwa-352.allsky_'+freqstr+'_10s.'+timestr_out + '.image_I.fits'
                        cpcommand = 'cp ' + allsky_fitsfile0 + ' ' + allsky_dir_fits_sub + allsky_fitsfile
                        os.system(cpcommand)
                        allsky_fitsfiles.append(allsky_dir_fits_sub + allsky_fitsfile)
                    else:
                        logging.info('All sky image {0:s} does not exist.'.format(\
                                        visdir_work + '/' + timestr1 + '_' + freqstr + '*_allsky-image.fits'))


                if delete_working_ms:
                    os.system('rm -rf '+ visdir_work + '/' + timestr1 + '_*'+freqstr+'*.ms*')
                    os.system('rm -rf '+ visdir_work + '/' + timestr1 + '_*'+freqstr+'*_self0*')
                    os.system('rm -rf '+ visdir_work + '/' + timestr1 + '_*'+freqstr+'*.cl')

            if len(allsky_fitsfiles) > 4:
                logging.info('Making allsky plots')
                timestr_ = os.path.basename(allsky_fitsfiles[0]).split('.')[2].split('-')
                allsky_dir_figs_sub = allsky_dir_figs + '/'+ timestr_[0] + '/' + timestr_[1] + '/' + timestr_[2][:2]+'/'
                if not os.path.exists(allsky_dir_figs_sub):
                    logging.debug('Putting allsky plots in '+allsky_dir_figs_sub)
                    os.makedirs(allsky_dir_figs_sub)
                fig, axes = ovis.make_allsky_image_plots(allsky_fitsfiles)
                figname_pre = os.path.basename(allsky_fitsfiles[0]).split('_')[0]
                figname_allsky = figname_pre + os.path.basename(allsky_fitsfiles[0])[len(figname_pre)+6:].replace('.fits', '.png')
                logging.info('Saving allsky plots to ' + allsky_dir_figs_sub + '/' + figname_allsky)
                fig.savefig(allsky_dir_figs_sub + '/' + figname_allsky)
                
        else:
            logging.debug('=====Selfcalibrated ms already exist for {0:s}. Proceed with imaging.========'.format(timestr)) 
            if delete_working_ms:
                for file1 in msfiles0_name:
                    timestr1 = utils.get_timestr_from_name(file1)
                    freqstr=utils.get_freqstr_from_name(file1)
                    os.system('rm -rf '+ visdir_work + '/' + timestr1 + '_*'+freqstr+'*.ms*')
                    os.system('rm -rf '+ visdir_work + '/' + timestr1 + '_*'+freqstr+'*_self0*')
                    os.system('rm -rf '+ visdir_work + '/' + timestr1 + '_*'+freqstr+'*.cl')

        # Do imaging
        print('<<',Time.now().isot,'>>','======= processed selfcaled ms files =====')
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
                            os.system('rm -rf '+ visdir_slfcaled + '/' + timestr1 + '_*'+freqstr+'*')
                            if do_selfcal:
                                os.system('rm -rf '+ gaintable_folder + '/' + timestr1 + '_*'+freqstr+'*')
                        return False
           
            fast_vis=check_fast_ms(msfiles_slfcaled[0])


            if do_imaging:
                #if fast_vis:
                #    nch_out=1
                fitsfiles=image_times(msfiles_slfcaled,imagedir_allch, nch_out=nch_out, \
                                   stokes=stokes, beam_fit_size=beam_fit_size, briggs=briggs)
            btime = Time(trange['begin']['m0']['value'], format='mjd')


        else:
            logging.error('For time {0:s}, less than 4 bands out of {1:d} bands were calibrated successfully. Abort....'.format(timestr, len(bands)))
            for file1 in msfiles0_name:
                timestr1 = utils.get_timestr_from_name(file1)
                freqstr=utils.get_freqstr_from_name(file1)
                os.system('rm -rf '+ visdir_slfcaled + '/' + timestr1 + '_*'+freqstr+'*.ms')
                if do_selfcal:
                    os.system('rm -rf '+ gaintable_folder + '/' + timestr1 + '_*'+freqstr+'*')
                os.system('rm -rf '+ visdir_work + '/' + timestr1 + '_*'+freqstr+'*.cl')
                os.system('rm -rf '+ visdir_work + '/' + timestr1 + '_*'+freqstr+'*.ms')
                os.system('rm -rf '+ visdir_work + '/' + timestr1 + '_*'+freqstr+'*.fits')
            return False

        new_caltables={}
        new_allsky_imgs={}
        for file1 in msfiles0_name:
            timestr1 = utils.get_timestr_from_name(file1)
            freqstr=utils.get_freqstr_from_name(file1)
            if delete_ms_slfcaled:
                os.system('rm -rf '+ visdir_slfcaled + '/' + timestr1 + '_*'+freqstr+'*.ms')
            new_caltables[freqstr]=len(glob.glob(os.path.join(gaintable_folder,timestr1+"_"+freqstr+"*.gcal"))) 
            new_allsky_imgs[freqstr]=len(glob.glob(os.path.join(visdir_work,timestr1+"_"+freqstr+"*allsky-image.fits")))
                           
        for calfile in prev_calfiles:
            freqstr=utils.get_freqstr_from_name(calfile)
            try:
                if new_caltables[freqstr]!=0 and do_selfcal:
                    os.system('rm -rf '+calfile) 
                    os.system('rm -rf '+calfile + '.fast') ### I am only deleting the previous calfile and
                                             ### and that also when this round has exited successfully.
            except KeyError:
                pass
        
        for imgfile in prev_allsky_imgs:
            freqstr=utils.get_freqstr_from_name(imgfile)
            try:
                if new_allsky_imgs[freqstr]!=0 and do_selfcal:
                    os.system('rm -rf '+imgfile.replace("-image.fits","*"))
            except KeyError:
                pass
        
        if 'fitsfiles' in locals():
            if len(fitsfiles) > 0:
                if len(fitsfiles[0]) > 1:
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
                    
                     
                    fits_images, plotted_image = compress_plot_images(fitsfiles, btime, datedir, imagedir_allch_combined, hdf_dir, \
                                            fig_mfs_dir, stokes, fast_vis=fast_vis)
                    
                    logging.info("Level 1 images plotted ok")
                    figname_to_copy=None
                    
                    # Do refraction corrections
                    if not fast_vis:    
                        if len(fits_images)==2:  #### 0  is mfs, 1 in fine channel
                            if do_refra:
                                refrafile = refradir + '/refra_coeff_' + datestr_synop + '.csv'
                                logging.info("Trying to do refraction correction")

                                refra_image, success=do_refraction_correction(fits_images, overbright, \
                                                            refrafile, datedir, imagedir_allch_combined, hdf_dir, \
                                                            fig_mfs_dir,btime)
                                if not success:
                                    logging.info('Refraction correction failed for '+ btime.isot)
                                    figname_to_copy=plotted_image
                                    figname_synop = os.path.basename(plotted_image).replace('.lev1_mfs', '.synop_mfs')    
                                else:
                                    figname_to_copy=refra_image
                                    figname_synop = os.path.basename(refra_image).replace('.lev1.5_mfs', '.synop_mfs')
                          
                        if not figname_to_copy:
                            figname_to_copy = plotted_image
                            figname_synop = os.path.basename(plotted_image).replace('.lev1', '.synop')             
                    else:
                        figname_to_copy = plotted_image
                        figname_synop = os.path.basename(plotted_image).replace('.lev1', '.synop')             

                    synoptic_image=os.path.join(fig_mfs_dir_sub_synop, figname_synop)    
                    os.system('cp '+ figname_to_copy + ' ' + synoptic_image)   
                    
                    if delete_working_fits:
                        os.system('rm -rf '+imagedir_allch + '*')
                    time_completed= timeit.default_timer() 
                    logging.debug('====All processing for time {0:s} is done in {1:.1f} minutes'.format(timestr, (time_completed-time_begin)/60.))
                    return True
                else:
                    if delete_working_fits:
                        os.system('rm -rf '+imagedir_allch + '*')
                    time_exit = timeit.default_timer()
                    logging.error('====Processing for time {0:s} failed in {1:.1f} minutes'.format(timestr, (time_exit-time_begin)/60.))
                    return False
            else:
                return False
        elif do_imaging:
            if delete_working_fits:
                os.system('rm -rf '+imagedir_allch + '*')
            time_exit = timeit.default_timer()
            logging.error('====Processing for time {0:s} failed in {1:.1f} minutes'.format(timestr, (time_exit-time_begin)/60.))
            return False
        else:
            if delete_working_fits:
                os.system('rm -rf '+imagedir_allch + '*')
            time_completed= timeit.default_timer() 
            logging.debug('====All processing for time {0:s} is done in {1:.1f} minutes'.format(timestr, (time_completed-time_begin)/60.))
            return True
    except Exception as e:
        logging.error(e)
        os.system('rm -rf '+imagedir_allch + '*')
        time_exit = timeit.default_timer()
        logging.error('====Processing for time {0:s} failed in {1:.1f} minutes'.format(image_time.isot, (time_exit-time_begin)/60.))
        return False

def parallel_task_runner(function_name,input_list,timeout=86400):
    '''
    Adapted from https://stackoverflow.com/questions/66051638/set-a-time-limit-on-the-pool-map-operation-when-using-multiprocessing
    This function will always use the same number of processes equal to the length of the list. Please see
    https://stackoverflow.com/questions/29494001/how-can-i-abort-a-task-in-a-multiprocessing-pool-after-a-timeout for implementing
    something when number of available threads is smaller than the length of the input_list. I have not taken a detailed look at this
    function, and putting this here for future reference if needed. Also note that there was some discussion about the second link, 
    in the discussion of the first link.
    '''
    pool = multiprocessing.pool.Pool(processes=len(input_list))
    results=[pool.apply_async(function_name,args=(item1,)) for item1 in input_list]
    starttime=timeit.default_timer()
    successful_results=[]
    time_to_wait=timeout
    for i,result in enumerate(results):
        try:
            return_value = result.get(time_to_wait)
        except multiprocessing.TimeoutError:
            pass
        else:
            successful_results.append(return_value)
        diff=timeit.default_timer()-starttime
        time_to_wait=timeout-diff
        if time_to_wait<0:
            time_to_wait=0
    pool.terminate()
    pool.close()
    pool.join()
    
    return successful_results
        
            

def image_times(msfiles_slfcaled, imagedir_allch, nch_out=12, stokes='I', beam_fit_size=2, briggs=-0.5):
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

    run_imager_partial = partial(run_imager, imagedir_allch=imagedir_allch, ephem=ephem, \
                nch_out=nch_out, stokes=stokes, beam_fit_size=beam_fit_size, briggs=briggs)

    timeout = 300.
    
    fitsfiles=parallel_task_runner(run_imager_partial,msfiles_slfcaled_success,timeout=timeout)
    if len(fitsfiles)==len(msfiles_slfcaled_success):
        time_img2 = timeit.default_timer()
        logging.debug('Imaging for all {0:d} bands is done in {1:.1f} s'.format(len(msfiles_slfcaled_success), time_img2-time_img1))
    else:
        logging.debug('Imaging for certain bands is incomplete in {0:.1f} s'.format(timeout))
        logging.debug('Proceed anyway')
                
    return fitsfiles

def compress_plot_images(fitsfiles, starttime, datedir, imagedir_allch_combined, \
                            hdf_dir, fig_mfs_dir, stokes, fast_vis=False):    

    if fast_vis:
        imagename_pre = 'ovro-lwa-48'
    else:
        imagename_pre = 'ovro-lwa-352'
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
    fits_mfs = imagedir_allch_combined_sub_lv10 + '/' + imagename_pre + '.lev1_mfs_10s.' + timestr_iso + '.image_'+stokes+'.fits' 
    #fitsfiles_mfs = glob.glob(imagedir_allch + '/' + timestr+ '*MFS-image.fits')
    fitsfiles_mfs = []
    for f in fitsfiles:
        if type(f) is list:
            #if 'MFS' in f[-1] and (not fast_vis):
            if 'MFS' in f[-1]:
                fitsfiles_mfs.append(f[-1])
            #elif fast_vis:
            #    fitsfiles_mfs+=f
             
        else:
            continue
    
    fitsfiles_mfs.sort()
    
    ndfits.wrap(fitsfiles_mfs, outfitsfile=fits_mfs)
    
    hdf_mfs = hdf_dir_sub_lv10 + os.path.basename(fits_mfs).replace('.fits', '.hdf')
    if fast_vis:
        utils.compress_fits_to_h5(fits_mfs, hdf_mfs, purge_corrupted=True)
    else:
        utils.compress_fits_to_h5(fits_mfs, hdf_mfs, purge_corrupted=False)
    
    #if not fast_vis:
    # fine channel spectral images
    fits_fch = imagedir_allch_combined_sub_lv10 + '/' + imagename_pre + '.lev1_fch_10s.' + timestr_iso + '.image_'+stokes+'.fits' 
    fitsfiles_fch = []
    for f in fitsfiles:
        if type(f) is list:
            fitsfiles_fch += f[:-1]
        else:
            continue
    fitsfiles_fch.sort()
    ndfits.wrap(fitsfiles_fch, outfitsfile=fits_fch)
    hdf_fch = hdf_dir_sub_lv10 + os.path.basename(fits_fch).replace('.fits', '.hdf')
    if fast_vis:
        utils.compress_fits_to_h5(fits_fch, hdf_fch, purge_corrupted=True)
    else:
        utils.compress_fits_to_h5(fits_fch, hdf_fch, purge_corrupted=False)
    
    fig, axes = ovis.slow_pipeline_default_plot(fits_mfs)
    figname_lv10 = os.path.basename(fits_mfs).replace('.fits', '.png')
    fig.savefig(fig_mfs_dir_sub_lv10 + '/' + figname_lv10)
    
    #if not fast_vis:
    return [fits_mfs, fits_fch], os.path.join(fig_mfs_dir_sub_lv10, figname_lv10)
    #else:
    #    return [fits_mfs], os.path.join(fig_mfs_dir_sub_lv10, figname_lv10)
    
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
   
    fits_mfs, fits_fch=fitsfiles

    px, py = orefr.refraction_fit_param(fits_fch, overbright=overbright)

    if (not np.isnan(px).any()) and (not np.isnan(py).any()): 
        df_new = pd.DataFrame({"Time":btime.isot, "px0":px[0], "px1":px[1], "py0":py[0], "py1":py[1]}, index=[0])
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


import os
import time
import shutil
from datetime import datetime

def remove_old_items(directory=".", minutes=45):
    """
    Remove files and folders older than specified minutes in the given directory
    """
    current_time = time.time()
    cutoff_time = current_time - (minutes * 60)
    
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)
        item_time = os.path.getctime(item_path)
        
        if item_time < cutoff_time:
            try:
                if os.path.isfile(item_path):
                    os.remove(item_path)
                else:
                    shutil.rmtree(item_path)
            except Exception as e:
                pass
        

def run_pipeline(time_start=Time.now(), time_end=None, time_interval=600., delay_from_now=180., do_selfcal=True, num_phase_cal=0, num_apcal=1, 
        server=None, lustre=True, file_path='slow', multinode=True, slurmmanaged=True, taskids='0123456789', delete_ms_slfcaled=True, slowfast='slow', 
        logger_dir = '/lustre/solarpipe/realtime_pipeline/logs/', logger_prefix='solar_realtime_pipeline', logger_level=20,
        #proc_dir = '/fast/solarpipe/realtime_pipeline/',
        proc_dir_mem = '/dev/shm/srtmp/', proc_dir = '/fast/solarpipe/realtime_pipeline/',
        proc_dir_isolation = True,
        save_dir = '/lustre/solarpipe/realtime_pipeline/',
        calib_dir = '/lustre/solarpipe/realtime_pipeline/caltables/',
        calib_file = '20240117_145752', altitude_limit=10., 
        beam_fit_size = 2,
        briggs=-0.5,
        delete_working_ms=True, do_refra=True, delete_working_fits=True,
        do_imaging=True, delete_allsky=False, save_allsky=False,
        save_selfcaltab=False,
        bands = ['32MHz', '36MHz', '41MHz', '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', '78MHz', '82MHz'],
        stop_at_sunset=True,
        do_daily_refracorr=True,
        slurm_kill_after_sunset=False,
        clear_old_files=True, clear_older_than=30, actively_rm_ms=True):
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
    :param slowfast: specify slow or fast visibility data to process
    :param delay_from_now: delay of the newest time to process compared to now.
    :param delete_ms_slfcaled: whether or not to delete the self-calibrated measurement sets.
    :param multinode: if True, will delay the start time by the node
    :param nodes: list of lwacalim nodes to be used. Default '0123456789'
    :param calib_file: calibration file to be used. Format yyyymmdd_hhmmss
    :param altitude_limit: lowest altitude to start the pipeline in degrees. Default to 15 deg.
    :param beam_fit_size: size of the beam area used for fitting to be passed to wsclean. See https://wsclean.readthedocs.io/en/v3.0/restoring_beam_size.html
    :param do_imaging: If True, imaging and other related steps will be done. Default: True
    :param briggs: briggs weighting parameter to be passed to wsclean. -1 is close to uniform and 1 is close to natural. See https://wsclean.readthedocs.io/en/latest/image_weighting.html
    :param do_imaging: If True (default), will do imaging. Otherwise just do all steps prior to imaging.
    :param bands: list of frequency bands to process. Can be a subset of the full list ['13MHz', '18MHz', '23MHz', '27MHz', '32MHz', '36MHz', '41MHz', '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', '78MHz', '82MHz']
    :param stop_at_sunset: if False, the program will attempt to continue across days untile being killed, otherwise it will exist after sunset.
    :param delete_allsky: if True, delete the allsky image after each run. Otherwise keep the latest frame. 
    :param save_allsky: if True, save the allsky image FITS file into save_dir + 'allsky/'
    '''


    if slurmmanaged:       
        task_id = int(os.environ.get('SLURM_PROCID', 0))
        task_count = int(os.environ.get('SLURM_NTASKS', 1))
        current_task_id = task_id    
    else:
        current_task_id = int(socket.gethostname()[-2:])

    # create proc_dir if not exist
    if not os.path.exists(proc_dir):
        os.makedirs(proc_dir)

    if not os.path.exists(proc_dir_mem):
        os.makedirs(proc_dir_mem)

    if proc_dir_isolation:
        proc_dir = proc_dir + 'task_' + str(current_task_id) + '/'
        proc_dir_mem = proc_dir_mem + 'task_' + str(current_task_id) + '/'

        if not os.path.exists(proc_dir):
            os.makedirs(proc_dir)
        if not os.path.exists(proc_dir_mem):
            os.makedirs(proc_dir_mem)
    

    print('<<',Time.now().isot,'>>','Starting solar real-time pipeline')
    try:
        time_start = Time(time_start)
    except Exception as e:
        logging.error(e)
        raise e

        
    (t_rise, t_set) = sun_riseset(time_start, altitude_limit=altitude_limit)
    # set up logging file
    server_runtime = socket.gethostname()
    datestr = Time(time_start.mjd, format='mjd').isot[:10].replace('-','')
    datedir = Time(time_start.mjd, format='mjd').isot[:10].replace('-','/') + '/'
    
    logger_file = logger_dir + datedir + logger_prefix + '_' + slowfast + '_'+ datestr + '_' + str(task_id) +'_' + server_runtime + '.log'  

    if not os.path.exists(os.path.dirname(logger_file)):
        print('Path to logger file {0:s} does not exist. Attempting to create the directory tree.'.format(logger_file))
        os.makedirs(os.path.dirname(logger_file))

    logging.basicConfig(filename=logger_file, filemode='at',
        format='%(asctime)s %(funcName)s %(lineno)d %(levelname)-8s %(message)s',
        level=logger_level,
        datefmt='%Y-%m-%d %H:%M:%S', force=True)

    logging.info('{0:s}: I am asked to start imaging for {1:s}'.format(socket.gethostname(), time_start.isot))
    if multinode:
        if slurmmanaged:
            # attribute the task_id as node number 
            nodenum = task_id
            nnode = task_count
            nodes_list=[int(n) for n in list(taskids)]

        else:
            # pdsh way, using hostname to manage
            nodenum = int(socket.gethostname()[-2:])
            nodes_list=[int(n) for n in list(taskids)]
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
        if slowfast.lower() == 'fast':
            twait += TimeDelta(600., format='sec') 
        logging.info('{0:s}: Start time {1:s} is before sunrise. Wait for {2:.1f} hours to start.'.format(socket.gethostname(), time_start.isot, twait.value * 24.))
        time_start += TimeDelta(twait.sec + 60., format='sec')
        sleep(twait.sec + 60.)
    else:
        logging.info("{0:s}: Start time {1:s} is after today's sunrise at {2:s}. Will try to proceed.".format(socket.gethostname(), time_start.isot, t_rise.isot))
    time_start += TimeDelta(delay_by_node, format='sec')
    logging.info('{0:s}: Delay {1:.1f} min to {2:s}'.format(socket.gethostname(), delay_by_node / 60., time_start.isot))
    sleep(delay_by_node)
    #loop_count=0
    while True:
        #loop_count+=1
        #if loop_count > 6:
            #gc.collect() # garbage collection every 6 frames
            #loop_count=0 
        
        if time_end:
            if time_start > Time(time_end):
                logging.info('The new imaging time now passes the provided end time. Ending the pipeline.'.format(Time(time_start).isot, Time(time_end).isot))
                print('<<',Time.now().isot,'>>','The new imaging time now passes the provided end time. Ending the pipeline.'.format(Time(time_start).isot, Time(time_end).isot))
                
                if slurm_kill_after_sunset:
                    # sleep for 15min for continuous imaging of all nodes to finish
                    sleep(900.)
                    # kill all the processes with scancel to jobname "solarpipedaily"
                    os.system('scancel -u solarpipe -n solarpipedaily')
                break

        time1 = timeit.default_timer()
        if time_start > Time.now() - TimeDelta(delay_from_now, format='sec'):
            twait = time_start - Time.now()
            logging.info('{0:s}: Start time {1:s} is too close to current time. Wait {2:.1f} m to start.'.format(socket.gethostname(), time_start.isot, (twait.sec + delay_from_now) / 60.))
            sleep(twait.sec + delay_from_now)
            time1 = timeit.default_timer()
        logging.info('{0:s}: Start processing {1:s}'.format(socket.gethostname(), time_start.isot))

        # do one round of cleaning up old files before pipeline_quick

        res = pipeline_quick(time_start, do_selfcal=do_selfcal, num_phase_cal=num_phase_cal, num_apcal=num_apcal, 
                            server=server, lustre=lustre, file_path=file_path, slowfast=slowfast, delete_ms_slfcaled=delete_ms_slfcaled,
                            logger_file=logger_file, proc_dir=proc_dir,  proc_dir_mem=proc_dir_mem, save_dir=save_dir, calib_dir=calib_dir, 
                            calib_file=calib_file, delete_working_ms=delete_working_ms,
                            delete_working_fits=delete_working_fits, do_refra=do_refra,
                            beam_fit_size=beam_fit_size, briggs=briggs, do_imaging=do_imaging, bands=bands, delete_allsky=delete_allsky, save_allsky=save_allsky,
                            clear_old_files=clear_old_files, clear_older_than=clear_older_than, save_selfcaltab=save_selfcaltab, actively_rm_ms=actively_rm_ms)

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
            
            # use last "worker" for daily refraction correction
            if do_daily_refracorr:
                if current_task_id == nodes_list[-1] and slowfast.lower()=='slow':
                    logging.info('{0:s}: Sun is setting. Done for the day. Doing refraction corrections for the full day.'.format(socket.gethostname())) 
                    daily_refra_correction(date_synop, save_dir=save_dir, overwrite=False, dointerp=True, interp_method='linear', max_dt=600.)
                    twait = t_rise_next - Time.now() 
                    logging.info('{0:s}: Refraction corrections done. Wait for {1:.1f} hours to start.'.format(socket.gethostname(), twait.value * 24.)) 
                else:
                    twait = t_rise_next - Time.now() 
                    logging.info('{0:s}: Sun is setting. Done for the day. Wait for {1:.1f} hours to start.'.format(socket.gethostname(), twait.value * 24.)) 

            if slowfast.lower() == 'fast':
                twait += TimeDelta(600., format='sec') 

            if stop_at_sunset:
                if slurm_kill_after_sunset:
                    # sleep for 15min for continuous imaging of all nodes to finish
                    sleep(900.)
                    logging.info('{0:s}: Sun is setting. Done for the day. Exiting.'.format(socket.gethostname()))
                    print('<<',Time.now().isot,'>>','Sun is setting. Done for the day. Exiting.')

                    # kill all the processes with scancel to jobname "solarpipedaily"
                    os.system('scancel -u solarpipe -n solarpipedaily')

                break
            else:
                time_start += TimeDelta(twait.sec + 60. + delay_by_node, format='sec')
                logging.info('{0:s}: Next time to image is {1:s} .'.format(socket.gethostname(), time_start.isot)) 
                t_rise = t_rise_next
                t_set = t_set_next
                # updating the logger file
                datestr = Time(t_rise.mjd, format='mjd').isot[:10].replace('-','')
                datedir = Time(t_rise.mjd, format='mjd').isot[:10].replace('-','/') + '/'
                logger_file = logger_dir + datedir + logger_prefix + '_' + slowfast + '_'+ datestr + '_' + str(current_task_id) + '.log'  

                if not os.path.exists(os.path.dirname(logger_file)):
                    print('<<',Time.now().isot,'>>','Path to logger file {0:s} does not exist. Attempting to create the directory tree.'.format(logger_file))
                    os.makedirs(os.path.dirname(logger_file))

                logging.basicConfig(filename=logger_file, filemode='at',
                    format='%(asctime)s %(funcName)s %(lineno)d %(levelname)-8s %(message)s',
                    level=logger_level,
                    datefmt='%Y-%m-%d %H:%M:%S', force=True)

                if twait.sec > 0:
                    sleep(twait.sec + 60. + delay_by_node) 

if __name__=='__main__':
    """
    Main routine of running the realtime pipeline. Example call
        Slow pipeline: pdsh -w lwacalim[00-09] 'conda activate suncasa && python /opt/devel/bin.chen/software/ovro-lwa-solar-ops/solar_realtime_pipeline.py 2024-05-15T14:10 --briggs -1.0 --slowfast slow --keep_allsky'
        Fast pipeline: pdsh -w lwacalim[00-09] 'conda activate suncasa && python /opt/devel/bin.chen/software/ovro-lwa-solar-ops/solar_realtime_pipeline.py 2024-05-15T14:10 --briggs -0.5 --slowfast fast --interval 100'
    Sometimes afer killing the pipeline (with ctrl c), one need to remove the temporary files and kill all the processes before restarting.
        pdsh -w lwacalim[00-09] 'rm -rf /fast/solarpipe/realtime_pipeline/slow_working/*'
        pdsh -w lwacalim[00-09] 'rm -rf /fast/solarpipe/realtime_pipeline/slow_slfcaled/*'
        pdsh -w lwacalim[00-09] 'rm -rf /fast/solarpipe/realtime_pipeline/fast_working/*'
        pdsh -w lwacalim[00-09] 'rm -rf /fast/solarpipe/realtime_pipeline/fast_slfcaled/*'
        pdsh -w lwacalim[00-09] 'rm -rf /dev/shm/srtmp/*' # clear ram disk
        pdsh -w lwacalim[00-09] 'pkill -u solarpipe -f wsclean'
        pdsh -w lwacalim[00-09] 'pkill -u solarpipe -f python'
    """
    parser = argparse.ArgumentParser(description='Solar realtime pipeline')
    parser.add_argument('--start_time', default=Time.now().isot, type=str, help='Timestamp for the start time. Format YYYY-MM-DDTHH:MM')
    parser.add_argument('--end_time', default='2030-01-01T00:00', help='End time in format YYYY-MM-DDTHH:MM')
    parser.add_argument('--interval', default=600., help='Time interval in seconds')
    parser.add_argument('--taskids', default='0123456789', help='List of taskids to use')
    parser.add_argument('--delay', default=60, help='Delay from current time in seconds')
    parser.add_argument('--server', default=None, help='Name of the server where the raw data is located. Must be defined in ~/.ssh/config.')
    parser.add_argument('--nolustre', default=False, help='If set, do NOT assume that the data are stored under /lustre/pipeline/ in the default tree', action='store_true')
    parser.add_argument('--file_path', default='slow/', help='Specify where the raw data is located')
    parser.add_argument('--proc_dir_mem', default='/dev/shm/srtmp/' )# default='/fast/solarpipe/realtime_pipeline/', help='Directory for processing')
    parser.add_argument('--proc_dir', default='/fast/solarpipe/realtime_pipeline/', help='Directory for processing')
    parser.add_argument('--save_dir', default='/lustre/solarpipe/realtime_pipeline/', help='Directory for saving fits files')
    parser.add_argument('--calib_dir', default='/lustre/solarpipe/realtime_pipeline/caltables/', help='Directory to calibration tables')
    parser.add_argument('--calib_file', default='', help='Calibration file to be used yyyymmdd_hhmmss')
    parser.add_argument('--alt_limit', default=15., help='Lowest solar altitude to start/end imaging')
    parser.add_argument('--bmfit_sz', default=2, help='Beam fitting size to be passed to wsclean')
    parser.add_argument('--briggs', default=-0.5, help='Briggs weighting parameter to be passed to wsclean')
    parser.add_argument('--do_refra', default=True, help='If True, do refraction correction', action='store_true')
    parser.add_argument('--singlenode', default=False, help='If True, delay the start time by the node', action='store_true')
    parser.add_argument('--logger_dir', default='/lustre/solarpipe/realtime_pipeline/logs/', help='Directory for logger files')
    parser.add_argument('--logger_prefix', default='solar_realtime_pipeline', help='Prefix for logger file')
    parser.add_argument('--logger_level', default=10, help='Specify logging level. Default to 10 (debug)')   
    parser.add_argument('--keep_working_ms', default=False, help='If True, keep the working ms files after imaging', action='store_true')
    parser.add_argument('--keep_working_fits', default=False, help='If True, keep the working fits files after imaging', action='store_true')
    parser.add_argument('--save_allsky', default=False, help='If True, save the band-averaged all sky images', action='store_true')
    parser.add_argument('--save_selfcaltab', default=True, help='If True, save the selfcalibration tables', action='store_true')
    parser.add_argument('--no_selfcal', default=False, help='If set, do not do selfcal regardless slow or fast', action='store_true')
    parser.add_argument('--no_imaging', default=False, help='If set, do not perform imaging', action='store_true')
    parser.add_argument('--nonstop', default=False, help='If set, the script will be run without stopping', action='store_true')
    parser.add_argument('--sleep_time', default=0.0, help='Process will sleep for these seconds before doing anything')
    parser.add_argument('--slowfast', default='slow', help='Specify slow or fast visibility data to be processed')
    parser.add_argument('--bands', '--item', action='store', dest='bands',
                    type=str, nargs='*', 
                    default=['32MHz', '36MHz', '41MHz', '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', '78MHz', '82MHz'],
                    help="Examples: --bands 32MHz 46MHz 64MHz")
    parser.add_argument('--no_actively_rm_ms', default=False, help='If set, actively remove the measurement sets after imaging', action='store_true')
    parser.add_argument('--no_refracorr', default=False, help='If set, do not do daily refraction correction', action='store_true')
    parser.add_argument('--slurm_kill_after_sunset', default=False, help='If set, kill all the processes with scancel after sunset', action='store_true')
    
    args = parser.parse_args()
    sleep(int(args.sleep_time))

    if args.slowfast.lower() == 'fast':
        do_selfcal=False
        delay=200.
    elif args.slowfast.lower() == 'slow' and args.no_selfcal:
        do_selfcal=False
    else:
        do_selfcal=True
    
    if len(args.calib_file) == 15:
        calib_file = args.calib_file
    else:
        logging.info('Calibration tables not provided or recognized. Attempting to find those from default location on lwacalim.')
        calib_tables = glob.glob('/lustre/solarpipe/realtime_pipeline/caltables_latest/*.bcal')
        if len(calib_tables) > 10:
            calib_file = os.path.basename(calib_tables[0])[:15]
            logging.info('Using calibration file {0:s}'.format(calib_file))
        else:
            logging.error('No calibration files found. Abort.')
            sys.exit(0)

    try:
        run_pipeline(time_start=args.start_time, time_end=Time(args.end_time), time_interval=float(args.interval), taskids=args.taskids, delay_from_now=float(args.delay),
                     server=args.server, lustre=(not args.nolustre), file_path=args.file_path,
                     proc_dir_mem=args.proc_dir_mem, proc_dir=args.proc_dir, save_dir=args.save_dir, calib_dir=args.calib_dir, calib_file=calib_file, 
                     altitude_limit=float(args.alt_limit), logger_dir = args.logger_dir, logger_prefix=args.logger_prefix, logger_level=int(args.logger_level), 
                     do_refra=args.do_refra, multinode= (not args.singlenode), delete_working_ms=(not args.keep_working_ms), 
                     delete_working_fits=(not args.keep_working_fits), save_allsky=args.save_allsky, beam_fit_size=args.bmfit_sz, briggs=args.briggs,
                     do_selfcal=do_selfcal, do_imaging=(not args.no_imaging), bands=args.bands, slowfast=args.slowfast, stop_at_sunset=(not args.nonstop),
                     do_daily_refracorr=(not args.no_refracorr), slurm_kill_after_sunset=args.slurm_kill_after_sunset, 
                     save_selfcaltab=args.save_selfcaltab, actively_rm_ms=(not args.no_actively_rm_ms))
    except Exception as e:
        logging.error(e)
        raise e



