import os
from astropy.time import Time,TimeDelta
import shlex, subprocess
from time import sleep
import socket,glob
import logging
import timeit
from ovrolwasolar import flagging
import numpy as np



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


def download_calibms(calib_time, download_fold = '/lustre/solarpipe/realtime_pipeline/ms_calib/',
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
    return ms_calib
