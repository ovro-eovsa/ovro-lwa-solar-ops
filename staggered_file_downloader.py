import os, sys, glob, getopt
from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, get_body, AltAz
from time import sleep
import argparse
import timeit
import shlex, subprocess
from casatasks import split
import numpy as np

def sun_riseset(date=Time.now(), observatory='ovro', altitude_limit=10):
    '''
    Given a date in Time object, determine the sun rise and set time as viewed from OVRO
    '''
    try:
        date_mjd = Time(date).mjd
    except Exception as e:
        logging.error(e)

    obs = EarthLocation.of_site(observatory)
    t0 = Time(int(date_mjd) + 13. / 24., format='mjd')
    sun_loc = get_body('sun', t0, location=obs)
    alt = sun_loc.transform_to(AltAz(obstime=t0, location=obs)).alt.degree
    while alt < 10.:
        t0 += TimeDelta(60., format='sec')
        alt = sun_loc.transform_to(AltAz(obstime=t0, location=obs)).alt.degree

    t1 = Time(int(date_mjd) + 22. / 24., format='mjd')
    sun_loc = get_body('sun', t1, location=obs)
    alt = sun_loc.transform_to(AltAz(obstime=t1, location=obs)).alt.degree
    while alt > altitude_limit:
        t1 += TimeDelta(60., format='sec')
        alt = sun_loc.transform_to(AltAz(obstime=t1, location=obs)).alt.degree

    return t0, t1

def download_msfiles(msfiles, destination='/fast/bin.chen/20231014_eclipse/slow_working/', bands=None, verbose=True, server=None, maxthread=5):
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
    
        
def list_msfiles(intime, lustre=True, file_path='slow', server=None, time_interval='10s', 
                # bands=['32MHz', '36MHz', '41MHz', '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', '78MHz', '82MHz']):
                bands=['46MHz','55MHz','78MHz']):
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
    print (file_path)
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
                pathstr = '{0:s}:{1:s}/{2:s}'.format(server, file_path, filename)
                tmpstr = filename[:15].replace('_', 'T')
                timestr = tmpstr[:4] + '-' + tmpstr[4:6] + '-' + tmpstr[6:11] + ':' + tmpstr[11:13] + ':' + tmpstr[13:]
                freqstr = filename[16:21]
                msfiles.append({'path': pathstr, 'name': filename, 'time': timestr, 'freq': freqstr})
    print (msfiles)
    return msfiles


def download_msfiles_cmd(msfile_path, server, destination):
    if server:
        p = subprocess.Popen(shlex.split('rsync -az --numeric-ids --info=progress2 --no-perms --no-owner --no-group {0:s}:{1:s} {2:s}'.format(server, msfile_path, destination)))
    else:
        p = subprocess.Popen(shlex.split('rsync -az --numeric-ids --info=progress2 --no-perms --no-owner --no-group {0:s} {1:s}'.format(msfile_path, destination)))
    std_out, std_err = p.communicate()
    if std_err:
        print(std_err)

def download_file(image_time,destination,min_nband=1,file_path='fast'):
    bands = ['32MHz', '36MHz', '41MHz', '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', '78MHz', '82MHz']
    
    msfiles0 = list_msfiles(image_time, file_path=file_path)
    if len(msfiles0) < min_nband:
        print('This time only has {0:d} subbands. Check nearby +-10s time.'.format(len(msfiles0)))
        image_time_before = image_time - TimeDelta(10., format='sec')
        msfiles0_before = list_msfiles(image_time_before, file_path=file_path)
        image_time_after = image_time + TimeDelta(10., format='sec')
        msfiles0_after = list_msfiles(image_time_after, file_path=file_path)
        if len(msfiles0_before) < min_nband and len(msfiles0_after) < min_nband:
            print('I cannot find a nearby time with at least {0:d} available subbands. Abort and wait for next time interval.'.format(min_nband))
            return False
        else:
            if len(msfiles0_before) > len(msfiles0_after):
                msfiles0 = msfiles0_before
                image_time = image_time_before
            else:
                msfiles0 = msfiles0_after
                image_time = image_time_after    
    msfiles = download_msfiles(msfiles0, destination=destination, bands=bands)
    return
    

            

def copy_files(time_start=Time.now(), time_end=None, time_interval=600., delay_from_now=180.,
        server=None, lustre=True, file_path='slow', multinode=True, nodes='0123456789', 
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
        
    if type(time_start) is str:
        time_start=Time(time_start)
    
    visdir_work = os.path.join(proc_dir , file_path+'_working/')
    if not os.path.exists(visdir_work):
        os.makedirs(visdir_work)
        
    (t_rise, t_set) = sun_riseset(time_start, altitude_limit=altitude_limit)
    # set up logging file
    server = socket.gethostname()
    datestr = Time(time_start.mjd, format='mjd').isot[:10].replace('-','')
    datedir = Time(time_start.mjd, format='mjd').isot[:10].replace('-','/') + '/'
    

    print('{0:s}: I am asked to start imaging for {1:s}'.format(socket.gethostname(), time_start.isot))
    if multinode:
        nodenum = int(socket.gethostname()[-2:])
        nodes_list=[int(n) for n in list(nodes)]
        nnode = len(nodes_list)
        delay_by_node = nodes_list.index(nodenum) * (time_interval/nnode) 
    else:
        print('{0:s}: I am running on a single node'.format(socket.gethostname()))
        nodenum = int(socket.gethostname()[-2:])
        nodes_list = [nodenum]
        delay_by_node = 0. 
    #while time_start > t_rise and time_start < Time.now() - TimeDelta(15.,format='sec'): 
    # find out when the Sun is high enough in the sky
    if time_start < t_rise:
        twait = t_rise - time_start
        print('{0:s}: Start time {1:s} is before sunrise. Wait for {2:.1f} hours to start.'.format(socket.gethostname(), \
                    time_start.isot, twait.value * 24.))
        time_start += TimeDelta(twait.sec + 60., format='sec')
        sleep(twait.sec + 60.)
    else:
        print("{0:s}: Start time {1:s} is after today's sunrise at {2:s}. Will try to proceed.".\
                        format(socket.gethostname(), time_start.isot, t_rise.isot))
    time_start += TimeDelta(delay_by_node, format='sec')
    print('{0:s}: Delay {1:.1f} min to {2:s}'.format(socket.gethostname(), delay_by_node / 60., time_start.isot))
    sleep(delay_by_node)
    while True:
        
        if time_end:
            if time_start > Time(time_end):
                print('The new imaging time now passes the provided end time. Ending the pipeline.'.\
                            format(Time(time_start).isot, Time(time_end).isot))
                break
        time1 = timeit.default_timer()
        if time_start > Time.now() - TimeDelta(delay_from_now, format='sec'):
            twait = time_start - Time.now()
            print('{0:s}: Start time {1:s} is too close to current time. Wait {2:.1f} m to start.'.\
                        format(socket.gethostname(), time_start.isot, (twait.sec + delay_from_now) / 60.))
            sleep(twait.sec + delay_from_now)
        print('{0:s}: Start processing {1:s}'.format(socket.gethostname(), time_start.isot))
        
        time2 = timeit.default_timer()
        
        files=glob.glob(os.path.join(visdir_work,"*.ms"))
        freqs=[]
        for file1 in files:
            cfreqidx = os.path.basename(file1).find('MHz') - 2
            cfreq = int(os.path.basename(file1)[cfreqidx:cfreqidx+2])
            freqs.append(cfreq)
        if len(freqs)!=0:
            freqs=np.array(freqs)
            unique_freq=np.unique(freqs)
            num_files=len(glob.glob(os.path.join(visdir_work,"*"+str(unique_freq[0])+"MHz*.ms")))
        else:
            num_files=0
        if num_files<max_files:
            download_file(time_start,visdir_work)
            time_start += TimeDelta(time_interval, format='sec')

            if time_start > t_set:
                (t_rise_next, t_set_next) = sun_riseset(t_set + TimeDelta(6./24., format='jd'))
                date_mjd = int(time_start.mjd)
                if time_start.mjd - date_mjd < 4./24.:
                    date_synop = Time(time_start.mjd - 1., format='mjd').isot[:10]
                else:
                    date_synop = Time(time_start.mjd, format='mjd').isot[:10]
                
                if int(socket.gethostname()[-2:]) == nodes_list[-1]:
                    print('{0:s}: Sun is setting. Done for the day.'.format(socket.gethostname())) 
                     
                else:
                    twait = t_rise_next - Time.now() 
                    print('{0:s}: Sun is setting. Done for the day. Wait for {1:.1f} hours to start.'.format(socket.gethostname(), twait.value * 24.)) 

                time_start += TimeDelta(twait.sec + 60. + delay_by_node, format='sec')
                t_rise = t_rise_next
                t_set = t_set_next
                # updating the logger file
                datestr = Time(t_rise.mjd, format='mjd').isot[:10].replace('-','')
                datedir = Time(t_rise.mjd, format='mjd').isot[:10].replace('-','/') + '/'
                logger_file = logger_dir + datedir + logger_prefix + '_' + datestr + '_' + server + '.log'  

                sleep(twait.sec + 60. + delay_by_node)
        else:
            sleep(5)




        
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Solar realtime pipeline')
    parser = argparse.ArgumentParser(description='Solar realtime pipeline')
    parser.add_argument('prefix', type=str, help='Timestamp for the start time. Format YYYY-MM-DDTHH:MM')
    parser.add_argument('--end_time', default='2030-01-01T00:00', help='End time in format YYYY-MM-DDTHH:MM')
    parser.add_argument('--interval', default=600., help='Time interval in seconds')
    parser.add_argument('--nodes', default='0123456789', help='List of nodes to use')
    parser.add_argument('--delay', default=60, help='Delay from current time in seconds')
    parser.add_argument('--server', default=None, help='Name of the server where the raw data is located. Must be defined in ~/.ssh/config.')
    parser.add_argument('--nolustre', default=False, help='If set, do NOT assume that the data are stored under /lustre/pipeline/ in the default tree', action='store_true')
    parser.add_argument('--file_path', default='slow/', help='Specify where the raw data is located')
    parser.add_argument('--proc_dir', default='/fast/bin.chen/realtime_pipeline/', help='Directory for processing')
    parser.add_argument('--save_dir', default='/lustre/bin.chen/realtime_pipeline/', help='Directory for saving fits files')
    parser.add_argument('--calib_dir', default='/lustre/bin.chen/realtime_pipeline/caltables/', help='Directory to calibration tables')
    parser.add_argument('--calib_file', default='', help='Calibration file to be used yyyymmdd_hhmmss')
    parser.add_argument('--alt_limit', default=15., help='Lowest solar altitude to start/end imaging')
    parser.add_argument('--bmfit_sz', default=2, help='Beam fitting size to be passed to wsclean')
    parser.add_argument('--do_refra', default=True, help='If True, do refraction correction', action='store_false')
    parser.add_argument('--singlenode', default=False, help='If True, delay the start time by the node', action='store_true')
    parser.add_argument('--logger_dir', default='/lustre/bin.chen/realtime_pipeline/logs/', help='Directory for logger files')
    parser.add_argument('--logger_prefix', default='solar_realtime_pipeline', help='Prefix for logger file')
    parser.add_argument('--logger_level', default=10, help='Specify logging level. Default to 10 (debug)')   
    parser.add_argument('--keep_working_ms', default=False, help='If True, keep the working ms files after imaging', action='store_true')
    parser.add_argument('--keep_working_fits', default=False, help='If True, keep the working fits files after imaging', action='store_true')  
    parser.add_argument('--no_selfcal', default=False, help='If True, perform selfcal', action='store_true')
    parser.add_argument('--no_imaging', default=False, help='If True, perform imaging', action='store_true')
    parser.add_argument('--bands', '--item', action='store', dest='bands',
                        type=str, nargs='*', 
                        default=['32MHz', '36MHz', '41MHz', '46MHz', '50MHz', '55MHz', '59MHz', '64MHz', '69MHz', '73MHz', '78MHz', '82MHz'],
                        help="Examples: --bands 32MHz 46MHz 64MHz")
    args = parser.parse_args()
    copy_files(args.prefix)        
    
