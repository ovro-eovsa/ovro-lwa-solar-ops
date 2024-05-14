import sys, os, glob
import shlex
import subprocess
from astropy.time import Time
import numpy as np
from time import sleep
import argparse
def get_filecopy_pid(mypid):
    # check whether a fits copy process is already running. Return PID if it is running, -1 if it is not.
    pidlist = subprocess.check_output(["pidof","python"]).split() # list of PIDs for all python processes
    for pid in pidlist:
        if mypid != pid:
            ps_out = str(subprocess.check_output(["ps", "-lfp" , pid]))
            ind = ps_out.find('daily_lwa_file_transfer.py') 
            if ind != -1:  # if this PID is running auto_beamcopy.py
                return str(pid)
    return -1


def copy_files(timestr, server='calim1', remote_dir='/lustre/bin.chen/realtime_pipeline/',
        local_dir='/nas6/ovro-lwa-data/', file_type='fits', file_ext='fits', file_level='slow/lev1/', ndays=1):
    """
    Purpose: rsync files from default directory of calim07:/lustre/bin.chen/realtime_pipeline/fits/*
        to pipeline:/nas6/ovro-lwa-data/fits/* and organize them according to date.
    timestr: input time '2024-01-02T18:00:00'
    """
    mjd = Time(timestr).mjd
    days = mjd - np.arange(ndays)
    for day in days:
        daystr = Time(day,format='mjd').iso
        yyyy = daystr[:4]
        mm = daystr[5:7]
        dd = daystr[8:10]
        if file_type.lower() == 'fits' or file_type.lower() == 'figs_mfs' or file_type.lower() == 'hdf':
            remotefolder = remote_dir + '/' + file_type + '/' + file_level + '/' + yyyy + '/' + mm + '/' + dd
        elif file_type.lower() == 'beam_plots':
            remotefolder = remote_dir + '/'
        else:
            print('file_type {0:s} not recognized! Abort...'.format(file_type)) 

        command = 'ssh {0:s} ls {1:s} | grep {2:s}'.format(server, remotefolder, file_ext)
        p = subprocess.run(shlex.split(command), capture_output=True)
        filenames = p.stdout.decode('utf-8').split('\n')[:-1]
        nfiles = len(filenames)
        print('{0:s}: Found {1:d} {2:s} files under {3:s}:{4:s}'.format(Time.now().iso[:19], nfiles, file_ext, server, remotefolder))
        if nfiles > 0:
            if file_type.lower() == 'fits' or file_type.lower() == 'figs_mfs' or file_type.lower() == 'hdf':
                localfolder = local_dir + '/' + file_level + '/' + yyyy + '/' + mm + '/' + dd + '/'
            elif file_type.lower() == 'beam_plots':
                localfolder = local_dir + '/'
            if not os.path.exists(localfolder):
                os.makedirs(localfolder)
                sleep(0.1)
            command = 'rsync -a --stats {0:s}:{1:s}/ {2:s}'.format(server, remotefolder, localfolder)
            print('{0:s}: Attemping to run {1:s}'.format(Time.now().iso[:19], command))
            res = subprocess.run(shlex.split(command), capture_output=True)
            output = res.stdout.decode('utf-8').split('\n')
            for outline in output:
                print(outline)
            if file_type.lower() == 'figs_mfs':
                lwa_html_movie(daystr, localfolder)
                files = glob.glob(localfolder + '*.png')
                files.sort()
                if len(files) > 0:
                    os.system('cp ' + files[-1] + ' ' + local_dir + '/../latest_' + file_level[:4] + '_image.png')
    return output

def copy_1min_plots(timestr, server='calim1', remote_img_folder='/lustre/bin.chen/realtime_pipeline/', 
        remote_beam_folder='/opt/devel/dgary/figs_beam/', 
        local_root_folder='/common/webplots/lwa-data/',
        local_img_folder='/common/webplots/lwa-data/qlook_images/',
        local_beam_staging_folder='/common/webplots/lwa-data/figs_beam/', 
        local_beam_final_folder='/common/webplots/lwa-data/qlook_spectra/', 
        ndays=2):
    ''' Automatically copy continually updating solar beam spectrogram files
        and solar images as png files, for display on
        the web.

        Since the files are now in subdirectories by date, we need to know how many such
        subdirectories to rsync, given by input: ndays.  The default is 3 days.  We also
        now need to know the current date, given by the iso time string tstr (e.g. '2023-12-10 20:23:12.000')

        Files are also saved in local folders local_img_final_folder and local_beam_final_folder
        with date subdirectories.
    '''
    ## Copy and reorganize beam spectrogram plots
    output = copy_files(timestr, server=server, remote_dir=remote_beam_folder, local_dir=local_beam_staging_folder, ndays=ndays, 
            file_type='beam_plots', file_ext='png', file_level='')

    files = glob.glob(local_beam_staging_folder + '*.png')
    files.sort()
    if len(files) > 0:
        os.system('cp ' + files[-1] + ' ' + local_root_folder + 'latest_spectrum.png')
    curdatestr = Time.now().iso.split(' ')[0].replace('-','')
    for f in files:
        filename = os.path.basename(f)
        dirname = local_beam_final_folder
        filetime = Time(filename[:4] + '-' + filename[4:6] + '-' + filename[6:8] + 'T' + filename[9:11] + ':' + filename[11:13])
        date_mjd = int(filetime.mjd)   
        if filetime.mjd - date_mjd < 4./24.:
            datestr_synop = Time(filetime.mjd - 1., format='mjd').isot[:10].replace('-','')
            datedir_synop = Time(filetime.mjd - 1., format='mjd').isot[:10].replace('-','/') + '/'
        else:
            datestr_synop = Time(filetime.mjd, format='mjd').isot[:10].replace('-','')
            datedir_synop = Time(filetime.mjd, format='mjd').isot[:10].replace('-','/') + '/'
    
        newpath = os.path.join(dirname, datedir_synop)
        if not os.path.exists(newpath):
            os.makedirs(newpath)
            sleep(0.1)
        outfile = os.path.join(newpath,filename)
        # Rsync file by file, logging only needed transfers
        blah = subprocess.run(["rsync", "-a", "--stats", f, newpath], capture_output=True)
        output = blah.stdout.decode('utf-8').split('\n')
        for outline in output:
            if outline.find('regular') != -1: 
                if int(outline.split(':')[1]) == 1: print(Time.now().iso[:19]+':', 'Transferred', outfile)

    # Copy synoptic slow and fast images
    output = copy_files(timestr, server=server, remote_dir=remote_img_folder, local_dir=local_img_folder, ndays=ndays, 
            file_type='figs_mfs', file_ext='png', file_level='slow/synop')

    output = copy_files(timestr, server=server, remote_dir=remote_img_folder, local_dir=local_img_folder, ndays=ndays, 
            file_type='figs_mfs', file_ext='png', file_level='fast/synop')

def copy_daily_fits(timestr, server='calim1', remote_data_folder='/lustre/bin.chen/realtime_pipeline/', 
        local_root_folder='/nas6/ovro-lwa-data/', ndays=2):
    ''' Automatically copy FITS and HDF data to local server
    '''
    # copy slow fits and hdf files
    local_fits_folder = local_root_folder + '/fits/'
    local_hdf_folder = local_root_folder + '/hdf/'

    output = copy_files(timestr, server=server, remote_dir=remote_data_folder, local_dir=local_fits_folder, ndays=ndays, 
            file_type='fits', file_ext='fits', file_level='slow/lev1/')
    output = copy_files(timestr, server=server, remote_dir=remote_data_folder, local_dir=local_fits_folder, ndays=ndays, 
            file_type='fits', file_ext='fits', file_level='slow/lev15/')
    output = copy_files(timestr, server=server, remote_dir=remote_data_folder, local_dir=local_fits_folder, ndays=ndays, 
            file_type='hdf', file_ext='hdf', file_level='slow/lev1/')
    output = copy_files(timestr, server=server, remote_dir=remote_data_folder, local_dir=local_fits_folder, ndays=ndays, 
            file_type='hdf', file_ext='hdf', file_level='slow/lev15/')

    # copy fast fits and hdf files
    output = copy_files(timestr, server=server, remote_dir=remote_data_folder, local_dir=local_hdf_folder, ndays=ndays, 
            file_type='fits', file_ext='fits', file_level='fast/lev1/')
    output = copy_files(timestr, server=server, remote_dir=remote_data_folder, local_dir=local_hdf_folder, ndays=ndays, 
            file_type='fits', file_ext='fits', file_level='fast/lev15/')
    output = copy_files(timestr, server=server, remote_dir=remote_data_folder, local_dir=local_hdf_folder, ndays=ndays, 
            file_type='hdf', file_ext='hdf', file_level='fast/lev1/')
    output = copy_files(timestr, server=server, remote_dir=remote_data_folder, local_dir=local_hdf_folder, ndays=ndays, 
            file_type='hdf', file_ext='hdf', file_level='fast/lev15/')

def copy_beam_data(timestr, server='calim1', remote_data_folder='/lustre/pipeline/beam02/',
        local_data_folder='/nas5/ovro-lwa-data/beam/beam-data/', ndays=2):
    """
    Automatically copy beamforming data to local server
    """
    blah = subprocess.run(["ssh", server, "ls " + remote_data_folder + " | tail -30"], capture_output=True)
    lines = blah.stdout.decode('utf-8').split('\n')
    # Example line is '060225_210100000000bf4b7d0'
    # Go through the list and find those within 24 hours of the time in tstr, or current time
    if timestr is None:
        mjd = Time.now().mjd
    else:
        mjd = Time(timestr).mjd
    mjdlim = min(mjd, Time.now().mjd - 1.05/24.)  # Most recent allowable file (63 min before current time)
    for line in lines:
        try:
            # Interpret the line as a filename encoding the mjd and clock time (see above example)
            linemjd = Time(line[:6], format='mjd').mjd + (float(line[7:9])+float(line[9:11])/60.+float(line[11:13])/3600.)/24.
            if linemjd >= mjd - ndays and linemjd <= mjdlim:
                # This is a line potentially needing to be transferred.  Check if it already exists.
                datstr = Time(linemjd - 3./24.,format='mjd').iso.replace('-','')[:8]   # Note, this places files 00:00-03:00 UTC in previous day's folder.
                monthdir = '/nas5/ovro-lwa-data/beam/beam-data/'+datstr[:6]
                if not os.path.exists(monthdir):
                    # The month directory doesn't exist, so create it
                    os.mkdir(monthdir)
                datedir = monthdir+'/beam'+datstr
                if not os.path.exists(datedir):
                    # The date directory doesn't exist, so create it.
                    os.mkdir(datedir)
                fname = datedir+'/'+line
                print('Transferring file', line, 'to', fname)
                blah = subprocess.run(["rsync", "-a", "--stats", server + ':' + remote_data_folder + line, fname], capture_output=True)
                output = blah.stdout.decode('utf-8').split('\n')
                print(Time.now().iso[:19]+':', 'RSYNC Summary for file', line)
                for outline in output:
                    if outline[:4] == 'sent': print('                    ',outline,'\n')
                return output
        except:
            return -1
            # Interpretation of line as mjd-coded filename failed
            if line != '':
                # If it is not just a blank line, emit a warning.
                print('***Warning: Line "'+line+'" is not a valid filename.  Skipping...')


def lwa_html_movie(t, image_dir=None):
    ''' This routine will be called after every update to the figs_mfs
        folder (in /common/webplots/lwa-data) to write the movie.html file that 
        allows them to be viewed as a movie.  Just call this with a Time() object 
        or iso time string containing the desired date.
    '''
    import glob, os
    import imageio.v2 as imageio
    try:
        datstr = t.iso[:10]
    except:
        datstr = t[:10]
    files = glob.glob(image_dir + '*.png')
    files.sort()
    nfiles = len(files)
    if nfiles == 0:
        print('No files (yet) in folder',image_dir)
        return
    # Read one file to determine its size
    img = imageio.imread(files[0])
    ysize, xsize, ncolors = img.shape
    f = open('/nas5/ovro-lwa-data/beam/software/html_movie_example.html', 'r')
    lines = f.readlines()
    nlines = len(lines)
    f.close()
    skiplines = []

    for i, line in enumerate(lines):
        k = line.find('var imax')
        j = line.find('var iwidth')
        l = line.find('NAME=animation')
        if k != -1:
            # Set number of frames
            lines[i] = line[:10] + '{:3d}'.format(nfiles) + line[13:]
        if j != -1:
            # Set width and height of images
            lines[i] = 'var iwidth = {:d}, iheight = {:d};\n'.format(xsize, ysize)
        if l != -1:
            # Set width and height of frame
            if xsize > 1125:
                # Reduce frame size for overly large images
                xfsize = xsize*3//4
                yfsize = ysize*3//4
            else:
                xfsize = xsize
                yfsize = ysize
            lines[i] = '<img NAME=animation ALT="FRAME" width='+str(xfsize)+' height='+str(yfsize)+'>'
        k = line.find('urls[')
        if k != -1:
            skiplines.append(i)
    #print skiplines
    htmlname = image_dir + '/movie_' + datstr + '.html'
    f = open(htmlname, 'w')
    for i in range(skiplines[1] - 1):
        f.write(lines[i])
    for i in range(nfiles):
        fname = os.path.basename(files[i])
        f.write('urls[{:d}]=url_path+"/'.format(i)+fname+'";\n')
    for i in range(skiplines[-1]+1,nlines):
        f.write(lines[i])
    f.close()
    print('html saved to {}'.format(htmlname))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Daily copy of OVRO-LWA data products')
    parser.add_argument('--plots', default=False, help='If set, copy and organize quicklook spectrogram and image plots', action='store_true')
    parser.add_argument('--beam', default=False, help='If set, copy and organize beam data', action='store_true')
    parser.add_argument('--fits', default=False, help='If set, copy and organize image fits and hdf data', action='store_true')
    timestr = Time.now().iso
    # ndays is the number of days to go back in time for checking the existence of data
    args = parser.parse_args()
    if args.plots:
        copy_1min_plots(timestr, ndays=1)
    if args.beam:
        copy_beam_data(timestr, ndays=2)
    if args.fits:
        copy_daily_fits(timestr, ndays=2)
    exit()
