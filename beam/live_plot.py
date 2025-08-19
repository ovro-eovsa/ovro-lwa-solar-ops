if __name__ == '__main__':
    import matplotlib
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from astropy.time import Time
import numpy as np
import glob, os
import h5py
from time import time, sleep
beam_data = "/lustre/pipeline/beam02/"

def plot_live_file(filename=None,timeskip=1,freqskip=1):
    if filename is None:
        files = glob.glob(beam_data+'*')
        files.sort()
        filename = os.path.basename(files[-1])

    f = h5py.File(beam_data+filename, 'r', libver='latest', swmr=True)
    freqs = f['Observation1']['Tuning1']['freq'][::freqskip]/1.e6
    df = freqs[1]-freqs[0]
    ts = f['Observation1']['time'][::timeskip]
    t0 = ts[1]
    t1 = ts[2]
    start_pd = Time(40587. + (t0[0] + t0[1])/86400.,format='mjd').plot_date
    dt = (t1[0]+t1[1]-t0[0]-t0[1])/86400.
    spec = f['Observation1']['Tuning1']['XX'][::timeskip,::freqskip]/1e4  # SFU
    nt, nf = spec.shape
    try:
        prev_idx = np.where(spec[:,0] == 0.)[0][0]
        live = True
    except:
        print('File',filename,'is not live.')
        live = False
    f.close()
    fig, ax = plt.subplots(1,1)
    im = ax.imshow(spec.T,vmin=0.4, vmax=1, aspect='auto',origin='lower',
                    extent=[start_pd,start_pd+dt*nt,freqs[0], freqs[0]+df*nf])
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))

    plt.show()
    if not live: return
    while 1:
        t1 = time()
        f = h5py.File(beamdata+filename, 'r', libver='latest', swmr=True)
        test = f['Observation1']['Tuning1']['XX'][::timeskip,0]
        try:
            last_idx = np.where(test == 0.)[0][0]
        except:
            print('Reached the end of the file.')
            return
        newdata = f['Observation1']['Tuning1']['XX'][prev_idx*timeskip:last_idx*timeskip:timeskip,::freqskip]/1e4
        f.close()
        print(prev_idx, last_idx)
        spec[prev_idx:last_idx] = newdata
        im.set_data(spec.T)
        prev_idx = last_idx
        plt.pause(0.01)
        print(str(time()-t1)[:3],'secs')

def format_axis(ax, t, live=False):
    if ax.get_xlabel() == '':
        # First setting up of axis
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
        ax.set_xlabel('Time [UT]')
        ax.set_ylabel('Frequency [MHz]')
        ax.text(0.025, 0.01, 'OVRO-LWA Solar Team (NJIT Center for Solar-Terrestrial Research)', fontsize=10, transform=plt.gcf().transFigure)
        ax.text(0.975, 0.01, 'OVRO Long Wavelength Array (Caltech)', ha='right',fontsize=10, transform=plt.gcf().transFigure)
    if live:
        # Live update of time
        ax.set_title('Live OVRO-LWA Spectrogram at '+t.iso[:19], fontsize=20)
    else:
        # Final archived plot at end of file
        ax.set_title('OVRO-LWA Spectrogram (Start time '+t.iso[:19]+')', fontsize=20)

def continuous_live_plot(timeskip=1,freqskip=1,sfufactor=20000):
    lastfile = ''
    plt.rcParams['font.size']=14
    while(1):
        files = glob.glob(beam_data+'*')
        files.sort()
        filename = os.path.basename(files[-1])
        if filename == lastfile:
            # No new file, so sleep for 60 s and try again
            print('Waiting for new file')
            sleep(60)
        else:
            # This is a new file, so do some initialization
            print('Detected a new file.  Initial setup for plotting.')
            f = h5py.File(beam_data+filename, 'r', libver='latest', swmr=True)
            ts = f['Observation1']['time'][::timeskip]
            t0 = ts[0]
            t1 = ts[1]
            if t0[0] == 0:
                # Rare case of first time being zero(!)
                t0 = ts[1]
                t1 = ts[2]
                if t0[0] == 0:
                    # Second time is also zero, so this may be a newly opened file with no data yet.
                    f.close()
                    # Sleep another 60 s and retry
                    sleep(60)
                    f = h5py.File(beam_data+filename, 'r', libver='latest', swmr=True)
                    ts = f['Observation1']['time'][::timeskip]
                    t0 = ts[0]
                    t1 = ts[1]
                    if t0[0] == 0:
                        # Rare case of first time being zero(!)
                        t0 = ts[1]
                        t1 = ts[2]
            lastfile = filename
            start_time = Time(40587. + (t0[0] + t0[1])/86400.,format='mjd')
            start_pd = start_time.plot_date
            datstr = start_time.iso[:19].replace('-','').replace(':','')
            dt = (t1[0]+t1[1]-t0[0]-t0[1])/86400.
            freqs = f['Observation1']['Tuning1']['freq'][::freqskip]/1.e6
            good_freqs, = np.where(freqs > 30)
            df = freqs[1]-freqs[0]
            spec = f['Observation1']['Tuning1']['XX'][::timeskip,::freqskip]/1e4/sfufactor  # SFU
            nt, nf = spec.shape
            nonzeroidx = np.where(ts['int'] != 0)[0]
            if len(nonzeroidx) == 0:
                prev_idx = 0   # File exists but has no data yet
            else:
                prev_idx = nonzeroidx[-1]+1    # Start of zero data
            if prev_idx < len(spec[:,0]):
                live = True
            else:
                prev_idx = len(spec[:,0])
                print('File',filename,'is not live.')
                live = False
            f.close()
            fig, ax = plt.subplots(1,1,figsize=(12,5))
            if prev_idx > 0:
                vec = spec[:prev_idx,good_freqs].flatten()
                vec.sort()
                vmax = vec[int(len(vec)*0.97)]
                vmin = vec[int(len(vec)*0.03)]
            else:
                vmax = 1.
                vmax = 0.1
            im = ax.imshow(spec.T,vmin=vmin, vmax=vmax, aspect='auto',origin='lower',
                            extent=[start_pd,start_pd+dt*nt,freqs[0], freqs[0]+df*nf])
            clb = fig.colorbar(im, ax=im.axes)
            clb.set_label('Flux [approx. sfu]')
            format_axis(ax, start_time, live)
            pngname = datstr.replace(' ','-')+'-lwa_beam.png'
            fig.subplots_adjust(right=1)
            figsave(pngname)
            print('New plot file written.')
            if live:
                sleep(10.)
                # Read continuously once per minute until the file ends
                while 1:
                    print('Attempt to read a live file')
                    f = h5py.File(beam_data+filename, 'r', libver='latest', swmr=True)
                    # Read times to check for zeros
                    test = f['Observation1']['time'][::timeskip]
                    # Find where the data ends, signified by all zeros at the end of the data
                    nonzeroidx = np.where(test['int'] != 0)[0]
                    last_idx = nonzeroidx[-1]+1
                    if last_idx > prev_idx:
                        print('Updating and saving latest data from the file.')
                        newdata = f['Observation1']['Tuning1']['XX'][prev_idx*timeskip:last_idx*timeskip:timeskip,::freqskip]/1e4/sfufactor  # SFU
                        ts = f['Observation1']['time'][prev_idx*timeskip:last_idx*timeskip:timeskip]
                        f.close()
                        while ts[-1][0] == 0.0: ts = ts[:-1]
                        t0 = ts[-1]
                        latest_time = Time(40587. + (t0[0] + t0[1])/86400.,format='mjd')
                        spec[prev_idx:last_idx] = newdata
                        im.set_data(spec.T)
                        vec = spec[:last_idx,good_freqs].flatten()
                        vec.sort()
                        vmax = vec[int(len(vec)*0.97)]
                        vmin = vec[int(len(vec)*0.03)]
                        im.norm.autoscale([vmin,vmax])
                        clb.update_normal(im)
                        format_axis(ax, latest_time, live)
                        prev_idx = last_idx
                        if matplotlib.get_backend() != 'Agg': plt.show()
                        figsave(pngname)
                        sleep(60.0)
                    else:
                        f.close()
                        print('Reached the end of the file.')
                        live = False
                        format_axis(ax, start_time, live)
                        figsave(pngname)
                        break

def figsave(pngname):
    datstr = pngname[:8]
    if pngname[9] == '0':
        # For times less than 10 UT, set datstr to previous day
        datstr = Time(Time(datstr[:4]+'-'+datstr[4:6]+'-'+datstr[6:8]).mjd-1,format='mjd').iso[:10].replace('-','')
    subyear = os.path.join('/opt/devel/solarpipe/figs_beam',datstr[:4])
    if not os.path.exists(subyear):
        os.makedirs(subyear)
    submonth = os.path.join(subyear,datstr[4:6])
    if not os.path.exists(submonth):
        os.makedirs(submonth)
    subday = os.path.join(submonth,datstr[6:])
    if not os.path.exists(subday):
        os.makedirs(subday)
    plt.savefig(subday+'/'+pngname)

if __name__ == '__main__':
    continuous_live_plot(timeskip=10, freqskip=4, sfufactor=1)