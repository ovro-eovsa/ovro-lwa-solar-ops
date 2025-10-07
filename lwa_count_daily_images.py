import sys, os
import fnmatch
from datetime import datetime, timedelta
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import h5py
from astropy.coordinates import SkyCoord, EarthLocation, get_body, AltAz
from astropy.time import Time, TimeDelta
import argparse

CSV_FILE = './lwa_daily_image_counts_db.fch.csv'
CSV_FILE_HOURLY = './lwa_hourly_image_counts_db.fch.csv'

def sun_riseset(date=Time.now(), observatory='ovro', altitude_limit=15.):
    '''
    Given a date in Time object, determine the sun rise and set time as viewed from OVRO
    :param date: input time in astropy.time.Time format
    :param observatory: name of the observatory recognized by astropy
    :param altitude_limit: lower limit of altitude to consider. Default to 15 degrees.


    :return trise, tset
    '''
    try:
        date_mjd = Time(date).mjd
    except Exception as e:
        print(e)

    t_now = Time(date_mjd, format='mjd')

    obs = EarthLocation.of_site(observatory)
    t0 = Time(int(date_mjd) + 12. / 24., format='mjd')

    sun_loc = get_body('sun', t0, location=obs)

    alt = sun_loc.transform_to(AltAz(obstime=t0, location=obs)).alt.degree
    while alt < altitude_limit:
        t0 += TimeDelta(60., format='sec')
        alt = sun_loc.transform_to(AltAz(obstime=t0, location=obs)).alt.degree

    t1 = Time(int(date_mjd) + 23. / 24., format='mjd')
    sun_loc = get_body('sun', t1, location=obs)
    alt = sun_loc.transform_to(AltAz(obstime=t1, location=obs)).alt.degree
    while alt > altitude_limit:
        t1 += TimeDelta(60., format='sec')
        alt = sun_loc.transform_to(AltAz(obstime=t1, location=obs)).alt.degree

    if t_now < t1-TimeDelta(1., format='jd'):
        # still before previous day sunset
        t1 = t1 - TimeDelta(1., format='jd')
        t0 = t0 - TimeDelta(1., format='jd')

    return t0, t1

def scan_all_dir(base_dir='/nas7/ovro-lwa-data/hdf/slow/lev1/', filetype='fch'):
    year_dirs= [f.path for f in os.scandir(base_dir) if f.is_dir()]
    year_dirs.sort()
    pattern = '*' + filetype + '*.hdf'
    tot_count = 0 
    data = []
    for year_dir in year_dirs:
        month_dirs = [f.path for f in os.scandir(year_dir) if f.is_dir()]
        month_dirs.sort()
        for month_dir in month_dirs:
            month_count = 0
            for entry in os.scandir(month_dir):
                # Only process subdirectories
                if entry.is_dir():
                    count = 0  # Initialize counter for this subdirectory
                    datestr = month_dir[-7:] + '/' + entry.name
                    date = datetime.strptime(datestr, "%Y/%m/%d")

                    # Walk through the current subdirectory
                    for root, dirs, files in os.walk(entry.path):
                        # Check each file against the pattern
                        for f in files:
                            if fnmatch.fnmatch(f, pattern):
                                count += 1

                    data.append((date, count))
                    month_count += count
            tot_count += month_count
            print(f"{month_dir}: {month_count}")

    print(f'Total: {tot_count}')
    df = pd.DataFrame(data, columns=["Date", "Daily Files"])
    df = df.sort_values("Date")
    return df

def scan_pastdays(ndays=7, startdate=Time.now(), base_dir='/nas7/ovro-lwa-data/hdf/slow/lev1/', filetype='fch'):
    from astropy.io import fits
    tot_file_count = 0 
    tot_image_count = 0
    data = []
    for i in range(ndays):
        date = startdate.datetime.date() - timedelta(days=i)
        subdir = date.strftime("%Y/%m/%d/")
        pattern = '*' + filetype + '*.hdf'
        # Only process subdirectories
        filedir = base_dir + subdir
        if os.path.isdir(filedir):
            # calculate the number of files and images if fully operational
            t0, t1 = sun_riseset(Time(date.strftime("%Y-%m-%d")), altitude_limit=15)
            if date < Time('2024-11-01').datetime.date():
                cadence = 60.
            else:
                cadence = 20.
            nfile_th = int((t1-t0).sec/cadence)
            nimage_th = int((t1-t0).sec/cadence*144.)
            file_count = 0  # Initialize counter for this subdirectory
            image_count = 0
            # Walk through the current subdirectory
            for root, dirs, files in os.walk(filedir):
                # Check each file against the pattern
                for f in files:
                    if fnmatch.fnmatch(f, pattern):
                        try:
                            with h5py.File(filedir + f, 'r') as f0:
                                header = dict(f0['ch_vals'].attrs)
                                datashape = header['original_shape']
                                nimages = datashape[1]
                            file_count += 1
                            image_count += nimages
                        except Exception as e:
                            print(f"An unexpected error occurred for file {f}: {e}")
                            continue

            datestr = date.strftime("%Y-%m-%d")
            data.append((datestr, file_count, image_count, nfile_th, nimage_th))
            print(f"{date}: {file_count} files, {image_count} images")
            tot_file_count += file_count
            tot_image_count += image_count

    if len(data) > 0:
        print(f'Total: {tot_file_count} files, {tot_image_count} images')
        df = pd.DataFrame(data, columns=["Date", "Daily Files", "Daily Images", "Max Daily Files", "Max Daily Images"])
        df = df.sort_values("Date", ascending=False)
        return df
    else:
        print('No new data generated')
        return pd.DataFrame()
        


def scan_pasthours(nhours=24, startdate=Time.now(), base_dir='/lustre/solarpipe/realtime_pipeline/', 
        filetype='fch', fileformat='fits'):
    from astropy.io import fits
    if fileformat == 'fits':
        base_dir += 'fits/slow/lev1/'
    elif fileformat == 'hdf':
        base_dir += 'hdf/slow/lev1/'
    else:
        print(f'{fileformat} not recognized')
    print(f'Scanning {base_dir}')
    tot_file_count = 0 
    tot_image_count = 0
    data = []
    if startdate.datetime.hour > 12:
        date = startdate.datetime.date()
    else:
        date = startdate.datetime.date() - timedelta(days=1)
    if date < Time('2024-11-01').datetime.date():
        cadence = 60.
    else:
        cadence = 20.
    for i in range(nhours):
        timeplot = (startdate.datetime - timedelta(hours=i)).replace(minute=0, second=0, microsecond=0)
        # trim to start of the hour
        t0, t1 = sun_riseset(timeplot, altitude_limit=15)
        nfile_hourly = int(3600./cadence)
        #print(f'timeplot:{timeplot}, t0:{t0.isot}, t1:{t1.isot}')
        if (timeplot < t1.datetime) and ((Time(timeplot) - t0).sec > -3600):
            nfile_pre = int((Time(timeplot) - t0).sec/cadence)
            nfile_end = int((t1 - Time(timeplot)).sec/cadence)
            if (nfile_pre >= nfile_hourly) and (nfile_end > nfile_hourly):
                nfile_hourly_th = nfile_hourly
            elif (nfile_pre >= nfile_hourly) and (nfile_end < nfile_hourly):
                nfile_hourly_th = nfile_end
            #elif (nfile_pre < nfile_hourly) and (nfile_end > nfile_hourly):
            #    nfile_hourly_th = nfile_pre
            elif nfile_pre < 0:
                nfile_hourly_th = nfile_pre + nfile_hourly
            #print(f'nfile_pre:{nfile_pre}, nfile_end:{nfile_end}, nfile_hourly:{nfile_hourly}, nfile_hourly_th:{nfile_hourly_th}')
            else:
                nfile_hourly_th = nfile_hourly

        else:
            nfile_hourly_th = 0
            nimage_hourly_th = 0
            #print(f'sunset, expect {nfile_hourly_th} hourly files')
        nimage_hourly_th = nfile_hourly_th * 144

        date = timeplot.date()
        hour = timeplot.hour
        subdir = date.strftime("%Y/%m/%d/")
        filedir = base_dir + subdir
        if os.path.isdir(filedir):
            datestr = date.strftime("%Y-%m-%d")
            datehrstr = '{0:s}T{1:02d}:00:00'.format(datestr,hour)
            pattern = '*{0:s}*{1:s}T{2:02d}*.{3:s}'.format(filetype, datestr, hour, fileformat)
            # calculate the number of files and images if fully operational
            file_count = 0  # Initialize counter for this subdirectory
            image_count = 0
            # Walk through the current subdirectory
            for root, dirs, files in os.walk(filedir):
                # Check each file against the pattern
                for f in files:
                    if fnmatch.fnmatch(f, pattern):
                        try:
                            if fileformat.lower() == 'hdf':
                                with h5py.File(filedir + f, 'r') as f0:
                                    header = dict(f0['ch_vals'].attrs)
                                    datashape = header['original_shape']
                                    nimages = datashape[1]
                                file_count += 1
                                image_count += nimages
                            elif fileformat.lower() == 'fits':
                                with fits.open(filedir + f) as hdu:
                                    nimages = hdu[0].data.shape[1]
                                file_count += 1
                                image_count += nimages
                            else:
                                print(f'File format {fileformat} not recognized')

                        except Exception as e:
                            print(f"An unexpected error occurred for file {f}: {e}")
                            continue

            data.append((datehrstr, file_count, image_count, nfile_hourly_th, nimage_hourly_th))
            print(f"{datehrstr}: {file_count} files ({nfile_hourly_th} expected), {image_count} images ({nimage_hourly_th} expected)")
            tot_file_count += file_count
            tot_image_count += image_count

    print(f'Total: {tot_file_count} files, {tot_image_count} images')
    if len(data) > 0:
        df = pd.DataFrame(data, columns=["Hour", "Hourly Files", "Hourly Images", "Max Hourly Files", "Max Hourly Images"])
        df = df.sort_values("Hour", ascending=False)
        return df
    else:
        print('No new data generated')
        return pd.DataFrame() 

def update_daily_database(ndays=7, startdate=Time.now(), filetype='fch', base_dir='/nas7/ovro-lwa-data/hdf/slow/lev1/'):
    new_df = scan_pastdays(ndays, startdate=startdate, filetype=filetype, base_dir=base_dir)

    if not new_df.empty:
        # Load existing CSV if exists
        if os.path.exists(CSV_FILE):
            df = pd.read_csv(CSV_FILE)
        else:
            df = new_df

        # Merge: keep all old rows, update rows that overlap
        df = pd.concat([df[~df["Date"].isin(new_df["Date"])], new_df], ignore_index=True)

        # Sort by date (optional)
        df = df.sort_values("Date", ascending=False).reset_index(drop=True)

        # Write back
        df.to_csv(CSV_FILE, index=False)
    else:
        print('Nothing to update')

def update_hourly_database(nhours=3, startdate=Time.now(), filetype='fch', base_dir='/lustre/solarpipe/realtime_pipeline/'):
    new_df = scan_pasthours(nhours, startdate=startdate, filetype=filetype, base_dir=base_dir)

    if not new_df.empty:
        # Load existing CSV if exists
        if os.path.exists(CSV_FILE_HOURLY):
            df = pd.read_csv(CSV_FILE_HOURLY)
        else:
            df = new_df

        # Merge: keep all old rows, update rows that overlap
        df = pd.concat([df[~df["Hour"].isin(new_df["Hour"])], new_df], ignore_index=True)

        # Sort by date (optional)
        df = df.sort_values("Hour", ascending=False).reset_index(drop=True)

        # Write back
        df.to_csv(CSV_FILE_HOURLY, index=False)
    else:
        print('Nothing to update')

def make_plot(ndays=10, alltime=True, figdir='/common/webplots/lwa-data/'):
    df = pd.read_csv(CSV_FILE)
    tot_files = df['Daily Files'].sum()
    tot_images = df['Daily Images'].sum()
    dates = [Time(date).datetime for date in df['Date']]

    df_hr = pd.read_csv(CSV_FILE_HOURLY)
    hours = [Time(t).datetime for t in df_hr['Hour']]
    today = datetime.now() 
    yesterday = datetime.now() - timedelta(days=1)

    # Plot for the past 24 hours
    sunrise, sunset = sun_riseset(today, altitude_limit=15)
    if sunset.datetime > today:
        junk, sunset_yesterday = sun_riseset(yesterday, altitude_limit=15)
        t0 = sunset_yesterday
        t1 = sunrise
    else:
        t0 = sunrise
        t1 = sunset

    if alltime:
        figname = figdir + 'slow_image_number_all'
        # This is for all time
        fig, axs = plt.subplots(2, 2, figsize=(12, 8))
        timestart = today - timedelta(hours=24)
        axs[0,0].step(hours, df_hr["Hourly Files"], where='mid', label='Actual') 
        #axs[0,0].plot(hours, df_hr["Max Hourly Files"], '--k', label='Theoretical')
        axs[0,0].step(hours, df_hr["Max Hourly Files"], where='mid', color='k', ls=':', label='Theoretical')
        axs[0,0].axvline(today.replace(minute=0, second=0, microsecond=0), color='k', ls='--')
        axs[0,0].axvspan(t0.datetime, t1.datetime, alpha=0.5, facecolor='gray')
        axs[0,0].text(today.replace(minute=5, second=0, microsecond=0), 2, 'Hour Now', ha='left', va='bottom') 
        hr_formatter = mdates.DateFormatter('%H')
        axs[0,0].set_ylabel("Hourly # of FITS Files")
        axs[0,0].set_xlabel("Time")
        axs[0,0].set_xlim([timestart, today.replace(minute=0, second=0, microsecond=0) + timedelta(hours=4)])
        axs[0,0].set_title(f"OVRO-LWA # of Fine-Channel FITS Files (Last 24 hrs)")
        axs[0,0].xaxis.set_major_formatter(hr_formatter)
        axs[0,0].grid(True)
        axs[0,0].legend()
        axs[0,0].set_ylim([0, 250])
        axs[0,0].legend(loc='upper left')
        ax2 = axs[0,0].twinx()
        ax2.step(hours, (1.-df_hr["Hourly Images"]/df_hr["Hourly Files"]/144.)*100, where='mid', 
                color='r', label='Chan Loss %') 
        ax2.set_ylim([0, 100])
        ax2.set_ylabel('Channel Loss Percentage')
        ax2.legend(loc='upper right')

        # This is for the past 10 days
        datestart = today - timedelta(days=ndays)
        datestartstr = datestart.strftime('%Y%m%d')
        my_formatter = mdates.DateFormatter('%m/%d')
        axs[0,1].step(dates, df["Daily Files"], where='mid', label='Actual') 
        axs[0,1].step(dates, df["Max Daily Files"], where='mid', color='k', ls=':', label='Theoretical')
        axs[0,1].axvline(today.date(), color='k', ls='--')
        axs[0,1].text(today.replace(hour=1), 20, 'Day Now', ha='left', va='bottom') 
        axs[0,1].set_ylabel("Daily # of HDF Files")
        axs[0,1].set_xlabel("Date")
        axs[0,1].set_xlim([datestart, today.date() + timedelta(days=1)])
        axs[0,1].set_title(f"OVRO-LWA Daily # of Fine-Channel HDF Files (Last {ndays} days)")
        axs[0,1].xaxis.set_major_formatter(my_formatter)
        axs[0,1].grid(True)
        axs[0,1].set_ylim([0, 2500])
        axs[0,1].legend(loc='upper left')

        ax2 = axs[0,1].twinx()
        ax2.step(dates, (1.-df["Daily Images"]/df["Daily Files"]/144.)*100, where='mid', 
                color='r', label='Chan Loss %') 
        ax2.set_ylim([0, 100])
        ax2.set_ylabel('Channel Loss Percentage')
        ax2.legend(loc='upper right')

        #axs[1,1].step(dates, df["Daily Images"], where='mid', label='Actual') 
        #axs[1,1].plot(df["Date"], df["Daily Images"], marker="o", fillstyle='none', linestyle='None')
        #axs[1,1].plot(dates, df["Max Daily Images"], '--k', label='Theoretical')
        #axs[0,1].axhline(y=60000, color='b', linestyle='-', label='SDO/AIA (all bands)')
        #plt.ylim(0, 13.9)
        #axs[1,1].set_ylabel("Daily # of Images")
        #axs[1,1].set_xlabel("Date")
        #axs[1,1].set_xlim([datestart, today])
        #axs[1,1].set_title(f"OVRO-LWA Daily # of Fine-Channel Images (Last {ndays} days)")
        #axs[1,1].xaxis.set_major_formatter(my_formatter)
        #axs[1,1].grid(True)
        #axs[1,1].legend()
        #axs[1,1].set_ylim([0, 350000])


        # This is for all times
        my_formatter = mdates.DateFormatter('%y/%m')
        axs[1,0].step(dates, df["Daily Files"], where='mid', label='Actual') 
        #axs[0,0].plot(dates, df["Daily Files"], marker="o", fillstyle='none', linestyle='None') 
        axs[1,0].plot(dates, df["Max Daily Files"], ':k', label='Theoretical')
        props = dict(boxstyle='square', facecolor='lightblue', alpha=0.5)
        axs[1,0].text(0.98, 0.02, '{0:d} in total'.format(tot_files), bbox=props, 
                transform=axs[1,0].transAxes, ha='right', va='bottom', fontsize=12)
        #plt.ylim(0, 13.9)
        axs[1,0].set_ylabel("Daily # of Files")
        axs[1,0].set_xlabel("Date")
        axs[1,0].set_ylim([0, 2500])
        axs[1,0].set_title(f"OVRO-LWA Daily # of Fine-Channel Files (All Time)")
        axs[1,0].xaxis.set_major_formatter(my_formatter)
        axs[1,0].legend()
        axs[1,0].grid(True)

        axs[1,1].step(dates, df["Daily Images"], where='mid', label='Actual') 
        #axs[0,1].plot(dates, df["Daily Images"], marker="o", fillstyle='none', linestyle='None')
        axs[1,1].plot(dates, df["Max Daily Images"], ':k', label='Theoretical')
        axs[1,1].axhline(y=60000, color='orange', linestyle='-', lw=1., label='SDO/AIA (all bands)')
        axs[1,1].text(0.98, 0.02, '{0:.1f} M in total'.format(tot_images/1e6), bbox=props,
                transform=axs[1,1].transAxes, ha='right', va='bottom', fontsize=12)
        #plt.ylim(0, 13.9)
        axs[1,1].set_ylabel("Daily # of Images")
        axs[1,1].set_xlabel("Date")
        axs[1,1].set_title(f"OVRO-LWA Daily # of Fine-Channel Images (All Time)")
        axs[1,1].set_ylim([0, 350000])
        axs[1,1].xaxis.set_major_formatter(my_formatter)
        axs[1,1].legend()
        axs[1,1].grid(True)
        plt.tight_layout()
    else:
        figname = figdir + 'slow_image_number'
        fig, ax = plt.subplots(1, 1, figsize=(6.5, 4))
        # This is for the past 24 hours 
        timestart = today - timedelta(hours=24)
        ax.step(hours, df_hr["Hourly Files"], where='mid', label='Actual') 
        ax.step(hours, df_hr["Max Hourly Files"], where='mid', color='k', ls=':', label='Theoretical')
        ax.axvline(today.replace(minute=0, second=0, microsecond=0), color='k', ls='--')
        ax.axvspan(t0.datetime, t1.datetime, alpha=0.5, facecolor='gray')
        ax.text(today.replace(minute=5, second=0, microsecond=0), 2, 'Hour Now', ha='left', va='bottom') 
        hr_formatter = mdates.DateFormatter('%H')
        ax.set_ylabel("Hourly # of FITS Files")
        ax.set_xlabel("Time")
        ax.set_xlim([timestart, today + timedelta(hours=5)])
        ax.set_title(f"OVRO-LWA # of Fine-Channel FITS Files (Last 24 hrs)")
        ax.xaxis.set_major_formatter(hr_formatter)
        ax.grid(True)
        ax.legend()
        ax.set_ylim([0, 250])
        ax.legend(loc='upper left')
        ax2 = ax.twinx()
        ax2.step(hours, (1.-df_hr["Hourly Images"]/df_hr["Hourly Files"]/144.)*100, where='mid', 
                color='r', label='Chan Loss %') 
        ax2.set_ylim([0, 100])
        ax2.set_ylabel('Channel Loss Percentage')
        ax2.legend(loc='upper right')
        plt.tight_layout()

    if os.path.exists(f'{figname}_latest.png'):
        os.system(f'cp {figname}_latest.png {figname}_later.png')
    plt.savefig(f'{figname}_latest.png')
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Scan LWA level 1 image files and update database')
    parser.add_argument('--scandays', default=False, action='store_true', help='If set, scan files daily')
    parser.add_argument('--scanhours', default=False, action='store_true', help='If set, scan files hourly')
    parser.add_argument('--ndays', default=3, help='Number of days to scan')
    parser.add_argument('--nhours', default=2, help='Number of hours to scan')
    parser.add_argument('--doplot', default=False, action='store_true', help='Whether or not to make plot')
    parser.add_argument('--alltime', default=False, action='store_true', help='Whether or not to make plot')
    args = parser.parse_args()
    if args.scandays:
        update_daily_database(ndays=int(args.ndays))
    if args.scanhours:
        update_hourly_database(nhours=int(args.nhours))
    if args.doplot:
        make_plot(alltime=False)
        make_plot(alltime=True)
