import solar_realtime_pipeline as srp
import argparse
import datetime

yesterday_date_str_in_utc = (datetime.datetime.now(datetime.timezone.utc)
                     - datetime.timedelta(days=1)).strftime('%Y-%m-%d')

def main():
    parser = argparse.ArgumentParser(description='Calculate the daily refraction correction for the LWA')

    # data not provided, default to yesterday
    # not mandatory to provide date
    parser.add_argument('--date', type=str, default=yesterday_date_str_in_utc, help='Date in UTC format (YYYY-MM-DD)')

    args = parser.parse_args()

    srp.daily_refra_correction(args.date)

if __name__ == '__main__':
    main()
